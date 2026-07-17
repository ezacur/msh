/* IntersectSurfaceRay_mx -- fast drop-in core of IntersectSurfaceRay.
 *
 *   [xyz,P_id,cell_id,t,ray_id] = IntersectSurfaceRay_mx( P , ray , MODE )
 *
 *   SAME SYNTAX as IntersectSurfaceRay.m; the ONLY extension is that several
 *   rays can be passed at once:
 *
 *     P      a mesh struct with .vertices/.faces or .xyz/.tri (double), a
 *            struct array, or a cell array of such structs. (Graphics handles
 *            and the interactive empty-ray form are NOT supported here: use
 *            IntersectSurfaceRay.m, which resolves them and delegates.)
 *     ray    2 x 3       [p0;p1], one LINE through p0 -> p1
 *            2 x 3 x N   N rays, one [p0;p1] page each          (EXTENSION)
 *            N x 6       N rays, one [p0,p1] row each           (EXTENSION)
 *            Parametrization: hit = p0 + t*(p1-p0); t UNBOUNDED except 'any'.
 *     MODE   'first' (default) | 'last' | 'all' | 'any'
 *
 *            'any' is an OCCLUSION (shadow-segment) query: per ray it reports
 *            SOME hit with  1e-9 < t < 1-1e-5  (guard bands sized so that p1
 *            may lie exactly ON the surface without self-reporting) and stops
 *            at the first one found -- unordered traversal with early exit,
 *            substantially cheaper than 'first'. Miss conventions as 'first'.
 *
 *     xyz,P_id,cell_id : one row per ray ('first'/'last'/'any'; misses are
 *                        NaN/0/0) or one row per hit ('all', sorted by ray,t).
 *     t,ray_id         : optional extras aligned with the rows above (t = ray
 *                        parameter of each hit; ray_id = which ray it belongs
 *                        to -- needed to read 'all' with several rays).
 *
 *   INTERNALS (single-threaded on purpose, no OpenMP; SIMD is used -- data
 *   parallelism inside the one thread):
 *     - per-mesh acceleration CACHED across calls in a small LRU, keyed by a
 *       4-lane 64-bit fingerprint of the raw vertex/face bytes: a cache hit
 *       reuses everything (no re-copy, no re-validation of the mesh);
 *     - binned SAH build (16 bins on the largest centroid-extent axis) with
 *       ADAPTIVE leaves: a range becomes a leaf when the SAH says splitting
 *       does not pay (min 4, max 16 triangles) -- larger leaves where the
 *       geometry is dense, so the 4-wide triangle kernel stays busy;
 *     - 4-ary BVH (BVH4): each 128-byte node holds the conservatively-rounded
 *       FLOAT bounds of its up-to-4 children in SoA layout, so ONE SSE pass
 *       slab-tests all four children (near/far bound picked per axis by the
 *       ray's direction sign -- branchless); rays with an exactly-zero
 *       component take a careful scalar containment variant instead (the
 *       arithmetic path would produce 0*Inf = NaN);
 *     - ORDERED traversal (children sorted by entry distance, nearest first)
 *       pruned by the current best t ('first'/'last'); 'any' traverses
 *       unordered, clips nodes against the [0,1] window and EXITS on the
 *       first accepting hit;
 *     - triangles PREPROCESSED as (v0,e1,e2) in 4-wide SoA blocks (PreTri4),
 *       leaf ranges contiguous and padded with degenerate lanes: the
 *       Moller-Trumbore kernel tests 4 triangles per iteration in AVX DOUBLE
 *       precision (no accuracy loss; tiny inclusive tolerance so rays through
 *       shared vertices/edges cannot slip between the incident triangles).
 *       AVX use is decided ONCE at runtime via CPUID (+OS ymm support); a
 *       scalar lane-loop fallback keeps pre-2011 CPUs working;
 *     - small workloads (few rays, or small meshes) take the brute path: the
 *       PreTri4 blocks are gathered ONCE per call and scanned linearly --
 *       never any hashing nor tree build;
 *     - PACKETS: for 'first'/'any' tree queries with k >= 16 rays (and AVX),
 *       rays are traced in bundles of 8 that SHARE the tree walk: each child
 *       box is slab-tested against the 8 rays in one AVX pass and carried with
 *       per-lane masks, amortizing the node fetches across the packet
 *       (coherent batches -- e.g. silhouette occlusion -- benefit the most;
 *       results are identical to single-ray tracing);
 *     - ROOT-PCA: at build time the mesh's principal axes are found (Jacobi)
 *       and, when the anisotropy is >= 2, vertices are rotated ONCE into that
 *       frame and each ray once at query time -- every node of the tree gets
 *       tighter for elongated DIAGONAL geometry at zero per-node cost
 *       (blob-like meshes keep the identity frame and today's behavior).
 *
 *   Compile:  mex IntersectSurfaceRay_mx.cpp
 */

#include "mex.h"
#include <math.h>
#include <float.h>
#include <string.h>
#include <ctype.h>
#include <vector>
#include <algorithm>
#include <immintrin.h>
#if defined(_MSC_VER)
#include <intrin.h>
#endif

#ifndef INFINITY
#define INFINITY (mxGetInf())
#endif

typedef unsigned long long u64;

#define MIN_LEAF    4         /* n <= MIN_LEAF is always a leaf                */
#define MAX_LEAF    16        /* SAH may grow a leaf up to this many tris      */
#define SAH_BINS    16
#define SAH_CTRAV   1.0       /* node-visit cost in SAH units                  */
#define SAH_CTRI    0.4       /* per-triangle cost (cheap: 4-wide MT kernel)   */
#define CACHE_SLOTS 4
#define UV_EPS      1e-9      /* inclusive barycentric tolerance (see header)  */
#define USE_PACKETS 0         /* 8-ray packet traversal: implemented and kept
                               * for reference, but MEASURED SLOWER than the
                               * single-ray path on this codebase (v2/v3 = 0.82
                               * even on coherent occlusion rays): the
                               * single-ray kernel is ALREADY 4-wide SIMD per
                               * node (BVH4) and sequential rays reuse the tree
                               * from L2, so packets only add mask bookkeeping
                               * and a weaker per-ray best-t prune. Flip to 1
                               * to re-evaluate on other hardware.             */
#define ANY_TLO     1e-9      /* 'any' accepts hits with ANY_TLO < t < ANY_THI */
#define ANY_THI     (1.0 - 1e-5)

/* ----------------------------- fast fingerprint ---------------------------- */
static u64 hashBytes( const void *data , size_t nbytes , u64 h )
{
  const u64 *w = (const u64 *)data;
  size_t nw = nbytes >> 3;
  u64 h1 = h ^ 0x9E3779B97F4A7C15ULL, h2 = h ^ 0xC2B2AE3D27D4EB4FULL,
      h3 = h ^ 0x165667B19E3779F9ULL, h4 = h ^ 0x27D4EB2F165667C5ULL;
  size_t i = 0;
  for ( ; i + 4 <= nw ; i += 4 ) {                      /* 4 independent lanes */
    h1 = ( h1 ^ w[i+0] ) * 0x9E3779B97F4A7C15ULL;
    h2 = ( h2 ^ w[i+1] ) * 0xC2B2AE3D27D4EB4FULL;
    h3 = ( h3 ^ w[i+2] ) * 0x165667B19E3779F9ULL;
    h4 = ( h4 ^ w[i+3] ) * 0x27D4EB2F165667C5ULL;
  }
  u64 hh = ( ( h1 ^ ( h2 >> 31 ) ) * 0x9E3779B97F4A7C15ULL )
         ^ ( ( h3 ^ ( h4 >> 29 ) ) * 0xC2B2AE3D27D4EB4FULL );
  for ( ; i < nw ; i++ ) hh = ( hh ^ w[i] ) * 0x100000001B3ULL;
  const unsigned char *b = (const unsigned char *)data + ( nw << 3 );
  for ( size_t j = 0 ; j < ( nbytes & 7 ) ; j++ ) hh = ( hh ^ b[j] ) * 0x100000001B3ULL;
  hh ^= hh >> 33;  hh *= 0xFF51AFD7ED558CCDULL;  hh ^= hh >> 33;
  return hh;
}

/* ---------------------- runtime AVX detection (once) ------------------------ */
static int g_avx = -1;
static bool useAVX( void )
{
  if ( g_avx < 0 ) {
#if defined(_MSC_VER)
    int ci[4];  __cpuid( ci , 1 );
    bool osxsave = ( ci[2] >> 27 ) & 1;
    bool avx     = ( ci[2] >> 28 ) & 1;
    g_avx = ( osxsave && avx && ( ( _xgetbv(0) & 6 ) == 6 ) ) ? 1 : 0;
#else
    g_avx = __builtin_cpu_supports("avx") ? 1 : 0;
#endif
  }
  return g_avx == 1;
}

/* --------------- 4-wide preprocessed triangle block (SoA) ------------------- */
struct alignas(32) PreTri4 {           /* Moller-Trumbore-ready, leaf-contiguous */
  double v0x[4], v0y[4], v0z[4];
  double e1x[4], e1y[4], e1z[4];
  double e2x[4], e2y[4], e2z[4];
  int    id[4];                        /* ORIGINAL global tri row (0-based); -1 = padding */
};

/* ------- BVH4 node: SoA float bounds of the 4 children, one struct ---------- */
/* conservative rounding; 96 + 32 = 128 bytes = two cache lines, four children  */
struct alignas(64) Node4 {
  float bminx[4], bmaxx[4];
  float bminy[4], bmaxy[4];
  float bminz[4], bmaxz[4];
  int   child[4];                      /* node index, or first PreTri4 block    */
  int   count[4];                      /* -1 empty ; 0 internal ; >0 leaf tris  */
};

struct ChildRef { int idx; int count; double bmin[3], bmax[3]; };

/* ------------------------------- ray precompute ----------------------------- */
struct Ray {
  double o[3], d[3];
  float  of[3], invf[3];
  int    sgn[3];                       /* d < 0                                 */
  bool   anyZero;                      /* some d[a] == 0: use the careful slab  */
};

static void makeRay( const double *p0 , const double *p1 , Ray &r )
{
  r.anyZero = false;
  for ( int a = 0 ; a < 3 ; a++ ) {
    r.o[a] = p0[a];  r.d[a] = p1[a] - p0[a];
    r.of[a]   = (float)r.o[a];
    r.invf[a] = (float)( 1.0 / r.d[a] );
    r.sgn[a]  = ( r.d[a] < 0.0 );
    if ( r.d[a] == 0.0 ) r.anyZero = true;
  }
}

/* ray into the tree's frame: applies the cached root-PCA rotation when active.
 * t is invariant under the rigid transform, so outputs need no un-mapping
 * (hit xyz is reconstructed from the ORIGINAL ray in the gateway).             */
struct TreeCache;
static void makeRayTC( const TreeCache *tc , const double *p0 , const double *p1 , Ray &r );

/* --------------------- 4-wide slab test of one Node4 ------------------------ */
/* branchless: near/far bound array picked per axis by the ray's sign (the sign
 * is a RAY property, so it is the same for all four children -> plain loads). */
static inline void slab4Fast( const Ray &r , const Node4 &nd , float *tn4 , float *tf4 )
{
  __m128 tn = _mm_set1_ps( -FLT_MAX ), tf = _mm_set1_ps( FLT_MAX );
  {
    const float *lo = r.sgn[0] ? nd.bmaxx : nd.bminx;
    const float *hi = r.sgn[0] ? nd.bminx : nd.bmaxx;
    __m128 o = _mm_set1_ps( r.of[0] ), iv = _mm_set1_ps( r.invf[0] );
    tn = _mm_max_ps( tn , _mm_mul_ps( _mm_sub_ps( _mm_loadu_ps(lo) , o ) , iv ) );
    tf = _mm_min_ps( tf , _mm_mul_ps( _mm_sub_ps( _mm_loadu_ps(hi) , o ) , iv ) );
  }
  {
    const float *lo = r.sgn[1] ? nd.bmaxy : nd.bminy;
    const float *hi = r.sgn[1] ? nd.bminy : nd.bmaxy;
    __m128 o = _mm_set1_ps( r.of[1] ), iv = _mm_set1_ps( r.invf[1] );
    tn = _mm_max_ps( tn , _mm_mul_ps( _mm_sub_ps( _mm_loadu_ps(lo) , o ) , iv ) );
    tf = _mm_min_ps( tf , _mm_mul_ps( _mm_sub_ps( _mm_loadu_ps(hi) , o ) , iv ) );
  }
  {
    const float *lo = r.sgn[2] ? nd.bmaxz : nd.bminz;
    const float *hi = r.sgn[2] ? nd.bminz : nd.bmaxz;
    __m128 o = _mm_set1_ps( r.of[2] ), iv = _mm_set1_ps( r.invf[2] );
    tn = _mm_max_ps( tn , _mm_mul_ps( _mm_sub_ps( _mm_loadu_ps(lo) , o ) , iv ) );
    tf = _mm_min_ps( tf , _mm_mul_ps( _mm_sub_ps( _mm_loadu_ps(hi) , o ) , iv ) );
  }
  _mm_storeu_ps( tn4 , tn );
  _mm_storeu_ps( tf4 , tf );
}

/* careful scalar variant for rays with zero direction components               */
static inline void slab4Zero( const Ray &r , const Node4 &nd , float *tn4 , float *tf4 )
{
  for ( int i = 0 ; i < 4 ; i++ ) {
    float tn = -FLT_MAX, tf = FLT_MAX;
    const float *bmin[3] = { nd.bminx , nd.bminy , nd.bminz };
    const float *bmax[3] = { nd.bmaxx , nd.bmaxy , nd.bmaxz };
    bool out = false;
    for ( int a = 0 ; a < 3 ; a++ ) {
      if ( r.d[a] == 0.0 ) {
        if ( r.of[a] < bmin[a][i] || r.of[a] > bmax[a][i] ) { out = true; break; }
      } else {
        float b0 = r.sgn[a] ? bmax[a][i] : bmin[a][i];
        float b1 = r.sgn[a] ? bmin[a][i] : bmax[a][i];
        float t0 = ( b0 - r.of[a] ) * r.invf[a];
        float t1 = ( b1 - r.of[a] ) * r.invf[a];
        tn = ( t0 > tn ) ? t0 : tn;
        tf = ( t1 < tf ) ? t1 : tf;
      }
    }
    if ( out ) { tn4[i] = 1.0f; tf4[i] = -1.0f; }        /* forced miss          */
    else       { tn4[i] = tn;   tf4[i] = tf;   }
  }
}

static inline void slab4( const Ray &r , const Node4 &nd , float *tn4 , float *tf4 )
{
  if ( r.anyZero ) slab4Zero( r , nd , tn4 , tf4 );
  else             slab4Fast( r , nd , tn4 , tf4 );
}

/* --------------- Moller-Trumbore, one lane of a PreTri4 (scalar) ------------ */
static inline bool mtLane( const Ray &r , const PreTri4 &T , int l , double *tOut )
{
  const double *d = r.d, *o = r.o;
  double e1[3] = { T.e1x[l] , T.e1y[l] , T.e1z[l] };
  double e2[3] = { T.e2x[l] , T.e2y[l] , T.e2z[l] };
  double p[3] = { d[1]*e2[2]-d[2]*e2[1] , d[2]*e2[0]-d[0]*e2[2] , d[0]*e2[1]-d[1]*e2[0] };
  double det  = e1[0]*p[0] + e1[1]*p[1] + e1[2]*p[2];
  if ( fabs(det) < 1e-300 ) return false;
  double inv = 1.0/det;
  double s[3] = { o[0]-T.v0x[l] , o[1]-T.v0y[l] , o[2]-T.v0z[l] };
  double u = ( s[0]*p[0] + s[1]*p[1] + s[2]*p[2] )*inv;
  if ( u < -UV_EPS || u > 1.0+UV_EPS ) return false;
  double q[3] = { s[1]*e1[2]-s[2]*e1[1] , s[2]*e1[0]-s[0]*e1[2] , s[0]*e1[1]-s[1]*e1[0] };
  double v = ( d[0]*q[0] + d[1]*q[1] + d[2]*q[2] )*inv;
  if ( v < -UV_EPS || u+v > 1.0+UV_EPS ) return false;
  *tOut = ( e2[0]*q[0] + e2[1]*q[1] + e2[2]*q[2] )*inv;
  return true;
}

/* --------------- Moller-Trumbore, 4 lanes at once (AVX double) --------------
 * returns a 4-bit mask of valid lanes and fills t4; padding lanes (e1=e2=0)
 * yield det = 0 and reject themselves. Same math, same UV_EPS: bit-identical
 * results to the scalar lane (AVX double = plain double arithmetic).          */
struct RayA {                          /* per-ray broadcasts, made once         */
  __m256d dx, dy, dz, ox, oy, oz;
};
static inline void makeRayA( const Ray &r , RayA &a )
{
  a.dx = _mm256_set1_pd( r.d[0] );  a.dy = _mm256_set1_pd( r.d[1] );  a.dz = _mm256_set1_pd( r.d[2] );
  a.ox = _mm256_set1_pd( r.o[0] );  a.oy = _mm256_set1_pd( r.o[1] );  a.oz = _mm256_set1_pd( r.o[2] );
}

static inline int mt4( const RayA &r , const PreTri4 &T , double *t4 )
{
  __m256d e1x = _mm256_loadu_pd( T.e1x ), e1y = _mm256_loadu_pd( T.e1y ), e1z = _mm256_loadu_pd( T.e1z );
  __m256d e2x = _mm256_loadu_pd( T.e2x ), e2y = _mm256_loadu_pd( T.e2y ), e2z = _mm256_loadu_pd( T.e2z );

  __m256d px = _mm256_sub_pd( _mm256_mul_pd( r.dy , e2z ) , _mm256_mul_pd( r.dz , e2y ) );
  __m256d py = _mm256_sub_pd( _mm256_mul_pd( r.dz , e2x ) , _mm256_mul_pd( r.dx , e2z ) );
  __m256d pz = _mm256_sub_pd( _mm256_mul_pd( r.dx , e2y ) , _mm256_mul_pd( r.dy , e2x ) );

  __m256d det = _mm256_add_pd( _mm256_add_pd(
                  _mm256_mul_pd( e1x , px ) , _mm256_mul_pd( e1y , py ) ) ,
                  _mm256_mul_pd( e1z , pz ) );

  const __m256d absmask = _mm256_castsi256_pd( _mm256_set1_epi64x( 0x7FFFFFFFFFFFFFFFLL ) );
  __m256d ok = _mm256_cmp_pd( _mm256_and_pd( det , absmask ) ,
                              _mm256_set1_pd( 1e-300 ) , _CMP_GE_OQ );
  if ( !_mm256_movemask_pd( ok ) ) return 0;

  __m256d inv = _mm256_div_pd( _mm256_set1_pd( 1.0 ) , det );

  __m256d sx = _mm256_sub_pd( r.ox , _mm256_loadu_pd( T.v0x ) );
  __m256d sy = _mm256_sub_pd( r.oy , _mm256_loadu_pd( T.v0y ) );
  __m256d sz = _mm256_sub_pd( r.oz , _mm256_loadu_pd( T.v0z ) );

  __m256d u = _mm256_mul_pd( _mm256_add_pd( _mm256_add_pd(
                _mm256_mul_pd( sx , px ) , _mm256_mul_pd( sy , py ) ) ,
                _mm256_mul_pd( sz , pz ) ) , inv );

  const __m256d lo = _mm256_set1_pd( -UV_EPS ), hi = _mm256_set1_pd( 1.0+UV_EPS );
  ok = _mm256_and_pd( ok , _mm256_cmp_pd( u , lo , _CMP_GE_OQ ) );
  ok = _mm256_and_pd( ok , _mm256_cmp_pd( u , hi , _CMP_LE_OQ ) );
  if ( !_mm256_movemask_pd( ok ) ) return 0;

  __m256d qx = _mm256_sub_pd( _mm256_mul_pd( sy , e1z ) , _mm256_mul_pd( sz , e1y ) );
  __m256d qy = _mm256_sub_pd( _mm256_mul_pd( sz , e1x ) , _mm256_mul_pd( sx , e1z ) );
  __m256d qz = _mm256_sub_pd( _mm256_mul_pd( sx , e1y ) , _mm256_mul_pd( sy , e1x ) );

  __m256d v = _mm256_mul_pd( _mm256_add_pd( _mm256_add_pd(
                _mm256_mul_pd( r.dx , qx ) , _mm256_mul_pd( r.dy , qy ) ) ,
                _mm256_mul_pd( r.dz , qz ) ) , inv );

  ok = _mm256_and_pd( ok , _mm256_cmp_pd( v , lo , _CMP_GE_OQ ) );
  ok = _mm256_and_pd( ok , _mm256_cmp_pd( _mm256_add_pd( u , v ) , hi , _CMP_LE_OQ ) );
  int m = _mm256_movemask_pd( ok );
  if ( !m ) return 0;

  __m256d t = _mm256_mul_pd( _mm256_add_pd( _mm256_add_pd(
                _mm256_mul_pd( e2x , qx ) , _mm256_mul_pd( e2y , qy ) ) ,
                _mm256_mul_pd( e2z , qz ) ) , inv );
  _mm256_storeu_pd( t4 , t );
  return m;
}

/* ------------------------------ tree cache ---------------------------------- */
struct TreeCache {
  u64                  key;
  mwSize               nV, nF;
  unsigned             age;
  ChildRef             root;
  bool                 pca;            /* root-PCA frame active (anisotropy>=2) */
  double               Rr[9];          /* world->local rows (v' = Rr*(v-Cc))    */
  double               Cc[3];
  std::vector<Node4>   nodes;
  std::vector<PreTri4> tris4;          /* leaf order, 4-wide padded blocks      */
  std::vector<int>     surfId;         /* per ORIGINAL tri: 1-based surface     */
  std::vector<int>     locId;          /* per ORIGINAL tri: 1-based row         */
};
static TreeCache g_cache[ CACHE_SLOTS ];
static unsigned  g_clock = 0;

static void makeRayTC( const TreeCache *tc , const double *p0 , const double *p1 , Ray &r )
{
  if ( tc && tc->pca ) {
    double q0[3], q1[3];
    for ( int a = 0 ; a < 3 ; a++ ) {
      double d00 = p0[0]-tc->Cc[0], d01 = p0[1]-tc->Cc[1], d02 = p0[2]-tc->Cc[2];
      double d10 = p1[0]-tc->Cc[0], d11 = p1[1]-tc->Cc[1], d12 = p1[2]-tc->Cc[2];
      q0[a] = tc->Rr[3*a]*d00 + tc->Rr[3*a+1]*d01 + tc->Rr[3*a+2]*d02;
      q1[a] = tc->Rr[3*a]*d10 + tc->Rr[3*a+1]*d11 + tc->Rr[3*a+2]*d12;
    }
    makeRay( q0 , q1 , r );
  } else {
    makeRay( p0 , p1 , r );
  }
}

/* ------------------------------- BVH build ---------------------------------- */
struct BuildCtx {
  const double *V;  const int *T;
  std::vector<double>   bb, cc;        /* per-tri AABB (6) and centroid (3)     */
  std::vector<int>      ord;
  std::vector<Node4>    *nodes;
  std::vector<PreTri4>  *tris4;
};

static void rangeBounds( BuildCtx &B , int lo , int hi , double bmin[3] , double bmax[3] )
{
  for ( int c = 0 ; c < 3 ; c++ ) { bmin[c] = INFINITY; bmax[c] = -INFINITY; }
  for ( int i = lo ; i < hi ; i++ ) {
    const double *a = &B.bb[ 6*B.ord[i] ];
    for ( int c = 0 ; c < 3 ; c++ ) {
      bmin[c] = fmin( bmin[c] , a[c]   );
      bmax[c] = fmax( bmax[c] , a[3+c] );
    }
  }
}

static double halfArea( const double bmin[3] , const double bmax[3] )
{
  double e0 = bmax[0]-bmin[0], e1 = bmax[1]-bmin[1], e2 = bmax[2]-bmin[2];
  if ( e0 < 0 || e1 < 0 || e2 < 0 ) return 0.0;
  return e0*e1 + e1*e2 + e0*e2;
}

/* leaf-vs-split decision for the range [lo,hi):
 *   true  -> make it a LEAF;
 *   false -> *mid receives the split (binned SAH on the largest centroid-extent
 *            axis; median fallback when the centroids are degenerate).         */
static bool wantsLeaf( BuildCtx &B , int lo , int hi , int *mid )
{
  int n = hi - lo;
  if ( n <= MIN_LEAF ) return true;

  double cmin[3] = {  INFINITY,  INFINITY,  INFINITY };
  double cmax[3] = { -INFINITY, -INFINITY, -INFINITY };
  for ( int i = lo ; i < hi ; i++ ) {
    const double *ctr = &B.cc[ 3*B.ord[i] ];
    for ( int c = 0 ; c < 3 ; c++ ) {
      cmin[c] = fmin( cmin[c] , ctr[c] );
      cmax[c] = fmax( cmax[c] , ctr[c] );
    }
  }
  int axis = 0; double ext = cmax[0]-cmin[0];
  if ( cmax[1]-cmin[1] > ext ) { axis = 1; ext = cmax[1]-cmin[1]; }
  if ( cmax[2]-cmin[2] > ext ) { axis = 2; ext = cmax[2]-cmin[2]; }

  if ( ext <= 0.0 ) {                                   /* coincident centroids  */
    if ( n <= MAX_LEAF ) return true;
    *mid = lo + n/2;                                    /* arbitrary halves      */
    return false;
  }

  /* ---- binned SAH on `axis` ---- */
  int    binN [SAH_BINS];
  double binMn[SAH_BINS][3], binMx[SAH_BINS][3];
  for ( int b = 0 ; b < SAH_BINS ; b++ ) {
    binN[b] = 0;
    for ( int c = 0 ; c < 3 ; c++ ) { binMn[b][c] = INFINITY; binMx[b][c] = -INFINITY; }
  }
  double scale = SAH_BINS / ext;
  for ( int i = lo ; i < hi ; i++ ) {
    int f = B.ord[i];
    int b = (int)( ( B.cc[3*f+axis] - cmin[axis] ) * scale );
    if ( b < 0 ) b = 0;  if ( b >= SAH_BINS ) b = SAH_BINS-1;
    binN[b]++;
    const double *a = &B.bb[6*f];
    for ( int c = 0 ; c < 3 ; c++ ) {
      binMn[b][c] = fmin( binMn[b][c] , a[c]   );
      binMx[b][c] = fmax( binMx[b][c] , a[3+c] );
    }
  }

  /* suffix sweep: area/count of everything right of each boundary              */
  double rArea[SAH_BINS];  int rCnt[SAH_BINS];
  {
    double m[3] = {  INFINITY,  INFINITY,  INFINITY };
    double M[3] = { -INFINITY, -INFINITY, -INFINITY };
    int cnt = 0;
    for ( int b = SAH_BINS-1 ; b >= 1 ; b-- ) {
      for ( int c = 0 ; c < 3 ; c++ ) {
        m[c] = fmin( m[c] , binMn[b][c] );  M[c] = fmax( M[c] , binMx[b][c] );
      }
      cnt += binN[b];
      rArea[b] = ( cnt > 0 ) ? halfArea( m , M ) : 0.0;
      rCnt [b] = cnt;
    }
  }
  /* prefix sweep + best boundary                                                */
  double bestCost = INFINITY;  int bestB = -1;
  {
    double m[3] = {  INFINITY,  INFINITY,  INFINITY };
    double M[3] = { -INFINITY, -INFINITY, -INFINITY };
    int cnt = 0;
    for ( int b = 0 ; b < SAH_BINS-1 ; b++ ) {
      for ( int c = 0 ; c < 3 ; c++ ) {
        m[c] = fmin( m[c] , binMn[b][c] );  M[c] = fmax( M[c] , binMx[b][c] );
      }
      cnt += binN[b];
      if ( cnt == 0 || rCnt[b+1] == 0 ) continue;
      double cost = halfArea( m , M )*cnt + rArea[b+1]*rCnt[b+1];
      if ( cost < bestCost ) { bestCost = cost; bestB = b; }
    }
  }
  if ( bestB < 0 ) {                                    /* everything in one bin */
    if ( n <= MAX_LEAF ) return true;
    *mid = lo + n/2;
    std::nth_element( B.ord.begin()+lo , B.ord.begin()+*mid , B.ord.begin()+hi ,
      [&]( int a , int b ){ return B.cc[3*a+axis] < B.cc[3*b+axis]; } );
    return false;
  }

  /* SAH says leaf? cost in per-tri units, normalized by this node's area       */
  double bmin[3], bmax[3];
  rangeBounds( B , lo , hi , bmin , bmax );
  double A = halfArea( bmin , bmax );
  if ( A > 0.0 && n <= MAX_LEAF ) {
    double splitCost = SAH_CTRAV + SAH_CTRI * bestCost / A;
    if ( splitCost >= SAH_CTRI * n ) return true;
  }

  /* partition by the chosen bin boundary                                        */
  double cut = cmin[axis] + ( bestB + 1 ) / scale;
  int *first = &B.ord[lo], *last = &B.ord[hi-1] + 1;
  int *pm = std::partition( first , last ,
              [&]( int f ){ return B.cc[3*f+axis] < cut; } );
  *mid = lo + (int)( pm - first );
  if ( *mid == lo || *mid == hi ) {                     /* numeric edge: median  */
    *mid = lo + n/2;
    std::nth_element( B.ord.begin()+lo , B.ord.begin()+*mid , B.ord.begin()+hi ,
      [&]( int a , int b ){ return B.cc[3*a+axis] < B.cc[3*b+axis]; } );
  }
  return false;
}

static void emitLeaf( BuildCtx &B , int lo , int hi , ChildRef &ref )
{
  int n = hi - lo;
  ref.idx   = (int)B.tris4->size();
  ref.count = n;
  for ( int base = 0 ; base < n ; base += 4 ) {
    PreTri4 blk;
    for ( int l = 0 ; l < 4 ; l++ ) {
      if ( base + l < n ) {
        int f = B.ord[ lo + base + l ];
        double x0[3];
        for ( int c = 0 ; c < 3 ; c++ ) x0[c] = B.V[ 3*B.T[3*f+0] + c ];
        blk.v0x[l] = x0[0];  blk.v0y[l] = x0[1];  blk.v0z[l] = x0[2];
        blk.e1x[l] = B.V[ 3*B.T[3*f+1] + 0 ] - x0[0];
        blk.e1y[l] = B.V[ 3*B.T[3*f+1] + 1 ] - x0[1];
        blk.e1z[l] = B.V[ 3*B.T[3*f+1] + 2 ] - x0[2];
        blk.e2x[l] = B.V[ 3*B.T[3*f+2] + 0 ] - x0[0];
        blk.e2y[l] = B.V[ 3*B.T[3*f+2] + 1 ] - x0[1];
        blk.e2z[l] = B.V[ 3*B.T[3*f+2] + 2 ] - x0[2];
        blk.id[l]  = f;
      } else {                                          /* degenerate padding    */
        blk.v0x[l]=blk.v0y[l]=blk.v0z[l] = 0.0;
        blk.e1x[l]=blk.e1y[l]=blk.e1z[l] = 0.0;         /* det = 0 -> never hits */
        blk.e2x[l]=blk.e2y[l]=blk.e2z[l] = 0.0;
        blk.id[l]  = -1;
      }
    }
    B.tris4->push_back( blk );
  }
}

/* 4-ary build: split the range once, then split each half again (when they do
 * not want to be leaves) -> up to four grandchildren become this node's slots. */
static ChildRef build4( BuildCtx &B , int lo , int hi )
{
  ChildRef ref;
  rangeBounds( B , lo , hi , ref.bmin , ref.bmax );

  int mid;
  if ( wantsLeaf( B , lo , hi , &mid ) ) {
    emitLeaf( B , lo , hi , ref );
    return ref;
  }

  int R[5];  int nR = 0;                                /* up to 4 sub-ranges    */
  int half[3] = { lo , mid , hi };
  R[nR++] = lo;
  for ( int h = 0 ; h < 2 ; h++ ) {
    int a = half[h], b = half[h+1], m2;
    if ( !wantsLeaf( B , a , b , &m2 ) ) R[nR++] = m2;  /* expand this half      */
    R[nR++] = b;
  }
  /* R holds nR+... boundaries: R[0..nR], ranges are (R[i],R[i+1])              */

  Node4 nd;
  int slotChild[4], slotCount[4];
  ChildRef cr[4];
  int nc = nR - 1;                                      /* 2..4 children         */
  for ( int i = 0 ; i < 4 ; i++ ) {
    if ( i < nc ) {
      cr[i] = build4( B , R[i] , R[i+1] );
      slotChild[i] = cr[i].idx;
      slotCount[i] = cr[i].count;
    } else {
      slotChild[i] = -1;  slotCount[i] = -1;            /* empty slot            */
      for ( int c = 0 ; c < 3 ; c++ ) { cr[i].bmin[c] = INFINITY; cr[i].bmax[c] = -INFINITY; }
    }
  }
  for ( int i = 0 ; i < 4 ; i++ ) {                     /* conservative floats   */
    float mnx, mny, mnz, mxx, mxy, mxz;
    if ( slotCount[i] < 0 ) {                           /* empty: self-missing   */
      mnx = mny = mnz =  FLT_MAX;
      mxx = mxy = mxz = -FLT_MAX;
    } else {
      mnx = nextafterf( (float)cr[i].bmin[0] , -FLT_MAX );
      mny = nextafterf( (float)cr[i].bmin[1] , -FLT_MAX );
      mnz = nextafterf( (float)cr[i].bmin[2] , -FLT_MAX );
      mxx = nextafterf( (float)cr[i].bmax[0] ,  FLT_MAX );
      mxy = nextafterf( (float)cr[i].bmax[1] ,  FLT_MAX );
      mxz = nextafterf( (float)cr[i].bmax[2] ,  FLT_MAX );
    }
    nd.bminx[i] = mnx;  nd.bmaxx[i] = mxx;
    nd.bminy[i] = mny;  nd.bmaxy[i] = mxy;
    nd.bminz[i] = mnz;  nd.bmaxz[i] = mxz;
    nd.child[i] = slotChild[i];
    nd.count[i] = slotCount[i];
  }

  ref.idx   = (int)B.nodes->size();
  ref.count = 0;
  B.nodes->push_back( nd );
  return ref;
}

/* symmetric 3x3 eigendecomposition by cyclic Jacobi (A = Q diag(w) Q^T)        */
static void jacobi3( double A[3][3] , double w[3] , double Q[3][3] )
{
  for ( int i = 0 ; i < 3 ; i++ )
    for ( int j = 0 ; j < 3 ; j++ ) Q[i][j] = ( i == j ) ? 1.0 : 0.0;
  for ( int sweep = 0 ; sweep < 32 ; sweep++ ) {
    double off = fabs(A[0][1]) + fabs(A[0][2]) + fabs(A[1][2]);
    if ( off < 1e-300 ) break;
    for ( int p = 0 ; p < 2 ; p++ ) for ( int q = p+1 ; q < 3 ; q++ ) {
      if ( fabs( A[p][q] ) < 1e-300 ) continue;
      double th = 0.5*( A[q][q] - A[p][p] )/A[p][q];
      double t  = ( th >= 0 ? 1.0 : -1.0 )/( fabs(th) + sqrt(th*th+1.0) );
      double c  = 1.0/sqrt(t*t+1.0), s = t*c;
      for ( int r = 0 ; r < 3 ; r++ ) {
        double arp = A[r][p], arq = A[r][q];
        A[r][p] = c*arp - s*arq;  A[r][q] = s*arp + c*arq;
      }
      for ( int r = 0 ; r < 3 ; r++ ) {
        double apr = A[p][r], aqr = A[q][r];
        A[p][r] = c*apr - s*aqr;  A[q][r] = s*apr + c*aqr;
      }
      for ( int r = 0 ; r < 3 ; r++ ) {
        double qrp = Q[r][p], qrq = Q[r][q];
        Q[r][p] = c*qrp - s*qrq;  Q[r][q] = s*qrp + c*qrq;
      }
    }
  }
  for ( int i = 0 ; i < 3 ; i++ ) w[i] = A[i][i];
}

/* root-PCA frame: one global principal-axes rotation stored in the cache and
 * applied ONCE per ray -- tightens EVERY node of the tree for elongated
 * DIAGONAL geometry at zero per-node cost. Gated on anisotropy >= 2 so
 * blob-like meshes keep today's exact behavior (identity frame).              */
static void computePCA( const double *V , mwSize nV , TreeCache &tc )
{
  tc.pca = false;
  for ( int i = 0 ; i < 9 ; i++ ) tc.Rr[i] = ( i % 4 == 0 ) ? 1.0 : 0.0;
  tc.Cc[0] = tc.Cc[1] = tc.Cc[2] = 0.0;
  if ( nV < 4 ) return;
  double c[3] = {0,0,0};
  for ( mwSize i = 0 ; i < nV ; i++ )
    for ( int a = 0 ; a < 3 ; a++ ) c[a] += V[3*i+a];
  for ( int a = 0 ; a < 3 ; a++ ) c[a] /= (double)nV;
  double C[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
  for ( mwSize i = 0 ; i < nV ; i++ ) {
    double d0 = V[3*i]-c[0], d1 = V[3*i+1]-c[1], d2 = V[3*i+2]-c[2];
    C[0][0]+=d0*d0; C[0][1]+=d0*d1; C[0][2]+=d0*d2;
    C[1][1]+=d1*d1; C[1][2]+=d1*d2; C[2][2]+=d2*d2;
  }
  C[1][0]=C[0][1]; C[2][0]=C[0][2]; C[2][1]=C[1][2];
  double w[3], Q[3][3];
  jacobi3( C , w , Q );
  double wmax = fmax(w[0],fmax(w[1],w[2])), wmin = fmin(w[0],fmin(w[1],w[2]));
  if ( !(wmin > 0.0) || sqrt( wmax/wmin ) < 2.0 ) return;   /* blob: identity  */
  tc.pca = true;
  for ( int a = 0 ; a < 3 ; a++ ) {
    tc.Cc[a] = c[a];
    for ( int b = 0 ; b < 3 ; b++ ) tc.Rr[3*a+b] = Q[b][a]; /* rows = axes     */
  }
}

static void buildTree( const double *V , const int *T , mwSize nV , mwSize nF , TreeCache &tc )
{
  computePCA( V , nV , tc );
  std::vector<double> Vr;
  if ( tc.pca ) {                                    /* rotate vertices ONCE    */
    Vr.resize( 3*nV );
    for ( mwSize i = 0 ; i < nV ; i++ ) {
      double d0 = V[3*i]-tc.Cc[0], d1 = V[3*i+1]-tc.Cc[1], d2 = V[3*i+2]-tc.Cc[2];
      for ( int a = 0 ; a < 3 ; a++ )
        Vr[3*i+a] = tc.Rr[3*a]*d0 + tc.Rr[3*a+1]*d1 + tc.Rr[3*a+2]*d2;
    }
    V = &Vr[0];
  }

  BuildCtx B;
  B.V = V;  B.T = T;
  B.bb.resize( 6*nF );  B.cc.resize( 3*nF );  B.ord.resize( nF );
  for ( mwSize f = 0 ; f < nF ; f++ ) {
    double *a = &B.bb[6*f];
    for ( int c = 0 ; c < 3 ; c++ ) {
      double x0 = V[ 3*T[3*f+0] + c ];
      double x1 = V[ 3*T[3*f+1] + c ];
      double x2 = V[ 3*T[3*f+2] + c ];
      a[c]   = fmin( x0 , fmin( x1 , x2 ) );
      a[3+c] = fmax( x0 , fmax( x1 , x2 ) );
      B.cc[3*f+c] = ( a[c] + a[3+c] )*0.5;
    }
    B.ord[f] = (int)f;
  }
  tc.nodes.clear();  tc.tris4.clear();
  tc.nodes.reserve( nF ? nF/4 : 1 );
  tc.tris4.reserve( nF/4 + 8 );
  B.nodes = &tc.nodes;  B.tris4 = &tc.tris4;
  tc.root = build4( B , 0 , (int)nF );
}

/* ------------------------------ tracing ------------------------------------- */
/* leaf kernel over blocks [first, first + ceil(count/4)): updates best/all;
 * returns TRUE only in 'any' mode when an accepting hit was found (early exit) */
static inline bool mtBlocks( const std::vector<PreTri4> &tris4 , int firstBlock , int count ,
                             const Ray &ry , bool avx , int mode ,
                             double &best , int &hb ,
                             std::vector<double> *tAll , std::vector<int> *hAll )
{
  RayA ra;
  if ( avx ) makeRayA( ry , ra );
  int nBlk = ( count + 3 ) >> 2;
  for ( int b = 0 ; b < nBlk ; b++ ) {
    const PreTri4 &T = tris4[ firstBlock + b ];
    double t4[4];
    int m;
    if ( avx ) {
      m = mt4( ra , T , t4 );
    } else {
      m = 0;
      for ( int l = 0 ; l < 4 ; l++ ) {
        double t;
        if ( mtLane( ry , T , l , &t ) ) { t4[l] = t; m |= ( 1 << l ); }
      }
    }
    while ( m ) {
      int l = m & (-m);                                  /* lowest set bit       */
      int lane = ( l == 1 ) ? 0 : ( l == 2 ) ? 1 : ( l == 4 ) ? 2 : 3;
      m &= m - 1;
      if ( T.id[lane] < 0 ) continue;                    /* padding              */
      double t = t4[lane];
      if      ( mode == 0 ) { if ( t < best ) { best = t; hb = T.id[lane] + 1; } }
      else if ( mode == 1 ) { if ( t > best ) { best = t; hb = T.id[lane] + 1; } }
      else if ( mode == 2 ) { tAll->push_back( t ); hAll->push_back( T.id[lane] + 1 ); }
      else {                                             /* 'any': early exit    */
        if ( t > ANY_TLO && t < ANY_THI ) { best = t; hb = T.id[lane] + 1; return true; }
      }
    }
  }
  return false;
}

static void trace4( const TreeCache &tc , const Ray &ry , int mode ,
                    double *tBest , int *hBest ,
                    std::vector<double> *tAll , std::vector<int> *hAll )
{
  double best = ( mode == 1 ) ? -INFINITY : INFINITY;
  int    hb   = 0;
  bool   avx  = useAVX();

  if ( tc.root.count > 0 ) {                             /* whole tree is a leaf */
    mtBlocks( tc.tris4 , tc.root.idx , tc.root.count , ry , avx , mode , best , hb , tAll , hAll );
    if ( mode != 2 ) { *tBest = hb ? best : mxGetNaN(); *hBest = hb; }
    return;
  }

  int stack[128];  int sp = 0;
  stack[sp++] = tc.root.idx;

  while ( sp ) {
    const Node4 &nd = tc.nodes[ stack[--sp] ];

    float tn4[4], tf4[4];
    slab4( ry , nd , tn4 , tf4 );

    int   cand[4];  float ckey[4];  int ncand = 0;
    for ( int i = 0 ; i < 4 ; i++ ) {
      if ( nd.count[i] < 0 ) continue;                   /* empty slot           */
      if ( tn4[i] > tf4[i] ) continue;                   /* slab miss            */
      if      ( mode == 0 ) { if ( tn4[i] > (float)best ) continue; }
      else if ( mode == 1 ) { if ( tf4[i] < (float)best ) continue; }
      else if ( mode == 3 ) { if ( tf4[i] < 0.0f || tn4[i] > 1.0f ) continue; }

      if ( nd.count[i] > 0 ) {                           /* leaf: resolve now    */
        if ( mtBlocks( tc.tris4 , nd.child[i] , nd.count[i] , ry , avx , mode , best , hb , tAll , hAll ) ) {
          *tBest = best;  *hBest = hb;  return;          /* 'any' found          */
        }
      } else {                                           /* internal: candidate  */
        cand[ncand] = nd.child[i];
        ckey[ncand] = ( mode == 1 ) ? -tf4[i] : tn4[i];  /* order key            */
        ncand++;
      }
    }

    if ( sp + ncand > 124 )
      mexErrMsgIdAndTxt( "IntersectSurfaceRay_mx:stack" , "traversal stack overflow." );

    if ( mode == 0 || mode == 1 ) {
      /* push FAR first so the NEAR one pops first (tightens `best` earliest)    */
      for ( int i = 1 ; i < ncand ; i++ ) {              /* insertion sort asc   */
        int  ci = cand[i];  float ki = ckey[i];  int j = i-1;
        while ( j >= 0 && ckey[j] > ki ) { cand[j+1] = cand[j]; ckey[j+1] = ckey[j]; j--; }
        cand[j+1] = ci;  ckey[j+1] = ki;
      }
      for ( int i = ncand-1 ; i >= 0 ; i-- ) stack[sp++] = cand[i];
    } else {
      for ( int i = 0 ; i < ncand ; i++ ) stack[sp++] = cand[i];
    }
  }

  if ( mode != 2 ) { *tBest = hb ? best : mxGetNaN(); *hBest = hb; }
}

/* ------------------ brute force over PreTri4 blocks ------------------------- */
static void buildBlocks( const double *V , const int *T , mwSize nF ,
                         std::vector<PreTri4> &blocks )
{
  mwSize nB = ( nF + 3 ) >> 2;
  blocks.resize( nB );
  for ( mwSize b = 0 ; b < nB ; b++ ) {
    PreTri4 &blk = blocks[b];
    for ( int l = 0 ; l < 4 ; l++ ) {
      mwSize f = 4*b + l;
      if ( f < nF ) {
        double x0[3];
        for ( int c = 0 ; c < 3 ; c++ ) x0[c] = V[ 3*T[3*f+0] + c ];
        blk.v0x[l] = x0[0];  blk.v0y[l] = x0[1];  blk.v0z[l] = x0[2];
        blk.e1x[l] = V[ 3*T[3*f+1] + 0 ] - x0[0];
        blk.e1y[l] = V[ 3*T[3*f+1] + 1 ] - x0[1];
        blk.e1z[l] = V[ 3*T[3*f+1] + 2 ] - x0[2];
        blk.e2x[l] = V[ 3*T[3*f+2] + 0 ] - x0[0];
        blk.e2y[l] = V[ 3*T[3*f+2] + 1 ] - x0[1];
        blk.e2z[l] = V[ 3*T[3*f+2] + 2 ] - x0[2];
        blk.id[l]  = (int)f;
      } else {
        blk.v0x[l]=blk.v0y[l]=blk.v0z[l] = 0.0;
        blk.e1x[l]=blk.e1y[l]=blk.e1z[l] = 0.0;
        blk.e2x[l]=blk.e2y[l]=blk.e2z[l] = 0.0;
        blk.id[l]  = -1;
      }
    }
  }
}

static void traceBrute( const std::vector<PreTri4> &blocks , mwSize nF ,
                        const Ray &ry , int mode ,
                        double *tBest , int *hBest ,
                        std::vector<double> *tAll , std::vector<int> *hAll )
{
  double best = ( mode == 1 ) ? -INFINITY : INFINITY;
  int    hb   = 0;
  bool   avx  = useAVX();
  mtBlocks( blocks , 0 , (int)nF , ry , avx , mode , best , hb , tAll , hAll );
  if ( mode != 2 ) { *tBest = hb ? best : mxGetNaN(); *hBest = hb; }
}

/* --------------------- packet traversal (8 rays, modes first/any) ------------
 * Coherent rays SHARE the tree walk: each Node4 child box is slab-tested
 * against all 8 rays in one AVX pass (sign-agnostic min/max form, so mixed
 * direction signs inside the packet are fine) and carried down with per-lane
 * masks -- the node fetches (the actual bottleneck) are amortized across the
 * packet. Per-ray best-t pruning ('first') / window clip + retirement ('any')
 * work exactly as in the single-ray path; leaves run the same 4-wide MT per
 * active lane, so results are IDENTICAL to single-ray tracing. Rays with a
 * zero direction component are routed to the single-ray path by the caller
 * (their slab needs the careful containment variant).                          */
struct Pkt8 {
  float  ox[8], oy[8], oz[8], ivx[8], ivy[8], ivz[8];
  double best[8];
  float  bestF[8];
  int    hb[8];
  int    ridx[8];                      /* index into the caller's Ray array     */
  int    n;
};

static void tracePkt8( const TreeCache &tc , const Ray *rays , Pkt8 &P , int mode )
{
  bool avx = useAVX();
  unsigned act = ( P.n >= 8 ) ? 255u : ( ( 1u << P.n ) - 1u );

  if ( tc.root.count > 0 ) {                             /* degenerate root leaf */
    for ( int l = 0 ; l < P.n ; l++ ) {
      const Ray &ry = rays[ P.ridx[l] ];
      double best = INFINITY; int hb = 0;
      mtBlocks( tc.tris4 , tc.root.idx , tc.root.count , ry , avx , mode , best , hb , NULL , NULL );
      if ( hb ) { P.best[l] = best; P.bestF[l] = (float)best; P.hb[l] = hb; }
    }
    return;
  }

  __m256 ox = _mm256_loadu_ps( P.ox ), oy = _mm256_loadu_ps( P.oy ), oz = _mm256_loadu_ps( P.oz );
  __m256 ix = _mm256_loadu_ps( P.ivx ), iy = _mm256_loadu_ps( P.ivy ), iz = _mm256_loadu_ps( P.ivz );
  const __m256 zero = _mm256_setzero_ps(), one = _mm256_set1_ps( 1.0f );

  struct SE { int node; unsigned mask; };
  SE stack[192];  int sp = 0;
  stack[sp].node = tc.root.idx;  stack[sp].mask = act;  sp++;

  while ( sp ) {
    SE e = stack[--sp];
    unsigned m0 = e.mask & act;
    if ( !m0 ) continue;
    const Node4 &nd = tc.nodes[ e.node ];

    int cand[4]; unsigned cmask[4]; float ckey[4]; int nc = 0;

    for ( int i = 0 ; i < 4 ; i++ ) {
      if ( nd.count[i] < 0 ) continue;
      /* sign-agnostic 8-wide slab vs the (scalar per child) box                 */
      __m256 t0, t1, a, b, tn, tf;
      t0 = _mm256_mul_ps( _mm256_sub_ps( _mm256_set1_ps( nd.bminx[i] ) , ox ) , ix );
      t1 = _mm256_mul_ps( _mm256_sub_ps( _mm256_set1_ps( nd.bmaxx[i] ) , ox ) , ix );
      tn = _mm256_min_ps( t0 , t1 );  tf = _mm256_max_ps( t0 , t1 );
      t0 = _mm256_mul_ps( _mm256_sub_ps( _mm256_set1_ps( nd.bminy[i] ) , oy ) , iy );
      t1 = _mm256_mul_ps( _mm256_sub_ps( _mm256_set1_ps( nd.bmaxy[i] ) , oy ) , iy );
      a  = _mm256_min_ps( t0 , t1 );  b  = _mm256_max_ps( t0 , t1 );
      tn = _mm256_max_ps( tn , a );   tf = _mm256_min_ps( tf , b );
      t0 = _mm256_mul_ps( _mm256_sub_ps( _mm256_set1_ps( nd.bminz[i] ) , oz ) , iz );
      t1 = _mm256_mul_ps( _mm256_sub_ps( _mm256_set1_ps( nd.bmaxz[i] ) , oz ) , iz );
      a  = _mm256_min_ps( t0 , t1 );  b  = _mm256_max_ps( t0 , t1 );
      tn = _mm256_max_ps( tn , a );   tf = _mm256_min_ps( tf , b );

      __m256 hit = _mm256_cmp_ps( tn , tf , _CMP_LE_OQ );
      if ( mode == 0 )
        hit = _mm256_and_ps( hit , _mm256_cmp_ps( tn , _mm256_loadu_ps( P.bestF ) , _CMP_LE_OQ ) );
      else
        hit = _mm256_and_ps( hit , _mm256_and_ps(
                _mm256_cmp_ps( tf , zero , _CMP_GE_OQ ) ,
                _mm256_cmp_ps( tn , one  , _CMP_LE_OQ ) ) );
      unsigned m = (unsigned)_mm256_movemask_ps( hit ) & m0;
      if ( !m ) continue;

      if ( nd.count[i] > 0 ) {                           /* leaf: per active lane */
        for ( int l = 0 ; l < 8 ; l++ ) {
          if ( !( m & ( 1u << l ) ) ) continue;
          const Ray &ry = rays[ P.ridx[l] ];
          if ( mode == 0 ) {
            double best = P.best[l];  int hb = P.hb[l];
            mtBlocks( tc.tris4 , nd.child[i] , nd.count[i] , ry , avx , 0 , best , hb , NULL , NULL );
            P.best[l] = best;  P.bestF[l] = (float)best;  P.hb[l] = hb;
          } else {
            double best = INFINITY;  int hb = 0;
            if ( mtBlocks( tc.tris4 , nd.child[i] , nd.count[i] , ry , avx , 3 , best , hb , NULL , NULL ) ) {
              P.best[l] = best;  P.hb[l] = hb;
              act &= ~( 1u << l );                       /* resolved: retire ray  */
              if ( !act ) return;
            }
          }
        }
      } else {
        float key = FLT_MAX;                             /* masked min entry      */
        if ( mode == 0 ) {
          float tnv[8];  _mm256_storeu_ps( tnv , tn );
          for ( int l = 0 ; l < 8 ; l++ ) if ( m & (1u<<l) ) key = ( tnv[l] < key ) ? tnv[l] : key;
        }
        cand[nc] = nd.child[i];  cmask[nc] = m;  ckey[nc] = key;  nc++;
      }
    }

    if ( sp + nc > 188 )
      mexErrMsgIdAndTxt( "IntersectSurfaceRay_mx:stack" , "packet traversal stack overflow." );

    if ( mode == 0 ) {                                   /* near pops first       */
      for ( int i = 1 ; i < nc ; i++ ) {
        int ci = cand[i]; unsigned mi = cmask[i]; float ki = ckey[i]; int j = i-1;
        while ( j >= 0 && ckey[j] > ki ) { cand[j+1]=cand[j]; cmask[j+1]=cmask[j]; ckey[j+1]=ckey[j]; j--; }
        cand[j+1]=ci; cmask[j+1]=mi; ckey[j+1]=ki;
      }
      for ( int i = nc-1 ; i >= 0 ; i-- ) { stack[sp].node = cand[i]; stack[sp].mask = cmask[i]; sp++; }
    } else {
      for ( int i = 0 ; i < nc ; i++ ) { stack[sp].node = cand[i]; stack[sp].mask = cmask[i]; sp++; }
    }
  }
}

/* --------------- surface scanning / gathering ------------------------------- */
struct SurfRef { const mxArray *mv, *mt; mwSize nV, nF; };

static void scanSurface( const mxArray *S , mwSize elem , int sid , std::vector<SurfRef> &list )
{
  const mxArray *mv = mxGetField( S , elem , "vertices" );
  const mxArray *mt = mxGetField( S , elem , "faces" );
  if ( !mv || !mt ) {
    mv = mxGetField( S , elem , "xyz" );
    mt = mxGetField( S , elem , "tri" );
  }
  if ( !mv || !mt )
    mexErrMsgIdAndTxt( "IntersectSurfaceRay_mx:P" ,
      "P(%d) must have fields .vertices/.faces or .xyz/.tri." , sid );
  if ( !mxIsDouble(mv) || ( mxGetN(mv) != 3 && !mxIsEmpty(mv) ) )
    mexErrMsgIdAndTxt( "IntersectSurfaceRay_mx:P" , "P(%d) vertices must be nV x 3 double." , sid );
  if ( !mxIsDouble(mt) || ( mxGetN(mt) != 3 && !mxIsEmpty(mt) ) )
    mexErrMsgIdAndTxt( "IntersectSurfaceRay_mx:P" , "P(%d) faces must be nF x 3 double (triangles only)." , sid );
  SurfRef r;  r.mv = mv;  r.mt = mt;  r.nV = mxGetM( mv );  r.nF = mxGetM( mt );
  list.push_back( r );
}

static void gatherSurfaces( const std::vector<SurfRef> &list ,
                            std::vector<double> &V , std::vector<int> &T ,
                            std::vector<int> &surfId , std::vector<int> &locId )
{
  for ( size_t s = 0 ; s < list.size() ; s++ ) {
    const SurfRef &r = list[s];
    const double *v = mxGetPr( r.mv );
    const double *t = mxGetPr( r.mt );
    int base = (int)( V.size()/3 );
    for ( mwSize i = 0 ; i < r.nV ; i++ ) {
      V.push_back( v[ i          ] );
      V.push_back( v[ i +   r.nV ] );
      V.push_back( v[ i + 2*r.nV ] );
    }
    for ( mwSize f = 0 ; f < r.nF ; f++ ) {
      for ( int c = 0 ; c < 3 ; c++ ) {
        double x = t[ f + c*r.nF ];
        if ( !(x >= 1) || x > (double)r.nV || x != floor(x) )
          mexErrMsgIdAndTxt( "IntersectSurfaceRay_mx:P" ,
            "P(%d) faces must hold integer indices in 1..nV." , (int)s+1 );
        T.push_back( (int)x - 1 + base );
      }
      surfId.push_back( (int)s + 1 );
      locId.push_back( (int)f + 1 );
    }
  }
}

/* ------------------------------- gateway ------------------------------------ */
void mexFunction( int nlhs , mxArray *plhs[] , int nrhs , const mxArray *prhs[] )
{
  if ( nrhs < 2 || nrhs > 3 )
    mexErrMsgIdAndTxt( "IntersectSurfaceRay_mx:nrhs" ,
      "usage: [xyz,P_id,cell_id,t,ray_id] = IntersectSurfaceRay_mx(P,ray[,MODE])" );

  /* ---- P: scan surfaces WITHOUT copying ------------------------------------ */
  std::vector<SurfRef> surfs;
  const mxArray *P = prhs[0];
  if ( mxIsStruct( P ) ) {
    for ( mwSize e = 0 ; e < mxGetNumberOfElements( P ) ; e++ )
      scanSurface( P , e , (int)e+1 , surfs );
  } else if ( mxIsCell( P ) ) {
    for ( mwSize e = 0 ; e < mxGetNumberOfElements( P ) ; e++ ) {
      const mxArray *S = mxGetCell( P , e );
      if ( !S || !mxIsStruct( S ) )
        mexErrMsgIdAndTxt( "IntersectSurfaceRay_mx:P" ,
          "P{%d} must be a mesh struct (graphics handles only via IntersectSurfaceRay.m)." , (int)e+1 );
      scanSurface( S , 0 , (int)e+1 , surfs );
    }
  } else {
    mexErrMsgIdAndTxt( "IntersectSurfaceRay_mx:P" ,
      "P must be a mesh struct or a cell array of them (graphics handles only via IntersectSurfaceRay.m)." );
  }
  mwSize nV = 0, nF = 0;
  for ( size_t s = 0 ; s < surfs.size() ; s++ ) { nV += surfs[s].nV; nF += surfs[s].nF; }

  /* ---- rays: 2 x 3 (one), 2 x 3 x N (pages), or N x 6 (rows) ---------------- */
  const mxArray *mR = prhs[1];
  if ( !mxIsDouble( mR ) )
    mexErrMsgIdAndTxt( "IntersectSurfaceRay_mx:ray" , "ray must be 2x3, 2x3xN or N x 6 double." );
  mwSize k;
  std::vector<double> RR;                                /* row-major k x 6      */
  mwSize ndim = mxGetNumberOfDimensions( mR );
  const mwSize *dims = mxGetDimensions( mR );
  if ( ndim == 3 && dims[0] == 2 && dims[1] == 3 ) {     /* 2 x 3 x N pages      */
    k = dims[2]; RR.resize( 6*k );
    const double *r = mxGetPr( mR );
    for ( mwSize i = 0 ; i < k ; i++ )
      for ( int c = 0 ; c < 3 ; c++ ) {
        RR[ 6*i + c     ] = r[ 0 + 2*c + 6*i ];          /* p0 = page row 1      */
        RR[ 6*i + c + 3 ] = r[ 1 + 2*c + 6*i ];          /* p1 = page row 2      */
      }
  } else if ( ndim == 2 && mxGetM(mR) == 2 && mxGetN(mR) == 3 ) {
    k = 1; RR.resize( 6 );
    const double *r = mxGetPr( mR );
    RR[0]=r[0]; RR[1]=r[2]; RR[2]=r[4];                  /* p0 = row 1           */
    RR[3]=r[1]; RR[4]=r[3]; RR[5]=r[5];                  /* p1 = row 2           */
  } else if ( ndim == 2 && ( mxGetN(mR) == 6 || mxIsEmpty(mR) ) ) {
    k = mxGetM( mR ); RR.resize( 6*k );
    const double *r = mxGetPr( mR );
    for ( mwSize i = 0 ; i < k ; i++ )
      for ( int c = 0 ; c < 6 ; c++ ) RR[ 6*i + c ] = r[ i + c*k ];
  } else {
    mexErrMsgIdAndTxt( "IntersectSurfaceRay_mx:ray" , "ray must be 2x3, 2x3xN or N x 6 double." );
  }

  /* ---- MODE ------------------------------------------------------------------ */
  int mode = 0;
  if ( nrhs >= 3 && !mxIsEmpty( prhs[2] ) ) {
    if ( mxIsChar( prhs[2] ) ) {
      char ms[16] = {0};
      mxGetString( prhs[2] , ms , 15 );
      for ( char *c = ms ; *c ; c++ ) *c = (char)tolower( *c );
      if      ( !strcmp( ms , "first" ) ) mode = 0;
      else if ( !strcmp( ms , "last"  ) ) mode = 1;
      else if ( !strcmp( ms , "all"   ) ) mode = 2;
      else if ( !strcmp( ms , "any"   ) ) mode = 3;
      else mexErrMsgIdAndTxt( "IntersectSurfaceRay_mx:MODE" , "unknown mode '%s'." , ms );
    } else {
      mode = (int)mxGetScalar( prhs[2] );
      if ( mode < 0 || mode > 3 )
        mexErrMsgIdAndTxt( "IntersectSurfaceRay_mx:MODE" , "MODE must be 'first', 'last', 'all' or 'any' (0/1/2/3)." );
    }
  }

  /* ---- acceleration policy ----------------------------------------------------
   * nF >= 1024 and k >= 16     -> cached tree (hash raw bytes; hit = NO copies
   *                               at all; miss = gather + build once)
   * nF >= 1024 and k <  16     -> cache LOOKUP only (a hit still serves single
   *                               picks); miss = gather + brute (no build)
   * nF <  1024                 -> gather + brute (never hash)                  */
  const TreeCache *tc = NULL;
  std::vector<double> V;  std::vector<int> T;
  std::vector<int> surfIdRaw, locIdRaw;
  const std::vector<int> *surfId = &surfIdRaw, *locId = &locIdRaw;
  std::vector<PreTri4> bruteBlocks;

  if ( nF >= 1024 ) {
    u64 key = 1469598103934665603ULL;
    for ( size_t s = 0 ; s < surfs.size() ; s++ ) {
      key = hashBytes( mxGetPr(surfs[s].mv) , (size_t)(3*surfs[s].nV)*sizeof(double) , key );
      key = hashBytes( mxGetPr(surfs[s].mt) , (size_t)(3*surfs[s].nF)*sizeof(double) , key );
      key ^= ( (u64)surfs[s].nV << 32 ) ^ (u64)surfs[s].nF;  key *= 0x100000001B3ULL;
    }
    int slot = -1;
    for ( int i = 0 ; i < CACHE_SLOTS ; i++ )
      if ( !g_cache[i].tris4.empty() && g_cache[i].key == key &&
           g_cache[i].nV == nV && g_cache[i].nF == nF ) { slot = i; break; }
    if ( slot < 0 && k >= 16 ) {                        /* build once, cache    */
      slot = 0;
      for ( int i = 1 ; i < CACHE_SLOTS ; i++ )
        if ( g_cache[i].tris4.empty() || g_cache[i].age < g_cache[slot].age ) slot = i;
      TreeCache &c = g_cache[slot];
      c.key = key;  c.nV = nV;  c.nF = nF;
      c.surfId.clear();  c.locId.clear();
      V.reserve( 3*nV );  T.reserve( 3*nF );
      gatherSurfaces( surfs , V , T , c.surfId , c.locId );
      buildTree( &V[0] , &T[0] , nV , nF , c );
      V.clear(); V.shrink_to_fit();  T.clear(); T.shrink_to_fit();
    }
    if ( slot >= 0 ) {
      g_cache[slot].age = ++g_clock;
      tc = &g_cache[slot];
      surfId = &tc->surfId;  locId = &tc->locId;
    }
  }
  if ( !tc && nF > 0 ) {                                /* brute path            */
    V.reserve( 3*nV );  T.reserve( 3*nF );
    gatherSurfaces( surfs , V , T , surfIdRaw , locIdRaw );
    buildBlocks( &V[0] , &T[0] , nF , bruteBlocks );    /* ONCE per call         */
  }

  /* ---- trace ------------------------------------------------------------------ */
  if ( mode != 2 ) {                                    /* first / last / any    */
    plhs[0] = mxCreateDoubleMatrix( k , 3 , mxREAL );
    mxArray *mP = mxCreateDoubleMatrix( k , 1 , mxREAL );
    mxArray *mC = mxCreateDoubleMatrix( k , 1 , mxREAL );
    mxArray *mT = mxCreateDoubleMatrix( k , 1 , mxREAL );
    mxArray *mI = mxCreateDoubleMatrix( k , 1 , mxREAL );
    double *xyz = mxGetPr( plhs[0] );
    double *pO = mxGetPr( mP ), *cO = mxGetPr( mC ), *tO = mxGetPr( mT ), *rO = mxGetPr( mI );

    std::vector<double> tbv( k , 0.0 );
    std::vector<int>    hbv( k , 0 );

    if ( USE_PACKETS && tc && ( mode == 0 || mode == 3 ) && k >= 16 && useAVX() ) {
      /* PACKET path: batches of 8 rays share the tree walk                     */
      std::vector<Ray> rays( k );
      for ( mwSize r = 0 ; r < k ; r++ ) makeRayTC( tc , &RR[6*r] , &RR[6*r+3] , rays[r] );
      Pkt8 P;  P.n = 0;
      for ( mwSize r = 0 ; r <= k ; r++ ) {
        bool flush = ( r == k );
        if ( !flush ) {
          if ( rays[r].anyZero ) {                       /* zero-dir: single-ray */
            double tb; int hb = 0;
            trace4( *tc , rays[r] , mode , &tb , &hb , NULL , NULL );
            tbv[r] = tb;  hbv[r] = hb;
            continue;
          }
          int l = P.n;
          P.ox[l]  = rays[r].of[0];    P.oy[l]  = rays[r].of[1];    P.oz[l]  = rays[r].of[2];
          P.ivx[l] = rays[r].invf[0];  P.ivy[l] = rays[r].invf[1];  P.ivz[l] = rays[r].invf[2];
          P.best[l] = INFINITY;  P.bestF[l] = FLT_MAX;  P.hb[l] = 0;  P.ridx[l] = (int)r;
          P.n++;
        }
        if ( P.n == 8 || ( flush && P.n > 0 ) ) {
          for ( int l = P.n ; l < 8 ; l++ ) {            /* neutral unused lanes */
            P.ox[l]=P.oy[l]=P.oz[l] = 0.0f;  P.ivx[l]=P.ivy[l]=P.ivz[l] = 0.0f;
            P.best[l] = INFINITY;  P.bestF[l] = FLT_MAX;  P.hb[l] = 0;  P.ridx[l] = P.ridx[0];
          }
          tracePkt8( *tc , &rays[0] , P , mode );
          for ( int l = 0 ; l < P.n ; l++ )
            if ( P.hb[l] ) { tbv[ P.ridx[l] ] = P.best[l];  hbv[ P.ridx[l] ] = P.hb[l]; }
          P.n = 0;
        }
      }
    } else {
      for ( mwSize r = 0 ; r < k ; r++ ) {
        Ray ry;
        double tb; int hb = 0;
        if ( tc ) {
          makeRayTC( tc , &RR[6*r] , &RR[6*r+3] , ry );
          trace4( *tc , ry , mode , &tb , &hb , NULL , NULL );
        } else if ( nF > 0 ) {
          makeRay( &RR[6*r] , &RR[6*r+3] , ry );
          traceBrute( bruteBlocks , nF , ry , mode , &tb , &hb , NULL , NULL );
        }
        tbv[r] = tb;  hbv[r] = hb;
      }
    }

    for ( mwSize r = 0 ; r < k ; r++ ) {
      rO[r] = (double)(r+1);
      int hb = hbv[r];
      if ( hb ) {
        double tb = tbv[r];
        pO[r] = (*surfId)[hb-1];  cO[r] = (*locId)[hb-1];  tO[r] = tb;
        for ( int c = 0 ; c < 3 ; c++ )                  /* from the ORIGINAL ray */
          xyz[ r + c*k ] = RR[6*r+c] + tb*( RR[6*r+3+c] - RR[6*r+c] );
      } else {
        pO[r] = 0;  cO[r] = 0;  tO[r] = mxGetNaN();
        for ( int c = 0 ; c < 3 ; c++ ) xyz[ r + c*k ] = mxGetNaN();
      }
    }
    if ( nlhs >= 2 ) plhs[1] = mP; else mxDestroyArray( mP );
    if ( nlhs >= 3 ) plhs[2] = mC; else mxDestroyArray( mC );
    if ( nlhs >= 4 ) plhs[3] = mT; else mxDestroyArray( mT );
    if ( nlhs >= 5 ) plhs[4] = mI; else mxDestroyArray( mI );
  } else {                                              /* all                   */
    std::vector<double> tsAll;  std::vector<int> hsAll;  std::vector<int> rsAll;
    std::vector<double> ts;     std::vector<int> hs;
    for ( mwSize r = 0 ; r < k ; r++ ) {
      Ray ry;
      ts.clear(); hs.clear();
      if ( tc ) {
        makeRayTC( tc , &RR[6*r] , &RR[6*r+3] , ry );
        trace4( *tc , ry , 2 , NULL , NULL , &ts , &hs );
      } else if ( nF > 0 ) {
        makeRay( &RR[6*r] , &RR[6*r+3] , ry );
        traceBrute( bruteBlocks , nF , ry , 2 , NULL , NULL , &ts , &hs );
      }
      std::vector<int> ord( ts.size() );
      for ( size_t i = 0 ; i < ord.size() ; i++ ) ord[i] = (int)i;
      std::sort( ord.begin() , ord.end() , [&]( int a , int b ){ return ts[a] < ts[b]; } );
      for ( size_t i = 0 ; i < ord.size() ; i++ ) {
        tsAll.push_back( ts[ord[i]] );  hsAll.push_back( hs[ord[i]] );  rsAll.push_back( (int)r+1 );
      }
    }
    mwSize m = (mwSize)tsAll.size();
    plhs[0] = mxCreateDoubleMatrix( m , 3 , mxREAL );
    mxArray *mP = mxCreateDoubleMatrix( m , 1 , mxREAL );
    mxArray *mC = mxCreateDoubleMatrix( m , 1 , mxREAL );
    mxArray *mT = mxCreateDoubleMatrix( m , 1 , mxREAL );
    mxArray *mI = mxCreateDoubleMatrix( m , 1 , mxREAL );
    double *xyz = mxGetPr( plhs[0] );
    double *pO = mxGetPr( mP ), *cO = mxGetPr( mC ), *tO = mxGetPr( mT ), *rO = mxGetPr( mI );
    for ( mwSize i = 0 ; i < m ; i++ ) {
      int h = hsAll[i], r = rsAll[i]-1;
      const double *p0 = &RR[6*r], *p1 = &RR[6*r+3];
      pO[i] = (*surfId)[h-1];  cO[i] = (*locId)[h-1];  tO[i] = tsAll[i];  rO[i] = rsAll[i];
      for ( int c = 0 ; c < 3 ; c++ ) xyz[ i + c*m ] = p0[c] + tsAll[i]*( p1[c]-p0[c] );
    }
    if ( nlhs >= 2 ) plhs[1] = mP; else mxDestroyArray( mP );
    if ( nlhs >= 3 ) plhs[2] = mC; else mxDestroyArray( mC );
    if ( nlhs >= 4 ) plhs[3] = mT; else mxDestroyArray( mT );
    if ( nlhs >= 5 ) plhs[4] = mI; else mxDestroyArray( mI );
  }
}
