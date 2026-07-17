/* bvhIntersectRay_mx  --  ray/line vs mesh triangles on the BVH blob.
 *
 *   [xyz, cell, t, rid] = bvhIntersectRay_mx( rays , B , mode , nthreads )
 *
 *     rays  N x 6 double  [p0 , p1] per row; the LINE hit = p0 + t*(p1-p0),
 *                         t UNBOUNDED (negatives = behind p0), same convention
 *                         as tools/IntersectSurfaceRay.
 *     B     struct        BVH v2 blob WITH packed leaves (pkV/pkS/pkT/pkE)
 *     mode  1 first | 2 last | 3 all | 4 any
 *     nthreads            OpenMP threads over rays (not for 'all')
 *
 *   'first'/'last' run an ORDERED traversal pruned by the best t so far:
 *   slots are visited by entry (exit) parameter and culled when their whole
 *   interval lies beyond the current best -- correct with NEGATIVE t too,
 *   since the ordering key is the interval bound itself. 'any' is the
 *   occlusion query (1e-9 < t < 1-1e-5, early exit); 'all' visits everything.
 *   Leaves use a 4-wide branchless AVX Moller-Trumbore (runtime CPUID, scalar
 *   fallback), double precision, INCLUSIVE barycentric tolerance 1e-9 as in
 *   IntersectSurfaceRay_mx. Only TRIANGLE cells are tested.
 *
 *   Compile:  mex COMPFLAGS="$COMPFLAGS /openmp" bvhIntersectRay_mx.cpp
 *
 * See also bvhIntersectRay, BVH, bvhClosestElement_mx.
 */

#include "mex.h"
#include <cmath>
#include <cstdint>
#include <cstring>
#include <vector>
#include <algorithm>
#include <limits>
#include <immintrin.h>
#if defined(_MSC_VER)
#include <intrin.h>
#endif
#ifdef _OPENMP
#include <omp.h>
#endif

static const double INF = std::numeric_limits<double>::infinity();

#define UV_EPS  1e-9
#define ANY_TLO 1e-9
#define ANY_THI (1.0 - 1e-5)

struct Hit { double t; int32_t cell; };

static int g_avx = -1;
static bool useAVX( void )
{
  if( g_avx < 0 ) {
#if defined(_MSC_VER)
    int ci[4];  __cpuid( ci, 1 );
    const bool osxsave = ( ci[2] >> 27 ) & 1;
    const bool avx     = ( ci[2] >> 28 ) & 1;
    g_avx = ( osxsave && avx && ( ( _xgetbv(0) & 6 ) == 6 ) ) ? 1 : 0;
#else
    g_avx = __builtin_cpu_supports( "avx" ) ? 1 : 0;
#endif
  }
  return g_avx == 1;
}

/* scalar Moller-Trumbore, double, inclusive tolerance; q = 9 doubles (A,B,C) */
static inline bool mtHit( const double* q, const double* o, const double* d, double& t )
{
  const double e1x=q[3]-q[0], e1y=q[4]-q[1], e1z=q[5]-q[2];
  const double e2x=q[6]-q[0], e2y=q[7]-q[1], e2z=q[8]-q[2];
  const double pvx = d[1]*e2z - d[2]*e2y;
  const double pvy = d[2]*e2x - d[0]*e2z;
  const double pvz = d[0]*e2y - d[1]*e2x;
  const double det = e1x*pvx + e1y*pvy + e1z*pvz;
  if( det == 0.0 ) return false;
  const double inv = 1.0/det;
  const double tvx = o[0]-q[0], tvy = o[1]-q[1], tvz = o[2]-q[2];
  const double u = ( tvx*pvx + tvy*pvy + tvz*pvz ) * inv;
  if( u < -UV_EPS || u > 1.0+UV_EPS ) return false;
  const double qvx = tvy*e1z - tvz*e1y;
  const double qvy = tvz*e1x - tvx*e1z;
  const double qvz = tvx*e1y - tvy*e1x;
  const double v = ( d[0]*qvx + d[1]*qvy + d[2]*qvz ) * inv;
  if( v < -UV_EPS || u+v > 1.0+UV_EPS ) return false;
  t = ( e2x*qvx + e2y*qvy + e2z*qvz ) * inv;
  return true;
}

static inline __m256d mm_dot3( __m256d ax, __m256d ay, __m256d az,
                               __m256d bx, __m256d by, __m256d bz )
{
  return _mm256_add_pd( _mm256_add_pd( _mm256_mul_pd(ax,bx), _mm256_mul_pd(ay,by) ),
                        _mm256_mul_pd(az,bz) );
}

/* 4-wide branchless Moller-Trumbore over one PreTri4 block: (v0,e1,e2) SoA,
 * aligned loads, no marshalling. Null-padding lanes have e1=e2=0 -> det=0 ->
 * a natural miss; hits are additionally id-filtered by the caller. */
static inline int mt4blk( const double* blk, const double* o,
                          const double* d, double tout[4] )
{
  const __m256d Ax  = _mm256_loadu_pd( blk      ), Ay  = _mm256_loadu_pd( blk +  4 ), Az  = _mm256_loadu_pd( blk +  8 );
  const __m256d e1x = _mm256_loadu_pd( blk + 12 ), e1y = _mm256_loadu_pd( blk + 16 ), e1z = _mm256_loadu_pd( blk + 20 );
  const __m256d e2x = _mm256_loadu_pd( blk + 24 ), e2y = _mm256_loadu_pd( blk + 28 ), e2z = _mm256_loadu_pd( blk + 32 );
  const __m256d ox = _mm256_set1_pd(o[0]), oy = _mm256_set1_pd(o[1]), oz = _mm256_set1_pd(o[2]);
  const __m256d dx = _mm256_set1_pd(d[0]), dy = _mm256_set1_pd(d[1]), dz = _mm256_set1_pd(d[2]);
  const __m256d one = _mm256_set1_pd(1.0), zero = _mm256_setzero_pd();
  const __m256d epn = _mm256_set1_pd(-UV_EPS), epp = _mm256_set1_pd(1.0+UV_EPS);

  const __m256d pvx = _mm256_sub_pd( _mm256_mul_pd(dy,e2z), _mm256_mul_pd(dz,e2y) );
  const __m256d pvy = _mm256_sub_pd( _mm256_mul_pd(dz,e2x), _mm256_mul_pd(dx,e2z) );
  const __m256d pvz = _mm256_sub_pd( _mm256_mul_pd(dx,e2y), _mm256_mul_pd(dy,e2x) );
  const __m256d det = mm_dot3( e1x,e1y,e1z, pvx,pvy,pvz );

  __m256d ok = _mm256_cmp_pd( det, zero, _CMP_NEQ_OQ );
  const __m256d inv = _mm256_div_pd( one,
      _mm256_blendv_pd( det, one, _mm256_cmp_pd( det, zero, _CMP_EQ_OQ ) ) );

  const __m256d tvx = _mm256_sub_pd(ox,Ax), tvy = _mm256_sub_pd(oy,Ay), tvz = _mm256_sub_pd(oz,Az);
  const __m256d u = _mm256_mul_pd( mm_dot3( tvx,tvy,tvz, pvx,pvy,pvz ), inv );
  ok = _mm256_and_pd( ok, _mm256_cmp_pd( u, epn, _CMP_GE_OQ ) );
  ok = _mm256_and_pd( ok, _mm256_cmp_pd( u, epp, _CMP_LE_OQ ) );

  const __m256d qvx = _mm256_sub_pd( _mm256_mul_pd(tvy,e1z), _mm256_mul_pd(tvz,e1y) );
  const __m256d qvy = _mm256_sub_pd( _mm256_mul_pd(tvz,e1x), _mm256_mul_pd(tvx,e1z) );
  const __m256d qvz = _mm256_sub_pd( _mm256_mul_pd(tvx,e1y), _mm256_mul_pd(tvy,e1x) );
  const __m256d v = _mm256_mul_pd( mm_dot3( dx,dy,dz, qvx,qvy,qvz ), inv );
  ok = _mm256_and_pd( ok, _mm256_cmp_pd( v, epn, _CMP_GE_OQ ) );
  ok = _mm256_and_pd( ok, _mm256_cmp_pd( _mm256_add_pd(u,v), epp, _CMP_LE_OQ ) );

  _mm256_storeu_pd( tout, _mm256_mul_pd( mm_dot3( e2x,e2y,e2z, qvx,qvy,qvz ), inv ) );
  return _mm256_movemask_pd( ok );
}

/* all 4 AABB slots in PURE FLOAT with per-ray sign-selected near/far bounds
 * (branchless); intervals WIDENED conservatively (relative + origin/scale
 * slack) so the float arithmetic can never lose a true hit -- exactness is
 * restored by the double MT tests at the leaves. Rays with a zero direction
 * component take the double scalar path instead. */
static inline int slot4AABBf( const float* nb, const __m128 of[3], const __m128 idf[3],
                              const int nOff[3], const int fOff[3], float slackT,
                              float wlo, float whi, float t0o[4], float t1o[4] )
{
  __m128 t0 = _mm_set1_ps( wlo ), t1 = _mm_set1_ps( whi );
  for( int a = 0; a < 3; ++a ) {
    t0 = _mm_max_ps( t0, _mm_mul_ps( _mm_sub_ps( _mm_loadu_ps( nb + nOff[a] ), of[a] ), idf[a] ) );
    t1 = _mm_min_ps( t1, _mm_mul_ps( _mm_sub_ps( _mm_loadu_ps( nb + fOff[a] ), of[a] ), idf[a] ) );
  }
  const __m128 sl  = _mm_set1_ps( slackT );
  const __m128 rel = _mm_set1_ps( 4e-7f );
  const __m128 am  = _mm_castsi128_ps( _mm_set1_epi32( 0x7FFFFFFF ) );
  t0 = _mm_sub_ps( t0, _mm_add_ps( _mm_mul_ps( _mm_and_ps( t0, am ), rel ), sl ) );
  t1 = _mm_add_ps( t1, _mm_add_ps( _mm_mul_ps( _mm_and_ps( t1, am ), rel ), sl ) );
  _mm_storeu_ps( t0o, t0 );  _mm_storeu_ps( t1o, t1 );
  return _mm_movemask_ps( _mm_cmple_ps( t0, t1 ) );
}

/* line vs slot volume: interval [t0,t1] clipped to the window; false = miss.
 * invd = 1/d precomputed ONCE per ray (12 divisions per node otherwise). */
static inline bool slotIntervalAABB( const float* nb, int k, const double* o,
                                     const double* d, const double* invd,
                                     double wlo, double whi,
                                     double& t0, double& t1 )
{
  t0 = wlo;  t1 = whi;
  for( int a = 0; a < 3; ++a ) {
    const double lo = (double)nb[ 8*a + k ], hi = (double)nb[ 8*a + 4 + k ];
    if( d[a] != 0.0 ) {
      const double ia = invd[a];
      double u = ( lo - o[a] )*ia, v = ( hi - o[a] )*ia;
      if( u > v ) { const double w = u; u = v; v = w; }
      if( u > t0 ) t0 = u;
      if( v < t1 ) t1 = v;
      if( t0 > t1 ) return false;
    } else if( o[a] < lo || o[a] > hi ) return false;
  }
  return true;
}

/* all 4 AABB slots in one AVX pass (double arithmetic on the float bounds ->
 * identical intervals to the scalar test); returns the hit mask */
static inline int slot4AABB( const float* nb, const double* o, const double* d,
                             const double* invd, double wlo, double whi,
                             double t0o[4], double t1o[4] )
{
  __m256d t0  = _mm256_set1_pd( wlo ), t1 = _mm256_set1_pd( whi );
  __m256d okm = _mm256_castsi256_pd( _mm256_set1_epi64x( -1 ) );
  for( int a = 0; a < 3; ++a ) {
    const __m256d lo = _mm256_cvtps_pd( _mm_loadu_ps( nb + 8*a ) );
    const __m256d hi = _mm256_cvtps_pd( _mm_loadu_ps( nb + 8*a + 4 ) );
    const __m256d ov = _mm256_set1_pd( o[a] );
    if( d[a] != 0.0 ) {
      const __m256d ia = _mm256_set1_pd( invd[a] );
      const __m256d u  = _mm256_mul_pd( _mm256_sub_pd( lo, ov ), ia );
      const __m256d v  = _mm256_mul_pd( _mm256_sub_pd( hi, ov ), ia );
      t0 = _mm256_max_pd( t0, _mm256_min_pd( u, v ) );
      t1 = _mm256_min_pd( t1, _mm256_max_pd( u, v ) );
    } else {
      okm = _mm256_and_pd( okm, _mm256_cmp_pd( lo, ov, _CMP_LE_OQ ) );
      okm = _mm256_and_pd( okm, _mm256_cmp_pd( ov, hi, _CMP_LE_OQ ) );
    }
  }
  okm = _mm256_and_pd( okm, _mm256_cmp_pd( t0, t1, _CMP_LE_OQ ) );
  _mm256_storeu_pd( t0o, t0 );  _mm256_storeu_pd( t1o, t1 );
  return _mm256_movemask_pd( okm );
}

static inline bool slotIntervalOBB( const float* nb, int k, const double* o,
                                    const double* d, double wlo, double whi,
                                    double& t0, double& t1 )
{
  t0 = wlo;  t1 = whi;
  for( int a = 0; a < 3; ++a ) {
    const double ax = (double)nb[ 4*(3*a  ) + k ];
    const double ay = (double)nb[ 4*(3*a+1) + k ];
    const double az = (double)nb[ 4*(3*a+2) + k ];
    const double oq = o[0]*ax + o[1]*ay + o[2]*az;
    const double dq = d[0]*ax + d[1]*ay + d[2]*az;
    const double lo = (double)nb[ 4*( 9+a) + k ], hi = (double)nb[ 4*(12+a) + k ];
    if( dq != 0.0 ) {
      const double ia = 1.0/dq;
      double u = ( lo - oq )*ia, v = ( hi - oq )*ia;
      if( u > v ) { const double w = u; u = v; v = w; }
      if( u > t0 ) t0 = u;
      if( v < t1 ) t1 = v;
      if( t0 > t1 ) return false;
    } else if( oq < lo || oq > hi ) return false;
  }
  return true;
}

/* kdop = the AABB slabs (same first-24 layout) + 4 diagonal slabs; od/dd4 are
 * the per-RAY projections of o and d onto the diagonals (1,±1,±1) */
static inline bool slotIntervalKDOP( const float* nb, int k, const double* o,
                                     const double* d, const double* invd,
                                     const double od[4], const double dd4[4],
                                     double wlo, double whi,
                                     double& t0, double& t1 )
{
  if( !slotIntervalAABB( nb, k, o, d, invd, wlo, whi, t0, t1 ) ) return false;
  for( int a = 0; a < 4; ++a ) {
    const double lo = (double)nb[ 4*(6+2*a) + k ], hi = (double)nb[ 4*(7+2*a) + k ];
    if( dd4[a] != 0.0 ) {
      const double ia = 1.0/dd4[a];
      double u = ( lo - od[a] )*ia, v = ( hi - od[a] )*ia;
      if( u > v ) { const double w = u; u = v; v = w; }
      if( u > t0 ) t0 = u;
      if( v < t1 ) t1 = v;
      if( t0 > t1 ) return false;
    } else if( od[a] < lo || od[a] > hi ) return false;
  }
  return true;
}

/* RSS rays: CONSERVATIVE test against the rect grown by r along its own axes
 * (an OBB superset of the swept volume; never misses, prunes slightly less) */
static inline bool slotIntervalRSS( const float* nb, int k, const double* o,
                                    const double* d, double wlo, double whi,
                                    double& t0, double& t1 )
{
  const double r = (double)nb[ 4*14 + k ];
  if( r < 0.0 ) return false;                        /* empty-slot marker */
  t0 = wlo;  t1 = whi;
  const double w0 = (double)nb[ 4*13 + k ];
  for( int a = 0; a < 3; ++a ) {
    const double ax = (double)nb[ 4*(3*a  ) + k ];
    const double ay = (double)nb[ 4*(3*a+1) + k ];
    const double az = (double)nb[ 4*(3*a+2) + k ];
    const double oq = o[0]*ax + o[1]*ay + o[2]*az;
    const double dq = d[0]*ax + d[1]*ay + d[2]*az;
    double lo, hi;
    if( a == 0 )      { lo = (double)nb[4* 9+k] - r;  hi = (double)nb[4*10+k] + r; }
    else if( a == 1 ) { lo = (double)nb[4*11+k] - r;  hi = (double)nb[4*12+k] + r; }
    else              { lo = w0 - r;                  hi = w0 + r;                 }
    if( dq != 0.0 ) {
      const double ia = 1.0/dq;
      double u = ( lo - oq )*ia, v = ( hi - oq )*ia;
      if( u > v ) { const double w = u; u = v; v = w; }
      if( u > t0 ) t0 = u;
      if( v < t1 ) t1 = v;
      if( t0 > t1 ) return false;
    } else if( oq < lo || oq > hi ) return false;
  }
  return true;
}

/* LSS rays: CONSERVATIVE test against the capsule's world-axis box */
static inline bool slotIntervalLSS( const float* nb, int k, const double* o,
                                    const double* d, const double* invd,
                                    double wlo, double whi,
                                    double& t0, double& t1 )
{
  const double r = (double)nb[ 4*6 + k ];
  if( r < 0.0 ) return false;
  t0 = wlo;  t1 = whi;
  for( int a = 0; a < 3; ++a ) {
    const double p0 = (double)nb[ 4*a + k ], p1 = (double)nb[ 4*(3+a) + k ];
    const double lo = ( p0 < p1 ? p0 : p1 ) - r;
    const double hi = ( p0 > p1 ? p0 : p1 ) + r;
    if( d[a] != 0.0 ) {
      const double ia = invd[a];
      double u = ( lo - o[a] )*ia, v = ( hi - o[a] )*ia;
      if( u > v ) { const double w = u; u = v; v = w; }
      if( u > t0 ) t0 = u;
      if( v < t1 ) t1 = v;
      if( t0 > t1 ) return false;
    } else if( o[a] < lo || o[a] > hi ) return false;
  }
  return true;
}

static inline bool slotIntervalSphere( const float* nb, int k, const double* o,
                                       const double* d, double dd, double wlo, double whi,
                                       double& t0, double& t1 )
{
  const double cx = (double)nb[k]-o[0], cy = (double)nb[4+k]-o[1], cz = (double)nb[8+k]-o[2];
  const double r  = (double)nb[12+k];
  if( r < 0.0 ) return false;                        /* empty-slot marker */
  const double oc2 = cx*cx + cy*cy + cz*cz;
  if( dd == 0.0 ) { if( oc2 > r*r ) return false;  t0 = wlo;  t1 = whi;  return true; }
  const double w = cx*d[0] + cy*d[1] + cz*d[2];
  const double disc = w*w - dd*( oc2 - r*r );
  if( disc < 0.0 ) return false;
  const double sq = std::sqrt( disc );
  t0 = ( w - sq )/dd;  t1 = ( w + sq )/dd;
  if( t0 < wlo ) t0 = wlo;  if( t1 > whi ) t1 = whi;
  return t0 <= t1;
}

void mexFunction( int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[] )
{
  if( nrhs < 3 )
    mexErrMsgIdAndTxt( "bvhIntersectRay_mx:nrhs", "expected rays, B, mode [, nthreads]." );
  if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) || mxGetN(prhs[0]) != 6 )
    mexErrMsgIdAndTxt( "bvhIntersectRay_mx:rays", "rays must be N x 6 double [p0,p1]." );

  const mwSize nR   = mxGetM( prhs[0] );
  const double* RY  = mxGetPr( prhs[0] );
  const int     mode = (int)mxGetScalar( prhs[2] );
  if( mode < 1 || mode > 4 )
    mexErrMsgIdAndTxt( "bvhIntersectRay_mx:mode", "mode must be 1..4." );
  int nt = ( nrhs > 3 ) ? (int)mxGetScalar( prhs[3] ) : 1;
  if( nt < 1 ) nt = 1;  if( nt > 64 ) nt = 64;

  /* ---- blob fields (v2 + packed leaves required) ---- */
  const mxArray* B = prhs[1];
  auto fld = [&]( const char* n ) -> const mxArray* {
    const mxArray* f = mxGetField( B, 0, n );
    if( !f ) mexErrMsgIdAndTxt( "bvhIntersectRay_mx:B",
                 "B lacks field '%s' (build with BVH, v2 packed blob).", n );
    return f;
  };
  const mxArray *fB4 = fld("bounds4"), *fC4 = fld("child4"), *fR4 = fld("srange");
  const mxArray *fkV = fld("pkV"), *fkT = fld("pkT"), *fkE = fld("pkE"), *fV = fld("vol");
  const mxArray *fP4 = fld("pk4"), *fPI = fld("pk4id"), *fS4 = fld("s4");
  if( !mxIsSingle(fB4) || !mxIsInt32(fC4) || !mxIsInt32(fR4) ||
      !mxIsDouble(fkV) || !mxIsInt32(fkT) || !mxIsInt32(fkE) ||
      !mxIsDouble(fP4) || !mxIsInt32(fPI) || !mxIsInt32(fS4) )
    mexErrMsgIdAndTxt( "bvhIntersectRay_mx:B", "wrong blob field types." );
  const int    vol = (int)mxGetScalar( fV );
  const mwSize nN  = mxGetN( fB4 );
  const mwSize nE  = mxGetN( fkV );
  const mwSize S   = ( vol == 1 ) ? 16 : ( vol == 2 ) ? 24 : ( vol == 3 ) ? 60 : ( vol == 4 ) ? 56 : ( vol == 5 ) ? 60 : 28;
  if( vol < 1 || vol > 6 )
    mexErrMsgIdAndTxt( "bvhIntersectRay_mx:B", "vol must be 1..6." );
  if( nN == 0 || mxGetM(fB4) != S || mxGetM(fC4) != 4 || mxGetN(fC4) != nN ||
      mxGetM(fR4) != 8 || mxGetN(fR4) != nN || mxGetM(fkV) != 12 ||
      mxGetNumberOfElements(fkT) != nE || mxGetNumberOfElements(fkE) != nE )
    mexErrMsgIdAndTxt( "bvhIntersectRay_mx:B", "inconsistent blob sizes." );

  const float*   W4 = (const float*)  mxGetData( fB4 );
  const int32_t* Wc = (const int32_t*)mxGetData( fC4 );
  const int32_t* Wr = (const int32_t*)mxGetData( fR4 );
  const double*  vv = mxGetPr( fkV );
  const int32_t* ety= (const int32_t*)mxGetData( fkT );
  const int32_t* eii= (const int32_t*)mxGetData( fkE );
  const double*  PK4 = mxGetPr( fP4 );
  const int32_t* PKI = (const int32_t*)mxGetData( fPI );
  const int32_t* Ws4 = (const int32_t*)mxGetData( fS4 );
  const mwSize   nBk = mxGetN( fP4 );
  if( mxGetM(fP4) != 36 || mxGetM(fPI) != 4 || mxGetN(fPI) != nBk ||
      mxGetM(fS4) != 8 || mxGetN(fS4) != nN )
    mexErrMsgIdAndTxt( "bvhIntersectRay_mx:B", "inconsistent PreTri4 pool." );
  for( mwSize i = 0; i < 4*nBk; ++i )
    if( PKI[i] < 0 || PKI[i] > (int32_t)nE )
      mexErrMsgIdAndTxt( "bvhIntersectRay_mx:B", "corrupt PreTri4 lane id." );

  for( mwSize i = 0; i < nN; ++i )                       /* bounds-check      */
    for( int k = 0; k < 4; ++k ) {
      const int32_t c = Wc[ i*4 + k ];
      if( c < -1 || c > (int32_t)nN )
        mexErrMsgIdAndTxt( "bvhIntersectRay_mx:B", "corrupt child index." );
      if( c != 0 ) {
        const int32_t lo = Wr[ i*8 + 2*k ], hi = Wr[ i*8 + 2*k + 1 ];
        if( lo < 1 || hi < lo || hi > (int32_t)nE )
          mexErrMsgIdAndTxt( "bvhIntersectRay_mx:B", "corrupt slot range." );
        const int32_t s0 = Ws4[ i*8 + 2*k ], sn = Ws4[ i*8 + 2*k + 1 ];
        if( sn < 0 || ( sn > 0 && ( s0 < 1 || (mwSize)( s0 + sn - 1 ) > nBk ) ) )
          mexErrMsgIdAndTxt( "bvhIntersectRay_mx:B", "corrupt PreTri4 slot range." );
      }
    }

  /* fused per-call node pool: bounds + children in one contiguous street */
  const size_t stride = (size_t)S*4 + 16;
  std::vector<char> fused( stride * nN );
  for( mwSize i = 0; i < nN; ++i ) {
    memcpy( &fused[ i*stride ], W4 + (size_t)i*S, (size_t)S*4 );
    memcpy( &fused[ i*stride + (size_t)S*4 ], Wc + (size_t)i*4, 16 );
  }
  const char* FZ = fused.data();

  /* scene scale for the float-slab conservative slack */
  float bScale = 1.f;
  if( vol == 2 )
    for( int r = 0; r < 24; ++r ) {
      const float b = std::fabs( W4[r] );
      if( b < 1e30f && b > bScale ) bScale = b;
    }

  const double wlo = ( mode == 4 ) ? ANY_TLO : -INF;
  const double whi = ( mode == 4 ) ? ANY_THI :  INF;
  const bool   avx = useAVX();

  std::vector<double>  bT( nR, mxGetNaN() );
  std::vector<int32_t> bC( nR, 0 );
  std::vector<std::vector<Hit>> allHits;                 /* mode 3 only */
  if( mode == 3 ) allHits.resize( nR );

  auto runRange = [&]( mwSize i0, mwSize i1 )
  {
    int32_t stk[192];  double stkT[192];
    for( mwSize q = i0; q < i1; ++q ) {
      const double o[3] = { RY[q], RY[q+nR], RY[q+2*nR] };
      const double d[3] = { RY[q+3*nR]-o[0], RY[q+4*nR]-o[1], RY[q+5*nR]-o[2] };
      const double dd   = d[0]*d[0] + d[1]*d[1] + d[2]*d[2];
      if( !std::isfinite(o[0]) || !std::isfinite(o[1]) || !std::isfinite(o[2]) ||
          !std::isfinite(d[0]) || !std::isfinite(d[1]) || !std::isfinite(d[2]) )
        continue;
      const double invd[3] = { ( d[0] != 0.0 ) ? 1.0/d[0] : 0.0,
                               ( d[1] != 0.0 ) ? 1.0/d[1] : 0.0,
                               ( d[2] != 0.0 ) ? 1.0/d[2] : 0.0 };
      const double od[4]  = { o[0]+o[1]+o[2], o[0]-o[1]+o[2],
                              o[0]+o[1]-o[2], o[0]-o[1]-o[2] };   /* kdop diags */
      const double dd4[4] = { d[0]+d[1]+d[2], d[0]-d[1]+d[2],
                              d[0]+d[1]-d[2], d[0]-d[1]-d[2] };

      /* float fast-slab setup (aabb, all direction components nonzero) */
      const bool fslab = ( vol == 2 ) && avx &&
                         d[0] != 0.0 && d[1] != 0.0 && d[2] != 0.0;
      __m128 ofv[3], idv[3];
      int nOff[3], fOff[3];
      float slackT = 0.f;
      if( fslab ) {
        for( int a = 0; a < 3; ++a ) {
          ofv[a] = _mm_set1_ps( (float)o[a] );
          idv[a] = _mm_set1_ps( (float)invd[a] );
          if( d[a] >= 0.0 ) { nOff[a] = 8*a;     fOff[a] = 8*a + 4; }
          else              { nOff[a] = 8*a + 4; fOff[a] = 8*a;     }
          slackT += ( (float)std::fabs( o[a] ) + bScale ) * 4e-7f * (float)std::fabs( invd[a] );
        }
        slackT += 1e-30f;
      }

      double  hitT = ( mode == 2 ) ? -INF : INF;
      int32_t hitC = 0;
      bool    found = false, stop = false;

      int top = 1;  stk[0] = 0;  stkT[0] = ( mode == 2 ) ? INF : -INF;
      while( top && !stop ) {
        --top;
        if( mode == 1 && stkT[top] >= hitT ) continue;   /* stale-pop culls  */
        if( mode == 2 && stkT[top] <= hitT ) continue;
        const int32_t  ni = stk[top];
        const char*    nz = FZ + (size_t)ni * stride;
        const float*   nb = (const float*)nz;
        const int32_t* nc = (const int32_t*)( nz + (size_t)S*4 );
        const int32_t* nr = Wr + (size_t)ni * 8;
        const int32_t* n4 = Ws4 + (size_t)ni * 8;

        double key[4];  int act[4];  int na = 0;
        if( fslab ) {                    /* pure-float sign-selected slabs */
          float t0v[4], t1v[4];
          const int mm = slot4AABBf( nb, ofv, idv, nOff, fOff, slackT,
                                     (float)wlo, (float)whi, t0v, t1v );
          for( int k = 0; k < 4; ++k ) {
            if( !( mm & (1<<k) ) || nc[k] == 0 ) continue;
            if( mode == 1 && (double)t0v[k] >= hitT ) continue;
            if( mode == 2 && (double)t1v[k] <= hitT ) continue;
            key[na] = ( mode == 2 ) ? -(double)t1v[k] : (double)t0v[k];
            act[na] = k;  ++na;
          }
        } else if( vol == 2 && avx ) {                   /* 4 slots, one pass */
          double t0v[4], t1v[4];
          const int mm = slot4AABB( nb, o, d, invd, wlo, whi, t0v, t1v );
          for( int k = 0; k < 4; ++k ) {
            if( !( mm & (1<<k) ) || nc[k] == 0 ) continue;
            if( mode == 1 && t0v[k] >= hitT ) continue;
            if( mode == 2 && t1v[k] <= hitT ) continue;
            key[na] = ( mode == 2 ) ? -t1v[k] : t0v[k];
            act[na] = k;  ++na;
          }
        } else {
          for( int k = 0; k < 4; ++k ) {
            if( nc[k] == 0 ) continue;
            double t0, t1;
            bool h;
            switch( vol ) {
              case 1:  h = slotIntervalSphere( nb, k, o, d, dd, wlo, whi, t0, t1 );          break;
              case 3:  h = slotIntervalOBB( nb, k, o, d, wlo, whi, t0, t1 );                 break;
              case 4:  h = slotIntervalKDOP( nb, k, o, d, invd, od, dd4, wlo, whi, t0, t1 ); break;
              case 5:  h = slotIntervalRSS( nb, k, o, d, wlo, whi, t0, t1 );                 break;
              case 6:  h = slotIntervalLSS( nb, k, o, d, invd, wlo, whi, t0, t1 );           break;
              default: h = slotIntervalAABB( nb, k, o, d, invd, wlo, whi, t0, t1 );          break;
            }
            if( !h ) continue;
            if( mode == 1 && t0 >= hitT ) continue;      /* beyond the best  */
            if( mode == 2 && t1 <= hitT ) continue;
            key[na] = ( mode == 2 ) ? -t1 : t0;          /* ascending order  */
            act[na] = k;  ++na;
          }
        }
        if( ( mode == 1 || mode == 2 ) && na > 1 ) {     /* insertion sort   */
          for( int a = 1; a < na; ++a ) {
            const double kk = key[a];  const int ak = act[a];
            int b = a-1;
            while( b >= 0 && key[b] > kk ) { key[b+1]=key[b]; act[b+1]=act[b]; --b; }
            key[b+1]=kk; act[b+1]=ak;
          }
        }

        for( int a = 0; a < na && !stop; ++a ) {         /* leaves in order  */
          const int k = act[a];
          if( nc[k] != -1 ) continue;
          if( mode == 1 &&  key[a] >= hitT ) continue;   /* re-check         */
          if( mode == 2 && -key[a] <= hitT ) continue;
          const int32_t lo = nr[2*k]-1, hi = nr[2*k+1]-1;
          const int32_t s0 = n4[2*k]-1, sn = n4[2*k+1];

          /* (with PreTri4's aligned loads the 4-wide kernel also wins for
           * 'any': one cheap block test beats 1-4 scalar tests) */
          if( avx && sn > 0 ) {
            for( int32_t b = 0; b < sn && !stop; ++b ) {
              const double*  blk = PK4 + (size_t)( s0 + b )*36;
              const int32_t* ids = PKI + (size_t)( s0 + b )*4;
              double tv4[4];
              int mm = mt4blk( blk, o, d, tv4 );
              while( mm ) {
                const int l = ( mm & 1 ) ? 0 : ( ( mm & 2 ) ? 1 : ( ( mm & 4 ) ? 2 : 3 ) );
                mm &= mm - 1;
                if( ids[l] <= 0 ) continue;              /* padding lane */
                const double t = tv4[l];
                const int32_t ce = eii[ ids[l]-1 ];
                switch( mode ) {
                  case 1: if( t < hitT ) { hitT = t; hitC = ce; found = true; } break;
                  case 2: if( t > hitT ) { hitT = t; hitC = ce; found = true; } break;
                  case 3: allHits[q].push_back( { t, ce } );                    break;
                  case 4: if( t > ANY_TLO && t < ANY_THI ) {
                            hitT = t; hitC = ce; found = true; stop = true; }   break;
                }
                if( stop ) break;
              }
            }
          } else {
            for( int32_t j = lo; j <= hi && !stop; ++j ) {
              if( ety[j] != 3 ) continue;
              double t;
              if( !mtHit( vv + (size_t)j*12, o, d, t ) ) continue;
              switch( mode ) {
                case 1: if( t < hitT ) { hitT = t; hitC = eii[j]; }  found = true; break;
                case 2: if( t > hitT ) { hitT = t; hitC = eii[j]; }  found = true; break;
                case 3: allHits[q].push_back( { t, eii[j] } );                     break;
                case 4: if( t > ANY_TLO && t < ANY_THI ) {
                          hitT = t; hitC = eii[j]; found = true; stop = true; }    break;
              }
            }
          }
        }
        for( int a = na-1; a >= 0 && !stop; --a ) {      /* internals, far 1st */
          const int k = act[a];
          if( nc[k] <= 0 ) continue;
          if( top > 188 )
            mexErrMsgIdAndTxt( "bvhIntersectRay_mx:stack", "stack overflow (corrupt blob?)." );
          stk[top] = nc[k]-1;
          stkT[top] = ( mode == 2 ) ? -key[a] : key[a];
          ++top;
        }
      }
      if( found ) { bT[q] = hitT;  bC[q] = hitC; }
      if( mode == 3 )
        std::sort( allHits[q].begin(), allHits[q].end(),
                   []( const Hit& a, const Hit& b ){ return a.t < b.t; } );
    }
  };

#ifdef _OPENMP
  if( nt > 1 && mode != 3 && nR > 64 ) {
    #pragma omp parallel num_threads( nt )
    {
      const int tid = omp_get_thread_num();
      const int nth = omp_get_num_threads();
      const mwSize per = ( nR + nth - 1 ) / nth;
      const mwSize i0  = (mwSize)tid * per;
      const mwSize i1  = ( i0 + per < nR ) ? i0 + per : nR;
      if( i0 < i1 ) runRange( i0, i1 );
    }
  } else runRange( 0, nR );
#else
  runRange( 0, nR );
#endif

  if( mode == 3 ) {                                      /* one row per hit */
    mwSize nH = 0;
    for( mwSize q = 0; q < nR; ++q ) nH += allHits[q].size();
    plhs[0] = mxCreateDoubleMatrix( nH, 3, mxREAL );
    mxArray* mc = mxCreateDoubleMatrix( nH, 1, mxREAL );
    mxArray* mt = mxCreateDoubleMatrix( nH, 1, mxREAL );
    mxArray* mr = mxCreateDoubleMatrix( nH, 1, mxREAL );
    double *ox = mxGetPr(plhs[0]), *oc = mxGetPr(mc), *ot = mxGetPr(mt), *orr = mxGetPr(mr);
    mwSize h = 0;
    for( mwSize q = 0; q < nR; ++q ) {
      const double o[3] = { RY[q], RY[q+nR], RY[q+2*nR] };
      const double d[3] = { RY[q+3*nR]-o[0], RY[q+4*nR]-o[1], RY[q+5*nR]-o[2] };
      for( const Hit& hh : allHits[q] ) {
        ox[h] = o[0]+hh.t*d[0];  ox[h+nH] = o[1]+hh.t*d[1];  ox[h+2*nH] = o[2]+hh.t*d[2];
        oc[h] = hh.cell;  ot[h] = hh.t;  orr[h] = (double)(q+1);
        ++h;
      }
    }
    if( nlhs > 1 ) plhs[1] = mc; else mxDestroyArray( mc );
    if( nlhs > 2 ) plhs[2] = mt; else mxDestroyArray( mt );
    if( nlhs > 3 ) plhs[3] = mr; else mxDestroyArray( mr );
  } else {                                               /* one row per ray */
    plhs[0] = mxCreateDoubleMatrix( nR, 3, mxREAL );
    mxArray* mc = mxCreateDoubleMatrix( nR, 1, mxREAL );
    mxArray* mt = mxCreateDoubleMatrix( nR, 1, mxREAL );
    mxArray* mr = mxCreateDoubleMatrix( nR, 1, mxREAL );
    double *ox = mxGetPr(plhs[0]), *oc = mxGetPr(mc), *ot = mxGetPr(mt), *orr = mxGetPr(mr);
    for( mwSize q = 0; q < nR; ++q ) {
      const double t = bT[q];
      if( bC[q] > 0 ) {
        const double o0=RY[q], o1=RY[q+nR], o2=RY[q+2*nR];
        ox[q]      = o0 + t*( RY[q+3*nR]-o0 );
        ox[q+nR]   = o1 + t*( RY[q+4*nR]-o1 );
        ox[q+2*nR] = o2 + t*( RY[q+5*nR]-o2 );
      } else {
        ox[q] = ox[q+nR] = ox[q+2*nR] = mxGetNaN();
      }
      oc[q] = (double)bC[q];  ot[q] = t;  orr[q] = (double)(q+1);
    }
    if( nlhs > 1 ) plhs[1] = mc; else mxDestroyArray( mc );
    if( nlhs > 2 ) plhs[2] = mt; else mxDestroyArray( mt );
    if( nlhs > 3 ) plhs[3] = mr; else mxDestroyArray( mr );
  }
}
