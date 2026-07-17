/* bvhClosestElement_mx  --  compiled core of bvhClosestElement.
 *
 *   [e, cp, d] = bvhClosestElement_mx( P , B , nthreads , Dmax )
 *
 *     P    nP x 3  double  query points (padded to 3 columns)
 *     B    struct          BVH blob (version >= 3): wide 4-ary nodes
 *                          (bounds4/child4/srange; sphere or AABB slots in
 *                          conservative float, pruning arithmetic in double)
 *                          + packed leaf data (pkV/pkS/pkT/pkE). The blob is
 *                          SELF-CONTAINED: the mesh itself is not needed here.
 *     nthreads scalar      OpenMP threads over the query points (default 1)
 *     Dmax scalar          search radius (default Inf): the best-so-far bound
 *                          is SEEDED with Dmax, so everything farther prunes
 *                          from the very root -- a point beyond Dmax costs one
 *                          node visit. Elements at d < Dmax are returned; a
 *                          point with no element within Dmax gives e = 0,
 *                          d = Inf, cp = NaN (non-finite points: d = NaN).
 *
 *   Celltypes by nonzero node count: 1 vertex, 2 segment, 3 triangle,
 *   4 TETRAHEDRON (always; quads are not a thing here).
 *
 *   Performance notes:
 *     - sqrt-free pruning:  |p-c|^2 > (r+best)^2 (spheres) / squared box
 *       distance vs best2 (AABBs); best updated by sqrt only on improvement.
 *     - stack entries carry their pruning key: stale nodes are culled again
 *       at pop time against the CURRENT best.
 *     - nearer child first; element-sphere pre-prune before exact tests.
 *     - WARM START: the previous point's winner seeds best2 -- with
 *       Morton-ordered queries consecutive points are spatial neighbours.
 *     - 4-wide branchless AVX Ericson kernel on all-triangle leaves
 *       (runtime CPUID, scalar fallback; degenerate lanes redone scalar).
 *     - iterative traversal, fixed stacks, all indices bounds-checked at
 *       entry: a corrupt B errors out, it cannot access-violate.
 *
 *   Compile (MSVC):  mex COMPFLAGS="$COMPFLAGS /openmp" bvhClosestElement_mx.cpp
 *
 * See also bvhClosestElement, BVH, bvhIntersectRay_mx.
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

struct Elem { double cx, cy, cz, r; };     /* element sphere, pkS column */

static inline double sq( double v ) { return v*v; }

/* ---------------------------------------------------------------- primitives */

static inline double d2Point( const double* a, const double* p, double* cp )
{
  cp[0]=a[0]; cp[1]=a[1]; cp[2]=a[2];
  return sq(p[0]-a[0]) + sq(p[1]-a[1]) + sq(p[2]-a[2]);
}

static inline double d2Seg( const double* a, const double* b, const double* p, double* cp )
{
  const double v0=b[0]-a[0], v1=b[1]-a[1], v2=b[2]-a[2];
  const double vv=v0*v0+v1*v1+v2*v2;
  double t = (p[0]-a[0])*v0 + (p[1]-a[1])*v1 + (p[2]-a[2])*v2;
  t = ( vv > 0.0 ) ? t/vv : 0.0;
  if( t < 0.0 ) t = 0.0; else if( t > 1.0 ) t = 1.0;
  cp[0]=a[0]+t*v0; cp[1]=a[1]+t*v1; cp[2]=a[2]+t*v2;
  return sq(p[0]-cp[0]) + sq(p[1]-cp[1]) + sq(p[2]-cp[2]);
}

/* closest point on a triangle: Ericson, "Real-Time Collision Detection" 5.1.5 */
static double d2Tri( const double* A, const double* B, const double* C,
                     const double* p, double* cp )
{
  const double ab0=B[0]-A[0], ab1=B[1]-A[1], ab2=B[2]-A[2];
  const double ac0=C[0]-A[0], ac1=C[1]-A[1], ac2=C[2]-A[2];
  const double ap0=p[0]-A[0], ap1=p[1]-A[1], ap2=p[2]-A[2];

  const double d1 = ab0*ap0 + ab1*ap1 + ab2*ap2;
  const double d2 = ac0*ap0 + ac1*ap1 + ac2*ap2;
  if( d1 <= 0.0 && d2 <= 0.0 ) return d2Point( A, p, cp );          /* vertex A */

  const double bp0=p[0]-B[0], bp1=p[1]-B[1], bp2=p[2]-B[2];
  const double d3 = ab0*bp0 + ab1*bp1 + ab2*bp2;
  const double d4 = ac0*bp0 + ac1*bp1 + ac2*bp2;
  if( d3 >= 0.0 && d4 <= d3 ) return d2Point( B, p, cp );           /* vertex B */

  const double vc = d1*d4 - d3*d2;
  if( vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 ) {                       /* edge AB */
    const double den = d1 - d3;
    const double t = ( den != 0.0 ) ? d1/den : 0.0;
    cp[0]=A[0]+t*ab0; cp[1]=A[1]+t*ab1; cp[2]=A[2]+t*ab2;
    return sq(p[0]-cp[0]) + sq(p[1]-cp[1]) + sq(p[2]-cp[2]);
  }

  const double cq0=p[0]-C[0], cq1=p[1]-C[1], cq2=p[2]-C[2];
  const double d5 = ab0*cq0 + ab1*cq1 + ab2*cq2;
  const double d6 = ac0*cq0 + ac1*cq1 + ac2*cq2;
  if( d6 >= 0.0 && d5 <= d6 ) return d2Point( C, p, cp );           /* vertex C */

  const double vb = d5*d2 - d1*d6;
  if( vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 ) {                       /* edge AC */
    const double den = d2 - d6;
    const double t = ( den != 0.0 ) ? d2/den : 0.0;
    cp[0]=A[0]+t*ac0; cp[1]=A[1]+t*ac1; cp[2]=A[2]+t*ac2;
    return sq(p[0]-cp[0]) + sq(p[1]-cp[1]) + sq(p[2]-cp[2]);
  }

  const double va = d3*d6 - d5*d4;
  if( va <= 0.0 && (d4-d3) >= 0.0 && (d5-d6) >= 0.0 ) {             /* edge BC */
    const double den = (d4-d3) + (d5-d6);
    const double t = ( den != 0.0 ) ? (d4-d3)/den : 0.0;
    cp[0]=B[0]+t*(C[0]-B[0]); cp[1]=B[1]+t*(C[1]-B[1]); cp[2]=B[2]+t*(C[2]-B[2]);
    return sq(p[0]-cp[0]) + sq(p[1]-cp[1]) + sq(p[2]-cp[2]);
  }

  const double den = va + vb + vc;                                  /* interior */
  if( den != 0.0 && std::isfinite( den ) ) {
    const double v = vb/den, w = vc/den;
    cp[0]=A[0]+v*ab0+w*ac0; cp[1]=A[1]+v*ab1+w*ac1; cp[2]=A[2]+v*ab2+w*ac2;
    return sq(p[0]-cp[0]) + sq(p[1]-cp[1]) + sq(p[2]-cp[2]);
  }

  /* fully degenerate triangle: closest of the three edges */
  double c2[3], c3[3];
  double q1 = d2Seg( A, B, p, cp );
  double q2 = d2Seg( A, C, p, c2 );
  double q3 = d2Seg( B, C, p, c3 );
  if( q2 < q1 ) { q1=q2; cp[0]=c2[0]; cp[1]=c2[1]; cp[2]=c2[2]; }
  if( q3 < q1 ) { q1=q3; cp[0]=c3[0]; cp[1]=c3[1]; cp[2]=c3[2]; }
  return q1;
}

static inline double det3( const double* u, const double* v, const double* w )
{
  return u[0]*(v[1]*w[2]-v[2]*w[1])
       - u[1]*(v[0]*w[2]-v[2]*w[0])
       + u[2]*(v[0]*w[1]-v[1]*w[0]);
}

static double d2Tet( const double* A, const double* B, const double* C,
                     const double* D, const double* p, double* cp )
{
  const double BA[3]={B[0]-A[0],B[1]-A[1],B[2]-A[2]};
  const double CA[3]={C[0]-A[0],C[1]-A[1],C[2]-A[2]};
  const double DA[3]={D[0]-A[0],D[1]-A[1],D[2]-A[2]};
  const double d0 = det3( BA, CA, DA );
  if( d0 != 0.0 && std::isfinite( d0 ) ) {
    const double Bp[3]={B[0]-p[0],B[1]-p[1],B[2]-p[2]};
    const double Cp[3]={C[0]-p[0],C[1]-p[1],C[2]-p[2]};
    const double Dp[3]={D[0]-p[0],D[1]-p[1],D[2]-p[2]};
    const double pA[3]={p[0]-A[0],p[1]-A[1],p[2]-A[2]};
    const double l1 = det3( Bp, Cp, Dp ) / d0;
    const double l2 = det3( pA, CA, DA ) / d0;
    const double l3 = det3( BA, pA, DA ) / d0;
    const double l4 = 1.0 - l1 - l2 - l3;
    const double tol = -1e-12;
    if( l1 >= tol && l2 >= tol && l3 >= tol && l4 >= tol ) {        /* inside */
      cp[0]=p[0]; cp[1]=p[1]; cp[2]=p[2];
      return 0.0;
    }
  }
  double c2[3];
  double q  = d2Tri( A, B, C, p, cp );
  double q2 = d2Tri( A, B, D, p, c2 );
  if( q2 < q ) { q=q2; cp[0]=c2[0]; cp[1]=c2[1]; cp[2]=c2[2]; }
  q2 = d2Tri( A, C, D, p, c2 );
  if( q2 < q ) { q=q2; cp[0]=c2[0]; cp[1]=c2[1]; cp[2]=c2[2]; }
  q2 = d2Tri( B, C, D, p, c2 );
  if( q2 < q ) { q=q2; cp[0]=c2[0]; cp[1]=c2[1]; cp[2]=c2[2]; }
  return q;
}

/* exact test against packed element j (verts at v = pkV + 12*j).
 * 4 nonzero nodes ALWAYS mean a tetrahedron. */
static inline double d2Elem( const double* v, int k, const double* p, double* cp )
{
  switch( k ) {
    case 1: return d2Point( v, p, cp );
    case 2: return d2Seg( v, v+3, p, cp );
    case 3: return d2Tri( v, v+3, v+6, p, cp );
    case 4: return d2Tet( v, v+3, v+6, v+9, p, cp );
  }
  return INF;   /* k == 0: empty 0-padded row */
}

/* ------------------------------------------------- 4-wide triangle kernel */

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

static inline __m256d mm_dot3( __m256d ax, __m256d ay, __m256d az,
                               __m256d bx, __m256d by, __m256d bz )
{
  return _mm256_add_pd( _mm256_add_pd( _mm256_mul_pd(ax,bx), _mm256_mul_pd(ay,by) ),
                        _mm256_mul_pd(az,bz) );
}

/* branchless Ericson over one PreTri4 block: (v0,e1,e2) in 4-wide SoA --
 * ALIGNED sequential loads, no marshalling, edges precomputed at build.
 * Same region order and same per-region arithmetic as the scalar d2Tri, so
 * selected lanes match it to the last bit; genuinely degenerate lanes are
 * redone with the scalar fallback. Null-padding lanes (id 0) produce a finite
 * garbage result that the CALLER filters by id. */
static void tri4blk( const double* blk, const double* p,
                     double d2o[4], double cpo[12] )
{
  const __m256d Ax  = _mm256_loadu_pd( blk      ), Ay  = _mm256_loadu_pd( blk +  4 ), Az  = _mm256_loadu_pd( blk +  8 );
  const __m256d abx = _mm256_loadu_pd( blk + 12 ), aby = _mm256_loadu_pd( blk + 16 ), abz = _mm256_loadu_pd( blk + 20 );
  const __m256d acx = _mm256_loadu_pd( blk + 24 ), acy = _mm256_loadu_pd( blk + 28 ), acz = _mm256_loadu_pd( blk + 32 );
  const __m256d Bx = _mm256_add_pd( Ax, abx ), By = _mm256_add_pd( Ay, aby ), Bz = _mm256_add_pd( Az, abz );
  const __m256d Cx = _mm256_add_pd( Ax, acx ), Cy = _mm256_add_pd( Ay, acy ), Cz = _mm256_add_pd( Az, acz );
  const __m256d px = _mm256_set1_pd( p[0] ), py = _mm256_set1_pd( p[1] ), pz = _mm256_set1_pd( p[2] );
  const __m256d zero = _mm256_setzero_pd(), one = _mm256_set1_pd( 1.0 );

  const __m256d apx = _mm256_sub_pd(px,Ax), apy = _mm256_sub_pd(py,Ay), apz = _mm256_sub_pd(pz,Az);
  const __m256d d1 = mm_dot3( abx,aby,abz, apx,apy,apz );
  const __m256d d2 = mm_dot3( acx,acy,acz, apx,apy,apz );
  const __m256d bpx = _mm256_sub_pd(px,Bx), bpy = _mm256_sub_pd(py,By), bpz = _mm256_sub_pd(pz,Bz);
  const __m256d d3 = mm_dot3( abx,aby,abz, bpx,bpy,bpz );
  const __m256d d4 = mm_dot3( acx,acy,acz, bpx,bpy,bpz );
  const __m256d cqx = _mm256_sub_pd(px,Cx), cqy = _mm256_sub_pd(py,Cy), cqz = _mm256_sub_pd(pz,Cz);
  const __m256d d5 = mm_dot3( abx,aby,abz, cqx,cqy,cqz );
  const __m256d d6 = mm_dot3( acx,acy,acz, cqx,cqy,cqz );
  const __m256d vc = _mm256_sub_pd( _mm256_mul_pd(d1,d4), _mm256_mul_pd(d3,d2) );
  const __m256d vb = _mm256_sub_pd( _mm256_mul_pd(d5,d2), _mm256_mul_pd(d1,d6) );
  const __m256d va = _mm256_sub_pd( _mm256_mul_pd(d3,d6), _mm256_mul_pd(d5,d4) );

#define LE(a,b) _mm256_cmp_pd( a, b, _CMP_LE_OQ )
#define GE(a,b) _mm256_cmp_pd( a, b, _CMP_GE_OQ )
  const __m256d mA  = _mm256_and_pd( LE(d1,zero), LE(d2,zero) );
  const __m256d mB  = _mm256_and_pd( GE(d3,zero), LE(d4,d3) );
  const __m256d mAB = _mm256_and_pd( _mm256_and_pd( LE(vc,zero), GE(d1,zero) ), LE(d3,zero) );
  const __m256d mC  = _mm256_and_pd( GE(d6,zero), LE(d5,d6) );
  const __m256d mAC = _mm256_and_pd( _mm256_and_pd( LE(vb,zero), GE(d2,zero) ), LE(d6,zero) );
  const __m256d d43 = _mm256_sub_pd( d4, d3 ), d56 = _mm256_sub_pd( d5, d6 );
  const __m256d mBC = _mm256_and_pd( _mm256_and_pd( LE(va,zero), GE(d43,zero) ), GE(d56,zero) );
#undef LE
#undef GE

  /* safe divisions (unselected lanes discard the garbage) */
  const __m256d eqz1 = _mm256_cmp_pd( _mm256_sub_pd(d1,d3), zero, _CMP_EQ_OQ );
  const __m256d denAB = _mm256_blendv_pd( _mm256_sub_pd(d1,d3), one, eqz1 );
  const __m256d tAB = _mm256_div_pd( d1, denAB );
  const __m256d eqz2 = _mm256_cmp_pd( _mm256_sub_pd(d2,d6), zero, _CMP_EQ_OQ );
  const __m256d denAC = _mm256_blendv_pd( _mm256_sub_pd(d2,d6), one, eqz2 );
  const __m256d tAC = _mm256_div_pd( d2, denAC );
  const __m256d sBC = _mm256_add_pd( d43, d56 );
  const __m256d eqz3 = _mm256_cmp_pd( sBC, zero, _CMP_EQ_OQ );
  const __m256d tBC = _mm256_div_pd( d43, _mm256_blendv_pd( sBC, one, eqz3 ) );
  const __m256d denI = _mm256_add_pd( _mm256_add_pd( va, vb ), vc );
  const __m256d eqzI = _mm256_cmp_pd( denI, zero, _CMP_EQ_OQ );
  const __m256d rden = _mm256_div_pd( one, _mm256_blendv_pd( denI, one, eqzI ) );
  const __m256d vI = _mm256_mul_pd( vb, rden ), wI = _mm256_mul_pd( vc, rden );

  /* priority cascade: A, B, AB, C, AC, BC, interior (Ericson order) */
  __m256d cx = Ax, cy = Ay, cz = Az;
  __m256d done = mA;
  __m256d sel;
#define PICK(m,qx,qy,qz) \
  sel = _mm256_andnot_pd( done, m ); \
  cx = _mm256_blendv_pd( cx, qx, sel ); \
  cy = _mm256_blendv_pd( cy, qy, sel ); \
  cz = _mm256_blendv_pd( cz, qz, sel ); \
  done = _mm256_or_pd( done, m );
  PICK( mB, Bx, By, Bz );
  PICK( mAB, _mm256_add_pd(Ax,_mm256_mul_pd(tAB,abx)),
             _mm256_add_pd(Ay,_mm256_mul_pd(tAB,aby)),
             _mm256_add_pd(Az,_mm256_mul_pd(tAB,abz)) );
  PICK( mC, Cx, Cy, Cz );
  PICK( mAC, _mm256_add_pd(Ax,_mm256_mul_pd(tAC,acx)),
             _mm256_add_pd(Ay,_mm256_mul_pd(tAC,acy)),
             _mm256_add_pd(Az,_mm256_mul_pd(tAC,acz)) );
  PICK( mBC, _mm256_add_pd(Bx,_mm256_mul_pd(tBC,_mm256_sub_pd(Cx,Bx))),
             _mm256_add_pd(By,_mm256_mul_pd(tBC,_mm256_sub_pd(Cy,By))),
             _mm256_add_pd(Bz,_mm256_mul_pd(tBC,_mm256_sub_pd(Cz,Bz))) );
  /* interior: everything not yet done */
  sel = _mm256_andnot_pd( done, _mm256_castsi256_pd( _mm256_set1_epi64x( -1 ) ) );
  cx = _mm256_blendv_pd( cx, _mm256_add_pd(Ax,_mm256_add_pd(_mm256_mul_pd(vI,abx),_mm256_mul_pd(wI,acx))), sel );
  cy = _mm256_blendv_pd( cy, _mm256_add_pd(Ay,_mm256_add_pd(_mm256_mul_pd(vI,aby),_mm256_mul_pd(wI,acy))), sel );
  cz = _mm256_blendv_pd( cz, _mm256_add_pd(Az,_mm256_add_pd(_mm256_mul_pd(vI,abz),_mm256_mul_pd(wI,acz))), sel );

  const __m256d dx = _mm256_sub_pd( px, cx ), dy = _mm256_sub_pd( py, cy ), dz = _mm256_sub_pd( pz, cz );
  const __m256d dd = mm_dot3( dx,dy,dz, dx,dy,dz );

  _mm256_storeu_pd( d2o, dd );
  double cxx[4], cyy[4], czz[4];
  _mm256_storeu_pd( cxx, cx );  _mm256_storeu_pd( cyy, cy );  _mm256_storeu_pd( czz, cz );
  for( int l = 0; l < 4; ++l ) { cpo[3*l]=cxx[l]; cpo[3*l+1]=cyy[l]; cpo[3*l+2]=czz[l]; }

  /* degenerate lanes (interior selected with a zeroed denominator): redo scalar */
  const int degm = _mm256_movemask_pd( _mm256_and_pd( sel, eqzI ) );
  if( degm ) {
    for( int l = 0; l < 4; ++l )
      if( degm & (1<<l) ) {
        const double A3[3] = { blk[   l], blk[ 4+l], blk[ 8+l] };
        const double B3[3] = { A3[0]+blk[12+l], A3[1]+blk[16+l], A3[2]+blk[20+l] };
        const double C3[3] = { A3[0]+blk[24+l], A3[1]+blk[28+l], A3[2]+blk[32+l] };
        d2o[l] = d2Tri( A3, B3, C3, p, cpo+3*l );
      }
  }
#undef PICK
}

/* ---------------------------------------------------------------- Morton */

static inline uint32_t expandBits( uint32_t v )
{
  v = ( v * 0x00010001u ) & 0xFF0000FFu;
  v = ( v * 0x00000101u ) & 0x0F00F00Fu;
  v = ( v * 0x00000011u ) & 0xC30C30C3u;
  v = ( v * 0x00000005u ) & 0x49249249u;
  return v;
}

/* ---------------------------------------------------------------- gateway */

void mexFunction( int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[] )
{
  if( nrhs < 2 )
    mexErrMsgIdAndTxt( "bvhClosestElement_mx:nrhs",
                       "expected P, B [, nthreads, Dmax]." );
  if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) || mxIsSparse(prhs[0]) ||
      mxGetN(prhs[0]) != 3 )
    mexErrMsgIdAndTxt( "bvhClosestElement_mx:P", "P must be nP x 3 double (pad it first)." );
  if( !mxIsStruct( prhs[1] ) )
    mexErrMsgIdAndTxt( "bvhClosestElement_mx:B", "B must be a struct (BVH blob)." );

  const mwSize  nP = mxGetM( prhs[0] );
  const double* P  = mxGetPr( prhs[0] );

  int nt = ( nrhs > 2 ) ? (int)mxGetScalar( prhs[2] ) : 1;
  if( nt < 1 ) nt = 1;  if( nt > 64 ) nt = 64;
  const double Dmax = ( nrhs > 3 ) ? mxGetScalar( prhs[3] ) : INF;
  if( !( Dmax >= 0.0 ) )   /* also rejects NaN */
    mexErrMsgIdAndTxt( "bvhClosestElement_mx:Dmax", "Dmax must be a nonnegative scalar." );
  const double Dmax2 = Dmax * Dmax;

  /* ---- blob fields ---- */
  auto fld = [&]( const char* n ) -> const mxArray* {
    const mxArray* f = mxGetField( prhs[1], 0, n );
    if( !f ) mexErrMsgIdAndTxt( "bvhClosestElement_mx:B",
                 "B lacks field '%s' (rebuild with BVH).", n );
    return f;
  };
  const mxArray *fB4 = fld("bounds4"), *fC4 = fld("child4"), *fR4 = fld("srange");
  const mxArray *fkV = fld("pkV"), *fkS = fld("pkS"), *fkT = fld("pkT"), *fkE = fld("pkE");
  const mxArray *fP4 = fld("pk4"), *fPI = fld("pk4id"), *fS4 = fld("s4");
  const mxArray *fV  = fld("vol");
  if( !mxIsSingle(fB4) || !mxIsInt32(fC4) || !mxIsInt32(fR4) ||
      !mxIsDouble(fkV) || !mxIsDouble(fkS) || !mxIsInt32(fkT) || !mxIsInt32(fkE) ||
      !mxIsDouble(fP4) || !mxIsInt32(fPI) || !mxIsInt32(fS4) )
    mexErrMsgIdAndTxt( "bvhClosestElement_mx:B", "wrong blob field types (rebuild with BVH)." );
  const int    vol = (int)mxGetScalar( fV );
  const mwSize nN  = mxGetN( fB4 );
  const mwSize nE  = mxGetN( fkV );
  const mwSize S   = ( vol == 1 ) ? 16 : ( vol == 2 ) ? 24 : ( vol == 3 ) ? 60 :
                     ( vol == 4 ) ? 56 : ( vol == 5 ) ? 60 : 28;
  if( ( vol < 1 || vol > 6 ) || nN == 0 || nE == 0 ||
      mxGetM(fB4) != S || mxGetM(fC4) != 4 || mxGetN(fC4) != nN ||
      mxGetM(fR4) != 8 || mxGetN(fR4) != nN ||
      mxGetM(fkV) != 12 || mxGetM(fkS) != 4 || mxGetN(fkS) != nE ||
      mxGetNumberOfElements(fkT) != nE || mxGetNumberOfElements(fkE) != nE )
    mexErrMsgIdAndTxt( "bvhClosestElement_mx:B", "inconsistent blob sizes (rebuild with BVH)." );

  const float*   W4  = (const float*)  mxGetData( fB4 );
  const int32_t* Wc  = (const int32_t*)mxGetData( fC4 );
  const int32_t* Wr  = (const int32_t*)mxGetData( fR4 );
  const double*  vv  = mxGetPr( fkV );
  const Elem*    ee  = (const Elem*)mxGetPr( fkS );
  const int32_t* ety = (const int32_t*)mxGetData( fkT );
  const int32_t* eii = (const int32_t*)mxGetData( fkE );
  const double*  PK4 = mxGetPr( fP4 );
  const int32_t* PKI = (const int32_t*)mxGetData( fPI );
  const int32_t* Ws4 = (const int32_t*)mxGetData( fS4 );
  const mwSize   nB  = mxGetN( fP4 );
  if( mxGetM(fP4) != 36 || mxGetM(fPI) != 4 || mxGetN(fPI) != nB ||
      mxGetM(fS4) != 8 || mxGetN(fS4) != nN )
    mexErrMsgIdAndTxt( "bvhClosestElement_mx:B", "inconsistent PreTri4 pool (rebuild with BVH)." );
  for( mwSize i = 0; i < 4*nB; ++i )
    if( PKI[i] < 0 || PKI[i] > (int32_t)nE )
      mexErrMsgIdAndTxt( "bvhClosestElement_mx:B", "corrupt PreTri4 lane id." );

  for( mwSize i = 0; i < nN; ++i )                       /* bounds-check      */
    for( int k = 0; k < 4; ++k ) {
      const int32_t c = Wc[ i*4 + k ];
      if( c < -1 || c > (int32_t)nN )
        mexErrMsgIdAndTxt( "bvhClosestElement_mx:B", "corrupt child index." );
      if( c != 0 ) {
        const int32_t lo = Wr[ i*8 + 2*k ], hi = Wr[ i*8 + 2*k + 1 ];
        if( lo < 1 || hi < lo || hi > (int32_t)nE )
          mexErrMsgIdAndTxt( "bvhClosestElement_mx:B", "corrupt slot range." );
        const int32_t s0 = Ws4[ i*8 + 2*k ], sn = Ws4[ i*8 + 2*k + 1 ];
        if( sn < 0 || ( sn > 0 && ( s0 < 1 || (mwSize)( s0 + sn - 1 ) > nB ) ) )
          mexErrMsgIdAndTxt( "bvhClosestElement_mx:B", "corrupt PreTri4 slot range." );
      }
    }

  /* ---- fused per-call node pool: bounds + children CONTIGUOUS (one memory
   *      street per visit instead of two far-apart column reads) ---- */
  const size_t stride = (size_t)S*4 + 16;
  std::vector<char> fused( stride * nN );
  for( mwSize i = 0; i < nN; ++i ) {
    memcpy( &fused[ i*stride ], W4 + (size_t)i*S, (size_t)S*4 );
    memcpy( &fused[ i*stride + (size_t)S*4 ], Wc + (size_t)i*4, 16 );
  }
  const char* FZ = fused.data();

  /* ---- outputs ---- */
  plhs[0] = mxCreateDoubleMatrix( nP, 1, mxREAL );
  mxArray* mxCP = mxCreateDoubleMatrix( nP, 3, mxREAL );
  mxArray* mxD  = mxCreateDoubleMatrix( nP, 1, mxREAL );
  double* oE  = mxGetPr( plhs[0] );
  double* oCP = mxGetPr( mxCP );
  double* oD  = mxGetPr( mxD );

  /* ---- Morton-order the queries (cache-coherent walks + useful warm starts) */
  std::vector<int32_t> order( nP );
  for( mwSize i = 0; i < nP; ++i ) order[i] = (int32_t)i;
  if( nP >= 128 ) {
    double mn[3]={INF,INF,INF}, mx[3]={-INF,-INF,-INF};
    for( mwSize i = 0; i < nP; ++i ) {
      const double x=P[i], y=P[i+nP], z=P[i+2*nP];
      if( std::isfinite(x) && std::isfinite(y) && std::isfinite(z) ) {
        if(x<mn[0])mn[0]=x; if(x>mx[0])mx[0]=x;
        if(y<mn[1])mn[1]=y; if(y>mx[1])mx[1]=y;
        if(z<mn[2])mn[2]=z; if(z>mx[2])mx[2]=z;
      }
    }
    const double sx = ( mx[0]>mn[0] ) ? 1023.0/(mx[0]-mn[0]) : 0.0;
    const double sy = ( mx[1]>mn[1] ) ? 1023.0/(mx[1]-mn[1]) : 0.0;
    const double sz = ( mx[2]>mn[2] ) ? 1023.0/(mx[2]-mn[2]) : 0.0;
    std::vector<uint32_t> code( nP );
    for( mwSize i = 0; i < nP; ++i ) {
      const double x=P[i], y=P[i+nP], z=P[i+2*nP];
      if( std::isfinite(x) && std::isfinite(y) && std::isfinite(z) ) {
        const uint32_t cx = (uint32_t)( (x-mn[0])*sx );
        const uint32_t cy = (uint32_t)( (y-mn[1])*sy );
        const uint32_t cz = (uint32_t)( (z-mn[2])*sz );
        code[i] = ( expandBits(cx) << 2 ) | ( expandBits(cy) << 1 ) | expandBits(cz);
      } else code[i] = 0xFFFFFFFFu;
    }
    std::sort( order.begin(), order.end(),
               [&code]( int32_t a, int32_t b ){ return code[a] < code[b]; } );
  }

  const int32_t* ord = order.data();
  const bool     avx = useAVX();
  const int      WS  = (int)S;

  auto runRange = [&]( mwSize i0, mwSize i1 )
  {
    int32_t jwarm = -1;              /* previous winner (packed index), per thread */
    int32_t stkN[192];  double stkD[192];  double stkR2[192];

    for( mwSize qi = i0; qi < i1; ++qi ) {
      const mwSize q = (mwSize)ord[qi];
      const double p[3] = { P[q], P[q+nP], P[q+2*nP] };

      if( !std::isfinite(p[0]) || !std::isfinite(p[1]) || !std::isfinite(p[2]) ) {
        oE[q] = 0.0;  oD[q] = mxGetNaN();
        oCP[q] = oCP[q+nP] = oCP[q+2*nP] = mxGetNaN();
        continue;
      }

      double best2 = Dmax2, best = Dmax;   /* Dmax seeds the bound: farther prunes */
      int32_t bestJ = -1;
      double bcp[3] = { 0.0, 0.0, 0.0 };

      if( jwarm >= 0 ) {                                   /* warm-start seed */
        double c[3];
        const double q2 = d2Elem( vv + (size_t)jwarm*12, ety[jwarm], p, c );
        if( q2 < best2 ) {
          best2 = q2;  best = std::sqrt( q2 );  bestJ = jwarm;
          bcp[0]=c[0]; bcp[1]=c[1]; bcp[2]=c[2];
        }
      }

      int top = 1;
      stkN[0] = 0;  stkD[0] = 0.0;  stkR2[0] = INF;        /* root: never culled */
      while( top ) {
        --top;
        const int32_t ni = stkN[top];
        if( vol == 1 || vol >= 5 ) {                       /* swept family    */
          const double rb = stkR2[top] + best;
          if( stkD[top] > rb*rb ) continue;                /* stale-pop cull */
        } else {
          if( stkD[top] >= best2 ) continue;
        }
        const char*    nz = FZ + (size_t)ni * stride;
        const float*   nb = (const float*)nz;
        const int32_t* nc = (const int32_t*)( nz + (size_t)WS*4 );
        const int32_t* nr = Wr + (size_t)ni * 8;
        const int32_t* n4 = Ws4 + (size_t)ni * 8;

        double key[4], rs[4];  int act[4];  int na = 0;
        for( int k = 0; k < 4; ++k ) {
          if( nc[k] == 0 ) continue;
          if( vol == 1 || vol >= 5 ) {
            /* swept family: prune with  core2 > (r + best)^2 */
            double d2c, r;
            if( vol == 1 ) {
              const double dx = p[0]-(double)nb[k], dy = p[1]-(double)nb[4+k], dz = p[2]-(double)nb[8+k];
              d2c = dx*dx + dy*dy + dz*dz;
              r   = (double)nb[12+k];
            } else if( vol == 6 ) {
              /* LSS/capsule: point-to-segment core (float segment, exact) */
              const double q0x=(double)nb[4*0+k], q0y=(double)nb[4*1+k], q0z=(double)nb[4*2+k];
              const double svx=(double)nb[4*3+k]-q0x, svy=(double)nb[4*4+k]-q0y, svz=(double)nb[4*5+k]-q0z;
              const double svv = svx*svx + svy*svy + svz*svz;
              double t = ( p[0]-q0x )*svx + ( p[1]-q0y )*svy + ( p[2]-q0z )*svz;
              t = ( svv > 0.0 ) ? t/svv : 0.0;
              if( t < 0.0 ) t = 0.0;  else if( t > 1.0 ) t = 1.0;
              const double cx = q0x+t*svx, cy = q0y+t*svy, cz = q0z+t*svz;
              d2c = ( p[0]-cx )*( p[0]-cx ) + ( p[1]-cy )*( p[1]-cy ) + ( p[2]-cz )*( p[2]-cz );
              r   = (double)nb[4*6+k];
            } else {
              /* RSS: point-to-rectangle core in the (near-orthonormal) float
               * axes basis -- deflate like the OBB */
              double du[3];
              for( int a = 0; a < 3; ++a ) {
                const double ax = (double)nb[ 4*(3*a  ) + k ];
                const double ay = (double)nb[ 4*(3*a+1) + k ];
                const double az = (double)nb[ 4*(3*a+2) + k ];
                du[a] = p[0]*ax + p[1]*ay + p[2]*az;
              }
              double eu = (double)nb[4* 9+k] - du[0];  if( du[0]-(double)nb[4*10+k] > eu ) eu = du[0]-(double)nb[4*10+k];  if( eu < 0 ) eu = 0;
              double ev = (double)nb[4*11+k] - du[1];  if( du[1]-(double)nb[4*12+k] > ev ) ev = du[1]-(double)nb[4*12+k];  if( ev < 0 ) ev = 0;
              const double ew = du[2] - (double)nb[4*13+k];
              d2c = ( eu*eu + ev*ev + ew*ew ) * ( 1.0 - 1e-5 );
              r   = (double)nb[4*14+k];
            }
            const double rb = r + best;
            if( d2c > rb*rb ) continue;
            key[na] = d2c;  rs[na] = r;  act[na] = k;  ++na;
          } else {
            double lb2;
            if( vol == 2 || vol == 4 ) {     /* AABB part (kdop shares layout) */
              double dx = (double)nb[   k] - p[0];  if( p[0]-(double)nb[ 4+k] > dx ) dx = p[0]-(double)nb[ 4+k];  if( dx < 0 ) dx = 0;
              double dy = (double)nb[ 8+k] - p[1];  if( p[1]-(double)nb[12+k] > dy ) dy = p[1]-(double)nb[12+k];  if( dy < 0 ) dy = 0;
              double dz = (double)nb[16+k] - p[2];  if( p[2]-(double)nb[20+k] > dz ) dz = p[2]-(double)nb[20+k];  if( dz < 0 ) dz = 0;
              lb2 = dx*dx + dy*dy + dz*dz;
              if( vol == 4 ) {               /* diagonal slabs, |dir|^2 = 3   */
                const double pr[4] = { p[0]+p[1]+p[2], p[0]-p[1]+p[2],
                                       p[0]+p[1]-p[2], p[0]-p[1]-p[2] };
                for( int a = 0; a < 4; ++a ) {
                  const double lo = (double)nb[ 4*(6+2*a) + k ];
                  const double hi = (double)nb[ 4*(7+2*a) + k ];
                  double e = lo - pr[a];  if( pr[a] - hi > e ) e = pr[a] - hi;
                  if( e > 0 ) { const double c = e*e/3.0;  if( c > lb2 ) lb2 = c; }
                }
              }
            } else {                         /* OBB: distance in the axes basis;
                                              * float axes are near-orthonormal,
                                              * deflate to stay conservative   */
              double acc = 0.0;
              for( int a = 0; a < 3; ++a ) {
                const double ax = (double)nb[ 4*(3*a  ) + k ];
                const double ay = (double)nb[ 4*(3*a+1) + k ];
                const double az = (double)nb[ 4*(3*a+2) + k ];
                const double qv = p[0]*ax + p[1]*ay + p[2]*az;
                const double lo = (double)nb[ 4*( 9+a) + k ];
                const double hi = (double)nb[ 4*(12+a) + k ];
                double dxa = lo - qv;  if( qv - hi > dxa ) dxa = qv - hi;
                if( dxa > 0 ) acc += dxa*dxa;
              }
              lb2 = acc * ( 1.0 - 1e-5 );
            }
            if( lb2 >= best2 ) continue;
            key[na] = lb2;  rs[na] = 0;  act[na] = k;  ++na;
          }
        }
        for( int a = 1; a < na; ++a ) {                    /* insertion sort <=4 */
          const double kk = key[a], rk = rs[a];  const int ak = act[a];
          int b = a-1;
          while( b >= 0 && key[b] > kk ) { key[b+1]=key[b]; rs[b+1]=rs[b]; act[b+1]=act[b]; --b; }
          key[b+1]=kk; rs[b+1]=rk; act[b+1]=ak;
        }
        for( int a = 0; a < na; ++a ) {                    /* leaves, near->far */
          const int k = act[a];
          if( nc[k] != -1 ) continue;
          if( vol == 1 || vol >= 5 ) { const double rb = rs[a]+best; if( key[a] > rb*rb ) continue; }
          else if( key[a] >= best2 ) continue;
          const int32_t lo = nr[2*k]-1, hi = nr[2*k+1]-1;
          const int32_t s0 = n4[2*k]-1, sn = n4[2*k+1];

          if( avx && sn > 0 ) {
            /* PreTri4: 4-wide branchless Ericson over aligned SoA blocks;
             * null-padding lanes are filtered by their id */
            for( int32_t b = 0; b < sn; ++b ) {
              const double*  blk = PK4 + (size_t)( s0 + b )*36;
              const int32_t* ids = PKI + (size_t)( s0 + b )*4;
              double d2v[4], cpv[12];
              tri4blk( blk, p, d2v, cpv );
              for( int l = 0; l < 4; ++l )
                if( ids[l] > 0 && d2v[l] < best2 ) {
                  best2 = d2v[l];  best = std::sqrt( best2 );  bestJ = ids[l]-1;
                  bcp[0]=cpv[3*l]; bcp[1]=cpv[3*l+1]; bcp[2]=cpv[3*l+2];
                }
            }
          } else {
            for( int32_t j = lo; j <= hi; ++j ) {
              const Elem& E = ee[j];
              const double dx=p[0]-E.cx, dy=p[1]-E.cy, dz=p[2]-E.cz;
              const double erb = E.r + best;
              if( dx*dx + dy*dy + dz*dz > erb*erb ) continue;
              double c[3];
              const double q2 = d2Elem( vv + (size_t)j*12, ety[j], p, c );
              if( q2 < best2 ) {
                best2 = q2;  best = std::sqrt( q2 );  bestJ = j;
                bcp[0]=c[0]; bcp[1]=c[1]; bcp[2]=c[2];
              }
            }
          }
        }
        for( int a = na-1; a >= 0; --a ) {                 /* internals, far 1st */
          const int k = act[a];
          if( nc[k] <= 0 ) continue;
          if( top > 188 )
            mexErrMsgIdAndTxt( "bvhClosestElement_mx:stack", "traversal stack overflow (corrupt blob?)." );
          stkN[top] = nc[k]-1;  stkD[top] = key[a];  stkR2[top] = rs[a];  ++top;
        }
      }

      jwarm = bestJ;
      if( bestJ >= 0 ) {
        oE[q]      = (double)eii[bestJ];
        oD[q]      = best;
        oCP[q]     = bcp[0];  oCP[q+nP] = bcp[1];  oCP[q+2*nP] = bcp[2];
      } else {                             /* nothing within Dmax */
        oE[q] = 0.0;  oD[q] = INF;
        oCP[q] = oCP[q+nP] = oCP[q+2*nP] = mxGetNaN();
      }
    }
  };

#ifdef _OPENMP
  if( nt > 1 && nP > 256 ) {
    #pragma omp parallel num_threads( nt )
    {
      const int tid = omp_get_thread_num();
      const int nth = omp_get_num_threads();
      const mwSize per = ( nP + nth - 1 ) / nth;
      const mwSize i0  = (mwSize)tid * per;
      const mwSize i1  = ( i0 + per < nP ) ? i0 + per : nP;
      if( i0 < i1 ) runRange( i0, i1 );
    }
  } else runRange( 0, nP );
#else
  runRange( 0, nP );
#endif

  if( nlhs > 1 ) plhs[1] = mxCP; else mxDestroyArray( mxCP );
  if( nlhs > 2 ) plhs[2] = mxD;  else mxDestroyArray( mxD );
}
