/* approximateClosestElement_mx  --  compiled core of approximateClosestElement.
 *
 *   [e, cp, d] = approximateClosestElement_mx( P , B , nthreads , Dmax )
 *
 *   1-RING-OF-NEAREST-VERTEX approximate locator, two fused stages per point:
 *     1) nearest VERTEX via a BVH over the mesh vertices (point elements,
 *        AABB nodes; pkS centers with r = 0 are the vertices themselves)
 *     2) EXACT distance to the vertex's incident-element fan (EsuP, packed
 *        as CSR int32 fanStart/fanEl) -- elements tested with the same
 *        primitives as bvhClosestElement_mx (point/segment/triangle/TET;
 *        4 nonzero nodes ALWAYS a tetrahedron; mixed 0-padded fine).
 *
 *   The result is an UPPER BOUND of the true distance (the fan is a real
 *   subset of the mesh); the winning element is typically exact for 95-100%
 *   of the queries (see bench_approximate).
 *
 *   B is the blob built by approximateClosestElement( M ): a vertex AABB
 *   blob (vol == 2) over the USED vertices, extended with
 *     fanStart  int32 (nV+1)   CSR offsets, aligned to the blob's point rows
 *     fanEl     int32 (nnz)    1-based mesh element ids
 *     elV       12 x nEl double  packed element vertices (BUILD space)
 *     elT       int32 nEl        nonzero-node count per element (0..4)
 *
 *   Dmax scalar | nP-vector: bounds the VERTEX distance (stage 1); a point
 *   whose nearest vertex is beyond Dmax gives e = 0, d = Inf, cp = NaN
 *   (non-finite query points: d = NaN). Same optimizations as the exact mex:
 *   Morton-ordered queries, warm start, sqrt-free pruning, stale-pop culls,
 *   fused node pool, fixed stacks, fully bounds-checked blob.
 *
 *   Compile (MSVC):  mex COMPFLAGS="$COMPFLAGS /openmp" approximateClosestElement_mx.cpp
 *
 * See also approximateClosestElement, bvhClosestElement_mx, BVH.
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

/* ------------------------- exact primitives (same as bvhClosestElement_mx) */

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

static double d2Tri( const double* A, const double* B, const double* C,
                     const double* p, double* cp )
{
  const double ab0=B[0]-A[0], ab1=B[1]-A[1], ab2=B[2]-A[2];
  const double ac0=C[0]-A[0], ac1=C[1]-A[1], ac2=C[2]-A[2];
  const double ap0=p[0]-A[0], ap1=p[1]-A[1], ap2=p[2]-A[2];

  const double d1 = ab0*ap0 + ab1*ap1 + ab2*ap2;
  const double d2 = ac0*ap0 + ac1*ap1 + ac2*ap2;
  if( d1 <= 0.0 && d2 <= 0.0 ) return d2Point( A, p, cp );

  const double bp0=p[0]-B[0], bp1=p[1]-B[1], bp2=p[2]-B[2];
  const double d3 = ab0*bp0 + ab1*bp1 + ab2*bp2;
  const double d4 = ac0*bp0 + ac1*bp1 + ac2*bp2;
  if( d3 >= 0.0 && d4 <= d3 ) return d2Point( B, p, cp );

  const double vc = d1*d4 - d3*d2;
  if( vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 ) {
    const double den = d1 - d3;
    const double t = ( den != 0.0 ) ? d1/den : 0.0;
    cp[0]=A[0]+t*ab0; cp[1]=A[1]+t*ab1; cp[2]=A[2]+t*ab2;
    return sq(p[0]-cp[0]) + sq(p[1]-cp[1]) + sq(p[2]-cp[2]);
  }

  const double cq0=p[0]-C[0], cq1=p[1]-C[1], cq2=p[2]-C[2];
  const double d5 = ab0*cq0 + ab1*cq1 + ab2*cq2;
  const double d6 = ac0*cq0 + ac1*cq1 + ac2*cq2;
  if( d6 >= 0.0 && d5 <= d6 ) return d2Point( C, p, cp );

  const double vb = d5*d2 - d1*d6;
  if( vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 ) {
    const double den = d2 - d6;
    const double t = ( den != 0.0 ) ? d2/den : 0.0;
    cp[0]=A[0]+t*ac0; cp[1]=A[1]+t*ac1; cp[2]=A[2]+t*ac2;
    return sq(p[0]-cp[0]) + sq(p[1]-cp[1]) + sq(p[2]-cp[2]);
  }

  const double va = d3*d6 - d5*d4;
  if( va <= 0.0 && (d4-d3) >= 0.0 && (d5-d6) >= 0.0 ) {
    const double den = (d4-d3) + (d5-d6);
    const double t = ( den != 0.0 ) ? (d4-d3)/den : 0.0;
    cp[0]=B[0]+t*(C[0]-B[0]); cp[1]=B[1]+t*(C[1]-B[1]); cp[2]=B[2]+t*(C[2]-B[2]);
    return sq(p[0]-cp[0]) + sq(p[1]-cp[1]) + sq(p[2]-cp[2]);
  }

  const double den = va + vb + vc;
  if( den != 0.0 && std::isfinite( den ) ) {
    const double v = vb/den, w = vc/den;
    cp[0]=A[0]+v*ab0+w*ac0; cp[1]=A[1]+v*ab1+w*ac1; cp[2]=A[2]+v*ab2+w*ac2;
    return sq(p[0]-cp[0]) + sq(p[1]-cp[1]) + sq(p[2]-cp[2]);
  }

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
    if( l1 >= tol && l2 >= tol && l3 >= tol && l4 >= tol ) {
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

static inline double d2Elem( const double* v, int k, const double* p, double* cp )
{
  switch( k ) {
    case 1: return d2Point( v, p, cp );
    case 2: return d2Seg( v, v+3, p, cp );
    case 3: return d2Tri( v, v+3, v+6, p, cp );
    case 4: return d2Tet( v, v+3, v+6, v+9, p, cp );
  }
  return INF;
}

/* ------------------------------------------------------- 4-wide kernels */

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

/* 4 point-distances at once over a pt4 SoA block [x0..x3 y0..y3 z0..z3] */
static inline void pt4blk( const double* blk, const double* p, double d2o[4] )
{
  const __m256d dx = _mm256_sub_pd( _mm256_set1_pd(p[0]), _mm256_loadu_pd( blk     ) );
  const __m256d dy = _mm256_sub_pd( _mm256_set1_pd(p[1]), _mm256_loadu_pd( blk + 4 ) );
  const __m256d dz = _mm256_sub_pd( _mm256_set1_pd(p[2]), _mm256_loadu_pd( blk + 8 ) );
  _mm256_storeu_pd( d2o, mm_dot3( dx,dy,dz, dx,dy,dz ) );
}

/* branchless Ericson over one PreTri4 block (A, AB, AC in 4-wide SoA) --
 * copied VERBATIM from bvhClosestElement_mx (same region order = same bits;
 * degenerate lanes redone scalar; padding lanes filtered by id) */
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
    mexErrMsgIdAndTxt( "approximateClosestElement_mx:nrhs",
                       "expected P, B [, nthreads, Dmax]." );
  if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) || mxIsSparse(prhs[0]) ||
      mxGetN(prhs[0]) != 3 )
    mexErrMsgIdAndTxt( "approximateClosestElement_mx:P", "P must be nP x 3 double (pad it first)." );
  if( !mxIsStruct( prhs[1] ) )
    mexErrMsgIdAndTxt( "approximateClosestElement_mx:B", "B must be a struct (approximate blob)." );

  const mwSize  nP = mxGetM( prhs[0] );
  const double* P  = mxGetPr( prhs[0] );

  int nt = ( nrhs > 2 ) ? (int)mxGetScalar( prhs[2] ) : 1;
  if( nt < 1 ) nt = 1;  if( nt > 64 ) nt = 64;

  const double* DmaxV = NULL;
  double Dmax = INF;
  if( nrhs > 3 ) {
    if( !mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) || mxIsSparse(prhs[3]) )
      mexErrMsgIdAndTxt( "approximateClosestElement_mx:Dmax",
                         "Dmax must be a double scalar or an nP-vector." );
    const mwSize nD = mxGetNumberOfElements( prhs[3] );
    if( nD == 1 ) {
      Dmax = mxGetScalar( prhs[3] );
      if( !( Dmax >= 0.0 ) )
        mexErrMsgIdAndTxt( "approximateClosestElement_mx:Dmax", "Dmax must be nonnegative." );
    } else if( nD == nP ) {
      DmaxV = mxGetPr( prhs[3] );
      for( mwSize i = 0; i < nP; ++i )
        if( !( DmaxV[i] >= 0.0 ) )
          mexErrMsgIdAndTxt( "approximateClosestElement_mx:Dmax",
                             "per-point Dmax must be nonnegative (no NaN)." );
    } else
      mexErrMsgIdAndTxt( "approximateClosestElement_mx:Dmax",
                         "Dmax must be a scalar or an nP-vector." );
  }
  const double Dmax2 = Dmax * Dmax;

  /* ---- blob fields ---- */
  auto fld = [&]( const char* n ) -> const mxArray* {
    const mxArray* f = mxGetField( prhs[1], 0, n );
    if( !f ) mexErrMsgIdAndTxt( "approximateClosestElement_mx:B",
                 "B lacks field '%s' (rebuild with approximateClosestElement(M)).", n );
    return f;
  };
  const mxArray *fB4 = fld("bounds4"), *fC4 = fld("child4"), *fR4 = fld("srange");
  const mxArray *fkS = fld("pkS"), *fkE = fld("pkE"), *fV = fld("vol");
  const mxArray *fFS = fld("fanStart"), *fFE = fld("fanEl");
  const mxArray *fEV = fld("elV"), *fET = fld("elT");
  const mxArray *fP4 = fld("pt4");
  /* fan4 trio: OPTIONAL (only pure-triangle meshes pack it) */
  const mxArray *fF4 = mxGetField( prhs[1], 0, "fan4" );
  const mxArray *fFI = mxGetField( prhs[1], 0, "fan4id" );
  const mxArray *fF0 = mxGetField( prhs[1], 0, "fan4Start" );
  if( !mxIsSingle(fB4) || !mxIsInt32(fC4) || !mxIsInt32(fR4) ||
      !mxIsDouble(fkS) || !mxIsInt32(fkE) || !mxIsDouble(fP4) ||
      !mxIsInt32(fFS) || !mxIsInt32(fFE) || !mxIsDouble(fEV) || !mxIsInt32(fET) )
    mexErrMsgIdAndTxt( "approximateClosestElement_mx:B", "wrong blob field types." );
  const bool hasF4 = ( fF4 && fFI && fF0 );
  if( hasF4 && ( !mxIsDouble(fF4) || !mxIsInt32(fFI) || !mxIsInt32(fF0) ) )
    mexErrMsgIdAndTxt( "approximateClosestElement_mx:B", "wrong fan4 field types." );
  if( (int)mxGetScalar( fV ) != 2 )
    mexErrMsgIdAndTxt( "approximateClosestElement_mx:B",
                       "the approximate blob must be an AABB vertex blob (vol == 2)." );

  const mwSize nN  = mxGetN( fB4 );
  const mwSize nV  = mxGetN( fkS );                      /* point elements    */
  const mwSize nEl = mxGetN( fEV );                      /* mesh elements     */
  const mwSize nFn = mxGetNumberOfElements( fFE );
  const int    S   = 24;
  if( nN == 0 || nV == 0 || nEl == 0 ||
      mxGetM(fB4) != (mwSize)S || mxGetM(fC4) != 4 || mxGetN(fC4) != nN ||
      mxGetM(fR4) != 8 || mxGetN(fR4) != nN ||
      mxGetM(fkS) != 4 || mxGetNumberOfElements(fkE) != nV ||
      mxGetNumberOfElements(fFS) != nV + 1 || mxGetM(fEV) != 12 ||
      mxGetNumberOfElements(fET) != nEl ||
      mxGetM(fP4) != 12 || mxGetN(fP4) != ( nV + 3 ) / 4 )
    mexErrMsgIdAndTxt( "approximateClosestElement_mx:B", "inconsistent blob sizes." );
  const mwSize nFB = hasF4 ? mxGetN( fF4 ) : 0;
  if( hasF4 && ( mxGetM(fF4) != 36 || mxGetM(fFI) != 4 || mxGetN(fFI) != nFB ||
                 mxGetNumberOfElements(fF0) != nV + 1 ) )
    mexErrMsgIdAndTxt( "approximateClosestElement_mx:B", "inconsistent fan4 sizes." );

  const float*   W4  = (const float*)  mxGetData( fB4 );
  const int32_t* Wc  = (const int32_t*)mxGetData( fC4 );
  const int32_t* Wr  = (const int32_t*)mxGetData( fR4 );
  const Elem*    ee  = (const Elem*)mxGetPr( fkS );
  const int32_t* eii = (const int32_t*)mxGetData( fkE );
  const int32_t* fanS = (const int32_t*)mxGetData( fFS );
  const int32_t* fanE = (const int32_t*)mxGetData( fFE );
  const double*  elV = mxGetPr( fEV );
  const int32_t* elT = (const int32_t*)mxGetData( fET );
  const double*  PT4 = mxGetPr( fP4 );
  const double*  FN4 = hasF4 ? mxGetPr( fF4 ) : NULL;
  const int32_t* FNI = hasF4 ? (const int32_t*)mxGetData( fFI ) : NULL;
  const int32_t* FN0 = hasF4 ? (const int32_t*)mxGetData( fF0 ) : NULL;

  /* full bounds-check: a corrupt blob errors out, it cannot access-violate */
  for( mwSize i = 0; i < nN; ++i )
    for( int k = 0; k < 4; ++k ) {
      const int32_t c = Wc[ i*4 + k ];
      if( c < -1 || c > (int32_t)nN )
        mexErrMsgIdAndTxt( "approximateClosestElement_mx:B", "corrupt child index." );
      if( c != 0 ) {
        const int32_t lo = Wr[ i*8 + 2*k ], hi = Wr[ i*8 + 2*k + 1 ];
        if( lo < 1 || hi < lo || hi > (int32_t)nV )
          mexErrMsgIdAndTxt( "approximateClosestElement_mx:B", "corrupt slot range." );
      }
    }
  for( mwSize i = 0; i < nV; ++i ) {
    if( eii[i] < 1 || (mwSize)eii[i] > nV )
      mexErrMsgIdAndTxt( "approximateClosestElement_mx:B", "corrupt packed vertex id." );
    if( fanS[i] < 0 || fanS[i] > fanS[i+1] )
      mexErrMsgIdAndTxt( "approximateClosestElement_mx:B", "corrupt fan offsets." );
  }
  if( (mwSize)fanS[nV] != nFn )
    mexErrMsgIdAndTxt( "approximateClosestElement_mx:B", "fan offsets do not match fanEl." );
  for( mwSize i = 0; i < nFn; ++i )
    if( fanE[i] < 1 || (mwSize)fanE[i] > nEl )
      mexErrMsgIdAndTxt( "approximateClosestElement_mx:B", "corrupt fan element id." );
  for( mwSize i = 0; i < nEl; ++i )
    if( elT[i] < 0 || elT[i] > 4 )
      mexErrMsgIdAndTxt( "approximateClosestElement_mx:B", "corrupt element type." );
  if( hasF4 ) {
    for( mwSize i = 0; i < nV; ++i )
      if( FN0[i] < 0 || FN0[i] > FN0[i+1] )
        mexErrMsgIdAndTxt( "approximateClosestElement_mx:B", "corrupt fan4 offsets." );
    if( (mwSize)FN0[nV] != nFB )
      mexErrMsgIdAndTxt( "approximateClosestElement_mx:B", "fan4 offsets do not match blocks." );
    for( mwSize i = 0; i < 4*nFB; ++i )
      if( FNI[i] < 0 || (mwSize)FNI[i] > nEl )
        mexErrMsgIdAndTxt( "approximateClosestElement_mx:B", "corrupt fan4 lane id." );
  }

  /* ---- fused node pool: bounds + children contiguous ---- */
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

  /* ---- Morton order (cache-coherent walks + useful warm starts) ---- */
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

  auto runRange = [&]( mwSize i0, mwSize i1 )
  {
    int32_t jwarm = -1;              /* previous winning vertex (packed index) */
    int32_t stkN[192];  double stkD[192];

    for( mwSize qi = i0; qi < i1; ++qi ) {
      const mwSize q = (mwSize)ord[qi];
      const double p[3] = { P[q], P[q+nP], P[q+2*nP] };

      if( !std::isfinite(p[0]) || !std::isfinite(p[1]) || !std::isfinite(p[2]) ) {
        oE[q] = 0.0;  oD[q] = mxGetNaN();
        oCP[q] = oCP[q+nP] = oCP[q+2*nP] = mxGetNaN();
        continue;
      }

      /* ---- stage 1: nearest vertex (AABB traversal, point leaves) ---- */
      double best2;
      if( DmaxV ) { best2 = DmaxV[q]*DmaxV[q]; }
      else        { best2 = Dmax2; }
      int32_t bestJ = -1;

      if( jwarm >= 0 ) {                                   /* warm-start seed */
        const Elem& E = ee[jwarm];
        const double dx=p[0]-E.cx, dy=p[1]-E.cy, dz=p[2]-E.cz;
        const double q2 = dx*dx + dy*dy + dz*dz;
        if( q2 < best2 ) { best2 = q2;  bestJ = jwarm; }
      }

      int top = 1;
      stkN[0] = 0;  stkD[0] = 0.0;
      while( top ) {
        --top;
        if( stkD[top] >= best2 ) continue;                 /* stale-pop cull */
        const int32_t ni = stkN[top];
        const char*    nz = FZ + (size_t)ni * stride;
        const float*   nb = (const float*)nz;
        const int32_t* nc = (const int32_t*)( nz + (size_t)S*4 );
        const int32_t* nr = Wr + (size_t)ni * 8;

        double key[4];  int act[4];  int na = 0;
        for( int k = 0; k < 4; ++k ) {
          if( nc[k] == 0 ) continue;
          double dx = (double)nb[   k] - p[0];  if( p[0]-(double)nb[ 4+k] > dx ) dx = p[0]-(double)nb[ 4+k];  if( dx < 0 ) dx = 0;
          double dy = (double)nb[ 8+k] - p[1];  if( p[1]-(double)nb[12+k] > dy ) dy = p[1]-(double)nb[12+k];  if( dy < 0 ) dy = 0;
          double dz = (double)nb[16+k] - p[2];  if( p[2]-(double)nb[20+k] > dz ) dz = p[2]-(double)nb[20+k];  if( dz < 0 ) dz = 0;
          const double lb2 = dx*dx + dy*dy + dz*dz;
          if( lb2 >= best2 ) continue;
          key[na] = lb2;  act[na] = k;  ++na;
        }
        for( int a = 1; a < na; ++a ) {                    /* insertion sort <=4 */
          const double kk = key[a];  const int ak = act[a];
          int b = a-1;
          while( b >= 0 && key[b] > kk ) { key[b+1]=key[b]; act[b+1]=act[b]; --b; }
          key[b+1]=kk; act[b+1]=ak;
        }
        for( int a = 0; a < na; ++a ) {                    /* leaves, near->far */
          const int k = act[a];
          if( nc[k] != -1 ) continue;
          if( key[a] >= best2 ) continue;
          const int32_t lo = nr[2*k]-1, hi = nr[2*k+1]-1;
          if( avx ) {
            /* Pt4: 4 distancias por bloque SoA; lanes fuera de [lo,hi]
             * (bordes y padding) filtradas por INDICE, sin ids */
            const int32_t b0 = lo >> 2, b1 = hi >> 2;
            for( int32_t b = b0; b <= b1; ++b ) {
              double d2v[4];
              pt4blk( PT4 + (size_t)b*12, p, d2v );
              const int32_t base = b*4;
              const int32_t l0 = ( base < lo ) ? lo - base : 0;
              const int32_t l1 = ( base + 3 > hi ) ? hi - base : 3;
              for( int32_t l = l0; l <= l1; ++l )
                if( d2v[l] < best2 ) { best2 = d2v[l];  bestJ = base + l; }
            }
          } else {
            for( int32_t j = lo; j <= hi; ++j ) {
              const Elem& E = ee[j];
              const double dx=p[0]-E.cx, dy=p[1]-E.cy, dz=p[2]-E.cz;
              const double q2 = dx*dx + dy*dy + dz*dz;
              if( q2 < best2 ) { best2 = q2;  bestJ = j; }
            }
          }
        }
        for( int a = na-1; a >= 0; --a ) {                 /* internals, far 1st */
          const int k = act[a];
          if( nc[k] <= 0 ) continue;
          if( top > 188 )
            mexErrMsgIdAndTxt( "approximateClosestElement_mx:stack", "traversal stack overflow (corrupt blob?)." );
          stkN[top] = nc[k]-1;  stkD[top] = key[a];  ++top;
        }
      }

      jwarm = bestJ;
      if( bestJ < 0 ) {                                    /* nearest vertex beyond Dmax */
        oE[q] = 0.0;  oD[q] = INF;
        oCP[q] = oCP[q+nP] = oCP[q+2*nP] = mxGetNaN();
        continue;
      }

      /* ---- stage 2: exact sweep of the vertex's fan ---- */
      const int32_t v  = eii[bestJ];                       /* 1-based blob row */
      double fb2 = INF, fcp[3] = { 0.0, 0.0, 0.0 };
      int32_t fe = 0;
      if( hasF4 && avx ) {
        /* fan4: el abanico como bloques PreTri4, kernel 4-wide del motor
         * exacto; lanes de relleno (id 0) filtradas */
        const int32_t b0 = FN0[v-1], b1 = FN0[v];
        for( int32_t b = b0; b < b1; ++b ) {
          const double*  blk = FN4 + (size_t)b*36;
          const int32_t* ids = FNI + (size_t)b*4;
          double d2v[4], cpv[12];
          tri4blk( blk, p, d2v, cpv );
          for( int l = 0; l < 4; ++l )
            if( ids[l] > 0 && d2v[l] < fb2 ) {
              fb2 = d2v[l];  fe = ids[l];
              fcp[0]=cpv[3*l]; fcp[1]=cpv[3*l+1]; fcp[2]=cpv[3*l+2];
            }
        }
      } else {
        const int32_t f0 = fanS[v-1], f1 = fanS[v];
        for( int32_t t = f0; t < f1; ++t ) {
          const int32_t e = fanE[t];
          double c[3];
          const double q2 = d2Elem( elV + (size_t)(e-1)*12, elT[e-1], p, c );
          if( q2 < fb2 ) {
            fb2 = q2;  fe = e;
            fcp[0]=c[0]; fcp[1]=c[1]; fcp[2]=c[2];
          }
        }
      }
      if( fe > 0 ) {
        oE[q]  = (double)fe;
        oD[q]  = std::sqrt( fb2 );
        oCP[q] = fcp[0];  oCP[q+nP] = fcp[1];  oCP[q+2*nP] = fcp[2];
      } else {                                             /* empty fan (should not happen) */
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
