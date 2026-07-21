/* bvhKernels.h  --  shared closest-point kernels for the BVH query MEXes.
 *
 *   Single source of truth for the geometric primitives (point/segment/triangle
 *   /tetrahedron distance), the REGION-EXACT barycentric emitters, the 4-wide
 *   AVX Ericson kernel, and the Morton bit-spread -- included by
 *   bvhClosestElement_mx.cpp and fanClosestElement_mx.cpp so a fix or a
 *   precision change lives in ONE place (they used to be copied verbatim).
 *
 *   All symbols are `static`/`static inline` -> internal linkage, one private
 *   copy per translation unit (each MEX is its own TU), no ODR concerns.
 *
 * See also bvhClosestElement_mx.cpp, fanClosestElement_mx.cpp.
 */
#ifndef BVH_KERNELS_H
#define BVH_KERNELS_H

#include <cmath>
#include <cstdint>
#include <limits>
#include <immintrin.h>
#if defined(_MSC_VER)
#include <intrin.h>
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

/* --------------------------------------------- region-exact barycentrics
 * The closest-point search classifies each hit into a REGION (vertex / edge /
 * face / interior). Emitting the barycentric coordinates FROM that region --
 * instead of reverse-engineering them from the (rounded) closest point -- makes
 * edge and vertex hits give EXACT zeros for the off-components and keeps the
 * whole thing well-conditioned on slivers (needle triangles), where the
 * from-cp cross-product forms divide by a vanishing area and lose 5-15 digits.
 * bcTri3 mirrors d2Tri's cascade bit-for-bit so the region matches the cp the
 * traversal already found. */
static void bcTri3( const double* A, const double* B, const double* C,
                    const double* p, double bc[3] )
{
  const double ab0=B[0]-A[0], ab1=B[1]-A[1], ab2=B[2]-A[2];
  const double ac0=C[0]-A[0], ac1=C[1]-A[1], ac2=C[2]-A[2];
  const double ap0=p[0]-A[0], ap1=p[1]-A[1], ap2=p[2]-A[2];
  const double d1 = ab0*ap0 + ab1*ap1 + ab2*ap2;
  const double d2 = ac0*ap0 + ac1*ap1 + ac2*ap2;
  if( d1 <= 0.0 && d2 <= 0.0 ) { bc[0]=1.0; bc[1]=0.0; bc[2]=0.0; return; }   /* A */

  const double bp0=p[0]-B[0], bp1=p[1]-B[1], bp2=p[2]-B[2];
  const double d3 = ab0*bp0 + ab1*bp1 + ab2*bp2;
  const double d4 = ac0*bp0 + ac1*bp1 + ac2*bp2;
  if( d3 >= 0.0 && d4 <= d3 ) { bc[0]=0.0; bc[1]=1.0; bc[2]=0.0; return; }    /* B */

  const double vc = d1*d4 - d3*d2;
  if( vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 ) {                                 /* AB */
    const double den = d1 - d3;
    const double t = ( den != 0.0 ) ? d1/den : 0.0;
    bc[0]=1.0-t; bc[1]=t; bc[2]=0.0; return;
  }

  const double cq0=p[0]-C[0], cq1=p[1]-C[1], cq2=p[2]-C[2];
  const double d5 = ab0*cq0 + ab1*cq1 + ab2*cq2;
  const double d6 = ac0*cq0 + ac1*cq1 + ac2*cq2;
  if( d6 >= 0.0 && d5 <= d6 ) { bc[0]=0.0; bc[1]=0.0; bc[2]=1.0; return; }    /* C */

  const double vb = d5*d2 - d1*d6;
  if( vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 ) {                                 /* AC */
    const double den = d2 - d6;
    const double t = ( den != 0.0 ) ? d2/den : 0.0;
    bc[0]=1.0-t; bc[1]=0.0; bc[2]=t; return;
  }

  const double va = d3*d6 - d5*d4;
  if( va <= 0.0 && (d4-d3) >= 0.0 && (d5-d6) >= 0.0 ) {                       /* BC */
    const double den = (d4-d3) + (d5-d6);
    const double t = ( den != 0.0 ) ? (d4-d3)/den : 0.0;
    bc[0]=0.0; bc[1]=1.0-t; bc[2]=t; return;
  }

  const double den = va + vb + vc;                                           /* interior */
  if( den != 0.0 && std::isfinite( den ) ) {
    const double v = vb/den, w = vc/den;
    bc[0]=1.0-v-w; bc[1]=v; bc[2]=w; return;
  }

  /* fully degenerate: nearest of the three edges, emit its 1-D parametrization */
  double c1[3], c2[3], c3[3];
  const double q1 = d2Seg( A, B, p, c1 );
  const double q2 = d2Seg( A, C, p, c2 );
  const double q3 = d2Seg( B, C, p, c3 );
  if( q1 <= q2 && q1 <= q3 ) {
    const double L=(ab0*ab0+ab1*ab1+ab2*ab2); double t=(L>0)?(d1/L):0.0; if(t<0)t=0; else if(t>1)t=1;
    bc[0]=1.0-t; bc[1]=t; bc[2]=0.0;
  } else if( q2 <= q3 ) {
    const double L=(ac0*ac0+ac1*ac1+ac2*ac2); double t=(L>0)?(d2/L):0.0; if(t<0)t=0; else if(t>1)t=1;
    bc[0]=1.0-t; bc[1]=0.0; bc[2]=t;
  } else {
    const double bcx=C[0]-B[0], bcy=C[1]-B[1], bcz=C[2]-B[2];
    const double L=bcx*bcx+bcy*bcy+bcz*bcz;
    double t=(L>0)?(((p[0]-B[0])*bcx+(p[1]-B[1])*bcy+(p[2]-B[2])*bcz)/L):0.0; if(t<0)t=0; else if(t>1)t=1;
    bc[0]=0.0; bc[1]=1.0-t; bc[2]=t;
  }
}

/* enforce the barycentric invariants as hard as floating point allows, WITHOUT
 * disturbing the region-exact structure: clamp negatives to 0 (rounding at
 * region boundaries / the tet interior tolerance can dip a hair below 0), then
 * renormalize so the weights sum to EXACTLY 1. Exact zeros survive (0 clamps to
 * 0, scales to 0) and the closest point still reconstructs to ~eps. A
 * pathological all-nonpositive vector (never happens for a real closest point)
 * snaps to the first vertex. */
static inline void finalizeBC( double bc[4] )
{
  if( bc[0] < 0.0 ) bc[0] = 0.0;   if( bc[1] < 0.0 ) bc[1] = 0.0;
  if( bc[2] < 0.0 ) bc[2] = 0.0;   if( bc[3] < 0.0 ) bc[3] = 0.0;
  const double s = bc[0] + bc[1] + bc[2] + bc[3];
  if( s > 0.0 ) { const double r = 1.0/s; bc[0]*=r; bc[1]*=r; bc[2]*=r; bc[3]*=r; }
  else          { bc[0] = 1.0; }
}

/* region-exact barycentrics of p's closest point on packed element (type k),
 * padded to 4 components (unused = 0). Caller applies finalizeBC. */
static void bcElem( const double* v, int k, const double* p, double bc[4] )
{
  bc[0]=bc[1]=bc[2]=bc[3]=0.0;
  switch( k ) {
    case 1: bc[0]=1.0; return;
    case 2: {
      const double* A=v; const double* B=v+3;
      const double abx=B[0]-A[0], aby=B[1]-A[1], abz=B[2]-A[2];
      const double L2=abx*abx+aby*aby+abz*abz;
      double t = (p[0]-A[0])*abx + (p[1]-A[1])*aby + (p[2]-A[2])*abz;
      t = ( L2 > 0.0 ) ? t/L2 : 0.0;
      if( t<0.0 ) t=0.0; else if( t>1.0 ) t=1.0;
      bc[0]=1.0-t; bc[1]=t; return;
    }
    case 3: { double b3[3]; bcTri3( v, v+3, v+6, p, b3 );
              bc[0]=b3[0]; bc[1]=b3[1]; bc[2]=b3[2]; return; }
    case 4: {
      const double* A=v; const double* B=v+3; const double* C=v+6; const double* D=v+9;
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
          bc[0]=l1; bc[1]=l2; bc[2]=l3; bc[3]=l4; return;
        }
      }
      /* outside/degenerate: nearest face, map its triangle bc onto the 4 verts */
      double c[3], cb[3], b3[3];  int face = 0;
      double q  = d2Tri( A, B, C, p, c  );
      double q2 = d2Tri( A, B, D, p, cb ); if( q2 < q ) { q=q2; face=1; }
      q2        = d2Tri( A, C, D, p, cb ); if( q2 < q ) { q=q2; face=2; }
      q2        = d2Tri( B, C, D, p, cb ); if( q2 < q ) { q=q2; face=3; }
      switch( face ) {
        case 0: bcTri3(A,B,C,p,b3); bc[0]=b3[0]; bc[1]=b3[1]; bc[2]=b3[2]; break;
        case 1: bcTri3(A,B,D,p,b3); bc[0]=b3[0]; bc[1]=b3[1]; bc[3]=b3[2]; break;
        case 2: bcTri3(A,C,D,p,b3); bc[0]=b3[0]; bc[2]=b3[1]; bc[3]=b3[2]; break;
        case 3: bcTri3(B,C,D,p,b3); bc[1]=b3[0]; bc[2]=b3[1]; bc[3]=b3[2]; break;
      }
      return;
    }
  }
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

#endif /* BVH_KERNELS_H */
