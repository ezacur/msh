/* bvhClosestElement_mx  --  compiled core of bvhClosestElement.
 *
 *   [e, cp, d, bc] = bvhClosestElement_mx( P , B , nthreads , Dmax )
 *
 *     bc   nP x 4  (only if requested) REGION-EXACT barycentrics of the closest
 *                  point (vertex/edge/face/interior), padded to 4 columns.
 *                  Computed from the query point in the region the search chose,
 *                  NOT reverse-engineered from cp -> edge/vertex hits are exact
 *                  zeros and slivers stay well-conditioned. Misses -> NaN.
 *
 *     P    nP x 3  double  query points (padded to 3 columns)
 *     B    struct          BVH blob (version >= 3): wide 4-ary nodes
 *                          (bounds4/child4/srange; sphere or AABB slots in
 *                          conservative float, pruning arithmetic in double)
 *                          + packed leaf data (pkV/pkS/pkT/pkE). The blob is
 *                          SELF-CONTAINED: the mesh itself is not needed here.
 *     nthreads scalar      OpenMP threads over the query points (default 1)
 *     Dmax scalar|nP-vec   search radius (default Inf): the best-so-far bound
 *                          is SEEDED with Dmax, so everything farther prunes
 *                          from the very root -- a point beyond Dmax costs one
 *                          node visit. An nP-VECTOR seeds a PER-POINT upper
 *                          bound (heuristics: nearest-vertex distance etc.).
 *                          Elements at d < Dmax are returned; a
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
 *   Compile (MSVC):  mex COMPFLAGS="$COMPFLAGS /openmp" -lut bvhClosestElement_mx.cpp
 *                    ( -lut links libut for Ctrl-C support; see mexInterrupt.h )
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

#include "bvhKernels.h"
#include "mexInterrupt.h"

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
  /* Dmax: scalar (uniform search radius) or nP-vector (PER-POINT bound seed:
   * an upper bound on each point's distance -- e.g. from a nearest-vertex
   * heuristic -- prunes the traversal exactly like Dmax does) */
  const double* DmaxV = NULL;
  double Dmax = INF;
  if( nrhs > 3 ) {
    if( !mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) || mxIsSparse(prhs[3]) )
      mexErrMsgIdAndTxt( "bvhClosestElement_mx:Dmax",
                         "Dmax must be a double scalar or an nP-vector." );
    const mwSize nD = mxGetNumberOfElements( prhs[3] );
    if( nD == 1 ) {
      Dmax = mxGetScalar( prhs[3] );
      if( !( Dmax >= 0.0 ) )   /* also rejects NaN */
        mexErrMsgIdAndTxt( "bvhClosestElement_mx:Dmax", "Dmax must be nonnegative." );
    } else if( nD == nP ) {
      DmaxV = mxGetPr( prhs[3] );
      for( mwSize i = 0; i < nP; ++i )
        if( !( DmaxV[i] >= 0.0 ) )   /* also rejects NaN */
          mexErrMsgIdAndTxt( "bvhClosestElement_mx:Dmax",
                             "per-point Dmax must be nonnegative (no NaN)." );
    } else
      mexErrMsgIdAndTxt( "bvhClosestElement_mx:Dmax",
                         "Dmax must be a scalar or an nP-vector." );
  }
  const double Dmax2 = Dmax * Dmax;
  /* 5th arg (benchmarking only): nonzero disables the WARM START, to isolate
   * how much of the pruning comes from it vs from external bound seeding */
  const bool noWarm = ( nrhs > 4 ) && ( mxGetScalar( prhs[4] ) != 0.0 );

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
  const bool wantBC = ( nlhs > 3 );      /* region-exact barycentrics on demand */
  plhs[0] = mxCreateDoubleMatrix( nP, 1, mxREAL );
  mxArray* mxCP = mxCreateDoubleMatrix( nP, 3, mxREAL );
  mxArray* mxD  = mxCreateDoubleMatrix( nP, 1, mxREAL );
  mxArray* mxBC = mxCreateDoubleMatrix( nP, wantBC ? 4 : 0, mxREAL );
  double* oE  = mxGetPr( plhs[0] );
  double* oCP = mxGetPr( mxCP );
  double* oD  = mxGetPr( mxD );
  double* oBC = mxGetPr( mxBC );

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
      if( mexInterrupted( (long long)( qi - i0 ) ) ) break;   /* Ctrl-C: bail, let MATLAB abort */
      const mwSize q = (mwSize)ord[qi];
      const double p[3] = { P[q], P[q+nP], P[q+2*nP] };

      if( !std::isfinite(p[0]) || !std::isfinite(p[1]) || !std::isfinite(p[2]) ) {
        oE[q] = 0.0;  oD[q] = mxGetNaN();
        oCP[q] = oCP[q+nP] = oCP[q+2*nP] = mxGetNaN();
        if( wantBC ) oBC[q]=oBC[q+nP]=oBC[q+2*nP]=oBC[q+3*nP]=mxGetNaN();
        continue;
      }

      double best, best2;                  /* (per-point) Dmax seeds the bound */
      if( DmaxV ) { best = DmaxV[q];  best2 = best * best; }
      else        { best = Dmax;      best2 = Dmax2;       }
      int32_t bestJ = -1;
      double bcp[3] = { 0.0, 0.0, 0.0 };

      if( !noWarm && jwarm >= 0 ) {                        /* warm-start seed */
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

        /* 4-wide node test: las 4 cajas/esferas del nodo de una vez (los
         * bounds YA estan en SoA de 4); slots vacios filtrados por nc[k] */
        double lb2v[4];  bool have4 = false;
        if( avx && ( vol == 1 || vol == 2 || vol == 4 ) ) {
          const __m256d px4 = _mm256_set1_pd( p[0] );
          const __m256d py4 = _mm256_set1_pd( p[1] );
          const __m256d pz4 = _mm256_set1_pd( p[2] );
          if( vol == 1 ) {                     /* esferas: d2 a los 4 centros */
            const __m256d dx = _mm256_sub_pd( px4, _mm256_cvtps_pd( _mm_loadu_ps( nb     ) ) );
            const __m256d dy = _mm256_sub_pd( py4, _mm256_cvtps_pd( _mm_loadu_ps( nb + 4 ) ) );
            const __m256d dz = _mm256_sub_pd( pz4, _mm256_cvtps_pd( _mm_loadu_ps( nb + 8 ) ) );
            _mm256_storeu_pd( lb2v, _mm256_add_pd( _mm256_add_pd(
                _mm256_mul_pd(dx,dx), _mm256_mul_pd(dy,dy) ), _mm256_mul_pd(dz,dz) ) );
          } else {                             /* aabb (y la parte aabb del kdop) */
            const __m256d zero4 = _mm256_setzero_pd();
            const __m256d xlo = _mm256_cvtps_pd( _mm_loadu_ps( nb      ) );
            const __m256d xhi = _mm256_cvtps_pd( _mm_loadu_ps( nb +  4 ) );
            const __m256d ylo = _mm256_cvtps_pd( _mm_loadu_ps( nb +  8 ) );
            const __m256d yhi = _mm256_cvtps_pd( _mm_loadu_ps( nb + 12 ) );
            const __m256d zlo = _mm256_cvtps_pd( _mm_loadu_ps( nb + 16 ) );
            const __m256d zhi = _mm256_cvtps_pd( _mm_loadu_ps( nb + 20 ) );
            const __m256d dx = _mm256_max_pd( _mm256_max_pd( _mm256_sub_pd(xlo,px4), _mm256_sub_pd(px4,xhi) ), zero4 );
            const __m256d dy = _mm256_max_pd( _mm256_max_pd( _mm256_sub_pd(ylo,py4), _mm256_sub_pd(py4,yhi) ), zero4 );
            const __m256d dz = _mm256_max_pd( _mm256_max_pd( _mm256_sub_pd(zlo,pz4), _mm256_sub_pd(pz4,zhi) ), zero4 );
            _mm256_storeu_pd( lb2v, _mm256_add_pd( _mm256_add_pd(
                _mm256_mul_pd(dx,dx), _mm256_mul_pd(dy,dy) ), _mm256_mul_pd(dz,dz) ) );
          }
          have4 = true;
        }

        double key[4], rs[4];  int act[4];  int na = 0;
        for( int k = 0; k < 4; ++k ) {
          if( nc[k] == 0 ) continue;
          if( vol == 1 || vol >= 5 ) {
            /* swept family: prune with  core2 > (r + best)^2 */
            double d2c, r;
            if( vol == 1 ) {
              if( have4 ) d2c = lb2v[k];
              else {
                const double dx = p[0]-(double)nb[k], dy = p[1]-(double)nb[4+k], dz = p[2]-(double)nb[8+k];
                d2c = dx*dx + dy*dy + dz*dz;
              }
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
              if( have4 ) lb2 = lb2v[k];
              else {
                double dx = (double)nb[   k] - p[0];  if( p[0]-(double)nb[ 4+k] > dx ) dx = p[0]-(double)nb[ 4+k];  if( dx < 0 ) dx = 0;
                double dy = (double)nb[ 8+k] - p[1];  if( p[1]-(double)nb[12+k] > dy ) dy = p[1]-(double)nb[12+k];  if( dy < 0 ) dy = 0;
                double dz = (double)nb[16+k] - p[2];  if( p[2]-(double)nb[20+k] > dz ) dz = p[2]-(double)nb[20+k];  if( dz < 0 ) dz = 0;
                lb2 = dx*dx + dy*dy + dz*dz;
              }
              if( vol == 4 && lb2 < best2 ) { /* diagonal slabs, |dir|^2 = 3  */
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
          _mm_prefetch( FZ + (size_t)( nc[k]-1 )*stride, _MM_HINT_T0 );
          stkN[top] = nc[k]-1;  stkD[top] = key[a];  stkR2[top] = rs[a];  ++top;
        }
      }

      jwarm = bestJ;
      if( bestJ >= 0 ) {
        oE[q]      = (double)eii[bestJ];
        oD[q]      = best;
        oCP[q]     = bcp[0];  oCP[q+nP] = bcp[1];  oCP[q+2*nP] = bcp[2];
        if( wantBC ) {                     /* region-exact bc for the winner */
          double b4[4];
          bcElem( vv + (size_t)bestJ*12, ety[bestJ], p, b4 );
          finalizeBC( b4 );                /* clamp >=0 + renormalize to sum 1 */
          oBC[q] = b4[0];  oBC[q+nP] = b4[1];  oBC[q+2*nP] = b4[2];  oBC[q+3*nP] = b4[3];
        }
      } else {                             /* nothing within Dmax */
        oE[q] = 0.0;  oD[q] = INF;
        oCP[q] = oCP[q+nP] = oCP[q+2*nP] = mxGetNaN();
        if( wantBC ) oBC[q]=oBC[q+nP]=oBC[q+2*nP]=oBC[q+3*nP]=mxGetNaN();
      }
    }
  };

  mexClearInterrupt();                 /* arm Ctrl-C detection for this call */
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
  if( nlhs > 3 ) plhs[3] = mxBC; else mxDestroyArray( mxBC );
}
