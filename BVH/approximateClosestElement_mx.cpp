/* approximateClosestElement_mx  --  compiled core of approximateClosestElement.
 *
 *   [e, cp, d, bc] = approximateClosestElement_mx( P , B , nthreads , Dmax )
 *
 *     bc   nP x 4 (only if requested) REGION-EXACT barycentrics of the winning
 *          fan element (same as bvhClosestElement_mx): edge/vertex hits are
 *          exact zeros, clamped >=0 and renormalized to sum 1. Misses -> NaN.
 *
 *   1-RING-OF-NEAREST-VERTEX approximate locator, two fused stages per point:
 *     1) nearest VERTEX via a BVH over the mesh vertices (point elements,
 *        AABB nodes; pkS centers with r = 0 are the vertices themselves), with
 *        a 4-wide AVX node test + pt4 leaf kernel;
 *     2) EXACT distance to the vertex's incident-element fan (EsuP, packed as
 *        CSR int32 fanStart/fanEl; pure-triangle fans also pre-packed as
 *        PreTri4 blocks fan4/fan4id/fan4Start swept with the shared tri4blk).
 *
 *   The result is an UPPER BOUND of the true distance; the winning element is
 *   typically exact for 95-99% of queries. Dmax (scalar | nP-vector) bounds the
 *   VERTEX distance (stage 1): a point whose nearest vertex is beyond Dmax
 *   gives e = 0, d = Inf, cp = NaN (non-finite query points: d = NaN).
 *
 *   Geometric primitives, region-exact barycentrics and the 4-wide Ericson
 *   kernel live in the shared header bvhKernels.h.
 *
 *   Compile (MSVC):  mex COMPFLAGS="$COMPFLAGS /openmp" -lut approximateClosestElement_mx.cpp
 *                    ( -lut links libut for Ctrl-C support; see mexInterrupt.h )
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

#include "bvhKernels.h"
#include "mexInterrupt.h"

/* --- pt4blk: 4 point-distances per SoA block [x0..x3 y0..y3 z0..z3] (the
 * approx locator's point-leaf kernel; every other primitive/kernel lives in
 * bvhKernels.h) --- */
static inline void pt4blk( const double* blk, const double* p, double d2o[4] )
{
  const __m256d dx = _mm256_sub_pd( _mm256_set1_pd(p[0]), _mm256_loadu_pd( blk     ) );
  const __m256d dy = _mm256_sub_pd( _mm256_set1_pd(p[1]), _mm256_loadu_pd( blk + 4 ) );
  const __m256d dz = _mm256_sub_pd( _mm256_set1_pd(p[2]), _mm256_loadu_pd( blk + 8 ) );
  _mm256_storeu_pd( d2o, mm_dot3( dx,dy,dz, dx,dy,dz ) );
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

  /* Dmax: scalar (uniform vertex-search radius) or nP-vector (per-point) */
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
  const bool wantBC = ( nlhs > 3 );      /* region-exact barycentrics on demand */
  plhs[0] = mxCreateDoubleMatrix( nP, 1, mxREAL );
  mxArray* mxCP = mxCreateDoubleMatrix( nP, 3, mxREAL );
  mxArray* mxD  = mxCreateDoubleMatrix( nP, 1, mxREAL );
  mxArray* mxBC = mxCreateDoubleMatrix( nP, wantBC ? 4 : 0, mxREAL );
  double* oE  = mxGetPr( plhs[0] );
  double* oCP = mxGetPr( mxCP );
  double* oBC = mxGetPr( mxBC );
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
      if( mexInterrupted( (long long)( qi - i0 ) ) ) break;   /* Ctrl-C: bail, let MATLAB abort */
      const mwSize q = (mwSize)ord[qi];
      const double p[3] = { P[q], P[q+nP], P[q+2*nP] };

      if( !std::isfinite(p[0]) || !std::isfinite(p[1]) || !std::isfinite(p[2]) ) {
        oE[q] = 0.0;  oD[q] = mxGetNaN();
        oCP[q] = oCP[q+nP] = oCP[q+2*nP] = mxGetNaN();
        if( wantBC ) oBC[q]=oBC[q+nP]=oBC[q+2*nP]=oBC[q+3*nP]=mxGetNaN();
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
        double lb2v[4];  bool have4 = false;
        if( avx ) {                          /* las 4 cajas del nodo de una vez */
          const __m256d px4 = _mm256_set1_pd( p[0] );
          const __m256d py4 = _mm256_set1_pd( p[1] );
          const __m256d pz4 = _mm256_set1_pd( p[2] );
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
          have4 = true;
        }
        for( int k = 0; k < 4; ++k ) {
          if( nc[k] == 0 ) continue;
          double lb2;
          if( have4 ) lb2 = lb2v[k];
          else {
            double dx = (double)nb[   k] - p[0];  if( p[0]-(double)nb[ 4+k] > dx ) dx = p[0]-(double)nb[ 4+k];  if( dx < 0 ) dx = 0;
            double dy = (double)nb[ 8+k] - p[1];  if( p[1]-(double)nb[12+k] > dy ) dy = p[1]-(double)nb[12+k];  if( dy < 0 ) dy = 0;
            double dz = (double)nb[16+k] - p[2];  if( p[2]-(double)nb[20+k] > dz ) dz = p[2]-(double)nb[20+k];  if( dz < 0 ) dz = 0;
            lb2 = dx*dx + dy*dy + dz*dz;
          }
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
          _mm_prefetch( FZ + (size_t)( nc[k]-1 )*stride, _MM_HINT_T0 );
          stkN[top] = nc[k]-1;  stkD[top] = key[a];  ++top;
        }
      }

      jwarm = bestJ;
      if( bestJ < 0 ) {                                    /* nearest vertex beyond Dmax */
        oE[q] = 0.0;  oD[q] = INF;
        oCP[q] = oCP[q+nP] = oCP[q+2*nP] = mxGetNaN();
        if( wantBC ) oBC[q]=oBC[q+nP]=oBC[q+2*nP]=oBC[q+3*nP]=mxGetNaN();
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
        if( wantBC ) {                     /* region-exact bc for the winner */
          double b4[4];
          bcElem( elV + (size_t)(fe-1)*12, elT[fe-1], p, b4 );
          finalizeBC( b4 );
          oBC[q] = b4[0];  oBC[q+nP] = b4[1];  oBC[q+2*nP] = b4[2];  oBC[q+3*nP] = b4[3];
        }
      } else {                                             /* empty fan (should not happen) */
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
