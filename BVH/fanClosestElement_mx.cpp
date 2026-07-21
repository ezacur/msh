/* fanClosestElement_mx  --  compiled core of fanClosestElement (stage 2).
 *
 *   [e, cp, d, bc] = fanClosestElement_mx( P , V , TRI , nodes , ...
 *                                          fanStart , fanEl , nthreads , Dmax )
 *
 *   Pure FAN SWEEP: for each query point, the EXACT distance to every element
 *   incident to its seed node (the node's fan, CSR fanStart/fanEl), nothing
 *   else. NO tree, NO seeding of anything: the caller decides where the seed
 *   nodes come from (a point-BVH search, a heuristic, a constant...). The
 *   result is an UPPER BOUND of the true mesh distance (the fan is a real
 *   subset of the mesh).
 *
 *     P        nP x 3   double  query points
 *     V        nV x 3   double  mesh vertices
 *     TRI      nEl x kk int32   elements, 1-based, 0-padded (kk in 1..4);
 *                               celltype by nonzero count: 1 vertex, 2 segment,
 *                               3 triangle, 4 TETRAHEDRON
 *     nodes    nP x 1   int32   seed node per point, 1-based; 0 = no seed ->
 *                               miss (e = 0, d = Inf, cp/bc = NaN)
 *     fanStart nV+1 x 1 int32   CSR: node n's fan is fanEl( fanStart(n) :
 *                               fanStart(n+1)-1 ), 1-based; empty fan -> miss
 *     fanEl    nnz x 1  int32   element ids, 1-based
 *     nthreads scalar           OpenMP threads over the points
 *     Dmax     scalar | nP-vec  only elements at d < Dmax are returned
 *                               (default Inf); beyond -> miss
 *
 *     bc       nP x 4 (only if requested) REGION-EXACT barycentrics of the
 *              winner (same contract as bvhClosestElement_mx): exact zeros on
 *              edge/vertex hits, clamped >= 0, renormalized to sum 1.
 *
 *   Geometric primitives and barycentrics live in the shared bvhKernels.h;
 *   Ctrl-C support in mexInterrupt.h (link -lut).
 *
 *   Compile (MSVC):  mex COMPFLAGS="$COMPFLAGS /openmp" -lut fanClosestElement_mx.cpp
 *
 * See also fanClosestElement, bvhClosestElement_mx, BVH.
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

void mexFunction( int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[] )
{
  if( nrhs < 7 ) mexErrMsgIdAndTxt( "fanClosestElement_mx:nrhs",
    "usage: [e,cp,d,bc] = fanClosestElement_mx(P,V,TRI,nodes,fanStart,fanEl,nthreads,Dmax)" );

  /* ---- inputs, fully bounds-checked at entry ---------------------------- */
  if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) || mxGetN(prhs[0]) != 3 )
    mexErrMsgIdAndTxt( "fanClosestElement_mx:P", "P must be nP x 3 double." );
  if( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) || mxGetN(prhs[1]) != 3 )
    mexErrMsgIdAndTxt( "fanClosestElement_mx:V", "V must be nV x 3 double." );
  if( !mxIsInt32(prhs[2]) || mxGetN(prhs[2]) < 1 || mxGetN(prhs[2]) > 4 )
    mexErrMsgIdAndTxt( "fanClosestElement_mx:TRI", "TRI must be nEl x kk int32, kk in 1..4." );
  if( !mxIsInt32(prhs[3]) )
    mexErrMsgIdAndTxt( "fanClosestElement_mx:nodes", "nodes must be int32." );
  if( !mxIsInt32(prhs[4]) || !mxIsInt32(prhs[5]) )
    mexErrMsgIdAndTxt( "fanClosestElement_mx:fan", "fanStart/fanEl must be int32." );

  const mwSize  nP  = mxGetM( prhs[0] );
  const mwSize  nV  = mxGetM( prhs[1] );
  const mwSize  nEl = mxGetM( prhs[2] );
  const int     kk  = (int)mxGetN( prhs[2] );
  const mwSize  nnz = mxGetNumberOfElements( prhs[5] );

  const double*  P    = mxGetPr( prhs[0] );
  const double*  Vp   = mxGetPr( prhs[1] );
  const int32_t* TRIp = (const int32_t*)mxGetData( prhs[2] );
  const int32_t* nod  = (const int32_t*)mxGetData( prhs[3] );
  const int32_t* fS   = (const int32_t*)mxGetData( prhs[4] );
  const int32_t* fE   = (const int32_t*)mxGetData( prhs[5] );

  if( mxGetNumberOfElements( prhs[3] ) != nP )
    mexErrMsgIdAndTxt( "fanClosestElement_mx:nodes", "nodes must have nP elements." );
  if( mxGetNumberOfElements( prhs[4] ) != nV + 1 )
    mexErrMsgIdAndTxt( "fanClosestElement_mx:fan", "fanStart must have nV+1 elements." );
  if( fS[0] != 1 || (mwSize)fS[nV] != nnz + 1 )
    mexErrMsgIdAndTxt( "fanClosestElement_mx:fan", "fanStart is not a valid 1-based CSR." );
  for( mwSize i = 0; i < nV; ++i )
    if( fS[i+1] < fS[i] )
      mexErrMsgIdAndTxt( "fanClosestElement_mx:fan", "fanStart must be non-decreasing." );
  for( mwSize i = 0; i < nnz; ++i )
    if( fE[i] < 1 || (mwSize)fE[i] > nEl )
      mexErrMsgIdAndTxt( "fanClosestElement_mx:fan", "fanEl out of range." );
  for( mwSize i = 0; i < nP; ++i )
    if( nod[i] < 0 || (mwSize)nod[i] > nV )
      mexErrMsgIdAndTxt( "fanClosestElement_mx:nodes", "nodes out of range." );
  for( mwSize i = 0; i < (mwSize)nEl*kk; ++i )
    if( TRIp[i] < 0 || (mwSize)TRIp[i] > nV )
      mexErrMsgIdAndTxt( "fanClosestElement_mx:TRI", "TRI out of range." );

  int nt = ( nrhs > 6 ) ? (int)mxGetScalar( prhs[6] ) : 1;
  if( nt < 1 ) nt = 1;  if( nt > 64 ) nt = 64;

  /* Dmax: scalar (uniform) or nP-vector (per point) */
  const double* DmaxV = NULL;
  double Dmax = INF;
  if( nrhs > 7 ) {
    if( !mxIsDouble(prhs[7]) || mxIsComplex(prhs[7]) )
      mexErrMsgIdAndTxt( "fanClosestElement_mx:Dmax", "Dmax must be double." );
    const mwSize nD = mxGetNumberOfElements( prhs[7] );
    if( nD == 1 )        Dmax  = mxGetScalar( prhs[7] );
    else if( nD == nP )  DmaxV = mxGetPr( prhs[7] );
    else mexErrMsgIdAndTxt( "fanClosestElement_mx:Dmax", "Dmax must be a scalar or an nP-vector." );
  }

  /* ---- outputs ---------------------------------------------------------- */
  const bool wantBC = ( nlhs > 3 );
  plhs[0] = mxCreateDoubleMatrix( nP, 1, mxREAL );
  mxArray* mxCP = mxCreateDoubleMatrix( nP, 3, mxREAL );
  mxArray* mxD  = mxCreateDoubleMatrix( nP, 1, mxREAL );
  mxArray* mxBC = mxCreateDoubleMatrix( nP, wantBC ? 4 : 0, mxREAL );
  double* oE  = mxGetPr( plhs[0] );
  double* oCP = mxGetPr( mxCP );
  double* oD  = mxGetPr( mxD );
  double* oBC = wantBC ? mxGetPr( mxBC ) : NULL;
  const double NaN = mxGetNaN();

  auto runRange = [&]( mwSize i0, mwSize i1 )
  {
    for( mwSize q = i0; q < i1; ++q ) {
      if( mexInterrupted( (long long)( q - i0 ) ) ) break;   /* Ctrl-C: bail, let MATLAB abort */
      const double p[3] = { P[q], P[q+nP], P[q+2*nP] };

      oE[q] = 0.0;
      oCP[q] = oCP[q+nP] = oCP[q+2*nP] = NaN;
      if( wantBC ) oBC[q] = oBC[q+nP] = oBC[q+2*nP] = oBC[q+3*nP] = NaN;

      if( !std::isfinite(p[0]) || !std::isfinite(p[1]) || !std::isfinite(p[2]) ) {
        oD[q] = NaN;                                       /* non-finite point */
        continue;
      }
      oD[q] = INF;                                         /* miss until proven */

      const int32_t n = nod[q];
      if( n < 1 ) continue;                                /* no seed -> miss  */

      const double dmx  = DmaxV ? DmaxV[q] : Dmax;
      const double dmx2 = ( dmx < INF ) ? dmx*dmx : INF;

      double  best2 = INF, bcp[3] = { 0.0, 0.0, 0.0 };
      int32_t bestE = -1;
      for( int32_t f = fS[n-1] - 1; f < fS[n] - 1; ++f ) {
        const int32_t ei = fE[f] - 1;
        double vv[12];  int kt = 0;
        for( int c = 0; c < kk; ++c ) {
          const int32_t id = TRIp[ (size_t)ei + (size_t)c*nEl ];
          if( id == 0 ) break;
          vv[3*kt  ] = Vp[ (size_t)id-1            ];
          vv[3*kt+1] = Vp[ (size_t)id-1 +   (size_t)nV ];
          vv[3*kt+2] = Vp[ (size_t)id-1 + 2*(size_t)nV ];
          ++kt;
        }
        if( !kt ) continue;                                /* all-zero row     */
        double cpe[3];
        const double dd = d2Elem( vv, kt, p, cpe );
        if( dd < best2 ) {
          best2 = dd;  bestE = ei;
          bcp[0] = cpe[0];  bcp[1] = cpe[1];  bcp[2] = cpe[2];
        }
      }
      if( bestE < 0 || !( best2 < dmx2 ) ) continue;       /* empty fan / beyond Dmax */

      oE[q]      = (double)( bestE + 1 );
      oD[q]      = std::sqrt( best2 );
      oCP[q]     = bcp[0];  oCP[q+nP] = bcp[1];  oCP[q+2*nP] = bcp[2];
      if( wantBC ) {                       /* region-exact bc of the winner   */
        double vv[12];  int kt = 0;
        for( int c = 0; c < kk; ++c ) {
          const int32_t id = TRIp[ (size_t)bestE + (size_t)c*nEl ];
          if( id == 0 ) break;
          vv[3*kt  ] = Vp[ (size_t)id-1            ];
          vv[3*kt+1] = Vp[ (size_t)id-1 +   (size_t)nV ];
          vv[3*kt+2] = Vp[ (size_t)id-1 + 2*(size_t)nV ];
          ++kt;
        }
        double b4[4];
        bcElem( vv, kt, p, b4 );
        finalizeBC( b4 );
        oBC[q] = b4[0];  oBC[q+nP] = b4[1];  oBC[q+2*nP] = b4[2];  oBC[q+3*nP] = b4[3];
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
