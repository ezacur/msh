/* triTriPairs_mx  --  interseccion de una LISTA de pares de triangulos.
 *
 *   [ K , XYZ , S , PR , CODE ] = triTriPairs_mx( V1 , F1 , V2 , F2 , pairs )
 *
 *     V1 (nV1 x 3), F1 (nF1 x 3, 1-based)   malla A
 *     V2 (nV2 x 3), F2 (nF2 x 3, 1-based)   malla B
 *     pairs (nP x 2, 1-based)               [ caraA , caraB ] a evaluar
 *
 *     K    nK x 8   claves canonicas: [ kind onB ia ib ja jb t s ]
 *                     kind 1 = vertice (ia)
 *                     kind 2 = punto en la arista canonica (ia<ib), parametro t
 *                     kind 3 = cruce de la arista (ia<ib) de A con la (ja<jb) de B
 *                     onB: 0 = la arista/vertice es de A, 1 = de B (kinds 1 y 2)
 *     XYZ  nK x 3   coordenadas, SOLO para salida (ninguna decision las usa)
 *     S    nS x 2   segmentos, como indices de fila en K (1-based)
 *     PR   nS x 2   el par ( caraA , caraB ) que produjo cada segmento
 *     CODE nP x 1   por par: 0 NONE, 1 SEGMENT, 2 TOUCH, 3 COPLANAR
 *
 *   Este es el nucleo de la interseccion malla-malla. Genera las claves; NO
 *   suelda: soldar es un `unique` por filas sobre K (mas la igualdad EXACTA de
 *   coordenadas para las claves de VERTICE, que es lo que resuelve el hueco de
 *   identidad entre mallas) y va en MATLAB.
 *
 *   La lista de pares se pasa desde fuera A PROPOSITO: primero se valida el
 *   pipeline con todos los pares (fuerza bruta) y despues se acelera la
 *   generacion de candidatos, que entonces tiene contra que compararse.
 *
 *   Compile:  mex triTriPairs_mx.cpp
 */
#include "mex.h"
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>
#include <immintrin.h>
#include "bvhKernels.h"

static void getTri( const double* V, mwSize nV, const double* F, mwSize nF,
                    mwSize f, double T[9], int g[3], const char* who )
{
  for( int k = 0; k < 3; ++k ) {
    const double gi = F[ f + k*nF ];
    const mwSize i  = (mwSize)gi - 1;
    if( gi != floor(gi) || gi < 1 || i >= nV )
      mexErrMsgIdAndTxt("triTriPairs_mx:conn","%s: node index out of range.",who);
    g[k] = (int)gi;
    T[3*k+0] = V[i]; T[3*k+1] = V[i+nV]; T[3*k+2] = V[i+2*nV];
  }
}

void mexFunction( int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[] )
{
  if( nrhs != 5 )
    mexErrMsgIdAndTxt("triTriPairs_mx:nargin","[K,XYZ,S,PR,CODE] = triTriPairs_mx(V1,F1,V2,F2,pairs)");
  for( int i = 0; i < 5; ++i )
    if( !mxIsDouble(prhs[i]) || mxIsComplex(prhs[i]) || mxIsSparse(prhs[i]) )
      mexErrMsgIdAndTxt("triTriPairs_mx:type","all inputs must be real double.");
  const mwSize nV1 = mxGetM(prhs[0]), nF1 = mxGetM(prhs[1]);
  const mwSize nV2 = mxGetM(prhs[2]), nF2 = mxGetM(prhs[3]);
  const mwSize nP  = mxGetM(prhs[4]);
  if( mxGetN(prhs[0]) != 3 || mxGetN(prhs[2]) != 3 )
    mexErrMsgIdAndTxt("triTriPairs_mx:V","V1,V2 must be nV x 3.");
  if( mxGetN(prhs[1]) != 3 || mxGetN(prhs[3]) != 3 )
    mexErrMsgIdAndTxt("triTriPairs_mx:F","F1,F2 must be nF x 3 (TRIANGLES only).");
  if( nP && mxGetN(prhs[4]) != 2 )
    mexErrMsgIdAndTxt("triTriPairs_mx:pairs","pairs must be nP x 2.");
  const double *V1 = mxGetPr(prhs[0]), *F1 = mxGetPr(prhs[1]);
  const double *V2 = mxGetPr(prhs[2]), *F2 = mxGetPr(prhs[3]);
  const double *PP = nP ? mxGetPr(prhs[4]) : NULL;

  std::vector<double> Kacc;   Kacc.reserve( 8*2*nP );
  std::vector<double> Xacc;   Xacc.reserve( 3*2*nP );
  std::vector<double> Sacc;   Sacc.reserve( 2*nP );
  std::vector<double> Racc;   Racc.reserve( 2*nP );
  std::vector<double> Cacc( nP , 0.0 );

  double T1[9], T2[9];  int g1[3], g2[3];

  for( mwSize p = 0; p < nP; ++p ) {
    const double fa = PP[p], fb = PP[p+nP];
    if( fa != floor(fa) || fb != floor(fb) || fa < 1 || fb < 1 ||
        (mwSize)fa > nF1 || (mwSize)fb > nF2 )
      mexErrMsgIdAndTxt("triTriPairs_mx:pairs","pair %llu out of range.",(unsigned long long)(p+1));
    const mwSize ia = (mwSize)fa - 1, ib = (mwSize)fb - 1;
    getTri( V1, nV1, F1, nF1, ia, T1, g1, "A" );
    getTri( V2, nV2, F2, nF2, ib, T2, g2, "B" );

    IxKey K[6];  int nK = 0;
    int code = triTriSegment( T1, g1, T2, g2, K, &nK );

    double XYZ[18];
    if( code == IX_COPLANAR ) {
      nK = triTriCoplanar( T1, g1, T2, g2, K, XYZ );
      if( nK == 0 ) { Cacc[p] = IX_NONE; continue; }
    } else {
      if( code == IX_NONE ) { Cacc[p] = IX_NONE; continue; }
      /* coordenadas de las claves del caso no coplanar */
      for( int k = 0; k < nK; ++k ) {
        /* IXK_EDGEX nombra SIEMPRE (arista de A , arista de B): su coordenada
         * sale de la arista de A, no de onB */
        const bool ex = ( K[k].kind == IXK_EDGEX );
        const double* Vk = ( ex || K[k].onB == 0 ) ? V1 : V2;
        const mwSize nVk = ( ex || K[k].onB == 0 ) ? nV1 : nV2;
        if( K[k].kind == IXK_VERTEX ) {
          const mwSize i = (mwSize)K[k].ia - 1;
          XYZ[3*k+0]=Vk[i]; XYZ[3*k+1]=Vk[i+nVk]; XYZ[3*k+2]=Vk[i+2*nVk];
        } else {
          const mwSize i = (mwSize)K[k].ia - 1, j = (mwSize)K[k].ib - 1;
          const double t = K[k].t;
          XYZ[3*k+0] = Vk[i]        + t*( Vk[j]        - Vk[i]        );
          XYZ[3*k+1] = Vk[i+nVk]    + t*( Vk[j+nVk]    - Vk[i+nVk]    );
          XYZ[3*k+2] = Vk[i+2*nVk]  + t*( Vk[j+2*nVk]  - Vk[i+2*nVk]  );
        }
      }
    }
    Cacc[p] = (double)code;

    const mwSize base = Kacc.size()/8;             /* 0-based, se pasa a 1 al final */
    for( int k = 0; k < nK; ++k ) {
      Kacc.push_back( (double)K[k].kind );
      Kacc.push_back( (double)K[k].onB  );
      Kacc.push_back( (double)K[k].ia   );
      Kacc.push_back( (double)K[k].ib   );
      Kacc.push_back( (double)K[k].ja   );
      Kacc.push_back( (double)K[k].jb   );
      Kacc.push_back( K[k].t );
      Kacc.push_back( K[k].s );
      Xacc.push_back( XYZ[3*k] ); Xacc.push_back( XYZ[3*k+1] ); Xacc.push_back( XYZ[3*k+2] );
    }
    /* segmentos: 2 claves -> 1 segmento; poligono coplanar de n>=3 -> n (cerrado) */
    if( nK == 2 ) {
      Sacc.push_back( (double)(base+1) ); Sacc.push_back( (double)(base+2) );
      Racc.push_back( fa );               Racc.push_back( fb );
    } else if( nK >= 3 ) {
      for( int k = 0; k < nK; ++k ) {
        const int k2 = ( k+1 == nK ) ? 0 : k+1;
        Sacc.push_back( (double)(base+1+k) ); Sacc.push_back( (double)(base+1+k2) );
        Racc.push_back( fa );                 Racc.push_back( fb );
      }
    }
  }

  const mwSize nK = Kacc.size()/8;
  const mwSize nS = Sacc.size()/2;
  plhs[0] = mxCreateDoubleMatrix( nK , 8 , mxREAL );
  double* oK = mxGetPr(plhs[0]);
  for( mwSize i = 0; i < nK; ++i )
    for( int c = 0; c < 8; ++c ) oK[i+c*nK] = Kacc[8*i+c];
  if( nlhs > 1 ) {
    plhs[1] = mxCreateDoubleMatrix( nK , 3 , mxREAL );
    double* oX = mxGetPr(plhs[1]);
    for( mwSize i = 0; i < nK; ++i )
      for( int c = 0; c < 3; ++c ) oX[i+c*nK] = Xacc[3*i+c];
  }
  if( nlhs > 2 ) {
    plhs[2] = mxCreateDoubleMatrix( nS , 2 , mxREAL );
    double* oS = mxGetPr(plhs[2]);
    for( mwSize i = 0; i < nS; ++i ) { oS[i] = Sacc[2*i]; oS[i+nS] = Sacc[2*i+1]; }
  }
  if( nlhs > 3 ) {
    plhs[3] = mxCreateDoubleMatrix( nS , 2 , mxREAL );
    double* oR = mxGetPr(plhs[3]);
    for( mwSize i = 0; i < nS; ++i ) { oR[i] = Racc[2*i]; oR[i+nS] = Racc[2*i+1]; }
  }
  if( nlhs > 4 ) {
    plhs[4] = mxCreateDoubleMatrix( nP , 1 , mxREAL );
    double* oC = mxGetPr(plhs[4]);
    for( mwSize i = 0; i < nP; ++i ) oC[i] = Cacc[i];
  }
}
