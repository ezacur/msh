/* bvhTriPairs_mx  --  pares de caras candidatos, recorriendo el BVH de B con cada
 *                     triangulo de A.
 *
 *   pairs = bvhTriPairs_mx( VA , FA , B )
 *
 *     VA (nVA x 3)   vertices de A EN EL BUILD FRAME DE B (el llamante los
 *                    transforma: el blob guarda su geometria en ese marco)
 *     FA (nFA x 3)   conectividad de A, 1-based
 *     B              blob de BVH(MB). Solo cage 'aabb' (vol == 2): es la que
 *                    tiene el test de nodo trivial y es la de por defecto.
 *
 *     pairs (nP x 2) [ caraA , caraB ], 1-based
 *
 *   Sustituye al prefiltro O(nFA*nFB) por AABB que hacia bvhIntersectMesh en
 *   MATLAB. El test de nodo es la AABB del triangulo contra la AABB del slot: 6
 *   comparaciones, CONSERVADOR (nunca descarta una interseccion real, solo deja
 *   pasar candidatos de mas). Un test de ejes separadores completo (13 ejes)
 *   podaria mejor; se deja para cuando se mida que hace falta.
 *
 *   Los bounds del blob son float redondeados hacia FUERA por el builder, asi que
 *   compararlos en doble sigue siendo conservador.
 *
 *   Compile:  mex bvhTriPairs_mx.cpp
 */
#include "mex.h"
#include <cmath>
#include <cstdint>
#include <vector>

void mexFunction( int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[] )
{
  if( nrhs != 3 )
    mexErrMsgIdAndTxt("bvhTriPairs_mx:nargin","pairs = bvhTriPairs_mx(VA,FA,B)");
  if( !mxIsDouble(prhs[0]) || mxGetN(prhs[0]) != 3 )
    mexErrMsgIdAndTxt("bvhTriPairs_mx:VA","VA must be nVA x 3 double.");
  if( !mxIsDouble(prhs[1]) || mxGetN(prhs[1]) != 3 )
    mexErrMsgIdAndTxt("bvhTriPairs_mx:FA","FA must be nFA x 3 double (triangles).");
  if( !mxIsStruct(prhs[2]) )
    mexErrMsgIdAndTxt("bvhTriPairs_mx:B","B must be a BVH blob struct.");

  const mwSize nVA = mxGetM(prhs[0]), nFA = mxGetM(prhs[1]);
  const double* VA = mxGetPr(prhs[0]);
  const double* FA = mxGetPr(prhs[1]);

  const mxArray* Bs = prhs[2];
  auto fld = [&]( const char* n ) -> const mxArray* {
    const mxArray* f = mxGetField( Bs , 0 , n );
    if( !f ) mexErrMsgIdAndTxt("bvhTriPairs_mx:B","B lacks field '%s' (rebuild with BVH).",n);
    return f;
  };
  const mxArray *fN4 = fld("node4"), *fC4 = fld("child4"), *fR4 = fld("srange");
  const mxArray *fkE = fld("pkE"),   *fV  = fld("vol");
  if( !mxIsUint8(fN4) || !mxIsInt32(fC4) || !mxIsInt32(fR4) || !mxIsInt32(fkE) )
    mexErrMsgIdAndTxt("bvhTriPairs_mx:B","wrong blob field types (rebuild with BVH).");
  const int vol = (int)mxGetScalar(fV);
  if( vol != 2 )
    mexErrMsgIdAndTxt("bvhTriPairs_mx:vol",
                      "only the 'aabb' cage is supported (vol=2, got %d). Rebuild: BVH(M).", vol);

  const mwSize nN = mxGetN(fC4);
  if( mxGetM(fC4) != 4 || mxGetM(fR4) != 8 || mxGetN(fR4) != nN || nN == 0 )
    mexErrMsgIdAndTxt("bvhTriPairs_mx:B","inconsistent blob sizes (rebuild with BVH).");
  const mwSize nE     = mxGetNumberOfElements(fkE);
  const size_t stride = mxGetNumberOfElements(fN4) / (size_t)nN;
  if( stride * (size_t)nN != mxGetNumberOfElements(fN4) || stride < 24*4 + 16 )
    mexErrMsgIdAndTxt("bvhTriPairs_mx:B","bad node4 pool (rebuild with BVH).");

  const char*    FZ  = (const char*)mxGetData(fN4);
  const int32_t* Wr  = (const int32_t*)mxGetData(fR4);
  const int32_t* eii = (const int32_t*)mxGetData(fkE);

  std::vector<int32_t> pa, pb;   pa.reserve( 4*nFA );  pb.reserve( 4*nFA );

  for( mwSize f = 0; f < nFA; ++f ) {
    /* AABB del triangulo de A */
    double tlo[3], thi[3];
    for( int c = 0; c < 3; ++c ) { tlo[c] =  mxGetInf(); thi[c] = -mxGetInf(); }
    bool bad = false;
    for( int k = 0; k < 3; ++k ) {
      const double gi = FA[ f + k*nFA ];
      const mwSize i  = (mwSize)gi - 1;
      if( gi != floor(gi) || gi < 1 || i >= nVA )
        mexErrMsgIdAndTxt("bvhTriPairs_mx:conn","FA: node index out of range.");
      for( int c = 0; c < 3; ++c ) {
        const double v = VA[ i + c*nVA ];
        if( !mxIsFinite(v) ) { bad = true; break; }
        if( v < tlo[c] ) tlo[c] = v;
        if( v > thi[c] ) thi[c] = v;
      }
      if( bad ) break;
    }
    if( bad ) continue;                    /* triangulo con coordenada no finita */

    /* descenso */
    mwSize stk[256];  int top = 0;
    stk[top++] = 0;                        /* raiz */
    while( top ) {
      const mwSize ni = stk[--top];
      const char*    nz = FZ + (size_t)ni * stride;
      const float*   nb = (const float*)  nz;
      const int32_t* nc = (const int32_t*)( nz + 24*4 );
      const int32_t* nr = Wr + (size_t)ni * 8;
      for( int k = 0; k < 4; ++k ) {
        if( nc[k] == 0 ) continue;         /* slot vacio */
        /* AABB del slot: SoA por eje  ( xlo x4 , xhi x4 , ylo x4 , ... ) */
        if( thi[0] < (double)nb[     k] || tlo[0] > (double)nb[ 4 + k] ) continue;
        if( thi[1] < (double)nb[ 8 + k] || tlo[1] > (double)nb[12 + k] ) continue;
        if( thi[2] < (double)nb[16 + k] || tlo[2] > (double)nb[20 + k] ) continue;
        if( nc[k] > 0 ) {                  /* nodo interno */
          if( top >= 250 )
            mexErrMsgIdAndTxt("bvhTriPairs_mx:stack","traversal stack overflow (corrupt blob?).");
          stk[top++] = (mwSize)( nc[k] - 1 );
        } else {                           /* hoja: emitir sus elementos */
          const int32_t lo = nr[2*k] - 1, hi = nr[2*k+1] - 1;
          if( lo < 0 || hi < lo || (mwSize)hi >= nE )
            mexErrMsgIdAndTxt("bvhTriPairs_mx:B","corrupt slot range (rebuild with BVH).");
          for( int32_t j = lo; j <= hi; ++j ) {
            pa.push_back( (int32_t)( f + 1 ) );
            pb.push_back( eii[j] );
          }
        }
      }
    }
  }

  const mwSize nP = pa.size();
  plhs[0] = mxCreateDoubleMatrix( nP , 2 , mxREAL );
  double* o = mxGetPr( plhs[0] );
  for( mwSize i = 0; i < nP; ++i ) { o[i] = (double)pa[i]; o[i+nP] = (double)pb[i]; }
}
