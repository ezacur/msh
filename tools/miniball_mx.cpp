#include "mex.h"

#include "miniball_mx.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  int     nX  = mxGetM( prhs[0] );
  int     nsd = (int) mxGetN( prhs[0] );
  double  *X  = mxGetPr( prhs[0] );
  
  std::vector<Seb::Point<double>> POINTS;
  std::vector<double> x(nsd);
  for( int i = 0 ; i < nX ; ++i ){
    for( int j = 0 ; j < nsd ; ++j ){
      x[j] = X[ i + j*nX ];
    }
    POINTS.push_back( Seb::Point<double>( nsd , x.begin() ) );
  }

  //'Seb' stands for "smalles enclosing ball"
  Seb::Smallest_enclosing_ball<double> MINIBALL( nsd , POINTS );

  plhs[0] = mxCreateDoubleMatrix( 1 , nsd , mxREAL );
  double *Y = mxGetPr( plhs[0] );

  Seb::Smallest_enclosing_ball<double>::Coordinate_iterator C = MINIBALL.center_begin();
  for( int j = 0 ; j < nsd ; ++j ){
    Y[j] = C[j];
  }
  
  if( nlhs > 1 ){
    plhs[1] = mxCreateDoubleScalar( (double) MINIBALL.radius() );
  }

}
