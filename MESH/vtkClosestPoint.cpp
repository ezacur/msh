#include "mex.h"
#include <stdlib.h>
#include <math.h>
#include "MESH2vtkPolyData.h"

#include "vtkPolyData.h"
#include "vtkPointLocator.h"

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[]){

  double    *points;
  long      n_points, p;
  double    point[3];
  double    closest[3];
  int       n;
  double    dist;
  vtkPolyData     *MESH;
  vtkPointLocator *LOC  = vtkPointLocator::New();

  long      n_xyz;
  double    *xyz;

  n_xyz    =  mxGetM( mxGetField( prhs[0] , 0, "xyz") );
  xyz      = mxGetPr( mxGetField( prhs[0] , 0, "xyz") );
  
  n_points =  mxGetM( prhs[1] );
  points   = mxGetPr( prhs[1] );
  
  MESH= MESH2vtkPolyData( prhs[0] );
  
  LOC->SetDataSet( MESH );
  if( nrhs > 2 ) {
    LOC->SetNumberOfPointsPerBucket( *mxGetPr( prhs[2] ) );
  }
  LOC->BuildLocator();
  
  /* //create the outputs */
  plhs[0]= mxCreateDoubleMatrix( n_points,1,mxREAL );
  if( nlhs > 1) {
    plhs[1]= mxCreateDoubleMatrix( n_points,3,mxREAL );
  }
  if( nlhs > 2) {
    plhs[2]= mxCreateDoubleMatrix( n_points,1,mxREAL );
  }

  for( p=0 ; p<n_points ; p++ ) {
    point[0]= points[p];
    point[1]= points[p+n_points];
    point[2]= points[p+2*n_points];
    
    n= LOC->FindClosestPoint( point );

    *(mxGetPr(plhs[0])+p)= n+1;           /* //update the output id */
    if( nlhs > 1) {                       /* //update the output point */
      *(mxGetPr(plhs[1])+p           )= xyz[n];
      *(mxGetPr(plhs[1])+p+  n_points)= xyz[n+n_xyz];
      *(mxGetPr(plhs[1])+p+2*n_points)= xyz[n+2*n_xyz];
    }
    if( nlhs > 2) {                       /* //update the output distance */
      *(mxGetPr(plhs[2])+p)= sqrt(
              (point[0]-xyz[n])*(point[0]-xyz[n]) + 
              (point[1]-xyz[n+n_xyz])*(point[1]-xyz[n+n_xyz]) + 
              (point[2]-xyz[n+2*n_xyz])*(point[2]-xyz[n+2*n_xyz])
                                  );
    }
  }
  LOC->FreeSearchStructure();

  LOC->Delete();
  MESH->Delete();
}

