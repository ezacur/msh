//mmex dir__     = fileparts( which( 'vec.m' ) );
//mmex VTK_source_dir__ = fullfile( fileparts( dir__ ) , 'MESHES\vtk\sources' );
//mmex 
//mmex INCLUDES = { '' ; 'Common' ; 'Filtering' ; 'IO' ; 'Graphics' };
//mmex INCLUDES = cellfun(@(f)fullfile(VTK_source_dir__,f),INCLUDES,'UniformOutput',false);
//mmex INCLUDES = [ INCLUDES ; dir__ ; fullfile( fileparts( dir__ ) , 'MESHES' ) ];
//mmex 
//mmex VTK_lib_dir__ = fullfile( fileparts( dir__ ) , 'MESHES\vtk' );
//mmex switch computer
//mmex    case 'PCWIN64', VTK_lib_dir__ = fullfile( VTK_lib_dir__ , 'w64\' );
//mmex    case 'MACI64',  VTK_lib_dir__ = fullfile( VTK_lib_dir__ , 'maci64/' );
//mmex end
//mmex LIBS = { VTK_lib_dir__ ; 'vtkCommon' ; 'vtkFiltering' ; 'vtkGraphics' ; 'vtkIO' };
//mmex 

/*

[ element_id , xyz_closest_point , distance , barycentric_coordinates ] = vtkClosestElement( MESH , [ point1 ; point2 ; ... ] );


 
 vtkClosestElement( MESH )                                  crea el locator
 [outs] = vtkClosestElement( [ point1 ; point2 ; ...] )     sobre el ultimo locator creado, calcula los puntos
 vtkClosestElement( [] , [] )                               libera el locator


 ============================================================================
 KNOWN BUG (4th output, barycentric_coordinates) + HOW TO FIX WHEN RECOMPILING
 ============================================================================
 The barycentric block below computes the weights by solving the ORIGIN-based
 linear system   M * w = closest   with   M = [ v1 v2 v3 ]   (the three
 triangle vertices as COLUMNS, i.e. position vectors from the origin), via
 Cramer's rule, dividing by
        DET = det(M) = v1 . (v2 x v3) = 6 * signed volume of the
                       tetrahedron ( ORIGIN , v1 , v2 , v3 ).
 For a point ON the triangle plane this DOES yield the true barycentric
 coordinates -- BUT ONLY while DET is well conditioned. DET -> 0 whenever the
 triangle's supporting plane passes near the coordinate ORIGIN, and DET == 0
 EXACTLY for any triangle that has a vertex at the origin (a zero column).
 There the 1/DET division blows up: the returned weights stop summing to 1 and
 go negative (e.g. [1,-1,0.5] for a centroid whose triangle touches (0,0,0)).
 It is a TRANSLATION-VARIANCE bug (the opposite of ipd, which loses precision
 FAR from the origin; this one loses it NEAR the origin). It bites meshes that
 touch/straddle the origin (unit sphere, normalized meshes, a vertex at 0);
 meshes at large offsets (e.g. DICOM coordinates) are safe (huge DET).

 FIX (translation-invariant, via CROSS products). Use the triangle's OWN edges
 instead of the origin. With a=v1, b=v2, c=v3, p=closest:
        e0 = b - a;   e1 = c - a;   e2 = p - a;
        n  = cross(e0,e1);          // triangle normal, |n| = 2*area
        nn = n.n;                   // = 4*area^2  (>0 for any non-degenerate
                                    //   triangle, origin-free, NO cancellation)
        gamma = ( n . cross(e0,e2) ) / nn;   // weight of v3 (area a-b-p / a-b-c)
        beta  = ( n . cross(e2,e1) ) / nn;   // weight of v2 (area a-p-c / a-b-c)
        alpha = 1 - beta - gamma;            // weight of v1
 Replace the DET/Cramer block below with this and the bug disappears.

 IMPORTANT -- do NOT use the Gram dot form  den = d00*d11 - d01*d01  (the
 classic Ericson p.47 code). It is origin-free and fixes the translation bug,
 but for a THIN/SLIVER triangle that determinant is (|e0|^2*|e1|^2)*(1-cos^2)
 = e.g. (0.25+h^2) - 0.25 -> CATASTROPHIC CANCELLATION: the weights lose all
 precision as the triangle gets thin (measured through the .m shadow: ~5% error
 at aspect 1e8, NaN at 1e10). The cross form above computes the SAME 4*area^2
 as |e0 x e1|^2 DIRECTLY, so it stays at machine precision even at aspect 1e12.
 (The .m shadow uses this cross form; it was switched away from the dot form
 after slivers were found to break it.)

 UNTIL THEN: this MEX is renamed vtkClosestElement_ and shadowed by
 vtkClosestElement.m, which delegates id/closest/distance here but recomputes
 the barycentric weights with the cross form above (and clamps the tiny [0,1]
 spill on edge points, then renormalizes so each row is a valid convex
 combination summing to 1). Once this block is fixed and recompiled, delete
 vtkClosestElement.m (and rename this back).
 (Also, the double-cast / 3-column checks that vtkClosestElement.m does on
 MESH.xyz, MESH.tri and the query points should ideally be enforced here too.)

 */

#include "mex.h"
#include <stdlib.h>
#include <math.h>
#include "MESH2vtkPolyData.h"

#include "vtkPolyData.h"
#include "vtkCellLocator.h"
#include "vtkGenericCell.h"

struct LOCATOR {
    vtkCellLocator *LOC;
    vtkPolyData    *MESH;
    vtkGenericCell *CELL;
};
typedef struct LOCATOR LOCATOR;

static LOCATOR L;


void clean(){

  try{ 
    if( L.CELL != NULL ){
      L.CELL->Delete(); 
    }
    L.CELL = NULL;
/* //     mexPrintf("CELL  deleted\n"); */
  } catch( char * str ) {
    mexPrintf("error DELETING CELL\n");
  }

  try{ 
    if( L.LOC != NULL ){
      L.LOC->Delete(); 
    }
    L.LOC = NULL;
/* //     mexPrintf("LOC  deleted\n"); */
  } catch( char * str ) {
    mexPrintf("error DELETING LOC\n");
  }

  try{ 
    if( L.MESH != NULL ){
      L.MESH->Delete(); 
    }
    L.MESH = NULL;
/* //     mexPrintf("MESH  deleted\n"); */
  } catch( char * str ) {
    mexPrintf("error DELETING MESH\n");
  }

}


void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ){
  double    *points;
  long      n_points, p;
  double    point[3];
  double    closest[3];
  int       subId, t;
  double    dist;
  double    *O_id, *O_x, *O_y, *O_z, *O_distance, *O_barycentric_coordinates;
  double    xyz[3], M[9], DET;
  
  
  if( nrhs == 1  &&  mxIsStruct( prhs[0] ) ){

   /*  //     vtkClosestElement( MESH ) */
    
    clean();

    L.LOC  = vtkCellLocator::New();  
    L.CELL = vtkGenericCell::New();  
    L.MESH = MESH2vtkPolyData( prhs[0] );

    L.LOC->SetDataSet( L.MESH );
    L.LOC->SetNumberOfCellsPerBucket( 5 );
    L.LOC->BuildLocator();
    
/* //     mexPrintf("\n\nCREADO\n\n"); */
    
  } else if( nrhs == 1  &&  mxIsNumeric(prhs[0]) ){
   /*  //     vtkClosestElement( [ point1 ; point2 ] ) */
    

    if( L.LOC == NULL  ||  L.CELL == NULL  || L.MESH == NULL ){
      mexPrintf("No hay creado un LOCATOR.");
    }

    point[0] = 0; point[1] = 0; point[2] = 0;
    try {
      L.LOC->FindClosestPoint( point , closest , L.CELL , ((vtkIdType&)t) , subId , dist);
    } catch( char * str ) {
      mexPrintf("el LOCATOR no responde\n");
      clean();
      mexPrintf("el LOCATOR no responde.");
    }    
    
    n_points =  mxGetM( prhs[0] );
    points   = mxGetPr( prhs[0] );

   /*  //create the outputs */
    plhs[0]= mxCreateDoubleMatrix( n_points,1,mxREAL );
    O_id   = mxGetPr( plhs[0] );
    if( nlhs > 1) {
      plhs[1] = mxCreateDoubleMatrix( n_points,3,mxREAL );
      O_x   = mxGetPr( plhs[1] );
      O_y   = O_x + n_points;
      O_z   = O_y + n_points;
    }
    if( nlhs > 2) {
      plhs[2] = mxCreateDoubleMatrix( n_points,1,mxREAL );
      O_distance = mxGetPr( plhs[2] );
    }
    if( nlhs > 3) {
      plhs[3] = mxCreateDoubleMatrix( n_points,3,mxREAL );
      O_barycentric_coordinates = mxGetPr( plhs[3] );
    }

    for( p=0 ; p<n_points ; p++ ) {
      point[0]= points[p];
      point[1]= points[p+n_points];
      point[2]= points[p+2*n_points];

      L.LOC->FindClosestPoint( point , closest , L.CELL , ((vtkIdType&)t) , subId , dist);

      O_id[ p ]= t+1;         /*   //update the output id */

      if( nlhs > 1) {                   /*     //update the output point */
        O_x[ p ] = closest[0];
        O_y[ p ] = closest[1];
        O_z[ p ] = closest[2];
      }
      if( nlhs > 2) {                    /*    //update the output distance */
        O_distance[p] = sqrt( dist );
      } 
      if( nlhs > 3) {                     /*   //barycentric coordinates */
        #define M(i,j)  M[ i-1 + (j-1)*3 ]
        L.MESH->GetPoint( L.MESH->GetCell(t)->GetPointId(0), xyz );
        M(1,1) = xyz[0];
        M(2,1) = xyz[1];
        M(3,1) = xyz[2];
        
        L.MESH->GetPoint( L.MESH->GetCell(t)->GetPointId(1), xyz );
        M(1,2) = xyz[0];
        M(2,2) = xyz[1];
        M(3,2) = xyz[2];
        
        L.MESH->GetPoint( L.MESH->GetCell(t)->GetPointId(2), xyz );
        M(1,3) = xyz[0];
        M(2,3) = xyz[1];
        M(3,3) = xyz[2];
        
        DET = 1.0/( M(3,1)*( M(1,3)*M(2,2) - M(1,2)*M(2,3) ) + M(2,1)*( M(1,2)*M(3,3) - M(1,3)*M(3,2) ) + M(1,1)*( M(2,3)*M(3,2) - M(2,2)*M(3,3) ) );
        
        O_barycentric_coordinates[ p              ] = ( closest[2]*( M(1,3)*M(2,2) - M(1,2)*M(2,3) ) + closest[1]*( M(1,2)*M(3,3) - M(1,3)*M(3,2) ) + closest[0]*( M(2,3)*M(3,2) - M(2,2)*M(3,3) ) )*DET;
        O_barycentric_coordinates[ p +   n_points ] = ( closest[2]*( M(1,1)*M(2,3) - M(1,3)*M(2,1) ) + closest[1]*( M(1,3)*M(3,1) - M(1,1)*M(3,3) ) + closest[0]*( M(2,1)*M(3,3) - M(2,3)*M(3,1) ) )*DET;
        O_barycentric_coordinates[ p + 2*n_points ] = ( closest[2]*( M(1,2)*M(2,1) - M(1,1)*M(2,2) ) + closest[1]*( M(1,1)*M(3,2) - M(1,2)*M(3,1) ) + closest[0]*( M(2,2)*M(3,1) - M(2,1)*M(3,2) ) )*DET;
      }

    }
    
  } else if( nrhs == 2  &&  mxIsStruct(prhs[0])  &&  mxIsNumeric(prhs[1]) ){
  /*   //     vtkClosestElement( MESH , [ point1 ; point2 ] ) */
    
    clean();    
    
    L.LOC  = vtkCellLocator::New();
    L.CELL = vtkGenericCell::New();
    L.MESH = MESH2vtkPolyData( prhs[0] );

    L.LOC->SetDataSet( L.MESH );
   /*  //     L.LOC->SetNumberOfCellsPerBucket( 5 ); */
    L.LOC->BuildLocator();

    
    n_points =  mxGetM( prhs[1] );
    points   = mxGetPr( prhs[1] );

   /*  //create the outputs */
    plhs[0]= mxCreateDoubleMatrix( n_points,1,mxREAL );
    O_id   = mxGetPr( plhs[0] );
    if( nlhs > 1) {
      plhs[1] = mxCreateDoubleMatrix( n_points,3,mxREAL );
      O_x   = mxGetPr( plhs[1] );
      O_y   = O_x + n_points;
      O_z   = O_y + n_points;
    }
    if( nlhs > 2) {
      plhs[2] = mxCreateDoubleMatrix( n_points,1,mxREAL );
      O_distance = mxGetPr( plhs[2] );
    }
    if( nlhs > 3) {
      plhs[3] = mxCreateDoubleMatrix( n_points,3,mxREAL );
      O_barycentric_coordinates = mxGetPr( plhs[3] );
    }

    
    for( p=0 ; p<n_points ; p++ ) {
      point[0]= points[p];
      point[1]= points[p+n_points];
      point[2]= points[p+2*n_points];

      L.LOC->FindClosestPoint( point , closest , L.CELL , ((vtkIdType&)t) , subId , dist);

      O_id[ p ]= t+1;         /*   //update the output id */

      if( nlhs > 1) {                /*       //update the output point */
        O_x[ p ] = closest[0];
        O_y[ p ] = closest[1];
        O_z[ p ] = closest[2];
      }
      if( nlhs > 2) {                      /*  //update the output distance */
        O_distance[p] = sqrt( dist );
      }
      
      if( nlhs > 3) {                     /*   //barycentric coordinates */
        #define M(i,j)  M[ i-1 + (j-1)*3 ]
        L.MESH->GetPoint( L.MESH->GetCell(t)->GetPointId(0), xyz );
        M(1,1) = xyz[0];
        M(2,1) = xyz[1];
        M(3,1) = xyz[2];
        
        L.MESH->GetPoint( L.MESH->GetCell(t)->GetPointId(1), xyz );
        M(1,2) = xyz[0];
        M(2,2) = xyz[1];
        M(3,2) = xyz[2];
        
        L.MESH->GetPoint( L.MESH->GetCell(t)->GetPointId(2), xyz );
        M(1,3) = xyz[0];
        M(2,3) = xyz[1];
        M(3,3) = xyz[2];
        
        DET = 1.0/( M(3,1)*( M(1,3)*M(2,2) - M(1,2)*M(2,3) ) + M(2,1)*( M(1,2)*M(3,3) - M(1,3)*M(3,2) ) + M(1,1)*( M(2,3)*M(3,2) - M(2,2)*M(3,3) ) );
        
        O_barycentric_coordinates[ p              ] = ( closest[2]*( M(1,3)*M(2,2) - M(1,2)*M(2,3) ) + closest[1]*( M(1,2)*M(3,3) - M(1,3)*M(3,2) ) + closest[0]*( M(2,3)*M(3,2) - M(2,2)*M(3,3) ) )*DET;
        O_barycentric_coordinates[ p +   n_points ] = ( closest[2]*( M(1,1)*M(2,3) - M(1,3)*M(2,1) ) + closest[1]*( M(1,3)*M(3,1) - M(1,1)*M(3,3) ) + closest[0]*( M(2,1)*M(3,3) - M(2,3)*M(3,1) ) )*DET;
        O_barycentric_coordinates[ p + 2*n_points ] = ( closest[2]*( M(1,2)*M(2,1) - M(1,1)*M(2,2) ) + closest[1]*( M(1,1)*M(3,2) - M(1,2)*M(3,1) ) + closest[0]*( M(2,2)*M(3,1) - M(2,1)*M(3,2) ) )*DET;
      }
      
      
    }

    clean();

  } else if(  nrhs == 2  &&  mxIsNumeric(prhs[0])  &&  mxIsNumeric(prhs[1])  ){
 /*   // vtkClosestElement( [] , [] ) */

    clean();
/* //     mexPrintf("LOCATOR borrado !!\n" ); */

  } else {
    /* Unrecognized call (this includes nrhs==0): the MEX does NOTHING useful
       here -- it just prints and returns, leaving any locator untouched and no
       outputs set. In particular vtkClosestElement_() (no arguments) does NOT
       free the locator; use ( [] , [] ) for that. The shadow vtkClosestElement.m
       turns a no-argument call into an error (vtkClosestElement:signature)
       instead of silently doing nothing. When this MEX is fixed/recompiled,
       consider mexErrMsgTxt here (and add the input validation the shadow does:
       MESH.xyz / MESH.tri as double with 3 columns, query points as double). */
    mexPrintf("No entiendo la llamada.");
  }

}
