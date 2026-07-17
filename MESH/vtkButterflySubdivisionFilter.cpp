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

#include "mex.h"
#include "myMEX.h"

#define vtkOBJ_TYPE      vtkButterflySubdivisionFilter

#include "vtkButterflySubdivisionFilter.h"
#include "MESH2vtkPolyData.h"
#include "vtkPolyData2MESH.h"

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){

  if(!nrhs){
   	
    mexPrintf("vtkButterflySubdivisionFilter( MESH )\n");
    mexPrintf("\n");
    if( nlhs ){ for (int i=0; i<nlhs; i++) plhs[i]=mxCreateDoubleMatrix( 0 , 0 , mxREAL ); }
    return;
  }
  
  ALLOCATES();
  vtkPolyData         		*MESH;
  vtkOBJ_TYPE			*SUB;

  MESH = MESH2vtkPolyData( prhs[0] );
  
  SUB = vtkOBJ_TYPE::New();
  SUB->SetInput( MESH );
  
  /*Defaults*/
  /*END Defaults*/
  
  /*Parsing arguments*/
  /*END Parsing arguments*/
  
  SUB->Update();  
  plhs[0]= vtkPolyData2MESH( SUB->GetOutput() );

  EXIT:
    SUB->Delete();
    MESH->Delete();
    myFreeALLOCATES();
}

