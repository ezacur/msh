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

#define   real       double
#define   mxREAL_CLASS       mxDOUBLE_CLASS

#include "vtkPolyDataReader.h"
#include "vtkTriangleFilter.h"

#include "vtkPolyData2MESH.h"

void mexFunction( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[]){
  char                STR[2048];
  
  if(!nrhs){
    mexPrintf("vtkPolyDataReader( FileName )\n");
    mexPrintf("\n");
    return;
  }
  
  ALLOCATES();

  mxGetString( prhs[0], STR, 1999 );
  
  
  vtkPolyDataReader *R = vtkPolyDataReader::New();
  R->SetFileName( STR );

  R->Update();
  

  vtkTriangleFilter   *T = vtkTriangleFilter::New();
  T->SetInput( R->GetOutput() );
  T->PassVertsOff();
  T->PassLinesOff();
  T->Update();  
  
  plhs[0] = vtkPolyData2MESH( T->GetOutput() );
  
  EXIT:
    T->Delete();
    R->Delete();
    myFreeALLOCATES();

}

