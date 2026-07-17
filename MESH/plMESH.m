function plMESH( M , celltype )

  PARAVIEW_executable = findfilename( fullfile( '<temp>','..' ,'ParaView','bin','paraview.exe') ,...
                                      fullfile( '<temp>','.'  ,'ParaView','bin','paraview.exe') ,...
                                      fullfile( 'C:\'   ,'Win','ParaView','bin','paraview.exe') ,...
                                      false );
  if ~isfile( PARAVIEW_executable )
    try,
      PARAVIEW_executable = getoption( 'PARAVIEW' , 'executable' );
      if isempty( PARAVIEW_executable )
        setoption( 'PARAVIEW' , 'executable' , '"C:\Win\Paraview\bin\paraview.exe"' );
      end
    end
    if isempty( PARAVIEW_executable ) || ~isfile( PARAVIEW_executable )
      try, setoption(); end
      error('Check the [PARAVIEW] executable option' );
    end
  end
  if isempty( PARAVIEW_executable ) || ~isfile( PARAVIEW_executable )
    error('Set the PARAVIEW_executable variable.' );
  end


  fn = inputname(1);
  if isempty( fn ), fn = 'mesh'; end

  if iscell( M )
    M = MeshAppend( M{:} ,'kp');
    fn = ['appended_',fn];
  end

  [fn,CLEANOUT] = sandfile( [ fn , '_???.vtk' ] );

  M.celltype = meshCelltype( M );
  if nargin > 1, M.celltype = celltype; end
  if isscalar(M.celltype) &&  M.celltype == 5
    M = rmfield( M , 'celltype' );
  end
  
  M.xyz(:,end+1:3) = 0;

  write_VTK( M , fn , 'binary' );

  [a,b] = system( [ '"' , PARAVIEW_executable , '" "' , fn , '"' ] );

end
