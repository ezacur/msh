function W = MeshClip_tetrahedralize_helper( W )

%   if ispc
%     cmd = getoption( 'VTK','tcl_interpreter' , [] );
%     if isempty( cmd )
%       cmd = sprintf('set VTK_PAKAGE_PATH=%s' , fullfile( fileparts( mfilename('fullpath') ) , 'TCL' , 'vtk' ) );
%       cmd = [ cmd , '&' , sprintf('set PATH=%%PATH%%;%s' , fullfile( fileparts( mfilename('fullpath') ) , 'TCL' , 'bin' ) ) ];
%       cmd = [ cmd , '&' , '"' , fullfile( fileparts( mfilename('fullpath') ) , 'TCL' , 'vtk.exe' ) , '"' ];
%     end
%   else,    cmd = getoption( 'VTK','tcl_interpreter' , 'export TCL_LIBRARY=/usr/lib/tcl8.5; export TK_LIBRARY=/usr/lib/tk8.5; /usr/local/vtk' );
%   end

  [p,f,e] = fileparts( mfilename('fullpath') );
  PY_script = fullfile( p , [ f , '.py' ] );
  
  [ INfile  , CLEAN_IN  ] = sandfile( 'IN_*****.vtk'  );
  [ OUTfile , CLEAN_OUT ] = sandfile( 'OUT_*****.vtk' );
  
  [W,IDSname] = MeshGenerateIDs( W , 'xyz_' );

  V = W;
  for f = fieldnames( V ).', f = f{1};
    if strcmp( f , 'celltype' ), continue; end
    if strcmp( f , 'xyz' ), continue; end
    if strcmp( f , IDSname ), continue; end
    if strncmp( f , 'tri' , 3 ), continue; end
    V = rmfield( V , f );
  end
  write_VTK_UNSTRUCTURED_GRID( V , INfile , 'binary' );

  cmd = 'python ';
  cmd = [ cmd , ' "' , PY_script  , '"' ];
  cmd = [ cmd , ' "' , INfile      , '"' ];
  cmd = [ cmd , ' "' , OUTfile     , '"' ];
  
  [status,cmdout] = system( cmd );
  if(status),fprintf(2,'Some error executing:');fprintf(1,'\n');disp(cmdout);fprintf(1,'\n\n');end

  T = read_VTK( OUTfile );
  
  T.tri(:,[3,4]) = T.tri(:,[4,3]);
%   try
%     w = meshQuality( T , 'orientation' ) < 0;
%     T.tri(w,[3,4]) = T.tri(w,[4,3]);
%   end
  W.celltype = T.celltype;
  
  if isequal( W.(IDSname) , T.(IDSname) )
    for f = fieldnames( T ).', f = f{1};
      if ~strncmp( f , 'tri' , 3 ), continue; end
      W.(f) = T.(f);
    end
  else
    warning('nodes coordinates (or ordering) had changed !!!');
  end
  
  W = rmfield( W , IDSname );

end
