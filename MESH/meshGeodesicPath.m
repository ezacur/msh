function G = meshGeodesicPath( M , nodes )


  nodes = M.xyz( nodes ,:);
  [Tdir,CLEAN] = sanddir( 'meshGeodesicPath_***' );

  M = Mesh( M ,0);
  M = MeshTidy( M ,0,true,[1,1,1,0]);
  nodes = knnsearch( M.xyz , nodes );

  write_OBJ( M , fullfile( Tdir , 'M.obj' ) );

  cmd = [];
  cmd = [ cmd , 'import potpourri3d as pp3d;' ];
  cmd = [ cmd , 'from scipy.io import savemat;' ];
  cmd = [ cmd , 'V,F = pp3d.read_mesh(''M.obj'');' ];
  cmd = [ cmd , 'path_solver = pp3d.EdgeFlipGeodesicSolver(V,F);' ];
  cmd = [ cmd , 'G = path_solver.find_geodesic_path(' , uneval( nodes(1) - 1 ) , ',' , uneval( nodes(2) - 1 ) , ');' ];
  cmd = [ cmd , 'savemat(''G.mat'', {''G'': G})' ];
  cmd = sprintf( 'cd /d "%s" & python -c "%s"' , Tdir , cmd );
  [a,v] = system( cmd );

  G = load( fullfile( Tdir , 'G.mat' ) );
  G = G.G;

end
