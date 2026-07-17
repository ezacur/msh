function M = MeshDefragTexture( M )

  [dn,dn_] = sanddir( 'MeshDefragTexture_*****' );
  write_GLB( M , fullfile( dn , 'M.glb' ) );
  
  cmd = '';
  if isempty( cmd )
    cmd = findfilename( false , {'<hpp>'} , 'MeshDefragTexture.exe' );
    if ~isempty( cmd ), cmd = [ '"' , cmd , '"' ]; end
  end
  if isempty( cmd )
    cmd = findfilename( false , {'<hpp>'} , 'MeshDefragTexture.py' );
    cmd = sprintf( 'conda activate env4meshlab & python "%s"', cmd );
  end
  if isempty( cmd )
    cmd = 'conda activate env4meshlab & python -c "';
    cmd = [ cmd , 'import pymeshlab; '];
    cmd = [ cmd , 'ms = pymeshlab.MeshSet();' ];
    cmd = [ cmd , 'ms.load_new_mesh( ''M.glb'' );' ];
    cmd = [ cmd , 'ms.apply_texmap_defragmentation( matchingthreshold=2 , boundarytolerance=0.2 , distortiontolerance=0.01 , globaldistortiontolerance=0.01 , uvreductionlimit=0 , offsetfactor=5 , timelimit=0 );' ];
    cmd = [ cmd , 'ms.compute_texcoord_transfer_wedge_to_vertex();' ];
    cmd = [ cmd , 'ms.save_current_mesh( ''M.ply'' , save_vertex_quality=False , save_face_color=False , save_vertex_color=False , save_wedge_texcoord=False );"' ];
  end
  cmd = [ sprintf('cd /d "%s"', dn ) , ' & ' , cmd ];

  [a,b] = system( cmd );
  
  M = read_PLY( fullfile( dn , 'M.ply' ) );
  M = rmfields( M ,'tritexnumber');

%   T = meshSeparate( MeshGenerateIDs(M,'tri'));
%   for a = 1:numel(T)
%     for b = a+1:numel(T)
%       w = distanceFrom( meshFacesCenter( T{a} ) , T{b} ) < 1e-5;
%       T{a} = MeshRemoveFaces( T{a} , w );
%     end
%   end
%   T = MeshAppend( T );

end
