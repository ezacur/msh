function test_meshNormals
%TEST_MESHNORMALS  Full regression suite for meshNormals.m (and for meshF2V /
% meshSeparate / meshEsuE / meshEsuP as exercised through it): 185 checks
% born out of the 2026-07 review + fix sessions. Any failure raises an error.
%
%   test_meshNormals_seg        exhaustive celltype-3 behavior: 2D, coplanar,
%                               marks, scatter-back, smoothing, triangles
%                               regression, mode validation (47)
%   test_meshNormals_fix        binormal-sign propagation (F1) + straight runs
%                               inside 3D components (F2), against the FROZEN
%                               pre-fix reference in old_meshNormals.mat (21)
%   test_meshNormals_F3         vertex-mode analogs for segments + the new
%                               'reciprocal' / 'quadratic' modes (31)
%   test_meshNormals_recip5     Max-weight 'reciprocal' on triangles: exact on
%                               sphere-inscribed meshes (8)
%   test_meshNormals_workflow   the smooth-then-aggregate recipe, k-mode
%                               CONTINUATION from M.triNORMALS (bit-exact),
%                               direction-field diffusion (21)
%   test_meshNormals_empty      empty meshes x every mode (53)
%   plus, inline here: the stale-triNORMALS and celltype guards.
%
% See also meshNormals.

  HERE = fileparts( mfilename( 'fullpath' ) );
  if isempty( which( 'meshNormals'          ) ), addpath( fileparts( HERE ) ); end
  if isempty( which( 'test_meshNormals_seg' ) ), addpath( HERE );             end

  test_meshNormals_seg;
  test_meshNormals_fix;
  test_meshNormals_F3;
  test_meshNormals_recip5;
  test_meshNormals_workflow;
  test_meshNormals_empty;
  guards;

  fprintf( '\n######## meshNormals: ALL SUITES GREEN ########\n' );
end


function guards
  fprintf( '\n--- guards: stale triNORMALS + unsupported celltype\n' );
  np = 0; nf = 0;

  M = struct( 'xyz',[0 0 0;1 0 0;0 1 0;1 1 0] , 'tri',[1 2 3;2 4 3] , ...
              'triNORMALS', zeros(5,3) );                %WRONG row count (5 vs 2)
  ok = false;
  try, meshNormals( M , 'uniform' ); catch e, ok = strcmp( e.identifier , 'meshNormals:triNORMALS' ); end
  [np,nf] = chk( ok , 'stale triNORMALS -> meshNormals:triNORMALS (vertex mode)' , np,nf );
  ok = false;
  try, meshNormals( M , 3 ); catch e, ok = strcmp( e.identifier , 'meshNormals:triNORMALS' ); end
  [np,nf] = chk( ok , 'stale triNORMALS -> meshNormals:triNORMALS (smoothing mode)' , np,nf );
  M.triNORMALS = repmat( [0 0 1] , 2 , 1 );              %right size -> accepted
  ok = isequal( meshNormals( M , 'uniform' ) , repmat( [0 0 1] , 4 , 1 ) );
  [np,nf] = chk( ok , 'right-sized triNORMALS still honored' , np,nf );

  Mt = struct( 'xyz',[0 0 0;1 0 0;0 1 0;0 0 1] , 'tri',[1 2 3 4] );   %a tet
  ok = false;
  try, meshNormals( Mt ); catch e, ok = strcmp( e.identifier , 'meshNormals:celltype' ); end
  [np,nf] = chk( ok , 'celltype 10 -> meshNormals:celltype' , np,nf );

  fprintf( '\n==== %d PASS, %d FAIL ====\n' , np , nf );
  if nf, error( 'FAILURES' ); end
end

function [np,nf] = chk( c , s , np , nf )
  if c, fprintf( '     PASS  %s\n' , s ); np = np+1;
  else, fprintf( '  ** FAIL  %s\n' , s ); nf = nf+1;
  end
end
