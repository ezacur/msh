function test_meshNormals_empty
  warning( 'on' , 'all' );
  np = 0; nf = 0;

  CASES = { ...
    'seg 3D, 0 verts'  , struct( 'xyz',zeros(0,3) , 'tri',zeros(0,2) ) ;
    'seg 2D, 0 verts'  , struct( 'xyz',zeros(0,2) , 'tri',zeros(0,2) ) ;
    'tri 3D, 0 verts'  , struct( 'xyz',zeros(0,3) , 'tri',zeros(0,3) ) ;
    'seg 3D, 5 verts'  , struct( 'xyz',randn(5,3) , 'tri',zeros(0,2) ) ;
    'seg 2D, 5 verts'  , struct( 'xyz',randn(5,2) , 'tri',zeros(0,2) ) ;
    'tri 3D, 5 verts'  , struct( 'xyz',randn(5,3) , 'tri',zeros(0,3) ) };

  MODES = { false , 3 , true , 'uniform' , 'angle' , 'area' , 'best' , 'reciprocal' , 'quadratic' };

  for ic = 1:size( CASES ,1)
    M  = CASES{ic,2};
    nV = size( M.xyz ,1);
    c  = 3;  if size( M.tri ,2) == 2 && size( M.xyz ,2) < 3, c = 2; end
    fprintf( '\n--- %s (expect faces 0x%d, vertices %dx%d all-NaN)\n' , CASES{ic,1} , c , nV , c );
    for im = 1:numel( MODES )
      md = MODES{im};
      if ischar( md ), lbl = md; else, lbl = mat2str( md ); end
      isVERTEX = ischar( md ) || ( islogical( md ) && md );
      if isequal( md , 'quadratic' ) && size( M.tri ,2) ~= 2, continue; end  %segments-only
      try
        N = meshNormals( M , md );
        if isVERTEX
          ok = isequal( size(N) , [nV,c] ) && all( isnan( N(:) ) );
          msg = sprintf( '%-12s -> %dx%d all-NaN' , lbl , size(N) );
        else
          ok = isequal( size(N) , [0,c] );
          msg = sprintf( '%-12s -> %dx%d empty' , lbl , size(N) );
        end
        [np,nf] = chk( ok , msg , np,nf );
      catch e
        [np,nf] = chk( false , sprintf( '%-12s -> ERROR "%s" (%s)' , lbl , e.message , e.identifier ) , np,nf );
      end
    end
  end

  %empty mesh + bad mode must STILL error (strictness intact)
  ok = false;
  try, meshNormals( CASES{1,2} , -1 ); catch e, ok = strcmp( e.identifier , 'meshNormals:mode' ); end
  [np,nf] = chk( ok , 'empty mesh + bad numeric mode still errors' , np,nf );

  fprintf( '\n==== %d PASS, %d FAIL ====\n' , np , nf );
  if nf, error('FAILURES'); end
end

function [np,nf] = chk( c , s , np , nf )
  if c, fprintf( '     PASS  %s\n' , s ); np = np+1;
  else, fprintf( '  ** FAIL  %s\n' , s ); nf = nf+1;
  end
end
