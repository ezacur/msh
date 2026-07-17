%TEST_bvhIntersectRay  bvhIntersectRay (unified blob) vs tools/IntersectSurfaceRay.
%
%  IntersectSurfaceRay_mx (BVH4 + Moller-Trumbore, its own internal LRU cache)
%  is the ORACLE: same ray conventions, same t parametrization, same inclusive
%  tolerance, same 'any' window -- results must agree. The point of the new
%  kernel is ARCHITECTURE: it runs on the SAME cached/transformable blob as
%  bvhClosestElement (one structure per mesh for both query types).

function test_bvhIntersectRay
  rng(5);
  addpath( fullfile( fileparts( mfilename('fullpath') ) , '..' , 'tools' ) );
  addpath( fullfile( fileparts( mfilename('fullpath') ) , '..' , 'MESH' ) );

  V = randn( 2000 ,3);  V = V ./ sqrt( sum( V.^2 ,2) );
  M = struct( 'xyz' , V , 'tri' , convhulln( V ) );
  N = 600;
  rays = [ randn(N,3)*3 , randn(N,3)*0.3 ];     %origins mostly out, targets mostly in
  B = BVH( M );

  %% modes vs oracle
  for MODE = { 'first' , 'last' , 'any' }
    [ x1 , c1 , t1 ] = bvhIntersectRay( M , rays , B , MODE{1} );
    [ x2 , ~ , c2 , t2 ] = IntersectSurfaceRay( M , rays , MODE{1} );
    hit1 = c1 > 0;  hit2 = c2 > 0;
    assert( isequal( hit1 , hit2 ) , '%s: hit/miss disagreement' , MODE{1} );
    if strcmp( MODE{1} , 'any' )
      %'any' reports SOME window hit: t may differ; verify OUR hits are valid
      assert( all( t1(hit1) > 1e-9 & t1(hit1) < 1-1e-5 ) , 'any: t out of window' );
    else
      assert( max( abs( t1(hit1) - t2(hit2) ) ) < 1e-9 , '%s: t differ' , MODE{1} );
      assert( max(max( abs( x1(hit1,:) - x2(hit2,:) ) )) < 1e-8 , '%s: xyz differ' , MODE{1} );
      same = c1(hit1) == c2(hit2);          %shared-edge hits may pick either face
      assert( mean( same ) > 0.99 , '%s: cell disagreement above tie level' , MODE{1} );
    end
    fprintf( 'ray %-5s    ok  (%d/%d rays hit)\n' , MODE{1} , sum(hit1) , N );
  end

  %'all': per-ray hit multiset must match
  [ ~ , c1 , t1 , r1 ] = bvhIntersectRay( M , rays , B , 'all' );
  [ ~ , ~ , c2 , t2 , r2 ] = IntersectSurfaceRay( M , rays , 'all' );
  assert( numel( t1 ) == numel( t2 ) , 'all: different hit counts' );
  for q = 1:N
    a = sort( t1( r1 == q ) );  b = sort( t2( r2 == q ) );
    assert( numel(a) == numel(b) && ( isempty(a) || max(abs(a-b)) < 1e-9 ) , ...
            'all: ray %d hit set differs' , q );
  end
  fprintf( 'ray all      ok  (%d total hits match)\n' , numel(t1) );

  %% frame: fold a similarity, query the transformed mesh with the SAME blob
  ang = 0.6;  s = 1.7;
  R = [ cos(ang) -sin(ang) 0 ; sin(ang) cos(ang) 0 ; 0 0 1 ] * s;
  T = [ R , [2;-1;0.5] ; 0 0 0 1 ];
  B2 = BVH( B , T );                                  %O(1) fold
  M2 = transform( M , T );
  rays2 = [ ( rays(:,1:3) * R.' + [2 -1 0.5] ) , ( rays(:,4:6) * R.' + [2 -1 0.5] ) ];
  [ x1 , c1 , t1 ] = bvhIntersectRay( M2 , rays2 , B2 , 'first' );
  [ x3 , c3 , t3 ] = bvhIntersectRay( M2 , rays2 , [] , 'first' );   %fresh build
  w = c1 > 0;
  assert( isequal( c1>0 , c3>0 ) && max( abs( t1(w) - t3(w) ) ) < 1e-9 , 'frame: t differ' );
  assert( max(max( abs( x1(w,:) - x3(w,:) ) )) < 1e-7 , 'frame: xyz differ' );
  fprintf( 'ray frame    ok  (folded blob == fresh build on transformed mesh)\n' );

  %% every volume type agrees on rays
  V = randn( 2000 ,3);  V = V ./ sqrt( sum( V.^2 ,2) );
  M = struct( 'xyz' , V , 'tri' , convhulln( V ) );
  B = BVH( M );
  [ ~ , cR , tR ] = bvhIntersectRay( M , rays , B , 'first' );
  for vv = { 'sphere' , 'obb' , 'kdop' , 'rss' , 'lss' }
    Bv = BVH( M , [] , vv{1} );
    [ ~ , cV , tV ] = bvhIntersectRay( M , rays , Bv , 'first' );
    w = cR > 0;
    assert( isequal( cV > 0 , w ) && max( abs( tV(w) - tR(w) ) ) < 1e-9 , ...
            'ray volumes: %s differs' , vv{1} );
  end
  fprintf( 'ray volumes  ok  (sphere/obb/kdop/rss/lss == aabb)\n' );

  %% mixed mesh: non-triangle cells are skipped
  Mx = struct( 'xyz' , V , 'tri' , [ M.tri ; [ (1:10).' , (11:20).' , zeros(10,1) ] ] );
  [ ~ , cM , tM ] = bvhIntersectRay( Mx , rays , [] , 'first' );
  [ ~ , cT , tT ] = bvhIntersectRay( M  , rays , B  , 'first' );
  wM = cM > 0;
  assert( isequal( wM , cT > 0 ) && max( abs( tM(wM) - tT(wM) ) ) < 1e-12 , ...
          'mixed: segments must not affect ray hits' );
  fprintf( 'ray mixed    ok  (segments skipped)\n' );

  %% timing vs IntersectSurfaceRay_mx (its LRU warmed by a first call)
  V = randn( 26000 ,3);  V = V ./ sqrt( sum( V.^2 ,2) );
  M = struct( 'xyz' , V , 'tri' , convhulln( V ) );
  Nb = 20000;
  rays = [ randn(Nb,3)*3 , randn(Nb,3)*0.3 ];
  B = BVH( M );
  IntersectSurfaceRay( M , rays(1:100,:) , 'first' );     %warm the oracle cache
  tic;  bvhIntersectRay( M , rays , B , 'first' );          tOur1 = toc;
  tic;  IntersectSurfaceRay( M , rays , 'first' );           tRef1 = toc;
  tic;  bvhIntersectRay( M , rays , B , 'any' );            tOur2 = toc;
  tic;  IntersectSurfaceRay( M , rays , 'any' );             tRef2 = toc;
  fprintf( 'timing:  %d tris, %d rays:  first %.2f vs %.2f us/ray | any %.2f vs %.2f us/ray  (ours vs IntersectSurfaceRay)\n' , ...
           size(M.tri,1) , Nb , 1e6*tOur1/Nb , 1e6*tRef1/Nb , 1e6*tOur2/Nb , 1e6*tRef2/Nb );

  fprintf( 'ALL bvhIntersectRay tests passed.\n' );
end
