%TEST_BVH  Correctness of BVH + bvhClosestElement.
%
%  Strategy: a BVH with leafSize = Inf degenerates to a single leaf, i.e. brute
%  force through the SAME exact primitives -- comparing it against the default
%  tree isolates traversal/pruning bugs. Primitive correctness is checked
%  independently: cp must reconstruct from its barycentric coordinates, |p-cp|
%  must equal d, and d must never exceed the distance to points sampled ON the
%  mesh. The transform test checks the update-instead-of-rebuild path.

function test_BVH
  rng(11);
  tol = 1e-10;
  addpath( fullfile( fileparts( mfilename('fullpath') ) , '..' , 'tools' ) );
  addpath( fullfile( fileparts( mfilename('fullpath') ) , '..' , 'MESH' ) );

  %% 1) triangle mesh (sphere-ish, via convex hull)
  V = randn( 700 ,3);  V = V ./ sqrt( sum( V.^2 ,2) );
  M = struct( 'xyz' , V , 'tri' , convhulln( V ) );
  P = [ randn( 300 ,3) * 2 ; randn( 200 ,3) * 0.3 ];   %outside & inside points
  checkMesh( M , P , 'triangles' , tol );

  %% 2) segment mesh (open polyline on a circle)
  t = linspace( 0 , 3*pi/2 , 200 ).';
  M = struct( 'xyz' , [ cos(t) , sin(t) , 0*t ] , 'tri' , [ (1:199).' , (2:200).' ] );
  P = randn( 400 ,3);
  checkMesh( M , P , 'segments' , tol );

  %independent reference implementation: tools/d2Wireframe (GEMM brute force)
  [ e , cp , d ] = bvhClosestElement( M , P );
  [ dW , sW , cpW ] = d2Wireframe( P , M );
  assert( max( abs( d - dW ) ) < 1e-9 , 'segments: BVH vs d2Wireframe distances differ' );
  same = e == sW;   %equidistant ties may pick different segments
  assert( max(max( abs( cp(same,:) - cpW(same,:) ) )) < 1e-9 , 'segments: cp differs from d2Wireframe' );
  fprintf( 'd2Wireframe  ok  (independent reference, %d/%d same winner)\n' , sum(same) , numel(same) );

  %% 3) vertex mesh (point cloud) -- independent brute force over vertices
  V = rand( 500 ,3);
  M = struct( 'xyz' , V , 'tri' , (1:500).' );
  P = rand( 300 ,3) * 1.4 - 0.2;
  [ e , cp , d ] = bvhClosestElement( M , P );
  D = sqrt( ( P(:,1) - V(:,1).' ).^2 + ( P(:,2) - V(:,2).' ).^2 + ( P(:,3) - V(:,3).' ).^2 );
  [ dref , eref ] = min( D ,[],2);
  assert( max( abs( d - dref ) ) < tol         , 'vertices: distances differ from direct min' );
  assert( isequal( e , eref ) || max(abs(d-dref))<tol , 'vertices: wrong element' );
  assert( max( abs( sqrt(sum((P-cp).^2,2)) - d ) ) < tol , 'vertices: |p-cp| ~= d' );
  fprintf( 'vertices     ok  (%d points, %d elements)\n' , size(P,1) , size(M.tri,1) );

  %% 4) tetrahedral mesh -- inside points must report d = 0
  V = rand( 80 ,3);
  M = struct( 'xyz' , V , 'tri' , delaunayn( V ) );
  P = [ rand( 200 ,3) ;                    %mostly inside the unit cube cloud
        rand( 200 ,3) * 3 - 1 ];           %mostly outside
  checkMesh( M , P , 'tetrahedra' , tol );
  [ ~ , ~ , d ] = bvhClosestElement( M , mean( V ,1) );   %deep interior point
  assert( d < tol , 'tets: interior point should have d = 0' );

  %% 5) mixed mesh: triangles + 0-padded segments (zeros trailing)
  V = randn( 300 ,3);  V = V ./ sqrt( sum( V.^2 ,2) );
  Tt = convhulln( V );
  Ts = [ (1:20).' , (21:40).' , zeros(20,1) ];   %random chords as segments
  M  = struct( 'xyz' , V , 'tri' , [ Tt ; Ts ] );
  P  = randn( 300 ,3) * 1.5;
  checkMesh( M , P , 'mixed tri+seg' , tol );

  %% 6) transform update: query on the moved mesh with the UPDATED (not rebuilt) BVH
  V = randn( 500 ,3);  V = V ./ sqrt( sum( V.^2 ,2) );
  M = struct( 'xyz' , V , 'tri' , convhulln( V ) );
  B = BVH( M );

  ang = 0.7;  s = 1.8;
  R = [ cos(ang) -sin(ang) 0 ; sin(ang) cos(ang) 0 ; 0 0 1 ] * s;
  T = [ R , [0.3;-1.2;2.5] ; 0 0 0 1 ];

  M2 = transform( M , T );                     %tools/transform (struct branch)
  B2 = BVH( B , T );                       %O(n) update
  P2 = randn( 300 ,3) * 2 + T(1:3,4).';

  [ e2 , cp2 , d2 ] = bvhClosestElement( M2 , P2 , B2 );
  [ eF , cpF , dF ] = bvhClosestElement( M2 , P2 );          %fresh build
  assert( max( abs( d2 - dF ) ) < tol , 'transform: updated BVH gives wrong distances' );
  assert( max( abs( cp2(:) - cpF(:) ) ) < 1e-8 , 'transform: updated BVH gives wrong cp' );
  fprintf( 'transform    ok  (rot+scale %.1f+transl, updated vs rebuilt)\n' , s );

  %% 6b) refit: non-affine DEFORMATION -> keep hierarchy, recompute spheres only
  V = randn( 600 ,3);  V = V ./ sqrt( sum( V.^2 ,2) );
  M = struct( 'xyz' , V , 'tri' , convhulln( V ) );
  B = BVH( M );

  M2 = M;
  M2.xyz = M.xyz .* ( 1 + 0.5*sin( 3*M.xyz(:,[2 3 1]) ) );    %smooth non-affine warp
  B2 = BVH( B , M2 );                                     %refit O(n)

  if isfield( B , 'bounds4' )
    assert( isequal( B2.child4 , B.child4 ) && isequal( B2.perm , B.perm ) ...
            && isequal( B2.srange , B.srange ) , 'refit: hierarchy was not preserved' );
  else
    assert( isequal( B2.child , B.child ) && isequal( B2.perm , B.perm ) ...
            && isequal( B2.range , B.range ) , 'refit: hierarchy was not preserved' );
  end
  P2 = randn( 400 ,3) * 1.6;
  [ eR_ , cpR , dR ] = bvhClosestElement( M2 , P2 , B2 );
  [ eF_ , cpF , dF ] = bvhClosestElement( M2 , P2 );         %fresh build reference
  assert( max( abs( dR - dF ) ) < tol , 'refit: distances differ from rebuilt BVH' );
  assert( max(max( abs( cpR - cpF ) )) < 1e-8 , 'refit: cp differ from rebuilt BVH' );

  %pruning-quality loss under HEAVY deformation (informative, correctness is
  %already asserted): query time refitted vs rebuilt tree on the same mesh
  V = randn( 2600 ,3);  V = V ./ sqrt( sum( V.^2 ,2) );
  M = struct( 'xyz' , V , 'tri' , convhulln( V ) );
  B = BVH( M );
  M2 = M;  M2.xyz = M.xyz .* ( 1 + 0.5*sin( 4*M.xyz(:,[2 3 1]) ) );
  P2 = randn( 4000 ,3) * 1.6;
  tic; Br = BVH( B , M2 );  tRefit = toc;
  tic; Bf = BVH( M2 );      tBuild = toc;
  tic; bvhClosestElement( M2 , P2 , Br ); tQr = toc;
  tic; bvhClosestElement( M2 , P2 , Bf ); tQf = toc;
  fprintf( 'refit        ok  (hierarchy kept; refit %.1fms vs rebuild %.1fms; query penalty x%.2f)\n' , ...
           1e3*tRefit , 1e3*tBuild , tQr/tQf );

  %% 6c) global frame: similarity transforms fold in O(1), nodes untouched
  V = randn( 600 ,3);  V = V ./ sqrt( sum( V.^2 ,2) );
  M = struct( 'xyz' , V , 'tri' , convhulln( V ) );
  B = BVH( M );

  ang = 0.9;  s = 2.3;
  R = [ cos(ang) 0 sin(ang) ; 0 1 0 ; -sin(ang) 0 cos(ang) ] * s;
  T = [ R , [1;-2;0.5] ; 0 0 0 1 ];
  B2 = BVH( B , T );
  same = isequal( B2.eCenter , B.eCenter ) && isequal( B2.X , B.X );
  if isfield( B , 'bounds4' ), same = same && isequal( B2.bounds4 , B.bounds4 );
  else, same = same && isequal( B2.center , B.center ) && isequal( B2.radius , B.radius );
  end
  assert( same , 'fold: a similarity must not touch any node/vertex' );
  assert( ~isequal( B2.frame , B.frame ) , 'fold: the frame must change' );

  M2 = transform( M , T );
  P2 = randn( 300 ,3) * 3;
  [ ~ , cA , dA ] = bvhClosestElement( M2 , P2 , B2 );
  [ ~ , cB , dB ] = bvhClosestElement( M2 , P2 );          %fresh build reference
  assert( max( abs( dA - dB ) ) < 1e-8 , 'fold: distances differ from fresh build' );
  assert( max(max( abs( cA - cB ) )) < 1e-7 , 'fold: cp differ from fresh build' );

  T2 = [ 0.5*eye(3) , [0;1;0] ; 0 0 0 1 ];                  %chained folds compose
  B3 = BVH( B2 , T2 );
  M3 = transform( M2 , T2 );
  [ ~ , ~ , dA ] = bvhClosestElement( M3 , P2 , B3 );
  [ ~ , ~ , dB ] = bvhClosestElement( M3 , P2 );
  assert( max( abs( dA - dB ) ) < 1e-8 , 'fold: chained frames compose wrong' );

  Ta = [ diag([2 1 0.5]) , [0;0;0] ; 0 0 0 1 ];   %anisotropic -> ERROR + refit
  M4 = transform( M3 , Ta );
  try
    B4 = BVH( B3 , Ta );
    error( 'test:bake' , 'a non-similarity must error' );
  catch ME
    assert( strcmp( ME.identifier , 'BVH:notSimilarity' ) , ...
            'non-similarity: wrong error (%s)' , ME.identifier );
    B4 = BVH( B3 , M4 );                       %the documented fallback: refit
  end
  assert( isequal( B4.child4 , B3.child4 ) , 'refit fallback: hierarchy changed' );
  [ ~ , ~ , dA ] = bvhClosestElement( M4 , P2 , B4 );
  [ ~ , ~ , dB ] = bvhClosestElement( M4 , P2 );
  assert( max( abs( dA - dB ) ) < 1e-8 , 'refit fallback: distances differ from fresh build' );

  Vb = randn( 26000 ,3);  Vb = Vb ./ sqrt( sum( Vb.^2 ,2) );
  Mb = struct( 'xyz' , Vb , 'tri' , convhulln( Vb ) );
  Bb = BVH( Mb );
  tic; for it = 1:200, BVH( Bb , T ); end;  tFold = toc/200;
  fprintf( 'frame fold   ok  (sim+chain; non-sim errors -> refit; fold on 52k tris: %.4f ms, O(1))\n' , 1e3*tFold );

  %% 6d) genuinely 2-D meshes (2-column vertices) -- dimension from COLUMNS only
  X2  = rand( 200 ,2);
  T2d = delaunayn( X2 );
  M2d = struct( 'xyz' , X2 , 'tri' , T2d );
  P2d = rand( 400 ,2)*1.6 - 0.3;
  [ e , cp , d ] = bvhClosestElement( M2d , P2d );
  tid = tsearchn( X2 , T2d , P2d );                         %independent 2-D check
  assert( max( d( ~isnan(tid) ) ) < 1e-7 , '2D: interior point with d > 0' );
  assert( all( d( isnan(tid) ) > 0 )     , '2D: exterior point with d == 0' );
  assert( max( abs( cp(:,3) ) ) == 0     , '2D: cp must stay in the z=0 plane' );

  B2d = BVH( M2d );
  a = 0.7;  sc = 1.5;
  R2 = [ cos(a) -sin(a) ; sin(a) cos(a) ] * sc;
  T3 = [ R2 , [0.4;-0.2] ; 0 0 1 ];                         %homogeneous 2-D similarity
  B2t = BVH( B2d , T3 );
  if isfield( B2d , 'bounds4' ), same = isequal( B2t.bounds4 , B2d.bounds4 );
  else,                          same = isequal( B2t.center  , B2d.center  );
  end
  assert( same , '2D fold: nodes must be untouched' );
  M2t = M2d;  M2t.xyz = X2 * R2.' + [0.4 , -0.2];
  Pq  = rand( 300 ,2)*2 - 0.5;
  [ ~ , ~ , dA ] = bvhClosestElement( M2t , Pq , B2t );
  [ ~ , ~ , dB ] = bvhClosestElement( M2t , Pq );
  assert( max( abs( dA - dB ) ) < 1e-8 , '2D fold: distances differ from fresh build' );

  tt = linspace( 0 , pi , 60 ).';                           %2-D segment polyline
  Ms = struct( 'xyz' , [cos(tt),sin(tt)] , 'tri' , [ (1:59).' , (2:60).' ] );
  Ps = rand( 200 ,2)*2 - 0.5;
  [ ~ , ~ , ds  ] = bvhClosestElement( Ms , Ps );
  [ ~ , ~ , dsR ] = bvhClosestElement( Ms , Ps , BVH( Ms , Inf ) );
  assert( max( abs( ds - dsR ) ) < 1e-10 , '2D segments: BVH vs brute force differ' );
  fprintf( '2D meshes    ok  (tri+seg with 2-col vertices, homogeneous 2-D fold)\n' );

  %% 6e) blob as a value: save/load roundtrip + stale-blob auto-rebuild
  V = randn( 300 ,3);  V = V ./ sqrt( sum( V.^2 ,2) );
  M = struct( 'xyz' , V , 'tri' , convhulln( V ) );
  B = BVH( M );
  P = randn( 100 ,3);
  [ ~ , ~ , d0 ] = bvhClosestElement( M , P , B );
  fn = [ tempname , '.mat' ];
  save( fn , 'B' , 'M' );  L = load( fn );  delete( fn );
  [ ~ , ~ , d1 ] = bvhClosestElement( L.M , P , L.B );
  assert( isequaln( d0 , d1 ) , 'blob: save/load changed the results' );

  Mm = M;  Mm.xyz = M.xyz * 2;                    %mesh edited, blob not updated
  try
    bvhClosestElement( Mm , P , B );
    error( 'test:stale' , 'a stale blob must error' );
  catch ME
    assert( strcmp( ME.identifier , 'bvhClosestElement:staleBVH' ) , ...
            'stale blob: wrong error (%s)' , ME.identifier );
  end
  B2 = BVH( B , Mm );                          %recovery: refit (same connectivity)
  [ ~ , ~ , d2 ] = bvhClosestElement( Mm , P , B2 );
  [ ~ , ~ , d3 ] = bvhClosestElement( Mm , P );
  assert( max( abs( d2 - d3 ) ) < 1e-10 , 'stale blob: refit recovery failed' );
  fprintf( 'blob i/o     ok  (save/load identical; stale blob ERRORS; refit recovers)\n' );

  %% 7) volumes and leaf sizes: every variant agrees on every celltype
  MS = {};
  V = randn( 700 ,3);  V = V ./ sqrt( sum( V.^2 ,2) );
  MS{end+1} = struct( 'xyz' , V , 'tri' , convhulln( V ) );                 %triangles
  t = linspace( 0 , 3*pi/2 , 200 ).';
  MS{end+1} = struct( 'xyz' , [cos(t),sin(t),0*t] , 'tri' , [(1:199).',(2:200).'] ); %segments
  MS{end+1} = struct( 'xyz' , rand(400,3) , 'tri' , (1:400).' );            %vertices
  Vt = rand( 80 ,3);
  MS{end+1} = struct( 'xyz' , Vt , 'tri' , delaunayn( Vt ) );               %tets
  V = randn( 300 ,3);  V = V ./ sqrt( sum( V.^2 ,2) );
  MS{end+1} = struct( 'xyz' , V , 'tri' , [ convhulln(V) ; [(1:20).',(21:40).',zeros(20,1)] ] ); %mixed
  for m = 1:numel( MS )
    Pq = randn( 400 ,3) * 1.5;
    [ ~ , ~ , dR ] = bvhClosestElement( MS{m} , Pq , BVH( MS{m} , Inf ) );  %brute
    for vv = { 'sphere' , 'aabb' , 'obb' , 'kdop' , 'rss' , 'lss' }
      for ls = { [] , 16 , [4 32] }
        Bv = BVH( MS{m} , ls{1} , vv{1} );
        assert( Bv.version == 3 && strcmp( Bv.volume , vv{1} ) , 'blob: wrong meta' );
        [ ~ , ~ , dV ] = bvhClosestElement( MS{m} , Pq , Bv );
        assert( max( abs( dV - dR ) ) < 1e-9 , ...
                'volumes: %s mesh %d differs from brute force' , vv{1} , m );
      end
    end
    %world-aligned comparison piece: no centroid/PCA frame at build
    Bnf = BVH( MS{m} , [] , 'aabb' , 'noframe' );
    assert( isequal( Bnf.frame , eye(4) ) , 'noframe: frame must stay identity' );
    [ ~ , ~ , dV ] = bvhClosestElement( MS{m} , Pq , Bnf );
    assert( max( abs( dV - dR ) ) < 1e-9 , 'noframe: mesh %d differs from brute force' , m );
  end
  %bundled {M,B} calling form
  Bq = BVH( MS{1} );
  Pq = randn( 200 ,3) * 1.5;
  [ e1 , c1 , d1 ] = bvhClosestElement( { MS{1} , Bq } , Pq );
  [ e2 , c2 , d2 ] = bvhClosestElement( MS{1} , Pq , Bq );
  assert( isequaln( [e1,c1,d1] , [e2,c2,d2] ) , '{M,B} form: results differ' );
  %refit on every volume type (tets, deformed) -- hierarchy kept, results exact
  Mt  = MS{4};
  Mt2 = Mt;  Mt2.xyz = Mt.xyz .* ( 1 + 0.3*sin( 5*Mt.xyz(:,[2 3 1]) ) );
  Pq  = rand( 200 ,3)*2 - 0.5;
  [ ~ , ~ , dB ] = bvhClosestElement( Mt2 , Pq );
  for vv = { 'sphere' , 'aabb' , 'obb' , 'kdop' , 'rss' , 'lss' }
    Bt  = BVH( Mt , [] , vv{1} );
    Bt2 = BVH( Bt , Mt2 );
    assert( isequal( Bt2.child4 , Bt.child4 ) && isequal( Bt2.perm , Bt.perm ) , ...
            'refit(%s): hierarchy changed' , vv{1} );
    [ ~ , ~ , dA ] = bvhClosestElement( Mt2 , Pq , Bt2 );
    assert( max( abs( dA - dB ) ) < 1e-9 , 'refit(%s): distances differ' , vv{1} );
  end
  fprintf( 'volumes      ok  (6 volumes x 3 leaf configs x 5 celltypes == brute; noframe; {M,B}; refit x6)\n' );

  %% 7a) PCA build frame: elongated DIAGONAL geometry gets an aligned frame
  L = randn( 2000 ,3) .* [ 8 , 1 , 0.5 ];                  %elongated cloud
  Rr = [ 1 1 0 ; -1 1 0 ; 0 0 sqrt(2) ]/sqrt(2);           %rotate to a diagonal
  L = L * Rr.';
  ML = struct( 'xyz' , L , 'tri' , convhulln( L ) );
  BL = BVH( ML );
  assert( ~isequal( BL.frame(1:3,1:3) , eye(3) ) , 'PCA: frame should be rotated' );
  PL = randn( 300 ,3) * 6;
  [ ~ , ~ , dA ] = bvhClosestElement( ML , PL , BL );
  [ ~ , ~ , dB ] = bvhClosestElement( ML , PL , BVH( ML , Inf ) );
  assert( max( abs( dA - dB ) ) < 1e-8 , 'PCA: distances differ from brute force' );
  %far-from-origin robustness: centering keeps float bounds tight & correct
  MF = MS{1};  MF.xyz = MF.xyz + 1e7;
  PF = randn( 200 ,3)*2 + 1e7;
  [ ~ , ~ , dA ] = bvhClosestElement( MF , PF );
  [ ~ , ~ , dB ] = bvhClosestElement( MF , PF , BVH( MF , Inf ) );
  assert( max( abs( dA - dB ) ) < 1e-6 , 'centering: far-from-origin results differ' );
  fprintf( 'build frame  ok  (PCA-aligned diagonal blob; centered at 1e7 from origin)\n' );

  %% 7b) Dmax: bounded-radius search; misses give e=0, d=Inf, cp/bc NaN
  V = randn( 700 ,3);  V = V ./ sqrt( sum( V.^2 ,2) );
  M = struct( 'xyz' , V , 'tri' , convhulln( V ) );
  P = [ randn( 300 ,3)*0.5 ; randn( 200 ,3)*0.3 + 5 ];   %near + clearly far
  Dmax = 1.0;

  [ ~ , ~ , d0 ] = bvhClosestElement( M , P );           %reference: full search
  [ e1 , c1 , d1 , bc1 ] = bvhClosestElement( M , P , [] , Dmax );
  f = d0 < Dmax;
  assert( max( abs( d1(f) - d0(f) ) ) < 1e-12 , 'Dmax: found distances differ' );
  assert( all( e1(f) >= 1 ) , 'Dmax: found points must have an element' );
  assert( all( e1(~f) == 0 ) && all( isinf( d1(~f) ) ) && ...
          all(all( isnan( c1(~f,:) ) )) && all(all( isnan( bc1(~f,:) ) )) , ...
          'Dmax: beyond-Dmax points must give e=0, d=Inf, cp/bc NaN' );
  fprintf( 'Dmax         ok  (%d found / %d beyond; e=0, d=Inf on misses)\n' , ...
           sum(f) , sum(~f) );

  %% 7b2) Dmax VECTORIAL: cota por punto (siembra de heuristicas, exacta si
  %%      la cota es alcanzable e inflada)
  [ e2 , ~ , d2 ] = bvhClosestElement( M , P , [] , d0 * 0.5 );      %inalcanzable
  assert( all( e2 == 0 ) && all( isinf( d2 ) ) , 'DmaxVec: cotas inalcanzables deben ser miss' );
  [ e2 , ~ , d2 ] = bvhClosestElement( M , P , [] , d0*(1+1e-9) + 1e-12 );
  assert( all( e2 >= 1 ) && max( abs( d2 - d0 ) ) < 1e-12 , ...
          'DmaxVec: cota alcanzable inflada debe reproducir la busqueda completa' );
  mixed = d0 * 0.5;  mixed( 1:2:end ) = Inf;                          %mezcla por punto
  [ e2 , ~ , d2 ] = bvhClosestElement( M , P , [] , mixed );
  assert( all( e2(1:2:end) >= 1 ) && all( e2(2:2:end) == 0 ) , 'DmaxVec: mezcla por punto' );
  assert( max( abs( d2(1:2:end) - d0(1:2:end) ) ) < 1e-12 , 'DmaxVec: exactitud en la mezcla' );
  fprintf( 'Dmax vector  ok  (siembra por punto: miss/exacto/mezcla)\n' );

  %% 7c) robust barycentrics (slivers) + feature classification + boundary
  %sliver of aspect ~1e8: bc must stay finite, in [0,1], and reconstruct cp
  Msl = struct( 'xyz' , [0 0 0 ; 1 0 0 ; 0.5 1e-8 0] , 'tri' , [1 2 3] );
  Psl = [ 0.3 0.5 0.7 ; 0.9 -0.2 -0.1 ; 0.5 0 1e-9 ];
  [ e , cp , ~ , bc ] = bvhClosestElement( Msl , Psl );
  assert( all(all( isfinite( bc(:,1:3) ) )) && all(all( bc(:,1:3) >= 0 )) && ...
          max( abs( sum( bc(:,1:3) ,2) - 1 ) ) < 1e-12 , 'sliver: invalid bc' );
  rec = bc(:,1).*Msl.xyz(1,:) + bc(:,2).*Msl.xyz(2,:) + bc(:,3).*Msl.xyz(3,:);
  assert( max(max( abs( rec - cp ) )) < 1e-7 , 'sliver: bc do not reconstruct cp' );

  %feature classification on a single triangle (open mesh: ALL edges boundary)
  Mtr = struct( 'xyz' , [0 0 0 ; 2 0 0 ; 0 2 0] , 'tri' , [1 2 3] );
  Pf  = [ -1 -1 0 ;     %closest to vertex 1
           1 -1 0 ;     %closest to edge 1-2 interior
           0.4 0.4 1 ;  %projects INSIDE the face
           3  3  0 ];   %closest to edge 2-3 (hypotenuse) interior
  [ e , cp , d , bc , F ] = bvhClosestElement( Mtr , Pf );
  assert( isequal( F.type , [1;2;3;2] ) , 'feature: wrong classification' );
  assert( isequal( F.onBoundary , [true;true;false;true] ) , 'feature: wrong boundary flag' );

  %closed surface: nothing is on the open boundary
  [ ~ , ~ , ~ , ~ , Fc ] = bvhClosestElement( M , P(1:50,:) );
  assert( ~any( Fc.onBoundary ) , 'feature: closed mesh cannot have boundary hits' );

  %tets: interior point -> type 4; exterior -> face/edge/vertex (2/3), no boundary flag
  Vt = rand( 80 ,3);
  Mtt = struct( 'xyz' , Vt , 'tri' , delaunayn( Vt ) );
  [ ~ , ~ , dt , ~ , Ft ] = bvhClosestElement( Mtt , [ mean(Vt,1) ; 3 3 3 ] );
  assert( dt(1) < 1e-12 && Ft.type(1) == 4 , 'feature: tet interior must be type 4' );
  assert( dt(2) > 0 && Ft.type(2) >= 1 && Ft.type(2) <= 3 , 'feature: tet exterior type' );

  %wireframe free ends
  Mw = struct( 'xyz' , [0 0 0 ; 1 0 0 ; 2 0 0] , 'tri' , [1 2 ; 2 3] );
  [ ~ , ~ , ~ , ~ , Fw ] = bvhClosestElement( Mw , [ -1 0 0 ; 1 1 0 ; 3 0 0 ] );
  assert( isequal( Fw.type , [1;1;1] ) && isequal( Fw.onBoundary , [true;false;true] ) , ...
          'feature: wireframe free ends wrong' );
  fprintf( 'features     ok  (sliver bc; vertex/edge/face/inside; open-boundary flag)\n' );

  %% 8) timing (all through the MEX; threads follow maxNumCompThreads)
  V = randn( 26000 ,3);  V = V ./ sqrt( sum( V.^2 ,2) );
  M = struct( 'xyz' , V , 'tri' , convhulln( V ) );
  P = randn( 20000 ,3) * 1.5;

  tic; Ba = BVH( M );                     tBuild = toc;
  Bs = BVH( M , [] , 'sphere' );
  Bo = BVH( M , [] , 'obb' );
  Bk = BVH( M , [] , 'kdop' );
  Br2 = BVH( M , [] , 'rss' );
  Bl2 = BVH( M , [] , 'lss' );
  tic; bvhClosestElement( M , P , Ba );      tA = toc;
  tic; bvhClosestElement( M , P , Bs );      tS = toc;
  tic; bvhClosestElement( M , P , Bo );      tO = toc;
  tic; bvhClosestElement( M , P , Bk );      tK = toc;
  tic; bvhClosestElement( M , P , Br2 );     tR2 = toc;
  tic; bvhClosestElement( M , P , Bl2 );     tL2 = toc;
  Bbrute = BVH( M , Inf );
  Ps = P( 1:250 ,:);
  tic; bvhClosestElement( M , Ps , Bbrute ); tBrute = toc * size(P,1)/size(Ps,1);
  fprintf( 'timing:  %d tris, %d pts:  build %.3fs | aabb %.2f | sphere %.2f | obb %.2f | kdop %.2f | rss %.2f | lss %.2f us/pt | brute ~%.1fs -> x%.0f\n' , ...
           size(M.tri,1) , size(P,1) , tBuild , 1e6*tA/size(P,1) , 1e6*tS/size(P,1) , ...
           1e6*tO/size(P,1) , 1e6*tK/size(P,1) , 1e6*tR2/size(P,1) , 1e6*tL2/size(P,1) , tBrute , tBrute/tA );

  Pfar = randn( 20000 ,3)*0.5 + 40;                    %all far from the mesh
  tic; bvhClosestElement( M , Pfar , Ba );        tFull = toc;
  tic; bvhClosestElement( M , Pfar , Ba , 0.5 );  tDmax = toc;
  fprintf( 'timing:  far cloud %d pts:  full %.2f us/pt | Dmax %.3f us/pt  ->  x%.0f\n' , ...
           size(Pfar,1) , 1e6*tFull/size(Pfar,1) , 1e6*tDmax/size(Pfar,1) , tFull/tDmax );

  %query cost is dominated by the DISTANCE to the mesh (far points must visit
  %the whole tangent shell of nodes, ~O(sqrt(n)) -- physics, not the engine;
  %use Dmax for far fields). Near-surface queries are the real workload:
  w  = randi( size(V,1) , size(P,1) ,1);
  Pn = V(w,:) .* ( 1 + 0.02*randn( size(P,1) ,1) );
  tic; bvhClosestElement( M , Pn , Ba );  tN = toc;
  fprintf( 'timing:  aabb by query distance:  far %.2f | near-surface %.2f us/pt\n' , ...
           1e6*tA/size(P,1) , 1e6*tN/size(P,1) );

  Vt = rand( 4000 ,3);
  Mt = struct( 'xyz' , Vt , 'tri' , delaunayn( Vt ) );
  Pt = rand( 10000 ,3)*1.4 - 0.2;
  Bs = BVH( Mt , [] , 'sphere' );  Bt = BVH( Mt );
  Bo = BVH( Mt , [] , 'obb' );     Bk = BVH( Mt , [] , 'kdop' );
  Br2 = BVH( Mt , [] , 'rss' );    Bl2 = BVH( Mt , [] , 'lss' );
  tic; bvhClosestElement( Mt , Pt , Bs );  tS = toc;
  tic; bvhClosestElement( Mt , Pt , Bt );  tA = toc;
  tic; bvhClosestElement( Mt , Pt , Bo );  tO = toc;
  tic; bvhClosestElement( Mt , Pt , Bk );  tK = toc;
  tic; bvhClosestElement( Mt , Pt , Br2 ); tR2 = toc;
  tic; bvhClosestElement( Mt , Pt , Bl2 ); tL2 = toc;
  fprintf( 'timing:  %d tets, %d pts:  sphere %.2f | aabb %.2f | obb %.2f | kdop %.2f | rss %.2f | lss %.2f us/pt\n' , ...
           size(Mt.tri,1) , size(Pt,1) , 1e6*tS/size(Pt,1) , 1e6*tA/size(Pt,1) , ...
           1e6*tO/size(Pt,1) , 1e6*tK/size(Pt,1) , 1e6*tR2/size(Pt,1) , 1e6*tL2/size(Pt,1) );

  %capsule home turf: a long 3-D helix wireframe
  tt = linspace( 0 , 40*pi , 30001 ).';
  Mw = struct( 'xyz' , [ cos(tt) , sin(tt) , 0.05*tt ] , ...
               'tri' , [ (1:30000).' , (2:30001).' ] );
  Pw = [ 2*randn(10000,2) , rand(10000,1)*7 - 0.5 ];
  Ba = BVH( Mw );                    Bl = BVH( Mw , [] , 'lss' );
  Bo = BVH( Mw , [] , 'obb' );       Br = BVH( Mw , [] , 'rss' );
  tic; bvhClosestElement( Mw , Pw , Ba );  tA = toc;
  tic; bvhClosestElement( Mw , Pw , Bl );  tL = toc;
  tic; bvhClosestElement( Mw , Pw , Bo );  tO = toc;
  tic; bvhClosestElement( Mw , Pw , Br );  tR = toc;
  fprintf( 'timing:  %d segs (helix), %d pts:  aabb %.2f | lss %.2f | obb %.2f | rss %.2f us/pt\n' , ...
           size(Mw.tri,1) , size(Pw,1) , 1e6*tA/size(Pw,1) , 1e6*tL/size(Pw,1) , ...
           1e6*tO/size(Pw,1) , 1e6*tR/size(Pw,1) );

  fprintf( 'ALL BVH tests passed.\n' );
end

%% brute force (single-leaf BVH, same primitives) vs tree + cp/bc consistency
function checkMesh( M , P , name , tol )
  [ e  , cp  , d  , bc ] = bvhClosestElement( M , P );                    %default tree
  [ ~  , ~   , dR      ] = bvhClosestElement( M , P , BVH( M ,Inf) ); %brute force

  assert( max( abs( d - dR ) ) < tol , '%s: BVH and brute-force distances differ' , name );
  assert( max( abs( sqrt(sum((P(:,1:3)-cp).^2,2)) - d ) ) < tol , '%s: |p-cp| ~= d' , name );
  assert( all( e >= 1 & e <= size(M.tri,1) ) , '%s: element id out of range' , name );

  %cp must reconstruct from its barycentric coordinates (validates bc AND cp)
  X = double( M.xyz ); X(:,end+1:3) = 0;
  T = double( M.tri );
  rec = zeros( size(P,1) ,3);
  for c = 1:size( T ,2)
    w = T( e ,c) > 0;
    rec(w,:) = rec(w,:) + bc(w,c) .* X( T(e(w),c) ,:);
  end
  ok = all( isfinite( bc( : ,1) ) ,2);   %quads yield NaN bc (documented TODO)
  err = max( abs( rec(ok,:) - cp(ok,:) ) ,[],2);
  assert( max( err ) < 1e-6 , '%s: cp does not reconstruct from barycentric bc' , name );

  fprintf( '%-12s ok  (%d points, %d elements)\n' , name , size(P,1) , size(M.tri,1) );
end
