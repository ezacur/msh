%TEST_bvhClosestElement_TETS  bvhClosestElement vs tsearchn on tet meshes.
%
%  tsearchn(X,TES,P) does POINT LOCATION on a Delaunay triangulation: the
%  containing simplex (NaN outside) + its barycentric coordinates.
%  bvhClosestElement SUBSUMES that semantics on any tet mesh:
%      tid = e;  tid( d > 0 ) = NaN;      %tsearchn equivalence
%  with the same barycentric coordinates for interior points -- and it also
%  answers what tsearchn cannot: the NEAREST element / closest point /
%  distance for EXTERIOR points, and it keeps working on NON-DELAUNAY
%  (deformed) tet meshes, where tsearchn is out of contract.

function test_bvhClosestElement_tets
  addpath( fullfile( fileparts( mfilename('fullpath') ) , '..' , 'MESH' ) );
  rng(3);

  %% A) semantics vs tsearchn on a genuine Delaunay tet mesh
  X   = rand( 300 ,3);
  TES = delaunayn( X );
  M   = struct( 'xyz' , X , 'tri' , TES );
  P   = [ rand( 1500 ,3)*0.9 + 0.05 ;          %mostly interior
          rand(  500 ,3)*3 - 1     ];          %mostly clearly exterior

  [ tid , bcT ]     = tsearchn( X , TES , P );
  [ e , cp , d , bc ] = bvhClosestElement( M , P );

  in  = ~isnan( tid );
  out =  isnan( tid );

  %tsearchn-interior  ->  our distance must vanish
  assert( max( d(in) ) < 1e-7 , 'interior point with nonzero distance' );
  %our-interior (d==0) ->  tsearchn must also see it inside
  assert( all( ~isnan( tid( d < 1e-12 ) ) ) , 'we report inside, tsearchn says outside' );
  %same containing tet (ties on shared faces may legitimately differ)
  sameT = in & ( e == tid );
  assert( mean( e(in) == tid(in) ) > 0.99 , 'containing-tet disagreement above tie level' );
  %same barycentric coordinates where the tet agrees
  assert( max(max( abs( bc(sameT,:) - bcT(sameT,:) ) )) < 1e-6 , 'barycentric coords differ' );

  %exterior: we must return a strictly positive distance, a cp ON the winning
  %tet (it reconstructs from bc), and |p-cp| == d
  assert( all( d(out) > 0 ) , 'exterior point with zero distance' );
  To  = TES( e(out) ,:);
  rec = zeros( sum(out) ,3);
  for c = 1:4, rec = rec + bc(out,c) .* X( To(:,c) ,:); end
  assert( max(max( abs( rec - cp(out,:) ) )) < 1e-6 , 'exterior cp does not reconstruct from bc' );
  assert( max( abs( sqrt(sum((P(out,:)-cp(out,:)).^2,2)) - d(out) ) ) < 1e-9 , '|p-cp| ~= d' );

  fprintf( 'semantics    ok  (%d tets; %d interior, %d exterior; same tet %.2f%%)\n' , ...
           size(TES,1) , sum(in) , sum(out) , 100*mean( e(in) == tid(in) ) );

  %% B) NON-Delaunay tets (deformed mesh): tsearchn out of contract, BVH fine
  W  = X + 0.15*sin( 4*X(:,[2 3 1]) );         %smooth warp: same TES, no longer Delaunay
  Mw = struct( 'xyz' , W , 'tri' , TES );
  Pw = rand( 800 ,3)*1.2 - 0.1;

  %ground truth containment: brute-force barycentric over ALL tets
  nT   = size( TES ,1);
  inGT = false( size(Pw,1) ,1);
  for t = 1:nT
    A  = W(TES(t,1),:);
    BA = W(TES(t,2),:) - A;  CA = W(TES(t,3),:) - A;  DA = W(TES(t,4),:) - A;
    d0 = det( [ BA ; CA ; DA ] );
    if d0 == 0, continue; end
    pA = Pw - A;
    l2 = pA * cross( CA , DA ).' / d0;
    l3 = pA * cross( DA , BA ).' / d0;
    l4 = pA * cross( BA , CA ).' / d0;
    l1 = 1 - l2 - l3 - l4;
    inGT( l1>=-1e-12 & l2>=-1e-12 & l3>=-1e-12 & l4>=-1e-12 ) = true;
  end

  [ tidW , ~ ]        = tsearchn( W , TES , Pw );
  [ ~ , ~ , dW , ~ ]  = bvhClosestElement( Mw , Pw );
  ourIn = dW < 1e-12;

  assert( isequal( ourIn , inGT ) , 'BVH containment differs from brute-force ground truth' );
  tsMiss = sum( xor( ~isnan(tidW) , inGT ) );
  if tsMiss == 0, tsTXT = 'agrees'; else, tsTXT = 'OUT OF CONTRACT'; end
  fprintf( 'non-Delaunay ok  (BVH == ground truth on %d pts; tsearchn %s, %d/%d misclassified)\n' , ...
           numel(ourIn) , tsTXT , tsMiss , sum(inGT) );

  %% C) performance (single-thread; ours = MEX backend, cached BVH)
  fprintf( '%8s | %12s | %14s | %14s\n' , 'tets' , 'ours us/pt' , 'tsearchn us/pt' , 'pointLoc us/pt' );
  for n = [ 300 , 1200 , 4000 ]
    X   = rand( n ,3);
    TES = delaunayn( X );
    M   = struct( 'xyz' , X , 'tri' , TES );
    P   = rand( 5000 ,3)*1.2 - 0.1;
    B   = BVH( M );

    tic;  [ e , cp , d , bc ] = bvhClosestElement( {M,B} , P );  tOur = toc;

    ns = min( 5000 , max( 100 , ceil( 2e7 / size(TES,1) ) ) );   %cap tsearchn work
    tic;  tsearchn( X , TES , P(1:ns,:) );  tTS = toc * 5000/ns;

    tPL = NaN;
    try
      DT = delaunayTriangulation( X );
      tic;  pointLocation( DT , P );  tPL = toc;
    end

    fprintf( '%8d | %12.2f | %14.1f | %14.2f   (tsearchn/ours: x%.0f)\n' , ...
             size(TES,1) , 1e6*tOur/5e3 , 1e6*tTS/5e3 , 1e6*tPL/5e3 , tTS/tOur );
  end

  fprintf( 'ALL tet tests passed.\n' );
end
