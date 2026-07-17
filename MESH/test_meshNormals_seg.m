function test_meshNormals_seg
  warning( 'on' , 'all' );
  PF( 'reset' );

  %------------------------------------------------- T1 2D single segment
  T('T1  2D single segment');
  M = mst( [0 0;1 0] , [1 2] );
  N = meshNormals( M );
  chk( isequal( N , [0 -1] ) , '(1,0) -> [0,-1] exact, 2 columns' );

  %------------------------------------------------- T2 2D CCW square -> outward
  T('T2  2D CCW square');
  M = mst( [0 0;1 0;1 1;0 1] , [1 2;2 3;3 4;4 1] );
  N = meshNormals( M );
  chk( isequal( N , [0 -1;1 0;0 1;-1 0] ) , 'CCW -> outward, exact' );
  M.tri = fliplr( M.tri );
  chk( isequal( meshNormals( M ) , -N ) , 'reversing p1<->p2 flips exactly' );

  %------------------------------------------------- T3 2D zero-length segment
  T('T3  2D zero-length segment');
  M = mst( [0 0;1 0] , [1 1] );
  N = meshNormals( M );
  chk( all( isnan( N(:) ) ) , 'degenerate -> [NaN NaN] (full-NaN row)' );

  %------------------------------------------------- T4 3D single segment (marked)
  T('T4  3D single segment along x');
  M = mst( [0 0 0;1 0 0] , [1 2] );
  N = meshNormals( M );
  chk( isnan( N(3) ) && ~any( isnan( N(1:2) ) ) , 'z is NaN (marked), xy finite' );
  chk( abs( hypot( N(1) , N(2) ) - 1 ) < 1e-15 , 'xy part unit' );
  chk( abs( N(1) ) < 1e-15 , 'perpendicular to +x (N(1)=0)' );

  %------------------------------------------------- T5 3D vertical segment
  T('T5  3D vertical segment');
  M = mst( [0 0 0;0 0 1] , [1 2] );
  N = meshNormals( M );
  chk( isequaln( N , [1 0 NaN] ) , 'vertical -> [1 0 NaN] exact' );

  %------------------------------------------------- T6 3D oblique straight polyline
  T('T6  3D oblique straight polyline (marked)');
  d5 = [1 1 1]/sqrt(3);
  M  = mst( (0:4).' * d5 , [ (1:4).' , (2:5).' ] );
  N  = meshNormals( M );
  chk( all( isnan( N(:,3) ) ) && ~any( isnan( N(:,1:2) ) ,'all') , 'all rows marked' );
  chk( size( unique( N(:,1:2) ,'rows') ,1) == 1 , 'same normal whole piece' );
  chk( abs( [ N(1,1:2) , 0 ] * d5.' ) < 1e-15 , '[nx,ny,0] truly perp to line dir' );

  %------------------------------------------------- T7 3D coplanar circle in XY == 2D branch
  T('T7  3D circle in z=0 plane == 2D branch');
  t  = linspace( 0 , 2*pi , 41 ).';  t(end) = [];
  P2 = [ cos(t) , sin(t) ];
  tr = [ (1:40).' , [2:40,1].' ];
  M3 = mst( [ P2 , zeros(40,1) ] , tr );
  M2 = mst( P2 , tr );
  N3 = meshNormals( M3 );  N2 = meshNormals( M2 );
  chk( max(abs( N3(:,1:2) - N2 ),[],'all') < 1e-12 , 'xy equal to 2D branch' );
  chk( max(abs( N3(:,3) )) < 1e-12 , 'z component ~0' );
  mid = ( M3.xyz( tr(:,1) ,:) + M3.xyz( tr(:,2) ,:) )/2;
  chk( all( sum( N3 .* mid ,2) > 0 ) , 'CCW -> outward' );

  %------------------------------------------------- T8 tilted coplanar circle
  T('T8  tilted coplanar circle');
  ax = [1 2 3]; ax = ax/norm(ax);         %rotation, det +1, keeps normal z>0
  R0 = expm( skew( ax * 0.8 ) );
  M8 = mst( M3.xyz * R0.' + [5 -3 2] , tr );
  N8 = meshNormals( M8 );
  chk( max(abs( N8 - N3 * R0.' ),[],'all') < 1e-9 , 'N == R0 * (flat N)  (frame covariant)' );
  sd = M8.xyz( tr(:,2) ,:) - M8.xyz( tr(:,1) ,:);
  chk( max(abs( sum( N8.*sd ,2) ./ sqrt(sum(sd.^2,2)) )) < 1e-12 , 'perp to every segment' );
  pn = R0(:,3);
  chk( max(abs( N8 * pn )) < 1e-12 , 'in-plane (perp to plane normal)' );

  %--- same circle, rotation that flips the plane normal below the horizon
  R1 = expm( skew( [1 0 0] * pi*0.9 ) );          %normal ends with z<0
  M8b = mst( M3.xyz * R1.' , tr );
  N8b = meshNormals( M8b );
  gflip = max(abs( N8b - N3*R1.' ),[],'all') < 1e-9;
  gkeep = max(abs( N8b + N3*R1.' ),[],'all') < 1e-9;
  info( sprintf( 'normal-below-horizon rotation: covariant=%d globally-flipped=%d (''+z'' convention)' , gflip , gkeep ) );

  %------------------------------------------------- T9 helix (genuinely 3D)
  T('T9  helix vs analytic curvature direction');
  t  = ( 0 : 0.15 : 4*pi ).';
  H  = [ cos(t) , sin(t) , 0.2*t ];
  trh = [ (1:numel(t)-1).' , (2:numel(t)).' ];
  MH = mst( H , trh );
  lastwarn('');
  NH = meshNormals( MH );
  [~,wid] = lastwarn();
  chk( isempty( wid ) , 'no warnings on a clean helix' );
  chk( ~any( isnan( NH(:) ) ) , 'no NaN' );
  tm = ( t(1:end-1) + t(2:end) )/2;
  NA = [ cos(tm) , sin(tm) , zeros(size(tm)) ];     %+- curvature line (sign = component-global)
  dp = sum( NH .* NA ,2);
  info( sprintf( 'dot(N,radial): min %.4f  mean %.4f  max %.4f' , min(dp) , mean(dp) , max(dp) ) );
  chk( nnz( diff( sign(dp) ) ) == 0 , 'consistent sign along the curve (no flips)' );
  chk( min( abs(dp) ) > 0.95 , 'within ~18 deg of the curvature line' );
  sd = H( trh(:,2) ,:) - H( trh(:,1) ,:);  sd = sd ./ sqrt(sum(sd.^2,2));
  info( sprintf( 'max |dot(N,segdir)| = %.2e (projection approx.)' , max(abs(sum(NH.*sd,2))) ) );

  %--- exact flip under reversal
  MHr = MH;  MHr.tri = fliplr( MHr.tri );
  chk( isequaln( meshNormals( MHr ) , -NH ) , 'reversal flips exactly' );

  %------------------------------------------------- T10 straight run inside a 3D component
  T('T10 straight run inside a genuinely-3D component');
  X = [ 0 0 0 ; 1 0 0 ; 2 0 0 ; 3 0 0 ; 3.7 .7 .7 ; 4 1.5 .2 ];
  M10 = mst( X , [ (1:5).' , (2:6).' ] );
  lastwarn('');
  N10 = meshNormals( M10 );
  [~,wid] = lastwarn();
  chk( isempty( wid ) , 'no warnings (straight run inherits its neighbours'' plane)' );
  chk( ~any( isnan( N10(:) ) ) , 'no NaN anywhere' );
  chk( max(abs( sqrt(sum(N10.^2,2)) - 1 )) < 1e-12 , 'all unit' );
  ca10 = sum( N10(1:end-1,:) .* N10(2:end,:) ,2);
  chk( all( ca10 > 0 ) , 'continuous along the polyline (run follows the bend, no flips)' );

  %------------------------------------------------- T11 near-collinear (wiggle 1e-7 > TH)
  T('T11 non-coplanar wiggle 1e-7 (just above TH=1e-8)');
  rng(1);
  s  = (0:9).';
  X  = [ s , zeros(10,1) , zeros(10,1) ] + 1e-7*randn(10,3);
  M11 = mst( X , [ (1:9).' , (2:10).' ] );
  lastwarn('');
  N11 = meshNormals( M11 );
  [~,wid] = lastwarn();
  info( sprintf( 'lastwarn id: %s' , wid ) );
  info( sprintf( 'marked rows: %d/9 ; normals(1:3,:) =' , sum(isnan(N11(:,3))) ) );
  disp( N11(1:3,:) );
  %--- same but wiggle 1e-9 (below TH): must be marked-collinear
  X  = [ s , zeros(10,1) , zeros(10,1) ] + 1e-9*randn(10,3);
  N11b = meshNormals( mst( X , [ (1:9).' , (2:10).' ] ) );
  chk( all( isnan( N11b(:,3) ) ) , 'wiggle 1e-9 < TH -> whole piece marked collinear' );

  %------------------------------------------------- T12 zero-length segment inside a circle
  T('T12 zero-length segment inside a coplanar circle');
  tr12 = [ tr ; 5 5 ];
  M12 = mst( M3.xyz , tr12 );
  N12 = meshNormals( M12 );
  chk( all( isnan( N12(end,:) ) ) , 'degenerate row -> full NaN' );
  chk( max(abs( N12(1:40,:) - N3 ),[],'all') < 1e-12 , 'rest of the circle untouched' );

  %------------------------------------------------- T13 noise sweep on the circle
  T('T13 z-noise sweep on the circle (branch selection)');
  rng(2);  nz = randn(40,1);
  for lvl = [ 1e-9 , 1e-6 , 1e-3 ]
    Mz = mst( [ P2 , lvl*nz ] , tr );
    Nz = meshNormals( Mz );
    ang = acosd( min( 1 , abs( sum( Nz(:,1:2).*N2 ,2) ) ./ sqrt(sum(Nz(:,1:2).^2,2)) ) );
    info( sprintf( 'noise %g: marked=%d  maxNaN=%d  max in-plane deviation %.3g deg' , ...
          lvl , sum(isnan(Nz(:,3))) , sum(any(isnan(Nz),2)) , max(ang) ) );
  end

  %------------------------------------------------- T14 mixed components + scatter-back
  T('T14 mixed mesh (straight + circle + helix), shuffled rows');
  nA = size( M.tri ,1);
  XYZ = [ (0:4).'*d5 + [10 0 0] ; M3.xyz ; H ];
  TRI = [ [ (1:4).' , (2:5).' ] ; tr+5 ; trh+45 ];
  rng(3);  pp = randperm( size( TRI ,1) ).';
  Mmix = mst( XYZ , TRI(pp,:) );
  Nmix = meshNormals( Mmix );
  Nref = [ meshNormals( mst( (0:4).'*d5+[10 0 0] , [ (1:4).' , (2:5).' ] ) ) ; N3 ; NH ];
  chk( isequaln( Nmix , Nref(pp,:) ) , 'scatter-back == per-component run, row-exact' );

  %------------------------------------------------- T15 vertex modes on segments
  T('T15 vertex modes on segment meshes');
  Nv = meshNormals( M3 , 'uniform' );
  chk( max(abs( sqrt(sum(Nv.^2,2)) - 1 )) < 1e-15 , '''uniform'' on circle: unit vertex normals' );
  chk( all( sum( Nv .* M3.xyz ,2) > 0.99 ) , '''uniform'' on circle: outward radial' );
  chk( isequaln( meshNormals( M3 ,'angle') , Nv ) , '''angle'' == ''uniform'' for segments (pseudonormal = bisector)' );
  Na = meshNormals( M3 , 'area' );
  chk( max(abs( sqrt(sum(Na.^2,2)) - 1 )) < 1e-15 && all( sum( Na.*M3.xyz ,2) > 0.99 ) , ...
       '''area'' -> length-weighted: works, unit + outward' );
  chk( isequaln( meshNormals( M3 ,'best') , Nv ) , '''best'' == ''uniform'' on a clean circle (valence 2)' );
  %marked component -> vertex normals stay marked with UNIT xy (contract restored)
  Vv = meshNormals( mst( (0:4).'*d5 , [ (1:4).' , (2:5).' ] ) , 'uniform' );
  chk( all( isnan( Vv(:,3) ) ) && max(abs( sqrt(Vv(:,1).^2+Vv(:,2).^2) - 1 )) < 1e-15 , ...
       'marked piece ''uniform'': [n1 n2 NaN] with UNIT xy at every vertex' );

  %------------------------------------------------- T16 smoothing on segments
  T('T16 numeric smoothing on segments');
  Ns = meshNormals( M3 , 3 );
  chk( max(abs( sqrt(sum(Ns.^2,2)) - 1 )) < 1e-15 , 'k=3 on circle: unit' );
  chk( all( sum( Ns .* mid ,2) > 0 ) , 'k=3 on circle: still outward' );
  Nsm = meshNormals( Mmix , 2 );
  w = isnan( Nref(pp,3) );
  chk( isequaln( Nsm(w,:) , Nref(pp(w),:) ) , 'marked rows survive smoothing untouched' );
  chk( ~any( isnan( Nsm(~w,3) ) ) , 'no NaN leaked into finite rows' );

  %------------------------------------------------- T17 triangles quick regression
  T('T17 triangles (celltype 5) quick regression');
  Mt = mst( [0 0 0;1 0 0;0 1 0] , [1 2 3] );
  chk( isequal( meshNormals( Mt ) , [0 0 1] ) , 'CCW triangle -> +z' );
  M2t = mst( [0 0;1 0;0 1] , [1 2 3] );
  chk( isequal( meshNormals( M2t ) , [0 0 1] ) , '2D triangle -> z-padded, +z, 3 cols' );
  Tt = [ 0 0 0 ; 1 0 0 ; 0 1 0 ; 0 0 1 ];
  Ft = [ 1 3 2 ; 1 2 4 ; 2 3 4 ; 1 4 3 ];    %outward wound tetra
  Nt = meshNormals( mst( Tt , Ft ) );
  ct = ( Tt(Ft(:,1),:) + Tt(Ft(:,2),:) + Tt(Ft(:,3),:) )/3 - mean(Tt,1);
  chk( all( sum( Nt.*ct ,2) > 0 ) , 'tetra: all faces outward' );

  %------------------------------------------------- T18 mode validation
  T('T18 mode validation');
  ok = 0;
  for bad = { -1 , 2.5 , Inf , NaN , [1 2] }
    try, meshNormals( M3 , bad{1} ); catch e, ok = ok + strcmp( e.identifier , 'meshNormals:mode' ); end
  end
  chk( ok == 5 , 'bad numeric modes (-1, 2.5, Inf, NaN, [1 2]) all error' );
  chk( isequaln( meshNormals( M3 , 0 ) , meshNormals( M3 , false ) ) , '0 == false == face normals' );
  chk( isequaln( meshNormals( M3 , true ) , meshNormals( M3 , 'uniform' ) ) , 'true == ''uniform''' );
  try, meshNormals( M3 , 'bogus' ); catch e, info( [ 'bad char mode -> "' , e.message , '" (id: ' , e.identifier , ')' ] ); end

  %------------------------------------------------- T19 triNORMALS shortcut
  T('T19 triNORMALS shortcut in vertex modes');
  Mtn = M3;  Mtn.triNORMALS = repmat( [1 0 0] , 40 , 1 );
  chk( isequal( meshNormals( Mtn , 'uniform' ) , repmat( [1 0 0] , 40 , 1 ) ) , 'triNORMALS honored' );

  %------------------------------------------------- T20 timing
  T('T20 timing (osculating loop)');
  t  = ( 0 : 0.005 : 4*pi ).';
  MB = mst( [ cos(t) , sin(t) , 0.2*t ] , [ (1:numel(t)-1).' , (2:numel(t)).' ] );
  tic;  meshNormals( MB );  tt = toc;
  info( sprintf( '%d segments (one 3D component): %.2f s  (%.2f ms/segment)' , size(MB.tri,1) , tt , 1000*tt/size(MB.tri,1) ) );
  t2 = ( linspace( 0 , 2*pi , numel(t) ) ).';
  MB2 = mst( [ cos(t2) , sin(t2) , 0*t2 ] , [ (1:numel(t2)-1).' , (2:numel(t2)).' ] );
  tic;  meshNormals( MB2 );  tt2 = toc;
  info( sprintf( '%d segments (coplanar comp.): %.3f s' , size(MB2.tri,1) , tt2 ) );

  pf = PF( 'get' );
  fprintf( '\n==== %d PASS, %d FAIL ====\n' , pf(1) , pf(2) );
  if pf(2), error('FAILURES'); end
end

%-------------------------------------------------------------- helpers
function M = mst( xyz , tri ), M = struct( 'xyz',xyz , 'tri',tri ); end
function s = skew( v ), s = [ 0 -v(3) v(2) ; v(3) 0 -v(1) ; -v(2) v(1) 0 ]; end
function T( s ), fprintf( '\n--- %s\n' , s ); end
function info( s ), fprintf( '     INFO  %s\n' , s ); end
function chk( c , s )
  if c, fprintf( '     PASS  %s\n' , s ); PF(1,1);
  else, fprintf( '  ** FAIL  %s\n' , s ); PF(2,1);
  end
end
function chk2( s , c ), chk( c , s ); end
function out = PF( a , ~ )
  persistent p
  if isempty( p ), p = [0 0]; end
  if ischar( a )
    if strcmp( a , 'reset' ), p = [0 0]; end
  elseif nargin == 2
    p(a) = p(a) + 1;
  end
  out = p;
end
