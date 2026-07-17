function test_MeshSubdivide_schemes
%TEST_MESHSUBDIVIDE_SCHEMES  The 2026-07-11 segment schemes: 'cornercutting'
% (Chaikin), '4points' (Dyn-Levin-Gregory) and 'pn<k>' (1->k PN split).
% Run standalone or via the suite; any failure raises an error.
  HERE = fileparts( mfilename( 'fullpath' ) );
  if isempty( which( 'MeshSubdivide' ) ), addpath( fileparts( HERE ) ); end
  warning( 'on' , 'all' );
  np = 0; nf = 0;

  n  = 32;
  th = ( 0:n-1 ).' * 2*pi/n;
  C  = mkcircle( th );
  C.triID = ( 1:n ).';
  a  = pi/n;

  %================================ CHAIKIN
  T('C1  chaikin on the closed circle');
  Mc = MeshSubdivide( C , 'cornercutting' );
  [np,nf] = chk( size(Mc.xyz,1) == 2*n && size(Mc.tri,1) == 2*n , '2n nodes, 2n segments (originals cut away)' , np,nf );
  QA = 0.75*C.xyz( C.tri(:,1) ,:) + 0.25*C.xyz( C.tri(:,2) ,:);
  QB = 0.25*C.xyz( C.tri(:,1) ,:) + 0.75*C.xyz( C.tri(:,2) ,:);
  [np,nf] = chk( isequal( sortrows( Mc.xyz ) , sortrows( [QA;QB] ) ) , 'cut points bit-exact at 1/4 and 3/4' , np,nf );
  val = accumarray( Mc.tri(:) , 1 );
  [np,nf] = chk( all( val == 2 ) , 'still a closed loop (all valences 2)' , np,nf );
  cnt = accumarray( Mc.triID , 1 , [n,1] );
  [np,nf] = chk( all( cnt == 2 ) , 'tri* inheritance: every parent -> its middle + its corner' , np,nf );
  M4 = MeshSubdivide( { C , 4 } , 'chaikin' );
  info( sprintf( '4 passes: radius in [%.4f, %.4f] (B-spline limit shrink)' , ...
        min(sqrt(sum(M4.xyz.^2,2))) , max(sqrt(sum(M4.xyz.^2,2))) ) );

  T('C2  chaikin on an open arc / a straight line / a star');
  na = 10;
  A  = mkchain( [ cos(linspace(0,pi/2,na)).' , sin(linspace(0,pi/2,na)).' , zeros(na,1) ] );
  Ma = MeshSubdivide( A , 'cornercutting' );
  [np,nf] = chk( size(Ma.tri,1) == 2*(na-1)+1 && size(Ma.xyz,1) == 2*(na-1)+2 , 'open chain: 2n+1 segments, 2n+2 nodes' , np,nf );
  [np,nf] = chk( ismember( A.xyz(1,:) , Ma.xyz ,'rows') && ismember( A.xyz(end,:) , Ma.xyz ,'rows') , ...
                 'endpoints KEPT exactly' , np,nf );
  [np,nf] = chk( ~any( ismember( A.xyz(2:end-1,:) , Ma.xyz ,'rows') ) , 'interior originals removed' , np,nf );
  %straight line + a linear per-node field: both stay exactly affine
  L = mkchain( (0:9).' * [1 0.5 0.25] );
  L.xyzF = 3*L.xyz(:,1) + 1;
  Ml = MeshSubdivide( L , 'cornercutting' );
  [np,nf] = chk( offline( Ml.xyz , L.xyz ) < 1e-14 , 'straight line stays on the line' , np,nf );
  [np,nf] = chk( max(abs( Ml.xyzF - ( 3*Ml.xyz(:,1) + 1 ) )) < 1e-12 , 'xyz* fields carried with the same weights' , np,nf );
  %star: 3 legs of 2 segments; center and leaves kept
  Xs = [ 0 0 0 ; 1 0 0 ; 2 0 0 ; 0 1 0 ; 0 2 0 ; 0 0 1 ; 0 0 2 ];
  Ts = [ 1 2 ; 2 3 ; 1 4 ; 4 5 ; 1 6 ; 6 7 ];
  Ms = MeshSubdivide( struct('xyz',Xs,'tri',Ts) , 'cornercutting' );
  [np,nf] = chk( size(Ms.tri,1) == 15 && size(Ms.xyz,1) == 16 , 'star: 6 middles + 3 corners + 6 stubs; 12 cuts + 4 kept' , np,nf );
  [np,nf] = chk( ismember( Xs(1,:) , Ms.xyz ,'rows') , 'branch node kept (network stays attached)' , np,nf );
  %orientation invariance of the point SET + subset W errors
  rng(3);  Cr = C;  rev = rand(n,1) < 0.5;  Cr.tri(rev,:) = fliplr( Cr.tri(rev,:) );
  Mr = MeshSubdivide( Cr , 'cornercutting' );
  [np,nf] = chk( isequal( sortrows( Mr.xyz ) , sortrows( Mc.xyz ) ) , 'reversed segments: same cut-point set' , np,nf );
  ok = false; try, MeshSubdivide( C , 1:5 , 'cornercutting' ); catch, ok = true; end
  [np,nf] = chk( ok , 'subset W -> error (whole-mesh only)' , np,nf );
  %'kp'
  Mk = MeshSubdivide( C , 'chaikin' , 'kp' );
  [np,nf] = chk( isequal( Mk.xyzParentEdge , [ double(C.tri) , 0.25*ones(n,1) ; double(C.tri) , 0.75*ones(n,1) ] ) , ...
                 '''kp'': parent edges with t = 1/4 and 3/4' , np,nf );

  %================================ 4POINTS
  T('C3  4points on the circle / arc / straight / star');
  Mf = MeshSubdivide( C , '4points' );
  [np,nf] = chk( size(Mf.tri,1) == 2*n && isequal( Mf.xyz(1:n,:) , C.xyz ) , 'interpolating: originals kept, 2n segments' , np,nf );
  e = max(abs( sqrt(sum(Mf.xyz.^2,2)) - 1 ));
  [np,nf] = chk( e/(3*a^4/8) > 0.9 && e/(3*a^4/8) < 1.05 , sprintf('4th order: %.3g == 3*th^4/8 (ratio %.3f)',e,e/(3*a^4/8)) , np,nf );
  Mfa = MeshSubdivide( A , '4points' );
  ea = max(abs( sqrt(sum(Mfa.xyz.^2,2)) - 1 ));
  [np,nf] = chk( isequal( Mfa.xyz(1:na,:) , A.xyz ) && ea < 5e-4 , sprintf('open arc: kept + one-sided ends fine (%.2g)',ea) , np,nf );
  Mfl = MeshSubdivide( L , '4points' );
  [np,nf] = chk( offline( Mfl.xyz , L.xyz ) < 1e-14 , 'straight polyline reproduced exactly' , np,nf );
  Mfs = MeshSubdivide( struct('xyz',Xs,'tri',Ts) , '4points' );
  [np,nf] = chk( size(Mfs.tri,1) == 12 && ~any(isnan(Mfs.xyz(:))) , 'star: runs (one-sided at the branch), finite' , np,nf );
  %non-uniform circle: documented degradation vs pn
  rng(4);  thn = cumsum( 0.05 + 0.15*rand(n,1) );  thn = thn*2*pi/thn(end);
  Cn = mkcircle( thn );
  e4 = max(abs( sqrt(sum( MeshSubdivide(Cn,'4points').xyz.^2 ,2)) - 1 ));
  ep = max(abs( sqrt(sum( MeshSubdivide(Cn,'pn'     ).xyz.^2 ,2)) - 1 ));
  [np,nf] = chk( e4/ep > 5 , sprintf('non-uniform sampling: 4points %.3g vs pn %.3g (%.0fx, as documented)',e4,ep,e4/ep) , np,nf );
  %SELECTIVE W: interpolating, so partial split == the same nodes of the full run
  Mw = MeshSubdivide( C , [1 5] , '4points' );
  ok = size( Mw.tri ,1) == n+2 && isequal( Mw.xyz(1:n,:) , C.xyz ) && ...
       isequal( Mw.xyz( n+1:n+2 ,:) , Mf.xyz( n+[1 5] ,:) );
  [np,nf] = chk( ok , 'subset W: only [1 5] split, nodes BIT-equal to the full run''s' , np,nf );
  Me4 = MeshSubdivide( mkcircle( ( 0:7 ).' * 2*pi/8 ) , -0.3 , '4points' );
  [np,nf] = chk( all( meshQuality( Me4 , 'lengths' ) <= 0.3 ) && ~any(isnan(Me4.xyz(:))) , ...
                 '-EL now composes with ''4points'' too' , np,nf );

  %================================ PN<K>
  T('C4  pn3 / pn<k>');
  M3 = MeshSubdivide( C , 'pn3' );
  [np,nf] = chk( size(M3.tri,1) == 3*n && size(M3.xyz,1) == n + 2*n , '1->3: 3n segments, 2n new nodes' , np,nf );
  e2 = max(abs( sqrt(sum( MeshSubdivide(C,'pn').xyz.^2 ,2)) - 1 ));
  e3 = max(abs( sqrt(sum( M3.xyz.^2 ,2)) - 1 ));
  [np,nf] = chk( e3 <= 1.02*e2 && e3 > 0 , sprintf('all nodes ON the gen-1 cubic: error %.3g <= pn''s %.3g',e3,e2) , np,nf );
  [np,nf] = chk( isequaln( MeshSubdivide( C , 'pn2' ) , MeshSubdivide( C , 'pn' ) ) , '''pn2'' == ''pn'' bit-exact' , np,nf );
  val3 = accumarray( M3.tri(:) , 1 );
  [np,nf] = chk( all( val3 == 2 ) , 'still one closed chain' , np,nf );
  cnt = accumarray( M3.triID , 1 , [n,1] );
  [np,nf] = chk( all( cnt == 3 ) , 'tri* fields: 3 children per parent' , np,nf );
  %straight marked piece: Bezier collapses to the exact linear t = j/k points
  d5 = [1 1 1]/sqrt(3);
  S5 = struct( 'xyz',(0:4).'*d5 , 'tri',[ (1:4).' , (2:5).' ] );
  S3 = MeshSubdivide( S5 , 'pn3' );
  LIN = [ S5.xyz ; (2/3)*S5.xyz(1:4,:) + (1/3)*S5.xyz(2:5,:) ; (1/3)*S5.xyz(1:4,:) + (2/3)*S5.xyz(2:5,:) ];
  [np,nf] = chk( max(abs( S3.xyz(:) - LIN(:) )) < 1e-14 , 'marked straight piece -> exact linear thirds' , np,nf );
  %carried normals + fields + kp at t = 1/3, 2/3
  Cn3 = C;  Cn3.xyzNORMALS = meshNormals( C , 'quadratic' );  Cn3.xyzF = th;
  Mk3 = MeshSubdivide( Cn3 , 'pn3' , 'kp' );
  NN  = Mk3.xyzNORMALS( n+1:end ,:);
  [np,nf] = chk( max(abs( sqrt(sum(NN.^2,2)) - 1 )) < 1e-12 , 'carried xyzNORMALS: unit at the new nodes' , np,nf );
  tf = [ (2/3)*th(C.tri(:,1)) + (1/3)*th(C.tri(:,2)) ; (1/3)*th(C.tri(:,1)) + (2/3)*th(C.tri(:,2)) ];
  [np,nf] = chk( max(abs( Mk3.xyzF( n+1:end ) - tf )) < 1e-12 , 'xyz* fields linear at each t' , np,nf );
  [np,nf] = chk( isequal( Mk3.xyzParentEdge( n+1:end ,3) , [ ones(n,1)/3 ; 2*ones(n,1)/3 ] ) , ...
                 '''kp'': t = 1/3 then 2/3' , np,nf );
  %pn4, -EL forwarding, and the error cases
  M4p = MeshSubdivide( C , 'pn4' );
  [np,nf] = chk( size(M4p.tri,1) == 4*n && max(abs( sqrt(sum(M4p.xyz.^2,2)) - 1 )) <= 1.02*e2 , 'pn4: 4n segments, still on the cubic' , np,nf );
  th8 = ( 0:7 ).' * 2*pi/8;
  Me  = MeshSubdivide( mkcircle( th8 ) , -0.1 , 'pn3' );
  [np,nf] = chk( all( meshQuality( Me , 'lengths' ) <= 0.1 ) && min( sqrt(sum(Me.xyz.^2,2)) ) > 0.985 , ...
                 '-EL forwards ''pn3'' (curved adaptive)' , np,nf );
  ok = 0;
  try, MeshSubdivide( struct('xyz',[0 0 0;1 0 0;0 1 0],'tri',[1 2 3]) , 'pn3' ); catch err, ok = ok + strcmp(err.identifier,'MeshSubdivide:pnK'); end
  try, MeshSubdivide( C , 'pn1' ); catch err, ok = ok + strcmp(err.identifier,'MeshSubdivide:pnK'); end
  ok2 = false; try, MeshSubdivide( C , 'bogus' ); catch err, ok2 = strcmp(err.identifier,'MeshSubdivide:segScheme'); end
  [np,nf] = chk( ok == 2 && ok2 , 'pn3-on-triangles / pn1 / unknown scheme -> clear errors' , np,nf );

  %================================ SQRT3 (Kobbelt)
  T('C5  sqrt3 on a closed sphere');
  rng(11);
  Ps = randn( 200 ,3);  Ps = Ps ./ sqrt( sum( Ps.^2 ,2) );
  Ks = convhull( Ps );
  Msph = struct( 'xyz',Ps , 'tri',Ks );
  NF = meshNormals( Msph );
  Cc = ( Ps(Ks(:,1),:) + Ps(Ks(:,2),:) + Ps(Ks(:,3),:) )/3;
  if mean( sum( NF.*Cc ,2) ) < 0, Msph.tri = Msph.tri(:,[1 3 2]); end
  nT0 = size( Msph.tri ,1);
  Msph.triID = ( 1:nT0 ).';
  S  = MeshSubdivide( Msph , 'sqrt3' );
  [np,nf] = chk( size(S.tri,1) == 3*nT0 && size(S.xyz,1) == 200+nT0 , '1->3: 3nT faces, nT centroids added' , np,nf );
  Es = sort( [ S.tri(:,[1 2]) ; S.tri(:,[2 3]) ; S.tri(:,[3 1]) ] ,2);
  [ue,~,ic] = unique( Es ,'rows');
  [np,nf] = chk( all( accumarray(ic,1) == 2 ) , 'still a closed manifold (every edge in 2 faces)' , np,nf );
  NFs = meshNormals( S );
  Ccs = ( S.xyz(S.tri(:,1),:) + S.xyz(S.tri(:,2),:) + S.xyz(S.tri(:,3),:) )/3;
  [np,nf] = chk( all( sum( NFs.*Ccs ,2) > 0 ) , 'orientation preserved (all faces outward)' , np,nf );
  [np,nf] = chk( ~any( ue(:,1) <= 200 & ue(:,2) <= 200 ) , 'ALL original edges flipped away (closed mesh)' , np,nf );
  vs  = accumarray( [ue(:,1);ue(:,2)] , 1 , [200+nT0,1] );
  ue0 = unique( sort( [ Ks(:,[1 2]) ; Ks(:,[2 3]) ; Ks(:,[3 1]) ] ,2) ,'rows');
  v0  = accumarray( [ue0(:,1);ue0(:,2)] , 1 , [200,1] );
  [np,nf] = chk( all( vs(201:end) == 6 ) && isequal( vs(1:200) , v0 ) , ...
                 'centroids valence 6; original vertices KEEP their valence' , np,nf );
  cnt = accumarray( S.triID , 1 , [nT0,1] );
  [np,nf] = chk( all( cnt == 3 ) , 'tri* inheritance: 3 children per parent' , np,nf );
  S2 = MeshSubdivide( { Msph , 2 } , 'sqrt3' );
  [np,nf] = chk( size( S2.tri ,1) == 9*nT0 , 'two passes: the realigned 3-adic 1->9' , np,nf );
  info( sprintf( 'sphere radius after 1 pass: [%.4f, %.4f] (approximating shrink)' , ...
        min(sqrt(sum(S.xyz.^2,2))) , max(sqrt(sum(S.xyz.^2,2))) ) );

  T('C6  sqrt3 on open meshes + the relaxation mask + errors');
  %planar grid: boundary kept, interior == the (1-a)p + a*mean(neighbours) mask
  [gx,gy] = meshgrid( 0:4 , 0:4 );
  GP = [ gx(:) , gy(:) , zeros(25,1) ];
  GT = delaunay( gx(:) , gy(:) );
  Gs = MeshSubdivide( struct('xyz',GP,'tri',GT) , 'sqrt3' );
  [np,nf] = chk( max(abs( Gs.xyz(:,3) )) == 0 , 'planar mesh stays exactly planar' , np,nf );
  bnd = GP(:,1)==0 | GP(:,1)==4 | GP(:,2)==0 | GP(:,2)==4;
  [np,nf] = chk( isequal( Gs.xyz(bnd,:) , GP(bnd,:) ) , 'boundary vertices NOT moved' , np,nf );
  UE0 = unique( sort( [GT(:,[1 2]);GT(:,[2 3]);GT(:,[3 1])] ,2) ,'rows');
  A0  = sparse( [UE0(:,1);UE0(:,2)] , [UE0(:,2);UE0(:,1)] , 1 , 25 , 25 );
  nv  = full( sum( A0 ,2) );
  alv = ( 4 - 2*cos( 2*pi./nv ) )/9;
  EXP = (1-alv).*GP + alv.*( (A0*GP)./nv );
  [np,nf] = chk( max(abs( Gs.xyz(~bnd,:) - EXP(~bnd,:) ),[],'all') < 1e-12 , ...
                 'interior vertices == Kobbelt valence mask, recomputed independently' , np,nf );
  %single triangle reduces to linear3
  T1s = struct( 'xyz',[0 0 0;1 0 0;0 1 0] , 'tri',[1 2 3] );
  [np,nf] = chk( isequal( sortrows( MeshSubdivide(T1s,'sqrt3').xyz ) , sortrows( MeshSubdivide(T1s,'linear3').xyz ) ) , ...
                 'single triangle: sqrt3 == linear3 (nothing to flip or relax)' , np,nf );
  %2-triangle strip: the interior edge is flipped away, boundary survives
  T2s = struct( 'xyz',[0 0 0;1 0 0;1 1 0;0 1 0] , 'tri',[1 2 3;1 3 4] );
  St  = MeshSubdivide( T2s , 'sqrt3' );
  ues = unique( sort( [St.tri(:,[1 2]);St.tri(:,[2 3]);St.tri(:,[3 1])] ,2) ,'rows');
  [np,nf] = chk( size(St.tri,1) == 6 && ~ismember( [1 3] , ues ,'rows') && isequal( St.xyz(1:4,:) , T2s.xyz ) , ...
                 'strip: 6 faces, interior edge GONE, (all-boundary) originals kept' , np,nf );
  %guards
  ok = 0;
  Tnm = struct( 'xyz',randn(5,3) , 'tri',[1 2 3;1 2 4;1 2 5] );
  try, MeshSubdivide( Tnm , 'sqrt3' ); catch err, ok = ok + strcmp( err.identifier , 'MeshSubdivide:sqrt3NonManifold' ); end
  Tw  = struct( 'xyz',[0 0 0;1 0 0;1 1 0;0 1 0] , 'tri',[1 2 3;3 1 4] );
  try, MeshSubdivide( Tw , 'sqrt3' ); catch err, ok = ok + strcmp( err.identifier , 'MeshSubdivide:sqrt3Winding' ); end
  try, MeshSubdivide( Msph , 1:5 , 'sqrt3' ); catch, ok = ok + 1; end
  [np,nf] = chk( ok == 3 , 'non-manifold / inconsistent winding / subset W -> clear errors' , np,nf );

  %================================ LINEAR9 / PN9 (triadic 1->9, selective)
  T('C7  linear9: triadic split, selective, conforming');
  T1t = struct( 'xyz',[0 0 0;3 0 0;0 3 0] , 'tri',[1 2 3] , 'triID',1 );
  L9 = MeshSubdivide( T1t , 'linear9' );
  ar = meshQuality( L9 , 'area' );
  [np,nf] = chk( size(L9.tri,1) == 9 && size(L9.xyz,1) == 10 , 'single triangle: 9 children, 10 nodes' , np,nf );
  [np,nf] = chk( max(abs( ar - 4.5/9 )) < 1e-12 , 'the 9 triadic children have EQUAL area' , np,nf );
  N9 = meshNormals( L9 );
  [np,nf] = chk( all( N9(:,3) > 0 ) , 'all children keep the parent orientation' , np,nf );
  [np,nf] = chk( all( L9.triID == 1 ) , 'tri* inherited' , np,nf );
  %selective: strip, W=1 -> neighbour gets the 1-edge fan; fan of 3 faces,
  %W=[1 3] -> the middle face gets the 2-edge (5-children) case
  T2t = struct( 'xyz',[0 0 0;3 0 0;0 3 0;3 3 0] , 'tri',[1 2 3;2 4 3] );
  S1 = MeshSubdivide( T2t , 1 , 'linear9' );
  Es = sort( [ S1.tri(:,[1 2]) ; S1.tri(:,[2 3]) ; S1.tri(:,[3 1]) ] ,2);
  [~,~,ic] = unique( Es ,'rows');
  [np,nf] = chk( size(S1.tri,1) == 12 && max( accumarray(ic,1) ) == 2 , ...
                 'strip W=1: 9 + 3-fan neighbour, CONFORMING (no T-junctions)' , np,nf );
  [np,nf] = chk( abs( sum( meshQuality( S1 ,'area') ) - 9 ) < 1e-12 , 'area conserved' , np,nf );
  T3t = struct( 'xyz',[0 0 0;3 0 0;0 3 0;3 3 0;0 6 0] , 'tri',[1 2 3;2 4 3;4 5 3] , 'triID',(1:3).' );
  S2 = MeshSubdivide( T3t , [1 3] , 'linear9' );
  cnt = accumarray( S2.triID , 1 , [3,1] );
  [np,nf] = chk( isequal( cnt , [9;5;9] ) , 'W=[1 3]: middle face -> the 2-edge 5-children case' , np,nf );
  Es = sort( [ S2.tri(:,[1 2]) ; S2.tri(:,[2 3]) ; S2.tri(:,[3 1]) ] ,2);
  [~,~,ic] = unique( Es ,'rows');
  ok = max( accumarray(ic,1) ) == 2 && all( meshQuality( S2 ,'area') > 0 ) && ...
       abs( sum( meshQuality( S2 ,'area') ) - sum( meshQuality( T3t ,'area') ) ) < 1e-12;
  [np,nf] = chk( ok , 'conforming, all children positive, area conserved' , np,nf );
  %fields at the new nodes: 2/3-1/3 on edges, mean-of-3 at the centroid
  T1f = T1t;  T1f.xyzF = [ 10 ; 20 ; 40 ];
  L9f = MeshSubdivide( T1f , 'linear9' );
  [np,nf] = chk( max(abs( sort( L9f.xyzF ) - sort( [10;20;40; (2*10+20)/3;(10+2*20)/3; (2*20+40)/3;(20+2*40)/3; (2*10+40)/3;(10+2*40)/3; 70/3] ) )) < 1e-12 , ...
                 'xyz* fields: thirds on edges, mean at the centroid' , np,nf );

  T('C8  pn9: the triadic PN patch, one pass, no accumulation');
  O = struct( 'xyz',[ 1 0 0;-1 0 0;0 1 0;0 -1 0;0 0 1;0 0 -1 ] , ...
              'tri',[ 1 3 5;3 2 5;2 4 5;4 1 5;3 1 6;2 3 6;4 2 6;1 4 6 ] );
  O.xyzNORMALS = O.xyz;
  P9 = MeshSubdivide( O , 'pn9' );
  [np,nf] = chk( size(P9.tri,1) == 72 && size(P9.xyz,1) == 6 + 2*12 + 8 , 'octahedron: 8->72 faces, 24 edge + 8 centroid nodes' , np,nf );
  %independent recompute: every new node must lie ON the generation-1 patch
  EO = unique( sort( [O.tri(:,[1 2]);O.tri(:,[2 3]);O.tri(:,[1 3])] ,2) ,'rows');
  Pu = O.xyz(EO(:,1),:);  Pv = O.xyz(EO(:,2),:);  Nu = O.xyzNORMALS(EO(:,1),:);  Nv = O.xyzNORMALS(EO(:,2),:);
  dd = Pv-Pu;
  g1 = (2*Pu+Pv)/3 - sum(dd.*Nu,2).*Nu/3;
  g2 = (Pu+2*Pv)/3 + sum(dd.*Nv,2).*Nv/3;
  QA = (8*Pu+12*g1+6*g2+Pv)/27;   QB = (Pu+6*g1+12*g2+8*Pv)/27;
  D  = [ QA ; QB ];
  ok = true;
  for r = 1:size(D,1)
    ok = ok && any( all( abs( P9.xyz - D(r,:) ) < 1e-13 ,2) );
  end
  [np,nf] = chk( ok , 'all edge nodes == independent Bezier recompute at t=1/3, 2/3' , np,nf );
  rc = sqrt( sum( P9.xyz.^2 ,2) );
  info( sprintf( 'pn9 octahedron radii in [%.4f, %.4f] (linear9 floor would be 0.577 at centroids)' , min(rc) , max(rc) ) );
  NN = P9.xyzNORMALS( 7:end ,:);
  [np,nf] = chk( max(abs( sqrt(sum(NN.^2,2)) - 1 )) < 1e-12 , 'carried xyzNORMALS: unit at every new node' , np,nf );
  %selective pn9 keeps the sphere closed and conforming
  Pw = MeshSubdivide( O , [1 2] , 'pn9' );
  Es = sort( [ Pw.tri(:,[1 2]) ; Pw.tri(:,[2 3]) ; Pw.tri(:,[3 1]) ] ,2);
  [~,~,ic] = unique( Es ,'rows');
  [np,nf] = chk( all( accumarray(ic,1) == 2 ) && ~any( isnan( Pw.xyz(:) ) ) , ...
                 'selective pn9: still a closed conforming manifold' , np,nf );
  %'pn9' on segments must STILL be the 1->9 pn<k> split
  th9 = ( 0:15 ).' * 2*pi/16;
  C9  = mkcircle( th9 );
  M9s = MeshSubdivide( C9 , 'pn9' );
  [np,nf] = chk( size( M9s.tri ,1) == 9*16 , '''pn9'' on a segment mesh is still the 1->9 pn<k>' , np,nf );
  ok = false;
  try, MeshSubdivide( O , 'pn5' ); catch err, ok = strcmp( err.identifier , 'MeshSubdivide:pnK' ); end
  [np,nf] = chk( ok , '''pn5'' on triangles -> clear error' , np,nf );

  %================================ NON-MANIFOLD edges (the 3-page "book")
  T('C9  non-manifold edge: linear/PN family tolerant, stencil schemes error');
  Bk = struct( 'xyz',[ 0 0 0 ; 1 0 0 ; 0.5 1 0 ; 0.5 -1 0 ; 0.5 0 1 ] , ...
               'tri',[ 1 2 3 ; 2 1 4 ; 1 2 5 ] );
  for s = { 'default' , 'linear3' , 'linear9' , 'pn' , 'pn9' }
    Sb = MeshSubdivide( Bk , s{1} );
    Es = sort( [ Sb.tri(:,[1 2]) ; Sb.tri(:,[2 3]) ; Sb.tri(:,[3 1]) ] ,2);
    [~,~,ic] = unique( Es ,'rows');
    ok = max( accumarray(ic,1) ) == 3 && ~any( isnan( Sb.xyz(:) ) );
    [np,nf] = chk( ok , sprintf('%-8s: subdivides the book conforming (max faces/edge stays 3)' , s{1}) , np,nf );
  end
  Sb = MeshSubdivide( Bk , 1 , 'linear9' );
  onE = abs(Sb.xyz(:,2)) < 1e-12 & abs(Sb.xyz(:,3)) < 1e-12 & Sb.xyz(:,1) > 0.01 & Sb.xyz(:,1) < 0.99;
  [np,nf] = chk( size(Sb.tri,1) == 15 && nnz(onE) == 2 , ...
                 'selective linear9 on the book: 9+3+3 children, ONE shared node pair on the seam' , np,nf );
  ok = 0;
  for s = { 'loop' , 'loop_matrix' , 'butterfly' }
    try, MeshSubdivide( Bk , s{1} ); catch err, ok = ok + strcmp( err.identifier , 'MeshSubdivide:nonManifold' ); end
  end
  [np,nf] = chk( ok == 3 , 'loop / loop_matrix / butterfly -> MeshSubdivide:nonManifold' , np,nf );

  fprintf( '\n==== %d PASS, %d FAIL ====\n' , np , nf );
  if nf, error( 'FAILURES' ); end
end

function M = mkcircle( th )
  n = numel( th );
  M = struct( 'xyz',[ cos(th) , sin(th) , zeros(n,1) ] , 'tri',[ (1:n).' , [2:n,1].' ] );
end
function M = mkchain( X )
  n = size( X ,1);
  M = struct( 'xyz',X , 'tri',[ (1:n-1).' , (2:n).' ] );
end
function d = offline( Q , L )
  u = ( L(end,:) - L(1,:) );  u = u / norm( u );
  V = Q - L(1,:);
  d = max( sqrt( sum( ( V - (V*u.')*u ).^2 ,2) ) );
end
function T( s ), fprintf( '\n--- %s\n' , s ); end
function info( s ), fprintf( '     INFO  %s\n' , s ); end
function [np,nf] = chk( c , s , np , nf )
  if c, fprintf( '     PASS  %s\n' , s ); np = np+1;
  else, fprintf( '  ** FAIL  %s\n' , s ); nf = nf+1;
  end
end
