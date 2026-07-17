function test_MeshSubdivide_pn
%TEST_MESHSUBDIVIDE_PN  PN subdivision for SEGMENT meshes (+ W / -EL modes for
% both families): the 2026-07-10 extension of MeshSubdivide. Run standalone or
% via the suite; any failure raises an error.
  HERE = fileparts( mfilename( 'fullpath' ) );
  if isempty( which( 'MeshSubdivide' ) ), addpath( fileparts( HERE ) ); end
  warning( 'on' , 'all' );
  np = 0; nf = 0;

  %================================ S1 circle: 4th order + carried normals
  T('S1  full ''pn'' on a circle (3D coplanar)');
  n  = 32;
  th = linspace( 0 , 2*pi , n+1 ).';  th(end) = [];
  C0 = mkcircle( th );
  C0.triID = ( 1:n ).';
  Mp = MeshSubdivide( C0 , 'pn' );
  [np,nf] = chk( size( Mp.tri ,1) == 2*n && size( Mp.xyz ,1) == 2*n , 'connectivity: n segments -> 2n, n new nodes' , np,nf );
  e  = max(abs( sqrt( sum( Mp.xyz.^2 ,2) ) - 1 ));
  t2 = ( pi/n )^4 * 3/8;
  [np,nf] = chk( e/t2 > 0.9 && e/t2 < 1.05 , sprintf('4th order: radial error %.3g == 3*th^4/8 (ratio %.3f)',e,e/t2) , np,nf );
  Ml = MeshSubdivide( C0 , 'default' );
  [np,nf] = chk( max(abs( sqrt(sum(Ml.xyz.^2,2)) - 1 )) > 100*e , 'linear midpoint is orders worse' , np,nf );
  [np,nf] = chk( isequal( sort( Mp.triID ) , sort( [C0.triID;C0.triID] ) ) , 'tri* fields inherited by the children' , np,nf );
  %carried normal field: updated at the new nodes, unit, ~radial
  Cn = C0;  Cn.xyzNORMALS = meshNormals( C0 , 'quadratic' );
  Mn = MeshSubdivide( Cn , 'pn' );
  NN = Mn.xyzNORMALS( n+1:end ,:);
  ok = max(abs( sqrt(sum(NN.^2,2)) - 1 )) < 1e-12 && ...
       max( acosd( min(1,abs(sum( NN .* ( Mn.xyz(n+1:end,:)./sqrt(sum(Mn.xyz(n+1:end,:).^2,2)) ) ,2))) ) ) < 0.5;
  [np,nf] = chk( ok , 'carried xyzNORMALS: new nodes get unit ~radial PN normals' , np,nf );
  %iterated: PN is INTERPOLATING, so first-generation nodes (and their error)
  %persist forever; later passes must not GROW the error (they refine between)
  M2 = MeshSubdivide( { Cn , 2 } , 'pn' );
  e2 = max(abs( sqrt( sum( M2.xyz.^2 ,2) ) - 1 ));
  [np,nf] = chk( e2 <= e*1.05 , sprintf('two passes: error stays bounded by 1st generation (%.3g <= %.3g)',e2,e) , np,nf );
  %the ORDER statement is per-pass from a given sampling: one pass at 2n is 16x
  th2 = linspace( 0 , 2*pi , 2*n+1 ).';  th2(end) = [];
  M1f = MeshSubdivide( mkcircle( th2 ) , 'pn' );
  e1f = max(abs( sqrt( sum( M1f.xyz.^2 ,2) ) - 1 ));
  [np,nf] = chk( e/e1f > 14 && e/e1f < 18 , sprintf('one pass from 2n: error /%.1f (4th order: 16x)',e/e1f) , np,nf );

  %================================ S2 selective W
  T('S2  selective W on segments');
  Mw = MeshSubdivide( C0 , [1 5] , 'pn' );
  [np,nf] = chk( size( Mw.tri ,1) == n+2 && size( Mw.xyz ,1) == n+2 , 'only 2 segments split' , np,nf );
  [np,nf] = chk( isequal( Mw.xyz(1:n,:) , C0.xyz ) , 'original nodes untouched' , np,nf );
  ch = sqrt( sum( Mw.xyz( n+1:end ,:).^2 ,2) );
  [np,nf] = chk( all( abs( ch - 1 ) < 1e-3 ) && all( ch > cos(pi/n) + 1e-6 ) , 'the 2 new nodes are LIFTED off their chords' , np,nf );

  %================================ S3 adaptive -EL with scheme
  T('S3  -EL adaptive mode, with and without scheme');
  n8  = 8;
  th8 = linspace( 0 , 2*pi , n8+1 ).';  th8(end) = [];
  C8  = mkcircle( th8 );
  EL  = 0.1;
  Mpn = MeshSubdivide( C8 , -EL , 'pn' );
  [np,nf] = chk( all( meshQuality( Mpn , 'length' ) <= EL ) , '-EL,''pn'': all edges <= EL' , np,nf );
  rpn = sqrt( sum( Mpn.xyz.^2 ,2) );
  Mln = MeshSubdivide( C8 , -EL );
  rln = sqrt( sum( Mln.xyz.^2 ,2) );
  %interpolating: the 1st-generation error (3*(pi/8)^4/8 = 8.9e-3, coarse start) persists
  [np,nf] = chk( min( rpn ) > 1 - 1.2*(3*(pi/n8)^4/8) && min( rln ) < cos(pi/n8) + 1e-9 , ...
                 sprintf('-EL,''pn'' tracks the circle (min r %.4f, theory %.4f) vs linear on the chords (min r %.4f)',min(rpn),1-3*(pi/n8)^4/8,min(rln)) , np,nf );
  [np,nf] = chk( all( meshQuality( Mln , 'length' ) <= EL ) , '-EL alone: backward compatible (linear, edges <= EL)' , np,nf );
  Mkp = MeshSubdivide( C8 , -EL , 'pn' , 'kp' );
  [np,nf] = chk( isfield( Mkp , 'xyzParentEdge' ) , '-EL forwards ''kp'' too' , np,nf );

  %================================ S4 marked straight piece stays straight
  T('S4  straight (marked) piece');
  d5 = [1 1 1]/sqrt(3);
  S5 = struct( 'xyz',(0:4).'*d5 , 'tri',[ (1:4).' , (2:5).' ] );
  Sp = MeshSubdivide( S5 , 'pn' );
  Sl = MeshSubdivide( S5 , 'default' );
  [np,nf] = chk( isequal( Sp.xyz , Sl.xyz ) , '''pn'' == ''default'' bit-exact (marked normals -> zero lift)' , np,nf );
  S5n = S5;  S5n.xyzNORMALS = meshNormals( S5 , 'uniform' );      %[n1 n2 NaN] marks
  Spn = MeshSubdivide( S5n , 'pn' );
  NN  = Spn.xyzNORMALS( 6:end ,:);
  [np,nf] = chk( all( isnan( NN(:,3) ) ) && max(abs( sqrt(NN(:,1).^2+NN(:,2).^2) - 1 )) < 1e-12 , ...
                 'carried marks: new normals stay [n1 n2 NaN], unit xy' , np,nf );

  %================================ S5 2D circle
  T('S5  2D segment mesh');
  C2 = struct( 'xyz',[ cos(th) , sin(th) ] , 'tri',C0.tri );
  Mp2 = MeshSubdivide( C2 , 'pn' );
  e2d = max(abs( sqrt( sum( Mp2.xyz.^2 ,2) ) - 1 ));
  [np,nf] = chk( size( Mp2.xyz ,2) == 2 && e2d/t2 > 0.9 && e2d/t2 < 1.05 , ...
                 sprintf('2 columns, same 4th order (%.3g)',e2d) , np,nf );

  %================================ S6 'kp' for segments
  T('S6  ''kp'' on segments (was silently ignored)');
  Mk = MeshSubdivide( C0 , [1 5] , 'kp' );
  ok = isfield( Mk , 'xyzParentEdge' ) && ...
       isequal( Mk.xyzParentEdge( 1:n ,:) , [ (1:n).' , zeros(n,2) ] ) && ...
       isequal( Mk.xyzParentEdge( n+1:end ,:) , [ double( C0.tri([1 5],:) ) , [0.5;0.5] ] );
  [np,nf] = chk( ok , 'xyzParentEdge: originals [i 0 0], new nodes [e1 e2 0.5]' , np,nf );

  %================================ S7 stale xyzNORMALS guards
  T('S7  stale xyzNORMALS -> clear error');
  Cs = C0;  Cs.xyzNORMALS = zeros( 5 , 3 );
  ok = false;
  try, MeshSubdivide( Cs , 'pn' ); catch err, ok = strcmp( err.identifier , 'MeshSubdivide:xyzNORMALS' ); end
  [np,nf] = chk( ok , 'segments: MeshSubdivide:xyzNORMALS' , np,nf );
  O = octa();  O.xyzNORMALS = zeros( 5 , 3 );
  ok = false;
  try, MeshSubdivide( O , 'pn' ); catch err, ok = strcmp( err.identifier , 'MeshSubdivide:xyzNORMALS' ); end
  [np,nf] = chk( ok , 'triangles: MeshSubdivide:xyzNORMALS' , np,nf );

  %================================ S8 triangles: W and -EL with 'pn'
  T('S8  triangles: ''pn'' with W and -EL');
  O = octa();  O.xyzNORMALS = O.xyz;                     %exact radial normals
  Op = MeshSubdivide( O , 'pn' );
  rn = sqrt( sum( Op.xyz( 7:end ,:).^2 ,2) );            %the new (edge) nodes
  Ol = MeshSubdivide( rmfield( O ,'xyzNORMALS' ) , 'default' );
  rl = sqrt( sum( Ol.xyz( 7:end ,:).^2 ,2) );
  %closed form: the PN edge midpoint of the octahedron is (5/8)*sqrt(2) = 0.8839
  [np,nf] = chk( max(abs( rn - 5*sqrt(2)/8 )) < 1e-12 && max(abs( rl - sqrt(2)/2 )) < 1e-12 , ...
                 sprintf('octahedron: PN new nodes r == 5*sqrt(2)/8 = %.4f exactly (linear: 0.7071)',5*sqrt(2)/8) , np,nf );
  Ow = MeshSubdivide( O , [1 2] , 'pn' );
  [np,nf] = chk( size( Ow.tri ,1) > 8 && ~any( isnan( Ow.xyz(:) ) ) , 'selective W + pn: runs, conforming, finite' , np,nf );
  Oe = MeshSubdivide( O , -0.6 , 'pn' );
  rr = sqrt( sum( Oe.xyz.^2 ,2) );
  %the PN limit surface of a COARSE octahedron is a "pillow": its interior dips
  %somewhat below the edge-curve radius 5*sqrt(2)/8 (measured 0.861) -- still
  %far above the linear 0.707 -- and never overshoots the unit sphere.
  ok = all( all( meshQuality( Oe , 'lengths' ) <= 0.6 ) ) && ...
       min( rr ) > 0.85 && max( rr ) < 1 + 1e-9;
  [np,nf] = chk( ok , sprintf('-EL,''pn'' triangles: edges <= EL, r in [%.3f,%.3f] (PN pillow; linear floor: 0.707)',min(rr),max(rr)) , np,nf );

  %================================ S9 helix end-to-end
  T('S9  helix: pn vs default, distance to the true curve');
  t  = ( 0 : 0.3 : 4*pi ).';
  H  = struct( 'xyz',[ cos(t) , sin(t) , 0.2*t ] , 'tri',[ (1:numel(t)-1).' , (2:numel(t)).' ] );
  Hp = MeshSubdivide( H , 'pn' );
  Hl = MeshSubdivide( H , 'default' );
  tm = ( t(1:end-1) + t(2:end) )/2;
  TR = [ cos(tm) , sin(tm) , 0.2*tm ];
  ep = max( sqrt( sum( ( Hp.xyz( numel(t)+1:end ,:) - TR ).^2 ,2) ) );
  el = max( sqrt( sum( ( Hl.xyz( numel(t)+1:end ,:) - TR ).^2 ,2) ) );
  [np,nf] = chk( el/ep > 10 , sprintf('new nodes: %.4g (linear) -> %.4g (pn), %.0fx closer',el,ep,el/ep) , np,nf );

  %================================ S10 orientation / sign invariance
  T('S10 geometry invariant to connectivity reversal and normal signs');
  rng(7);
  Cr = C0;  rev = rand(n,1) < 0.4;  Cr.tri(rev,:) = fliplr( Cr.tri(rev,:) );
  R  = [ cos(th) , sin(th) , zeros(n,1) ];
  Cn0 = C0;  Cn0.xyzNORMALS = R;
  Cnr = Cr;  Cnr.xyzNORMALS = R .* ( 1 - 2*( rand(n,1) < 0.5 ) );
  A = MeshSubdivide( Cn0 , 'pn' );
  B = MeshSubdivide( Cnr , 'pn' );
  [np,nf] = chk( isequal( A.xyz , B.xyz ) , 'same new GEOMETRY, bit-exact, whatever the orientations/signs' , np,nf );

  fprintf( '\n==== %d PASS, %d FAIL ====\n' , np , nf );
  if nf, error( 'FAILURES' ); end
end

function M = mkcircle( th )
  n = numel( th );
  M = struct( 'xyz',[ cos(th) , sin(th) , zeros(n,1) ] , 'tri',[ (1:n).' , [2:n,1].' ] );
end
function M = octa()
  M = struct( 'xyz',[ 1 0 0;-1 0 0;0 1 0;0 -1 0;0 0 1;0 0 -1 ] , ...
              'tri',[ 1 3 5;3 2 5;2 4 5;4 1 5;3 1 6;2 3 6;4 2 6;1 4 6 ] );
end
function T( s ), fprintf( '\n--- %s\n' , s ); end
function [np,nf] = chk( c , s , np , nf )
  if c, fprintf( '     PASS  %s\n' , s ); np = np+1;
  else, fprintf( '  ** FAIL  %s\n' , s ); nf = nf+1;
  end
end
