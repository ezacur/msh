function test_meshNormals_fix
  warning( 'on' , 'all' );
  SP = fileparts( mfilename('fullpath') );
  old = load( fullfile( SP , 'old_meshNormals.mat' ) );
  np = 0; nf = 0;

  %---------------- R1 upright helix: regression vs OLD (old had 0 flips there)
  T('R1  upright helix == old (where old was already consistent)');
  MH = mst( old.H , old.trh );
  lastwarn('');
  NH = meshNormals( MH );
  [np,nf] = chk( max( 1 - sum( NH .* old.N_helix ,2) ) < 1e-9 , 'same normals as before (same sign, <1e-4 deg)' , np,nf );
  [np,nf] = chk( isempty( lastwarn() ) , 'no warnings' , np,nf );

  %---------------- R2 tilted helix: THE fix (old had 9 sign flips)
  T('R2  tilted helix (old: 9 flips)');
  MT = mst( old.H * old.R0.' , old.trh );
  NT = meshNormals( MT );
  tm = ( old.trh(:,1) + 0.5 - 1 ) * 0.15;   %midpoint parameter (t = 0:0.15:4pi)
  NAt = [ cos(tm) , sin(tm) , 0*tm ] * old.R0.';
  dp = sum( NT .* NAt ,2);
  info( sprintf( 'dot(N, radial-line): [%.4f, %.4f], flips: %d (old: 9)' , min(dp) , max(dp) , nnz(diff(sign(dp))) ) );
  [np,nf] = chk( nnz( diff( sign(dp) ) ) == 0 , 'ZERO sign flips along the tilted helix' , np,nf );
  [np,nf] = chk( min( abs(dp) ) > 0.999 , 'still on the curvature line everywhere' , np,nf );

  %---------------- R3 random smooth 3D curve (old: 11 jumps > 90 deg)
  T('R3  random smooth curve (old: 11 jumps)');
  rng(7);
  s = linspace( 0 , 1 , 400 ).';
  C = zeros( numel(s) , 3 );
  for h = 1:4
    C = C + ( randn(1,3) .* cos( 2*pi*h*s + 2*pi*rand ) + randn(1,3) .* sin( 2*pi*h*s ) ) / h^2;
  end
  trc = [ (1:numel(s)-1).' , (2:numel(s)).' ];
  lastwarn('');
  N3 = meshNormals( mst( C , trc ) );
  ca = sum( N3(1:end-1,:) .* N3(2:end,:) ,2);
  info( sprintf( 'consecutive-normal cos in [%.4f, %.4f]' , min(ca) , max(ca) ) );
  [np,nf] = chk( all( ca > 0 ) , 'no jumps > 90 deg anymore' , np,nf );
  [np,nf] = chk( isempty( lastwarn() ) , 'no warnings' , np,nf );

  %---------------- R4 S-curve: matches its planar twin (2D-rotated semantics)
  T('R4  gently-non-planar S-curve vs its planar twin');
  x  = linspace( -2*pi , 2*pi , 200 ).';
  trs = [ (1:numel(x)-1).' , (2:numel(x)).' ];
  N4  = meshNormals( mst( [ x , sin(x) , 0.05*x.^2/10 ] , trs ) );  %osculating branch
  N4p = meshNormals( mst( [ x , sin(x) , 0*x           ] , trs ) );  %coplanar branch (z=0)
  dp4 = sum( N4 .* N4p ,2);
  xc  = ( x(1:end-1) + x(2:end) )/2;
  %at the sin inflections the PARABOLIC z-curvature dominates: the true osculating
  %plane turns vertical there and the normal passes smoothly THROUGH +-z; away from
  %them it must track the twin's line. Side switches only via those smooth passes.
  ca4 = sum( N4(1:end-1,:) .* N4(2:end,:) ,2);
  [np,nf] = chk( all( ca4 > 0 ) , 'no reversals along the curve (vertical pass under-sampled: dx=0.063 > transition |x|<~0.01)' , np,nf );
  info( sprintf( 'min consecutive cos = %.3f (at the fast vertical pass near x~0)' , min(ca4) ) );
  md = abs( interp1( xc , dp4 , [-4.7 -1.6 1.6 4.7] ) );
  [np,nf] = chk( all( md > 0.999 ) , 'tracks the planar twin''s line away from inflections' , np,nf );
  [np,nf] = chk( nnz( diff( sign(dp4) ) ) <= 3 , 'side changes only through real vertical-curvature passes' , np,nf );
  info( sprintf( 'dot(osc,twin) at midpoints: %s; sign changes: %d (smooth, |Nz|->%.2f at x~0)' , ...
        sprintf('%.4f ',md) , nnz(diff(sign(dp4))) , max(abs(N4(abs(xc)<0.3,3))) ) );

  %---------------- R5 straight run inside 3D component (old: arbitrary + 2 warnings)
  T('R5  short straight run inside 3D component');
  M10 = mst( old.X10 , [ (1:5).' , (2:6).' ] );
  lastwarn('');
  N10 = meshNormals( M10 );
  [np,nf] = chk( isempty( lastwarn() ) , 'no getPlane:collinear warnings' , np,nf );
  [np,nf] = chk( all( sum( N10(1:end-1,:) .* N10(2:end,:) ,2) > 0 ) , 'run continuous with the bend' , np,nf );
  [np,nf] = chk( max(abs( sqrt(sum(N10.^2,2)) - 1 )) < 1e-12 && ~any(isnan(N10(:))) , 'unit + finite' , np,nf );

  %---------------- R6 LONG straight run (old: 98 warnings in one call)
  T('R6  100-segment straight run inside 3D component (old: 98 warnings)');
  X = [ (0:99).' , zeros(100,2) ; 100.7 .7 .7 ; 101 1.5 .2 ];
  M6 = mst( X , [ (1:101).' , (2:102).' ] );
  w0 = warning( 'off' , 'backtrace' );
  lastwarn('');
  out = evalc( 'N6 = meshNormals( M6 );' );
  warning( w0 );
  [np,nf] = chk( ~contains( out , 'collinear' ) && isempty( lastwarn() ) , 'ZERO warnings' , np,nf );
  ca6 = sum( N6(1:end-1,:) .* N6(2:end,:) ,2);
  [np,nf] = chk( all( ca6(1:98) > 0.9999 ) , 'run-interior normals constant (inherited plane; only the 2 bend pairs rotate)' , np,nf );
  [np,nf] = chk( all( ca6 > 0 ) , 'continuous through the bend too' , np,nf );

  %---------------- R7 CHIRALITY: random reversals flip EXACTLY those normals
  T('R7  chirality: reversed segments flip their own normal, exactly');
  rng(5);
  rev = rand( size( old.trh ,1) ,1) < 0.3;
  TRr = old.trh;  TRr(rev,:) = fliplr( TRr(rev,:) );
  NTr = meshNormals( mst( old.H * old.R0.' , TRr ) );
  sgn = 1 - 2*rev;
  [np,nf] = chk( isequaln( NTr , NT .* sgn ) , '3D branch: N(rev) == -N, bit-exact, others untouched' , np,nf );
  %same property on a coplanar circle (was already true; the user rule)
  t8 = linspace( 0 , 2*pi , 41 ).'; t8(end) = [];
  tr8 = [ (1:40).' , [2:40,1].' ];
  rev8 = rand(40,1) < 0.5;
  TR8 = tr8;  TR8(rev8,:) = fliplr( TR8(rev8,:) );
  Nc  = meshNormals( mst( [cos(t8),sin(t8),0*t8] , tr8 ) );
  Nc8 = meshNormals( mst( [cos(t8),sin(t8),0*t8] , TR8 ) );
  [np,nf] = chk( isequaln( Nc8 , Nc .* (1-2*rev8) ) , 'coplanar branch: same rule, bit-exact' , np,nf );

  %---------------- R8 continuity across the coplanar/osculating threshold
  T('R8  nearly-coplanar curve agrees with the coplanar branch (sign included)');
  rng(2);  nz = randn(40,1);
  Nflat = meshNormals( mst( [cos(t8),sin(t8),0*t8] , tr8 ) );
  Nnear = meshNormals( mst( [cos(t8),sin(t8),1e-6*nz] , tr8 ) );
  dpn = sum( Nnear .* Nflat ,2);
  [np,nf] = chk( all( dpn > 0.999 ) , 'same normals AND same sign (seed ''+z'' == coplanar ''+z'')' , np,nf );

  %---------------- R9 3D star (mostly-collinear neighbourhoods, hub seeds)
  T('R9  3D star: 3 straight legs sharing one node');
  Xs = [ 0 0 0 ; 1 0 0 ; 2 0 0 ; 0 1 0 ; 0 2 0 ; 0 0 1 ; 0 0 2 ];
  Ts = [ 1 2 ; 2 3 ; 1 4 ; 4 5 ; 1 6 ; 6 7 ];
  lastwarn('');
  Ns = meshNormals( mst( Xs , Ts ) );
  [np,nf] = chk( isempty( lastwarn() ) && ~any(isnan(Ns(:))) && max(abs(sqrt(sum(Ns.^2,2))-1)) < 1e-12 , ...
                 'no warnings, finite, unit' , np,nf );
  sdir = Xs( Ts(:,2) ,:) - Xs( Ts(:,1) ,:);
  [np,nf] = chk( max(abs( sum( Ns.*sdir ,2) )) < 1e-12 , 'all perpendicular to their segment' , np,nf );

  %---------------- R10 determinism
  T('R10 determinism');
  [np,nf] = chk( isequaln( meshNormals( MT ) , NT ) , 'same input -> bit-identical output' , np,nf );

  %---------------- R11 timing
  T('R11 timing (was 0.10 s / 2513 segs)');
  t  = ( 0 : 0.005 : 4*pi ).';
  MB = mst( [ cos(t) , sin(t) , 0.2*t ] , [ (1:numel(t)-1).' , (2:numel(t)).' ] );
  tic;  meshNormals( MB );  tt = toc;
  info( sprintf( '%d segments: %.3f s (%.3f ms/segment)' , size(MB.tri,1) , tt , 1000*tt/size(MB.tri,1) ) );

  fprintf( '\n==== %d PASS, %d FAIL ====\n' , np , nf );
  if nf, error('FAILURES'); end
end

function M = mst( xyz , tri ), M = struct( 'xyz',xyz , 'tri',tri ); end
function T( s ), fprintf( '\n--- %s\n' , s ); end
function info( s ), fprintf( '     INFO  %s\n' , s ); end
function [np,nf] = chk( c , s , np , nf )
  if c, fprintf( '     PASS  %s\n' , s ); np = np+1;
  else, fprintf( '  ** FAIL  %s\n' , s ); nf = nf+1;
  end
end
