function test_meshNormals_F3
  warning( 'on' , 'all' );
  np = 0; nf = 0;

  %================================================= 1) routing & aliases
  T('F3-1 routing / aliases (non-uniform 2D circle)');
  rng(4);
  th = cumsum( 0.05 + 0.15*rand( 60 ,1) );  th = th * 2*pi/th(end);
  n  = numel(th);
  P2 = [ cos(th) , sin(th) ];
  tr = [ (1:n).' , [2:n,1].' ];
  M2 = mst( P2 , tr );
  NF = meshNormals( M2 );
  Nu = meshNormals( M2 , 'uniform' );
  [np,nf] = chk( isequaln( meshNormals(M2,'angle') , Nu ) , '''angle'' == ''uniform''' , np,nf );
  Nr = meshNormals( M2 , 'reciprocal' );
  [np,nf] = chk( isequaln( meshNormals(M2,'reciproco') , Nr ) && isequaln( meshNormals(M2,'r') , Nr ) , ...
                 '''reciproco'' / ''r'' aliases' , np,nf );
  Nq = meshNormals( M2 , 'quadratic' );
  [np,nf] = chk( isequaln( meshNormals(M2,'q') , Nq ) , '''q'' alias' , np,nf );
  %manual check of 'area' -> length routing
  L  = sqrt( sum( ( P2(tr(:,2),:) - P2(tr(:,1),:) ).^2 ,2) );
  Sl = [ accumarray( tr(:), repmat(L,2,1).*[NF(:,1);NF(:,1)] ,[n,1]) , ...
         accumarray( tr(:), repmat(L,2,1).*[NF(:,2);NF(:,2)] ,[n,1]) ];
  Na = meshNormals( M2 , 'area' );
  [np,nf] = chk( max(abs( Na - Sl./sqrt(sum(Sl.^2,2)) ),[],'all') < 1e-15 , '''area'' == length-weighted sum' , np,nf );

  %================================================= 2) accuracy on the circle
  T('F3-2 accuracy: exact radial reference');
  e = @(N) max( acosd( min(1, sum(N.*P2,2) ) ) );
  [np,nf] = chk( e(Nr) < 1e-6 , sprintf('''reciprocal'' EXACT on the circle (max %.2g deg)',e(Nr)) , np,nf );
  [np,nf] = chk( max(abs( sum( Nq.*Nr ,2) - 1 )) < 1e-12 , '''quadratic'' == ''reciprocal'' in 2D' , np,nf );
  info( sprintf( 'max error: uniform %.3f, area(length) %.3f, reciprocal %.2g, quadratic %.2g deg' , ...
        e(Nu) , e(Na) , e(Nr) , e(Nq) ) );
  for X = { Nu , Na , Nr , Nq }
    if max(abs( sqrt(sum(X{1}.^2,2)) - 1 )) > 1e-15, nf = nf+1; fprintf('  ** FAIL unit norm\n'); end
  end
  [np,nf] = chk( true , 'all four outputs unit (norm-1 loop)' , np,nf );

  %================================================= 3) tilted 3D coplanar circle
  T('F3-3 tilted coplanar circle in 3D');
  ax = [1 2 3]; ax = ax/norm(ax);
  R0 = expm( skw( ax*0.8 ) );
  M3 = mst( [ P2 , zeros(n,1) ] * R0.' + [5 -3 2] , tr );
  R3 = P2 * R0(1:3,1:2).';                       %exact radial, rotated
  Nr3 = meshNormals( M3 , 'reciprocal' );
  Nq3 = meshNormals( M3 , 'quadratic' );
  [np,nf] = chk( max( acosd( min(1,abs(sum(Nr3.*R3,2))) ) ) < 1e-5 , '''reciprocal'' exact on the tilted circle' , np,nf );
  [np,nf] = chk( max(abs( sum( Nq3.*Nr3 ,2) - 1 )) < 1e-9 , '''quadratic'' == ''reciprocal'' on coplanar pieces' , np,nf );

  %================================================= 4) helix: quadratic vs analytic tangent
  T('F3-4 helix: quadratic perpendicular to the analytic tangent');
  t  = ( 0 : 0.15 : 4*pi ).';
  H  = [ cos(t) , sin(t) , 0.2*t ];
  trh = [ (1:numel(t)-1).' , (2:numel(t)).' ];
  MH = mst( H , trh );
  NqH = meshNormals( MH , 'quadratic' );
  TA  = [ -sin(t) , cos(t) , 0.2*ones(size(t)) ] / sqrt(1.04);
  dpt = abs( sum( NqH .* TA ,2) );
  [np,nf] = chk( max( dpt(2:end-1) ) < 5e-3 , sprintf('interior: |dot(N,tangent)| < 5e-3 (max %.2g)',max(dpt(2:end-1))) , np,nf );
  [np,nf] = chk( max( dpt([1 end]) ) < 5e-2 , sprintf('one-sided ends: |dot| < 5e-2 (max %.2g)',max(dpt([1 end]))) , np,nf );
  Sb = meshNormals( MH , 'uniform' );
  [np,nf] = chk( all( sum( NqH.*Sb ,2) > 0 ) , 'sign follows the incident face normals' , np,nf );

  %================================================= 5) marked pieces
  T('F3-5 marked straight piece: all vertex modes stay marked, unit xy');
  d5 = [1 1 1]/sqrt(3);
  Ms = mst( (0:4).'*d5 , [ (1:4).' , (2:5).' ] );
  for md = { 'uniform' , 'angle' , 'area' , 'reciprocal' , 'quadratic' }
    Vm = meshNormals( Ms , md{1} );
    ok = all( isnan( Vm(:,3) ) ) && max(abs( sqrt(Vm(:,1).^2+Vm(:,2).^2) - 1 )) < 1e-15;
    [np,nf] = chk( ok , sprintf('%-10s -> [n1 n2 NaN], unit xy' , md{1}) , np,nf );
  end

  %================================================= 6) open arc: one-sided ends
  T('F3-6 open quarter-circle arc, non-uniform');
  rng(9);
  ta = cumsum( 0.03 + 0.1*rand( 20 ,1) );  ta = ta * (pi/2)/ta(end);  ta = [0;ta];
  Pa = [ cos(ta) , sin(ta) ];
  na = numel(ta);
  Ma = mst( Pa , [ (1:na-1).' , (2:na).' ] );
  NuA = meshNormals( Ma , 'uniform'   );
  NqA = meshNormals( Ma , 'quadratic' );
  eu = acosd( min(1, sum( NuA([1 na],:) .* Pa([1 na],:) ,2) ) );
  eq = acosd( min(1, sum( NqA([1 na],:) .* Pa([1 na],:) ,2) ) );
  info( sprintf( 'endpoint error: uniform (face normal) [%.3f %.3f] deg, quadratic one-sided [%.4f %.4f] deg' , eu , eq ) );
  [np,nf] = chk( all( eq < eu/3 ) , 'one-sided parabola >> face normal at open ends' , np,nf );
  ei = max( acosd( min(1, sum( NqA(2:end-1,:).*Pa(2:end-1,:) ,2) ) ) );
  [np,nf] = chk( ei < 1e-6 , sprintf('interior still exact (max %.2g deg)',ei) , np,nf );
  %lone segment: quadratic falls back to its face normal
  M1 = mst( [0 0;2 0] , [1 2] );
  [np,nf] = chk( isequaln( meshNormals(M1,'quadratic') , meshNormals(M1,'uniform') ) , ...
                 'lone segment: quadratic == face normal (fallback)' , np,nf );

  %================================================= 7) 'best' on segments
  T('F3-7 ''best'' on segments');
  [np,nf] = chk( isequaln( meshNormals(M3,'best') , meshNormals(M3,'uniform') ) , ...
                 'clean circle: best == uniform (bisector already max-margin)' , np,nf );
  Xs = [ 0 0 0 ; 1 0 0 ; 0 1 0 ; 0 0 1 ; -1 -1 -1 ];   %branching star, 4 legs
  Ts = [ 1 2 ; 1 3 ; 1 4 ; 1 5 ];
  w0 = warning( 'off' , 'meshNormals:emptyNormalCone' );
  Nb = meshNormals( mst( Xs , Ts ) , 'best' );
  warning( w0 );
  [np,nf] = chk( max(abs( sqrt(sum(Nb(~any(isnan(Nb),2),:).^2,2)) - 1 )) < 1e-15 , ...
                 'branching star: runs, unit rows' , np,nf );

  %================================================= 8) triangles: byte-equal routing
  T('F3-8 triangles (celltype 5): unchanged');
  Tt = [ 0 0 0 ; 1 0 0 ; 0 1 0 ; 0 0 1 ];
  Ft = [ 1 3 2 ; 1 2 4 ; 2 3 4 ; 1 4 3 ];
  Mt = mst( Tt , Ft );
  NFt = meshNormals( Mt );
  for md = { {'uniform','sum'} , {'angle','angles'} , {'area','area'} }
    ref = meshF2V( Mt , NFt , md{1}{2} );
    ref = ref ./ sqrt( sum( ref.^2 ,2) );
    D = abs( meshNormals( Mt , md{1}{1} ) - ref );
    [np,nf] = chk( max( D(:) ) < 1e-15 , sprintf('''%s'' routes to meshF2V ''%s''' , md{1}{1} , md{1}{2}) , np,nf );
  end
  [np,nf] = chk( isequaln( meshNormals(Mt,'best') , meshNormals(Mt,'angle') ) , '''best'' == ''angle'' on the tetra (no bad cones)' , np,nf );
  ok = 0;
  try, meshNormals( Mt , 'quadratic'  ); catch e, ok = ok + strcmp( e.identifier , 'meshNormals:mode' ); end
  try, meshNormals( Mt , 'bogus'      ); catch e, ok = ok + strcmp( e.identifier , 'meshNormals:mode' ); end
  [np,nf] = chk( ok == 2 , 'quadratic on triangles + unknown mode -> meshNormals:mode errors' , np,nf );
  Nrt = meshNormals( Mt , 'reciprocal' );      %now valid on triangles (Max weights)
  cc = Tt - mean( Tt ,1);
  [np,nf] = chk( all( sum( Nrt.*cc ,2) > 0 ) && max(abs( sqrt(sum(Nrt.^2,2)) - 1 )) < 1e-15 , ...
                 '''reciprocal'' now works on triangles: outward + unit on the tetra' , np,nf );

  %================================================= 9) degenerate segment poisoning
  T('F3-9 zero-length segment in the circle');
  M9 = mst( [ P2 , zeros(n,1) ] , [ tr ; 5 5 ] );
  for md = { 'uniform' , 'reciprocal' , 'quadratic' }
    V9 = meshNormals( M9 , md{1} );
    okv = all( isnan( V9(5,:) ) );                          %vertex 5 poisoned by its degenerate cell
    others = true(n,1); others(5) = false;
    oku = max(abs( sqrt(sum(V9(others,:).^2,2)) - 1 )) < 1e-15 && ~any(any(isnan(V9(others,:))));
    [np,nf] = chk( okv && oku , sprintf('%-10s: only the degenerate''s vertex poisoned' , md{1}) , np,nf );
  end

  fprintf( '\n==== %d PASS, %d FAIL ====\n' , np , nf );
  if nf, error('FAILURES'); end
end

function M = mst( xyz , tri ), M = struct( 'xyz',xyz , 'tri',tri ); end
function s = skw( v ), s = [ 0 -v(3) v(2) ; v(3) 0 -v(1) ; -v(2) v(1) 0 ]; end
function T( s ), fprintf( '\n--- %s\n' , s ); end
function info( s ), fprintf( '     INFO  %s\n' , s ); end
function [np,nf] = chk( c , s , np , nf )
  if c, fprintf( '     PASS  %s\n' , s ); np = np+1;
  else, fprintf( '  ** FAIL  %s\n' , s ); nf = nf+1;
  end
end
