function test_meshNormals_workflow
% The user's usual recipe on noisy meshes:
%     M.triNORMALS = meshNormals( M , 5 );    %smoothed FACE normals
%     Nv           = meshNormals( M ,'area'); %vertex normals FROM them
% Verify it works and helps, for triangles AND segments (2D/3D), incl. marks.
  warning( 'on' , 'all' );
  np = 0; nf = 0;

  %================================ A) noisy sphere, triangles
  T('A  noisy sphere (radial noise 2%), triangles');
  rng(13);
  P = randn( 500 ,3);  P = P ./ sqrt( sum( P.^2 ,2) );
  K = convhull( P );
  M = struct( 'xyz', P .* ( 1 + 0.02*randn(500,1) ) , 'tri', K );
  NF = meshNormals( M );
  C  = ( M.xyz(K(:,1),:) + M.xyz(K(:,2),:) + M.xyz(K(:,3),:) )/3;
  if mean( sum( NF.*C ,2) ) < 0, M.tri = M.tri(:,[1 3 2]); end

  err   = @(N) mean( acosd( min(1, sum( N.*P ,2) ) ) );    %vs TRUE sphere normal
  Nraw  = meshNormals( M , 'area' );                       %no triNORMALS yet
  SM    = meshNormals( M , 5 );
  M.triNORMALS = SM;
  Nsm   = meshNormals( M , 'area' );
  info( sprintf( 'mean error vs true radial: raw ''area'' %.3f deg -> smoothed(5)+''area'' %.3f deg' , err(Nraw) , err(Nsm) ) );
  [np,nf] = chk( err(Nsm) < err(Nraw) , 'smoothing first improves the vertex normals' , np,nf );
  ref = meshF2V( M , SM , 'area' );  ref = ref ./ sqrt( sum( ref.^2 ,2) );
  [np,nf] = chk( max(abs( Nsm - ref ),[],'all') < 1e-15 , '''area'' really USED M.triNORMALS' , np,nf );
  [np,nf] = chk( ~isequaln( Nsm , Nraw ) , '...and it differs from the raw-normal version' , np,nf );
  for md = { 'uniform' , 'angle' , 'best' , 'reciprocal' }
    V = meshNormals( M , md{1} );
    ok = max(abs( sqrt(sum(V.^2,2)) - 1 )) < 1e-15 && ~any( isnan(V(:)) );
    [np,nf] = chk( ok , sprintf('%-10s with triNORMALS: unit + finite' , md{1}) , np,nf );
  end
  %numeric mode CONTINUES from triNORMALS: (5 stored) + 10 more == 15 from raw
  M0 = rmfield( M , 'triNORMALS' );
  [np,nf] = chk( isequaln( meshNormals( M , 10 ) , meshNormals( M0 , 15 ) ) , ...
                 'continuation: field(5) + 10 == raw + 15, BIT-exact' , np,nf );
  [np,nf] = chk( isequaln( meshNormals( M , false ) , meshNormals( M0 , false ) ) , ...
                 'meshNormals(M,false) ignores the field: the RESET idiom' , np,nf );
  %k diffuses WHATEVER lives in triNORMALS: a random unit direction field
  rng(3);
  Mdf = M0;
  D0  = randn( size(K,1) ,3);  D0 = D0 ./ sqrt( sum( D0.^2 ,2) );
  Mdf.triNORMALS = D0;
  D1 = meshNormals( Mdf , 20 );
  A  = double( meshEsuE( Mdf , false ,'n') );
  [ii,jj] = find( triu( A ,1) );
  sm0 = mean( sum( D0(ii,:).*D0(jj,:) ,2) );
  sm1 = mean( sum( D1(ii,:).*D1(jj,:) ,2) );
  info( sprintf( 'direction-field diffusion (20 passes): mean neighbour dot %.3f -> %.3f' , sm0 , sm1 ) );
  [np,nf] = chk( sm1 > 0.8 && max(abs( sqrt(sum(D1.^2,2)) - 1 )) < 1e-15 , 'custom per-face field diffused + unit' , np,nf );

  %================================ B) noisy 3D circle, segments
  T('B  noisy 3D circle, segments (noise BELOW the curvature deflection h^2k/2)');
  th = linspace( 0 , 2*pi , 101 ).';  th(end) = [];
  tr = [ (1:100).' , [2:100,1].' ];
  R  = [ cos(th) , sin(th) , zeros(100,1) ];
  Mc = struct( 'xyz', R + 5e-4*randn(100,3) , 'tri', tr );   %sigma << h^2/2 = 0.002
  errc = @(N) mean( acosd( min(1, abs(sum( N.*R ,2)) ) ) );  %vs true radial (sign-free)
  Vraw = meshNormals( Mc , 'uniform' );
  Mc.triNORMALS = meshNormals( Mc , 5 );
  Vsm  = meshNormals( Mc , 'uniform' );
  info( sprintf( 'mean error vs true radial: raw ''uniform'' %.3f deg -> smoothed(5)+''uniform'' %.3f deg' , errc(Vraw) , errc(Vsm) ) );
  [np,nf] = chk( errc(Vsm) < errc(Vraw)/2 , 'smoothing improves segment vertex normals too (>2x)' , np,nf );
  %the caveat regime, documented in the help: sigma >> h^2k/2 -> osculating
  %planes see noise, not curvature; no normal smoothing recovers it
  rng(21);
  Mx = struct( 'xyz', R + 0.01*randn(100,3) , 'tri', tr );
  ex0 = errc( meshNormals( Mx ) );
  Mx.triNORMALS = meshNormals( Mx , 25 );
  info( sprintf( 'CAVEAT sigma=h/6: faces %.1f deg, k=25+uniform %.1f deg (unrecoverable -> smooth the GEOMETRY)' , ...
        ex0 , errc( meshNormals( Mx ,'uniform') ) ) );
  for md = { 'uniform' , 'area' , 'angle' , 'reciprocal' , 'best' , 'quadratic' }
    V = meshNormals( Mc , md{1} );
    ok = max(abs( sqrt(sum(V.^2,2)) - 1 )) < 1e-15 && ~any( isnan(V(:)) );
    [np,nf] = chk( ok , sprintf('%-10s with triNORMALS: unit + finite' , md{1}) , np,nf );
  end
  %continuation on segments too
  [np,nf] = chk( isequaln( meshNormals( Mc , 7 ) , meshNormals( rmfield(Mc,'triNORMALS') , 12 ) ) , ...
                 'segments: field(5) + 7 == raw + 12, BIT-exact' , np,nf );

  %================================ C) 2D segments
  T('C  noisy 2D circle, segments (2-col all the way)');
  M2 = struct( 'xyz', R(:,1:2) + 0.01*randn(100,2) , 'tri', tr );
  V2r = meshNormals( M2 , 'area' );
  M2.triNORMALS = meshNormals( M2 , 5 );
  V2  = meshNormals( M2 , 'area' );
  err2 = @(N) mean( acosd( min(1, abs(sum( N.*R(:,1:2) ,2)) ) ) );
  info( sprintf( '2D: raw ''area'' %.3f deg -> smoothed+''area'' %.3f deg' , err2(V2r) , err2(V2) ) );
  ok = size(V2,2) == 2 && max(abs( sqrt(sum(V2.^2,2)) - 1 )) < 1e-15 && err2(V2) < err2(V2r);
  [np,nf] = chk( ok , '2 columns, unit, improved' , np,nf );

  %================================ D) marked component through the workflow
  T('D  mixed mesh with a marked straight piece');
  d5 = [1 1 1]/sqrt(3);
  Xm = [ (0:4).'*d5 + [10 0 0] ; R ];
  Tm = [ [ (1:4).' , (2:5).' ] ; tr+5 ];
  Mm = struct( 'xyz',Xm , 'tri',Tm );
  Mm.triNORMALS = meshNormals( Mm , 5 );
  Vm = meshNormals( Mm , 'area' );
  mk = 1:5;  rest = 6:105;
  okm = all( isnan( Vm(mk,3) ) ) && max(abs( sqrt(Vm(mk,1).^2+Vm(mk,2).^2) - 1 )) < 1e-15;
  [np,nf] = chk( okm , 'marked piece: vertex normals [n1 n2 NaN], unit xy' , np,nf );
  [np,nf] = chk( max(abs( sqrt(sum(Vm(rest,:).^2,2)) - 1 )) < 1e-15 && ~any(any(isnan(Vm(rest,:)))) , ...
                 'circle vertices: unit + finite' , np,nf );

  fprintf( '\n==== %d PASS, %d FAIL ====\n' , np , nf );
  if nf, error('FAILURES'); end
end

function T( s ), fprintf( '\n--- %s\n' , s ); end
function info( s ), fprintf( '     INFO  %s\n' , s ); end
function [np,nf] = chk( c , s , np , nf )
  if c, fprintf( '     PASS  %s\n' , s ); np = np+1;
  else, fprintf( '  ** FAIL  %s\n' , s ); nf = nf+1;
  end
end
