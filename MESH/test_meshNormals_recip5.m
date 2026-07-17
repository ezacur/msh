function test_meshNormals_recip5
  warning( 'on' , 'all' );
  np = 0; nf = 0;

  %---------------- sphere: exactness + comparison
  T('R5-1 random irregular sphere (400 pts, convhull)');
  rng(11);
  P = randn( 400 , 3 );  P = P ./ sqrt( sum( P.^2 ,2) );
  K = convhull( P );
  M = struct( 'xyz',P , 'tri',K );
  NF = meshNormals( M );
  C  = ( P(K(:,1),:) + P(K(:,2),:) + P(K(:,3),:) )/3;
  if mean( sum( NF.*C ,2) ) < 0, M.tri = M.tri(:,[1 3 2]); NF = meshNormals( M ); end

  err = @(N) max( acosd( min(1, sum( N.*P ,2) ) ) );
  Nr = meshNormals( M , 'reciprocal' );
  [np,nf] = chk( err(Nr) < 1e-4 , sprintf('''reciprocal'' EXACT on the sphere (max %.2g deg)',err(Nr)) , np,nf );
  info( sprintf( 'uniform %.2f, area %.2f, angle %.2f, reciprocal %.2g deg' , ...
        err(meshNormals(M,'uniform')) , err(meshNormals(M,'area')) , err(meshNormals(M,'angle')) , err(Nr) ) );
  [np,nf] = chk( max(abs( sqrt(sum(Nr.^2,2)) - 1 )) < 1e-15 , 'unit rows' , np,nf );

  %---------------- equals the raw Max formula
  Vm = zeros( 400 , 3 );
  Tt = double( M.tri );
  for r = 1:3
    v  = Tt(:,r);  a = Tt(:,mod(r,3)+1);  b = Tt(:,mod(r+1,3)+1);
    e1 = P(a,:) - P(v,:);   e2 = P(b,:) - P(v,:);
    w  = cross( e1 , e2 ,2) ./ ( sum(e1.^2,2) .* sum(e2.^2,2) );
    Vm = Vm + [ accumarray(v,w(:,1),[400,1]) , accumarray(v,w(:,2),[400,1]) , accumarray(v,w(:,3),[400,1]) ];
  end
  Vm = Vm ./ sqrt( sum( Vm.^2 ,2) );
  [np,nf] = chk( max(abs( Nr - Vm ),[],'all') < 1e-12 , '== raw cross/(|e1|^2|e2|^2) Max formula' , np,nf );

  %---------------- winding flip -> exact negation
  Mf = M;  Mf.tri = Mf.tri(:,[1 3 2]);
  Dw = abs( meshNormals( Mf ,'reciprocal') + Nr );   %corner roles swap r-passes -> sum
  [np,nf] = chk( max( Dw(:) ) < 1e-12 , 'flipped winding -> -N (to summation-order ulps)' , np,nf );

  %---------------- 2D triangle
  M2 = struct( 'xyz',[0 0;1 0;0 1] , 'tri',[1 2 3] );
  [np,nf] = chk( isequal( meshNormals( M2 ,'reciprocal') , repmat([0 0 1],3,1) ) , '2D triangle -> [0 0 1] exact' , np,nf );

  %---------------- triNORMALS honored
  Mtn = M;  Mtn.triNORMALS = repmat( [1 0 0] , size(K,1) , 1 );
  [np,nf] = chk( isequal( meshNormals( Mtn ,'reciprocal') , repmat([1 0 0],400,1) ) , 'triNORMALS honored' , np,nf );

  %---------------- degenerate face poisons only its vertices
  Md = M;  Md.tri = [ Md.tri ; 1 2 2 ];
  Nd = meshNormals( Md , 'reciprocal' );
  others = true(400,1);  others([1 2]) = false;
  [np,nf] = chk( all(all(isnan( Nd([1 2],:) ))) && isequaln( Nd(others,:) , Nr(others,:) ) , ...
                 'degenerate face: its 2 vertices NaN, all others bit-identical' , np,nf );

  %---------------- celltype 10 -> error
  Mt10 = struct( 'xyz',[0 0 0;1 0 0;0 1 0;0 0 1] , 'tri',[1 2 3 4] , 'triNORMALS',[0 0 1] );
  ok = false;
  try, meshNormals( Mt10 , 'reciprocal' ); catch e, ok = strcmp( e.identifier , 'meshNormals:mode' ); end
  [np,nf] = chk( ok , 'celltype 10 -> meshNormals:mode error' , np,nf );

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
