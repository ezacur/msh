function bench_BVH
%BENCH_BVH  Speed & accuracy shoot-out:
%   A) bvhClosestElement (all 6 cages) vs vtkClosestElement (VTK cell locator)
%   B) bvhIntersectRay vs IntersectSurfaceRay_mx (specialized ray MEX)
%
%  Ground truth for A: OUR brute force (single-leaf blob = exact double
%  Ericson over every element) -- so vtk is also being validated against an
%  independent implementation. Accuracy metrics:
%    max|dd|   worst distance deviation from the ground truth
%    bc resid  worst barycentric reconstruction residual |sum(bc*verts) - cp|
%    e~ %      element-id disagreements (equidistant ties are legitimate)
%  Speed: build (locator/blob) and per-query cost, SEPARATED.

  rng(21);
  addpath( fullfile( fileparts( mfilename('fullpath') ) , '..' , 'tools' ) );
  addpath( fullfile( fileparts( mfilename('fullpath') ) , '..' , 'MESH' ) );
  %the VTK 5 DLLs the old mex links against live inside the repo:
  vtkbin = fullfile( fileparts( mfilename('fullpath') ) , '..' , 'MESH' , 'vtk' );
  if ~contains( getenv('PATH') , vtkbin )
    setenv( 'PATH' , [ vtkbin , ';' , getenv('PATH') ] );
  end

  fprintf( '\n==================== A) CLOSEST ELEMENT ====================\n' );
  for nV = [ 2600 , 26000 ]
    V = randn( nV ,3);  V = V ./ sqrt( sum( V.^2 ,2) );
    M = struct( 'xyz' , V , 'tri' , convhulln( V ) );
    nP = 20000;
    w  = randi( nV , nP ,1);
    Pnear = V(w,:) .* ( 1 + 0.02*randn( nP ,1) );
    Pfar  = randn( nP ,3) * 1.5;

    for PP = { {Pnear,'near-surface'} , {Pfar,'far'} }
      P = PP{1}{1};  tag = PP{1}{2};
      [ ~ , ~ , dT ] = bvhClosestElement( M , P , BVH( M , Inf ) );   %ground truth
      [ eT , ~ , ~ ] = bvhClosestElement( M , P );                       %tie reference

      fprintf( '\n--- %d tris, %d pts (%s) ---\n' , size(M.tri,1) , nP , tag );
      fprintf( '%-10s %10s %10s %12s %12s %8s %8s\n' , ...
               'engine' , 'build ms' , 'us/pt' , 'max|dd|' , 'bc resid' , 'tie %' , 'BAD e %' );

      %vtk cell locator (persistent form: build once, query many)
      try
        vtkClosestElement( [] , [] );
        tic;  vtkClosestElement( M );                 tB = toc;
        tic;  [ e , cp , d ] = vtkClosestElement( P );  tQ = toc;
        [ ~ , cp4 , ~ , bc ] = vtkClosestElement( P );
        vtkClosestElement( [] , [] );
        idx = M.tri( e ,:);
        rec = bc(:,1).*V(idx(:,1),:) + bc(:,2).*V(idx(:,2),:) + bc(:,3).*V(idx(:,3),:);
        dif = e ~= eT;
        fprintf( '%-10s %10.1f %10.3f %12.2e %12.2e %8.2f %8.2f\n' , 'vtk' , ...
                 1e3*tB , 1e6*tQ/nP , max(abs(d-dT)) , max(max(abs(rec-cp4))) , ...
                 100*mean( dif & abs(d-dT) <= 1e-12 ) , 100*mean( dif & abs(d-dT) > 1e-12 ) );
      catch ME
        fprintf( '%-10s UNAVAILABLE: %s\n' , 'vtk' , ME.message );
      end

      %ours, all six cages
      for vv = { 'aabb' , 'sphere' , 'obb' , 'kdop' , 'rss' , 'lss' }
        tic;  B = BVH( M , [] , vv{1} );             tB = toc;
        tic;  [ e , cp , d ] = bvhClosestElement( {M,B} , P );  tQ = toc;
        [ ~ , ~ , ~ , bc ] = bvhClosestElement( {M,B} , P );
        idx = M.tri( e ,:);
        rec = bc(:,1).*V(idx(:,1),:) + bc(:,2).*V(idx(:,2),:) + bc(:,3).*V(idx(:,3),:);
        dif = e ~= eT;
        fprintf( '%-10s %10.1f %10.3f %12.2e %12.2e %8.2f %8.2f\n' , ['msh-',vv{1}] , ...
                 1e3*tB , 1e6*tQ/nP , max(abs(d-dT)) , max(max(abs(rec-cp))) , ...
                 100*mean( dif & abs(d-dT) <= 1e-12 ) , 100*mean( dif & abs(d-dT) > 1e-12 ) );
      end
    end
  end

  %sliver micro-case: needle triangle, analytic reference
  fprintf( '\n--- sliver accuracy (needle triangle, aspect 1e8) ---\n' );
  Msl = struct( 'xyz' , [0 0 0 ; 1 0 0 ; 0.5 1e-8 0] , 'tri' , [1 2 3] );
  Psl = [ 0.5 , 0.3 , 0.4 ];
  dRef = norm( Psl - [0.5 , 1e-8 , 0] );         %the apex vertex is the winner
  try
    [ ~ , ~ , d , bc ] = vtkClosestElement( Msl , Psl );
    fprintf( 'vtk       : d err %.2e | bc = [%g %g %g] (sum %g)\n' , ...
             abs(d-dRef) , bc , sum(bc) );
  catch ME
    fprintf( 'vtk       : UNAVAILABLE (%s)\n' , ME.message );
  end
  [ ~ , ~ , d , bc ] = bvhClosestElement( Msl , Psl );
  fprintf( 'msh       : d err %.2e | bc = [%g %g %g] (sum %g)\n' , ...
           abs(d-dRef) , bc(1:3) , sum(bc(1:3)) );

  fprintf( '\n==================== B) RAYS ================================\n' );
  V = randn( 26000 ,3);  V = V ./ sqrt( sum( V.^2 ,2) );
  M = struct( 'xyz' , V , 'tri' , convhulln( V ) );
  nR = 20000;
  rays = [ randn(nR,3)*3 , randn(nR,3)*0.3 ];
  B = BVH( M );
  IntersectSurfaceRay( M , rays(1:100,:) , 'first' );        %warm the oracle LRU

  for MODE = { 'first' , 'any' }
    tic;  [ ~ , c1 , t1 ] = bvhIntersectRay( {M,B} , rays , MODE{1} );  tO = toc;
    tic;  [ ~ , ~ , c2 , t2 ] = IntersectSurfaceRay( M , rays , MODE{1} );  tR = toc;
    h1 = c1 > 0;  h2 = c2 > 0;
    if strcmp( MODE{1} , 'any' )
      dt = NaN;                                  %'any' may pick different hits
    else
      dt = max( abs( t1(h1&h2) - t2(h1&h2) ) );
    end
    fprintf( '%-6s  msh %6.3f vs isr %6.3f us/ray | hit mismatch %d/%d | max|dt| %.2e | cell~ %.2f%%\n' , ...
             MODE{1} , 1e6*tO/nR , 1e6*tR/nR , sum( h1 ~= h2 ) , nR , dt , ...
             100*mean( c1(h1&h2) ~= c2(h1&h2) ) );
  end

  %'all': the full hit multiset must match
  tic;  [ ~ , ~ , t1 , r1 ] = bvhIntersectRay( {M,B} , rays , 'all' );  tO = toc;
  tic;  [ ~ , ~ , ~ , t2 , r2 ] = IntersectSurfaceRay( M , rays , 'all' );  tR = toc;
  n1 = accumarray( r1 , 1 , [nR,1] );  n2 = accumarray( r2 , 1 , [nR,1] );
  if isequal( n1 , n2 )
    dt = max( abs( sort(t1) - sort(t2) ) );      %same per-ray counts + both sorted
    fprintf( '%-6s  msh %6.3f vs isr %6.3f us/ray | %d hits both | max|dt| %.2e\n' , ...
             'all' , 1e6*tO/nR , 1e6*tR/nR , numel(t1) , dt );
  else
    fprintf( '%-6s  HIT-COUNT MISMATCH on %d rays\n' , 'all' , sum( n1 ~= n2 ) );
  end

  fprintf( '\ndone.\n' );
end
