function VP = meshPrctile( M , V , n )
if 0


M = Mesh( ndmat(linspace(0,2,250),linspace(-1,1.5,350)) , 'delaunay');

r = fro( M.xyz ,2);
P = meshPrctile( M , r , 101 );
P = griddedInterpolant( P(:,1) , P(:,2) , 'pchip','nearest');
D = @(x) arrayfun( @(x)NumericalDiff( @(z)P(z) ,x, 'c' ) , x );

x = linspace( -0.2 , max(r)*1.2 , 100);
hist( r , 100 );
hplot( x , 450*D(x), '-2r')


%%
end



  if nargin < 3, n = 25; end

  if ~ismatrix( V ) || size( V ,2) ~= 1
    error('A one-column array is expected as V.');
  end

  if isscalar( n )
    Vs = linspace( min( V(:) ) , max( V(:) ) , n );
  else
    Vs = n;
  end
  VP(:,1) = Vs(:);
  for v = 1:numel( Vs )
    VP(v,2) = meshSurface( MeshClip( M , V - Vs(v) ) );
  end




end