function d = meshDistanceToBoundary( M , alg )
if 0
  M = Mesh( [0,1,0;0,-2,0;5,0,0] , [1,2,3] );
  M = MeshSubdivide( {M,3} );
  d = meshDistanceToBoundary( M ) - 0.1;


end



  if nargin < 2, alg = 'boundary'; end

  switch lower(alg)
    case {'dijstra','dijstrageodesic'}
      M.xyz = double( M.xyz );
      d = min( meshGeodesicDistance( M , find( meshBoundaryNodes( M ) ) , 'dijstra' ) , [] , 2 );

    case {'poisson','poissongeodesic'}
      M.xyz = double( M.xyz );
      d = min( meshGeodesicDistance( M , find( meshBoundaryNodes( M ) ) , 'poisson' ) , [] , 2 );

    case {'boundary','distance2boundary'}
      d = distanceFrom( M.xyz , MeshBoundary( M ) );

    otherwise
      error('Unknown algorithm.');
  end

end
