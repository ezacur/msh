function X = meshFarthestPointSampling( M , IDS , minD , maxN , maxSubs , VERBOSE )

  if nargin < 2 || isempty( IDS  ), IDS = 1; end
  if nargin < 3 || isempty( minD ), minD = 0;   end
  if nargin < 4 || isempty( maxN ), maxN = Inf; end
  if nargin < 5 || isempty( maxSubs ), maxSubs = Inf; end
  if nargin < 6 || isempty( VERBOSE ), VERBOSE = false; end

  if iscell( M )
    X = cell( size( M ) );
    for m = 1:numel(M)
      X{m} = meshFarthestPointSampling( M{m} , IDS , minD , maxN , maxSubs , VERBOSE );
    end
    return;
  end


  if any( IDS < 1 )
    error('initial IDS should be all indexes (greater than zero)');
  end
  if any( IDS > size( M.xyz ,1) )
    error('initial IDS should be all valid indexes (smaller than number of nodes)');
  end
%   IDS = M.xyz( IDS(:) ,:);

  M = Mesh( M ,0);
  M = MeshTidy( M ,0,true);
  s = 0;
  while size( M.xyz ,1) < 1e6  &&  s < maxSubs
  	M = MeshSubdivide( M );
    s = s + 1;
  end
%   IDS = vtkClosestPoint( struct('xyz',double(M.xyz)) , double( IDS ) );


  X = FarthestPointSampling( M.xyz , IDS , minD , maxN , [] , VERBOSE );
  
end
