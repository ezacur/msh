function C = meshFacesCenter( M , d )
%MESHFACESCENTER  Center of every cell: the plain average of its vertices.
%
%   C = meshFacesCenter( M )       one row per M.tri row: the MEAN of the
%           cell's vertex coordinates. Works for ANY cell arity and keeps the
%           input's column count and class:
%             segments  (celltype  3) -> the midpoints            (exact)
%             triangles (celltype  5) -> the centroids            (exact)
%             tetras    (celltype 10) -> the barycenters          (exact)
%             quads     (celltype  9) -> the VERTEX average, which is NOT the
%                       area centroid unless the quad is a parallelogram
%                       (e.g. trapezoid (0,0),(4,0),(3,1),(1,1): vertex mean
%                       y = 0.5, area centroid y = 0.444).
%           (for simplices vertex mean == centroid; that is why the first
%           three are exact.)
%
%   C = meshFacesCenter( M , d )   return only the coordinate(s) d:
%           d = 1|2|3 or 'x'|'y'|'z' -> that single column;
%           any OTHER value is used verbatim as a column index, so d = [3 1]
%           returns the (z,x) columns -- but it is NOT validated: an unknown
%           char (say 'w') indexes column double('w') and dies with an
%           unrelated out-of-bounds error.
%
%   NaN vertices propagate to their cell's center. The accumulation loop runs
%   over the (few) columns of M.tri, not over cells: O(arity) vectorized
%   passes. Widely used as the per-face probe point (meshIsInterior,
%   MeshQuery, MeshAddTexture, plotMESH labels, g4Remesh, ...); note that
%   silhouette / backFaceCullingSplit carry their own LOCAL (V,F) variant,
%   which shadows this one inside those files only.
%
% See also meshF2V, meshNormals, meshQuality, meshVolume.

%   C = mean( permute( reshape( M.xyz( M.tri(:) , :) , [ size( M.tri , 1 ) , size(M.tri,2) , size( M.xyz,2) ] ) , [1 3 2] ) , 3 );

    C = 0;
    for c = 1:size( M.tri ,2)
        C = C + M.xyz( M.tri(:,c) ,:);
    end
    C = C / size( M.tri ,2);

    if nargin > 1
      if 0
      elseif numel(d) == 1 && ( d == 1 || ( ischar( d ) && strcmpi( d , 'x' ) ) ),  C = C(:,1);
      elseif numel(d) == 1 && ( d == 2 || ( ischar( d ) && strcmpi( d , 'y' ) ) ),  C = C(:,2);
      elseif numel(d) == 1 && ( d == 3 || ( ischar( d ) && strcmpi( d , 'z' ) ) ),  C = C(:,3);
      else, C = C(:,d);
%         error('Invalid coordinate.');
      end
    end


end
