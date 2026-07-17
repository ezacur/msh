function [E,C,A] = meshCellsContact( M )
%MESHCELLSCONTACT  Unique edges of a triangle mesh, the cells touching each
%                  edge, and the inter-face angle across it.
%
%   [E,C,A] = meshCellsContact( M )
%
%   E : (nE x 2) unique UNDIRECTED edges (sorted node pairs, sorted rows).
%   C : (nE x 1) cell; C{k} = the faces (rows of M.tri) containing edge k.
%   A : (nE x 1) angle in DEGREES between the normals of the two faces of
%       edge k, computed with the chord formula 2*asind(|N1-N2|/2) (exact and
%       stable near 0). A = 0 for coplanar faces.
%         boundary edge (1 face)      -> A = NaN
%         non-manifold edge (>2 faces) -> A = Inf
%       NOTE: A depends on a CONSISTENT face orientation (a flipped face reads
%       as 180-angle); run MeshFixCellOrientation first if unsure.
%
%   Only triangle meshes (celltype 5). Used by MeshSplit's crease-threshold
%   ( MeshSplit(M,-ang) ) and 'nonmanifold' forms.
%
% See also MeshSplit, meshEdges, meshNormals, MeshFixCellOrientation.

  if meshCelltype( M ) ~= 5
    error('meshCellsContact:celltype','only triangle meshes (celltype 5) are supported (got %d).', meshCelltype( M ) );
  end

  C = ( 1:size( M.tri ,1) ).';
  E = [ M.tri(:,[1,2]) ;...
        M.tri(:,[2,3]) ;...
        M.tri(:,[1,3]) ];
  C = repmat( C , 3,1);

  E = sort( E ,2);
  [E,~,c] = unique( E , 'rows' );

  if nargout > 1
    C = accumarray( c , C , [], @(x){x(:).'} );
  end

  if nargout > 2

    A = NaN( size(E,1) ,1);

    n = cellfun('prodofsize',C);
    A( n > 2 ) = Inf;

    w = n == 2;
    if any( w )                     %without this, a mesh with NO 2-manifold edge
      PC = cell2mat( C(w) );        %(single triangle, pure non-manifold fan)
                                    %crashed at PC(:,1) on the empty 0x0 PC
      N = meshNormals( M );
      A(w) = 2*asind( fro( N( PC(:,1) ,:) - N( PC(:,2) ,:) ,2)/2 );
    end

  end

end
