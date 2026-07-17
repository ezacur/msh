function M = MeshBoundary( M , tidyAfter )
%MESHBOUNDARY  The boundary of a mesh: the facets belonging to a single cell.
%
%   B = MeshBoundary( M )              the boundary of mesh M. A facet (a cell's
%       sub-entity one dimension lower) is on the boundary when it belongs to
%       EXACTLY ONE cell; the boundary is the set of those facets, with each
%       cell's ORIENTATION preserved. The boundary drops one dimension:
%           tetrahedra (celltype 10) -> boundary TRIANGLES (5)
%           triangles  (celltype  5) -> boundary EDGES     (3)
%           segments   (celltype  3) -> boundary END POINTS(1)  (the degree-1
%                                       nodes, as a 1-column M.tri)
%       A closed mesh has an empty boundary. Per-node data (M.xyz and any
%       'xyz*' field) is UNCHANGED -- the boundary keeps the full node list, so
%       nodes not on the boundary are left dangling; pass tidyAfter=true to
%       compact them (via MeshTidy). Per-cell 'tri*' fields are inherited by
%       each boundary facet from its parent cell.
%   B = MeshBoundary( M , true )       also MeshTidy the result (drop the
%                                       unused nodes, renumber).
%   C = MeshBoundary( TRI )            with a plain connectivity matrix TRI
%                                       (not a struct), returns the boundary
%                                       connectivity matrix directly.
%
%   The output is a mesh struct with no 'celltype' field (its type is inferred
%   from the facet width). Unsupported celltypes error (MeshBoundary:celltype).
%
%   LIMITATION: degenerate cells (a repeated vertex, e.g. tri = [1 1 2]) produce
%   spurious facets and a wrong boundary -- clean the mesh first.
%
%   Examples:
%       MeshBoundary( struct('xyz',[0 0;1 0;1 1],'tri',[1 2 3]) ).tri  % 3 edges
%       B = MeshBoundary( M , true );                                  % + tidy
%
%   See also meshCelltype, MeshTidy, meshEsuE.

  if nargin < 2, tidyAfter = false; end

  asMESH = true;
  if ~isstruct( M )
    asMESH = false;
    M = struct('tri',M);
  end
%   if ~isfield( M , 'xyz' )
%     M.xyz = zeros( 0 , 3 );
%   end
  
  M.celltype = meshCelltype( M );

  
  switch M.celltype
    case 3
      allF = [ M.tri(:,1) ; M.tri(:,2) ];
      IDXs = repmat( (1:size(M.tri,1)).' ,[1 2] ); IDXs = IDXs(:);
      F    = allF;
      [u,~,c] = unique( F );
      c = accumarray( c(:) , 1);
      u = u( c == 1 ,:);
      w = ismember( F , u );
      IDXs  = IDXs(w);
      M.tri = allF( w ,:);

      M.celltype = 1;          %endpoints of a polyline are POINTS (1 node/cell)

    case 5
      allF = [ M.tri(:,[1,2]) ; M.tri(:,[2,3]) ; M.tri(:,[3,1]) ];
      IDXs = repmat( (1:size(M.tri,1)).' ,[1 3] ); IDXs = IDXs(:);
      F    = sort( allF ,2);
      [u,~,c] = unique( F , 'rows' );
      c = accumarray( c(:) , 1);
      u = u( c == 1 ,:);
      w = ismember( F , u , 'rows' );
      IDXs  = IDXs(w);
      M.tri = allF( w ,:);

      M.celltype = 3;
      
      
    case 10
      allF = [ M.tri(:,[2,3,4]) ; M.tri(:,[4,3,1]) ; M.tri(:,[1,2,4]) ; M.tri(:,[1,3,2]) ];
      IDXs = repmat( (1:size(M.tri,1)).' ,[1 4] ); IDXs = IDXs(:);
      F    = sort( allF ,2);
      [u,~,c] = unique( F , 'rows' );
      c = accumarray( c(:) , 1);
      u = u( c == 1 ,:);
      w = ismember( F , u , 'rows' );
      IDXs  = IDXs(w);
      M.tri = allF( w ,:);

      M.celltype = 5;

    otherwise
      error( 'MeshBoundary:celltype' , ...
             'unsupported celltype %g (only segments=3, triangles=5, tetrahedra=10).' , M.celltype );
  end
  
  if asMESH

    for f = fieldnames(M).'
      if ~strncmp( f{1} , 'tri' , 3 ) || strcmp( f{1} , 'tri' ), continue; end
      M.( f{1} ) = M.( f{1} )( IDXs ,:,:,:,:,:,:);
    end
    M = rmfield( M ,'celltype' );   %let the type be re-inferred from the facet width

    if tidyAfter
      M = MeshTidy( M ,NaN);
    end
  
  else
    
    M = M.tri;
  
  end

end
