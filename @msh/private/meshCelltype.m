function c = meshCelltype( M )
%MESHCELLTYPE  VTK-style cell-type code of a mesh, from its connectivity.
%
%   c = meshCelltype( M )
%
%   Returns the (scalar) cell type of M in the VTK numbering:
%       0  EMPTY_CELL     1  VERTEX (point)     3  LINE (segment)
%       5  TRIANGLE       9  QUAD              10  TETRA
%   The whole toolbox dispatches on this: 3 = polyline, 5 = triangle surface,
%   10 = tetrahedral volume.
%
%   With NO 'celltype' field the type is INFERRED from nodes-per-cell
%   (size(M.tri,2)) and spatial dimensions (size(M.xyz,2), default 3):
%       1 node  -> 1 (vertex)      2 nodes -> 3 (line)      3 nodes -> 5 (triangle)
%       4 nodes -> 9 (QUAD) if the mesh is 2D, else 10 (TETRA) in 3D
%   anything else (including a point cloud with no .tri) -> 0.
%   NOTE the 4-node case is ambiguous: a QUAD SURFACE embedded in 3D is reported
%   as TETRA (10); set an explicit M.celltype to override.
%
%   With a 'celltype' field it is HONORED and CHECKED against the nodes-per-cell
%   (e.g. celltype 5 requires 3 nodes/cell). A uniform per-face celltype VECTOR
%   collapses to that scalar; a genuinely MIXED vector is returned unchanged
%   (callers that switch on a scalar must handle that case).
%
% See also meshQuality, meshEdges, meshNormals, meshF2V.

  % VTK cell-type codes: 0 EMPTY, 1 VERTEX, 2 POLY_VERTEX, 3 LINE, 4 POLY_LINE,
  % 5 TRIANGLE, 6 TRIANGLE_STRIP, 7 POLYGON, 8 PIXEL, 9 QUAD, 10 TETRA

  nsd = 3;  %default NumberOfSpatialDimensions
  if isfield( M , 'xyz' ), nsd = size( M.xyz ,2); end
  nnf = []; %default NumberOfNodesPerFace
  if isfield( M , 'tri' ), nnf = size( M.tri ,2); end
  
  if isfield( M , 'celltype' )
    c = M.celltype;
    if ~isscalar(c) 
      if numel(c) ~= size( M.tri ,1)
        error('incorrect number of celltypes specified.');
      end
      if all( c == c(1) )
        c = c(1);
      end
    end
    if isscalar( c )
      switch c
        case 0,  error('celltype = 0?');
        case 1,  if ~isempty( M.tri ) && ~isempty(nnf) &&  nnf ~= 1, error('celltype is 1 (POINT), then NumberOfNodesPerFace should be 1.'); end
        case 2,  error('celltype = 2?');
        case 3,  if ~isempty( M.tri ) && ~isempty(nnf) &&  nnf ~= 2, error('celltype is 3 (LINE), then NumberOfNodesPerFace should be 2.'); end
        case 4,  error('celltype = 4?');
        case 5,  if ~isempty( M.tri ) && ~isempty(nnf) &&  nnf ~= 3, error('celltype is 5 (TRIANGLE), then NumberOfNodesPerFace should be 3.'); end
        case 6,  error('celltype = 6?');
        case 7,  error('celltype = 7?');
        case 8,  error('celltype = 8?');
        case 9,  if ~isempty( M.tri ) && ~isempty(nnf) &&  nnf ~= 4, error('celltype is 9 (QUAD), then NumberOfNodesPerFace should be 4.'); end
        case 10, if ~isempty( M.tri ) && ~isempty(nnf) &&  nnf ~= 4, error('celltype is 10 (TETRA), then NumberOfNodesPerFace should be 4.'); end
      end
    end
  else
    Mtype = [ nnf , nsd ];
    if     0
    elseif isequal( Mtype , [1 3] ) || isequal( Mtype , [1 2] )
      c = 1;
    elseif isequal( Mtype , [2 3] ) || isequal( Mtype , [2 2] )
      c = 3;
    elseif isequal( Mtype , [3 3] ) || isequal( Mtype , [3 2] )
      c = 5;
    elseif isequal( Mtype , [4 2] )
      c = 9;
    elseif isequal( Mtype , [4 3] )
      c = 10;
    else
      c = 0;
      %error('unknown type of mesh');
    end
  end

end
