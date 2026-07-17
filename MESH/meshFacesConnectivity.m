function C = meshFacesConnectivity( F )
%MESHFACESCONNECTIVITY  Label the node-connected components of a mesh's cells.
%
%   C = meshFacesConnectivity( M )   for a mesh struct M (or a plain
%       connectivity matrix), returns C (nCells-by-1): a component label per
%       cell. Two cells are in the SAME component when they share a NODE, and
%       transitively so (connected components of the share-a-vertex graph).
%       Labels are 1..K, numbered by the FIRST cell of each component (so the
%       component containing cell 1 is label 1, etc.).
%
%   PER-CELL CONSTRAINT: any per-cell field of M named 'tri*' (except 'tri' and
%       'triID') further SPLITS the components -- two cells are only joined when
%       they also carry the SAME value in every such field. (Constant columns
%       are ignored.) This is what lets meshSeparate keep, say, differently
%       labelled regions apart even where they touch.
%
%   Implementation: connected components of the bipartite cell<->(vertex,key)
%   graph via GRAPH / CONNCOMP -- fast and dependency-free (no VTK).
%
%   See also meshSeparate, conncomp, graph, meshEsuE.

  %----- gather the per-cell 'tri*' constraint fields (varying columns only)
  X = [];
  if isstruct( F )
    for f = fieldnames( F ).', f = f{1};
      if strcmp( f , 'tri' ) || strcmp( f , 'triID' ) || ~strncmp( f , 'tri' , 3 ), continue; end
      try
        thisX = F.(f);  thisX = thisX(:,:);
        thisX( : , all( bsxfun( @eq , thisX , thisX(1,:) ) , 1 ) ) = [];   %drop constant columns
        X = [ X , thisX ];   %#ok<AGROW>
      end
    end
    F = F.tri;
  end

  F  = double( F );
  nF = size( F , 1 );
  if nF == 0, C = zeros( 0 , 1 ); return; end
  nS = size( F , 2 );          %nodes per cell

  %----- a per-cell key from the constraint fields (all-1 when there is none)
  if isempty( X ) || size( X , 2 ) == 0
    xk = ones( nF , 1 );
  else
    [ ~ , ~ , xk ] = unique( X , 'rows' );
  end

  %----- bipartite graph: every cell links to a (vertex,key) node, so two cells
  %      meet only when they share a vertex AND carry the same key
  fi = repmat( ( 1:nF ).' , nS , 1 );                 %the cell of each incidence
  vv = F(:);                                          %its vertex
  xx = repmat( xk , nS , 1 );                         %its cell's key
  [ ~ , ~ , vn ] = unique( [ vv , xx ] , 'rows' );    %compact (vertex,key) node id
  nVN = max( vn );
  vn  = vn + nF;                                       %place after the nF cell nodes

  G = graph( [ fi ; vn ] , [ vn ; fi ] , [] , nF + nVN );

  comp = conncomp( G );                                %component id per graph node
  [ ~ , ~ , C ] = unique( comp( 1:nF ).' , 'stable' ); %relabel by first cell
  C = C(:);

end
