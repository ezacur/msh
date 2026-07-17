function L = meshPsuP( T , mode )
%MESHPSUP  Points surrounding points: the vertex-vertex (edge) adjacency.
%
%   L = meshPsuP( M )            SPARSE nV x nV logical adjacency (symmetric, no
%                               diagonal): L(i,j) = true iff vertices i,j share an
%                               EDGE.
%   L = meshPsuP( M , 'cell' )  nV x 1 CELL: L{i} = the neighbours of vertex i
%   L = meshPsuP( M , 'index' ) nV x maxValence padded matrix (NaN-filled) of the
%                               neighbours per vertex
%   L = meshPsuP( M , mode )    mode is 'sparse'|'cell'|'index' (or, for backward
%                               compatibility, false=sparse / true=cell)
%
%   M is a mesh struct (M.xyz sets the vertex count, M.tri the cells) or a plain
%   connectivity matrix T (then nV = max(T(:))). Edges come from meshEdges (each
%   undirected edge ONCE), so the neighbour lists carry no duplicates. This is the
%   graph adjacency behind Laplacian smoothing and one of the three fundamental
%   mesh adjacencies (Loehner's "surrounding" family).
%
% See also meshEsuP, meshEsuE, meshEdges.

  if nargin < 2 || isempty( mode ), mode = 'sparse'; end
  mode = parseMode( mode );

  if isstruct( T )
    nV = size( T.xyz ,1);
    T = T.tri;
  else
    nV = double( max( T(:) ) );
  end

  E = meshEdges( T );

  switch mode
    case 'sparse'
      L = sparse( double(E(:,1)) , double(E(:,2)) , true , nV , nV );
      L = L | L.';
    case 'cell'
      L = accumarray( [ E(:,1) ; E(:,2) ] , [ E(:,2) ; E(:,1) ] , [ nV , 1 ] , @(x){x} );
    case 'index'
      L = cell2index( accumarray( [ E(:,1) ; E(:,2) ] , [ E(:,2) ; E(:,1) ] , [ nV , 1 ] , @(x){x} ) );
  end

end

function m = parseMode( mode )
  if ( isnumeric(mode) || islogical(mode) ) && isscalar( mode )
    if mode, m = 'cell'; else, m = 'sparse'; end
    return;
  end
  if ~ischar( mode ), error('meshPsuP:mode','invalid mode'); end
  switch lower( mode )
    case {'sparse','s','sp'}, m = 'sparse';
    case {'cell','c'},        m = 'cell';
    case {'index','i'},       m = 'index';
    otherwise, error('meshPsuP:mode','invalid mode: use ''sparse'', ''cell'' or ''index''.');
  end
end

function I = cell2index( C )
% pack a cell of neighbour lists into a fixed-width matrix, NaN-padded.
  d = cellfun( 'prodofsize' , C );
  I = NaN( numel(C) , max([ 0 ; d(:) ]) );
  for n = 1:numel( C )
    if d(n), I( n , 1:d(n) ) = C{n}(:).'; end
  end
end
