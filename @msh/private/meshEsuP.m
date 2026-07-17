function L = meshEsuP( T , mode )
%MESHESUP  Elements surrounding points: which cells are incident to each vertex.
%
%   L = meshEsuP( M )            SPARSE nT x nV logical incidence: L(t,v) = true
%                               iff cell t contains vertex v (so L(:,v) selects
%                               the cells around vertex v).
%   L = meshEsuP( M , 'cell' )  nV x 1 CELL: L{v} = the cells incident to vertex v
%   L = meshEsuP( M , 'index' ) nV x maxAdj padded matrix (NaN-filled) of the
%                               incident cells per vertex
%   L = meshEsuP( M , mode )    mode is 'sparse'|'cell'|'index' (or, for backward
%                               compatibility, false=sparse / true=cell)
%
%   M is a mesh struct (M.xyz sets the vertex count, M.tri the cells) or a plain
%   connectivity matrix T (then nV = max(T(:))). One of the three fundamental
%   mesh adjacencies (Loehner's "surrounding" family).
%
% See also meshPsuP, meshEsuE, meshEdges.

  if nargin < 2 || isempty( mode ), mode = 'sparse'; end
  mode = parseMode( mode );

  if isstruct( T )
    nV = size( T.xyz ,1);
    T = T.tri;
  else
    nV = double( max( T(:) ) );
  end

  nT   = size( T ,1);
  Tids = repmat( ( 1:nT ).' , [ size( T ,2) , 1 ] );

  switch mode
    case 'sparse'
      L = sparse( Tids , double( T(:) ) , true , nT , nV );
    case 'cell'
      L = accumarray( double( T(:) ) , Tids , [ nV , 1 ] , @(x){x} );
    case 'index'
      L = cell2index( accumarray( double( T(:) ) , Tids , [ nV , 1 ] , @(x){x} ) );
  end

end

function m = parseMode( mode )
  if ( isnumeric(mode) || islogical(mode) ) && isscalar( mode )
    if mode, m = 'cell'; else, m = 'sparse'; end
    return;
  end
  if ~ischar( mode ), error('meshEsuP:mode','invalid mode'); end
  switch lower( mode )
    case {'sparse','s','sp'}, m = 'sparse';
    case {'cell','c'},        m = 'cell';
    case {'index','i'},       m = 'index';
    otherwise, error('meshEsuP:mode','invalid mode: use ''sparse'', ''cell'' or ''index''.');
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
