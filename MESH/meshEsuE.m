function L = meshEsuE( T , mode , byMODE )
%MESHESUE  Elements surrounding elements: the cell-cell adjacency.
%
%   L = meshEsuE( M )                 SPARSE nT x nT logical adjacency: L(s,t) is
%                                     true iff cells s,t are neighbours.
%   L = meshEsuE( M , 'cell' )        nT x 1 CELL: L{t} = the neighbours of cell t
%   L = meshEsuE( M , 'index' )       nT x maxAdj padded matrix (NaN-filled)
%   L = meshEsuE( M , mode , byMODE ) choose the OUTPUT (mode = 'sparse'|'cell'|
%                                     'index', or false=sparse / true=cell) and
%                                     the NEIGHBOUR criterion byMODE:
%       'bynode' (default)  cells sharing a VERTEX are neighbours
%       'byedge'            cells sharing an EDGE (2 nodes) -- triangles & tets
%       'byface'            cells sharing a FACE (3 nodes) -- tetrahedra only
%     byMODE also accepts the shorthands 'n'/'v', 'e', 'f' or 1/2/3.
%
%   M is a mesh struct or a plain connectivity matrix T. One of the three
%   fundamental mesh adjacencies (Loehner's "surrounding" family). 'byedge' /
%   'byface' separate blobs that only touch at a single vertex / edge, which
%   'bynode' keeps joined (see meshSeparate).
%
%   'byedge' / 'byface' build the adjacency in O(n log n): the K-vertex
%   sub-entities are grouped and every within-group cell pair is emitted at once
%   (no per-entity loop).
%
% See also meshPsuP, meshEsuP, meshSeparate, meshEdges.

  if nargin < 2 || isempty( mode ),   mode   = 'sparse'; end
  if nargin < 3 || isempty( byMODE ), byMODE = 'bynode'; end
  mode = parseMode( mode );

  if isstruct( T ), T = T.tri; end
  nT = size( T ,1);

  if isnumeric( byMODE )
    switch byMODE
      case 1, byMODE = 'bynode';
      case 2, byMODE = 'byedge';
      case 3, byMODE = 'byface';
      otherwise, error('incorrect byMODE option');
    end
  elseif numel( byMODE ) == 1
    switch lower( byMODE )
      case {'n','v'}, byMODE = 'bynode';
      case 'e',       byMODE = 'byedge';
      case 'f',       byMODE = 'byface';
      otherwise, error('no well specified byMODE option');
    end
  end

  switch lower( byMODE )
    case {'bynode','node','bynodes','nodes','byvertices','vertices','byvertice','vertice'}
      ESUP = meshEsuP( T , 'sparse' );

      F = sparse( [] , [] , [] , 0 , nT );
      for c = 1:size(T,2)
        F = [ F ; ESUP( : , T(:,c) ) ];
      end

      [ I , J ] = find( F );
      I = rem( I-1 , nT ) + 1;

      if strcmp( mode , 'sparse' )
        w = J <= I;  J(w) = []; I(w) = [];
        IJ = unique( [ J , I ] , 'rows' );
        L = sparse( IJ(:,1) , IJ(:,2) , true , nT , nT );
        L = L | L.';
      else
        w = J == I;  J(w) = []; I(w) = [];
        IJ = unique( [ J , I ] , 'rows' );
        C = accumarray( IJ(:,1) , IJ(:,2) , [nT,1] , @(x){x} );
        if strcmp( mode , 'cell' ), L = C; else, L = cell2index( C ); end
      end

    case {'byedge','edge','byedges','edges'}
      if size( T , 2 ) < 2, error('meshEsuE:byedge','''byedge'' needs cells with at least 2 vertices.'); end
      [ c , fid ] = entityGroups( T , 2 );
      L = adjFromGroups( c , fid , nT , mode );

    case {'byface','face','byfaces','faces'}
      if size( T , 2 ) < 3
        error('meshEsuE:byface','''byface'' adjacency needs cells with at least 3 vertices (use it for tetrahedra).');
      end
      [ c , fid ] = entityGroups( T , 3 );
      L = adjFromGroups( c , fid , nT , mode );

    otherwise
      error('unknow byMODE');
  end

end

function [ c , fid ] = entityGroups( T , K )
% for every cell, its K-vertex sub-entities (K=2 edges, K=3 faces): a shared-
% entity id c per (cell,entity) row and the owning cell id fid.
  nT  = size( T ,1);
  cmb = nchoosek( 1:size(T,2) , K );
  E   = zeros( nT*size(cmb,1) , K );
  fid = zeros( nT*size(cmb,1) , 1 );
  for i = 1:size(cmb,1)
    r = (i-1)*nT + (1:nT);
    E( r ,:) = T( : , cmb(i,:) );
    fid( r ) = ( 1:nT ).';
  end
  [ ~ , ~ , c ] = unique( sort( E , 2 ) , 'rows' );
end

function L = adjFromGroups( c , fid , nT , mode )
% cell-cell adjacency from shared-entity groups: two cells are neighbours iff
% they carry the same entity id. All within-group pairs are produced by pairing
% rows d apart after sorting by entity id (d = 1..maxGroupSize-1) -- no per-entity
% loop, so it is O(n log n).
  [ cs , ord ] = sort( c(:) );
  fs = fid( ord );
  n  = numel( cs );
  src = zeros(0,1); dst = zeros(0,1);
  d = 1;
  while d < n
    lo   = ( 1 : n-d ).';
    same = cs(lo) == cs(lo+d);
    if ~any( same ), break; end
    lo  = lo( same );
    src = [ src ; fs(lo)   ];
    dst = [ dst ; fs(lo+d) ];
    d   = d + 1;
  end

  switch mode
    case 'sparse'
      if isempty( src ), L = logical( sparse( nT , nT ) );
      else,              L = logical( sparse( [src;dst] , [dst;src] , 1 , nT , nT ) );
      end
    otherwise
      if isempty( src ), C = cell( nT , 1 );
      else,              C = accumarray( [src;dst] , [dst;src] , [nT,1] , @(x){ unique(x(:)).' } );
      end
      if strcmp( mode , 'cell' ), L = C; else, L = cell2index( C ); end
  end
end

function m = parseMode( mode )
  if ( isnumeric(mode) || islogical(mode) ) && isscalar( mode )
    if mode, m = 'cell'; else, m = 'sparse'; end
    return;
  end
  if ~ischar( mode ), error('meshEsuE:mode','invalid mode'); end
  switch lower( mode )
    case {'sparse','s','sp'}, m = 'sparse';
    case {'cell','c'},        m = 'cell';
    case {'index','i'},       m = 'index';
    otherwise, error('meshEsuE:mode','invalid mode: use ''sparse'', ''cell'' or ''index''.');
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
