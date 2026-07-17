function MS = meshSeparate( M , varargin )
%MESHSEPARATE  Split a mesh into pieces (connected components or by a label),
%              with an optional chain of piece-selection operations.
%
%   MS = meshSeparate( M )               split M into its CONNECTED COMPONENTS
%                                        (default: cells sharing a VERTEX join)
%   MS = meshSeparate( M , 'byedge' )    join cells only across a shared EDGE
%                                        (triangles & tetrahedra): parts touching
%                                        at a single VERTEX are split apart
%   MS = meshSeparate( M , 'byface' )    join cells only across a shared FACE
%                                        (tetrahedra only)
%   MS = meshSeparate( M , 'triFIELD' )  split by the per-face field M.triFIELD
%   MS = meshSeparate( M , labels )      split by a per-face label vector
%                                        (numel(labels) == number of faces)
%   MS = meshSeparate( CELLofMESHES )    separate each mesh and concatenate
%   MS = meshSeparate( ... , 'KeepNodes' )
%   OUT = meshSeparate( ... , OP1 , OP2 , ... )   chain the operations below
%
%   Without a field / label, components are found by CONNECTIVITY. The default
%   ('bynode') uses meshFacesConnectivity: cells sharing a NODE are in the same
%   piece. 'byedge' / 'byface' instead demand a shared EDGE (2 nodes) / FACE
%   (3 nodes), so blobs that meet only at a single vertex (or edge) are separated
%   -- 'byedge' works for triangles (5) and tets (10), 'byface' for tets only.
%   MS is a column CELL array of mesh structs; every carried field (xyz*, tri*,
%   ...) is propagated to each piece. By default each piece is passed through a
%   tidy step so it only keeps its own nodes; 'KeepNodes' skips that and keeps the
%   FULL node array in every piece (so node indices stay comparable across
%   pieces).
%
%   OPERATIONS (applied left to right after the split). Two kinds:
%
%   REDUCERS -- return a SINGLE mesh struct and MUST be the last operation:
%     'largest' / 'smallest'   piece with the most / fewest faces
%     'longest' / 'shortest'   piece with the largest / smallest total element
%                              size (sum of meshQuality(.,'size'))
%     'maxx'|'minx'|'maxy'|'miny'|'maxz'|'minz'
%                              piece reaching the extreme coordinate
%     'combine' / 'append'     glue every piece back with MeshAppend
%
%   TRANSFORMERS -- keep (and return) the CELL, so they can be chained:
%     'sort' / 'order' , KEY   sort pieces ASCENDING by KEY: 'ntri', 'nxyz', or
%                              a handle @(M)->scalar
%     'select' , W             keep pieces W: numeric indices (W<=0 counts from
%                              the end, 0 = last) or a predicate @(M)->logical
%     'remove' / 'delete' , P  drop the pieces where the predicate P(M) is true
%
%   NOTES
%   - Splitting by a per-NODE field is not implemented (only per-face fields /
%     labels, whose name must start with 'tri').
%   - The split is done in a SINGLE pass (local buildPieces): each piece only
%     touches its own nodes, so cost is O(nFaces+nNodes) rather than the old
%     O(nComponents x nNodes) of a per-piece MeshRemoveFaces+MeshTidy loop. This
%     matters a lot for meshes that fragment into many small pieces.
%
% See also meshFacesConnectivity, MeshRemoveFaces, MeshTidy, MeshAppend,
%          meshQuality, MeshFixCellOrientation.

% TODO
%  - PICK PIECES BY NODE: keep the piece(s) that CONTAIN a given node, addressed
%    in the ORIGINAL M numbering -- via a NEW op (e.g. 'atnode' / 'containing'),
%    NOT the current 'select' whose numeric argument means a PIECE index. Build a
%    node->piece lookup from fc (before buildPieces renumbers the nodes) so it
%    works with or without 'KeepNodes'.
%  - PICK PIECES BY PROXIMITY to a query point P = [x y z]: a reducer in the
%    maxx/minx family (e.g. 'nearest' , P) returning the piece with the minimum
%    Euclidean distance from P to its nodes/surface; optionally a K-nearest
%    variant returning the K closest pieces.

  if nargin == 1 && iscell( M )
    MS = {};
    for m = 1:numel( M )
      MS = [ MS ; meshSeparate( M{m} ) ];
    end

    return;
  end




  KeepNodes = false;
  try,[varargin,KeepNodes] = parseargs(varargin,'KeepNodes','$FORCE$',{true,KeepNodes}); end

  % connectivity criterion for the default (no field / label) separation:
  %   'bynode' (default) cells sharing a VERTEX are connected
  %   'byedge'           cells sharing an EDGE  are connected (tri & tet)
  %   'byface'           cells sharing a FACE   are connected (tet only)
  MODE = 'bynode';
  MODEKW = {'bynode','bynodes','byvertex','byvertices','byedge','byedges','byface','byfaces'};
  im = find( cellfun( @(a) ischar(a) && any(strcmpi(a,MODEKW)) , varargin ) , 1 );
  if ~isempty( im )
    m = lower( varargin{im} ); varargin(im) = [];
    if     any(strcmp(m,{'byedge','byedges'})), MODE = 'byedge';
    elseif any(strcmp(m,{'byface','byfaces'})), MODE = 'byface';
    else,                                       MODE = 'bynode';
    end
  end

  % ---- per-face label vector fc: one integer per face; pieces = its unique values
  if numel( varargin ) && ischar( varargin{1} ) && ...
     strncmp( varargin{1} , 'tri' , 3 ) && ...
     isfield( M , varargin{1} )

    fc = M.(varargin{1}); varargin(1) = [];
    [~,~,fc] = unique( fc );

  elseif numel( varargin ) && isnumeric( varargin{1} ) && ...
         numel( varargin{1} ) == size( M.tri , 1 )

    fc = varargin{1}; varargin(1) = [];
    [~,~,fc] = unique( fc );

  else

    switch MODE
      case 'bynode'
        fc = cellConnectivity( double(M.tri) , 1 );   %cells sharing a VERTEX (K=1)
      case 'byedge'
        if ~any( meshCelltype( M ) == [5 10] )
          error('meshSeparate:byedge','''byedge'' connectivity needs triangles (celltype 5) or tetrahedra (celltype 10).');
        end
        fc = cellConnectivity( double(M.tri) , 2 );    %cells sharing an EDGE (2 nodes)
      case 'byface'
        if meshCelltype( M ) ~= 10
          error('meshSeparate:byface','''byface'' connectivity is only defined for tetrahedra (celltype 10).');
        end
        fc = cellConnectivity( double(M.tri) , 3 );    %cells sharing a FACE (3 nodes)
    end

  end

  MS = buildPieces( M , fc , KeepNodes );

  while numel(varargin)
    op = varargin{1}; varargin(1) = [];
    switch lower( op )
      case {'maxx','minx','maxy','miny','maxz','minz'}
        if numel( varargin )
          error('no further options are allowed after minC/maxC.');
        end
        ord = zeros(size(MS));
        for s = 1:numel(MS)
          x = MS{s}.xyz( MS{s}.tri(:) ,:);
          switch lower( op )
            case 'minx', ord(s) =  min(x(:,1));
            case 'maxx', ord(s) = -max(x(:,1));
            case 'miny', ord(s) =  min(x(:,2));
            case 'maxy', ord(s) = -max(x(:,2));
            case 'minz', ord(s) =  min(x(:,3));
            case 'maxz', ord(s) = -max(x(:,3));
          end
        end
        [~,ord] = min( ord );
        MS = MS{ ord };
        return;

      case {'largest'}
        if numel( varargin )
          error('no further options are allowed after this (terminal) selector.');
        end
        ord = zeros(size(MS));
        for s = 1:numel(MS)
          ord(s) = size( MS{s}.tri ,1);
        end
        [~,ord] = max( ord );
        MS = MS{ ord };
        return;

      case {'smallest'}
        if numel( varargin )
          error('no further options are allowed after this (terminal) selector.');
        end
        ord = zeros(size(MS));
        for s = 1:numel(MS)
          ord(s) = size( MS{s}.tri ,1);
        end
        [~,ord] = min( ord );
        MS = MS{ ord };
        return;

      case {'longest'}
        if numel( varargin )
          error('no further options are allowed after this (terminal) selector.');
        end
        ord = zeros(size(MS));
        for s = 1:numel(MS)
          ord(s) = sum( meshQuality( MS{s} ,'size') );
        end
        [~,ord] = max( ord );
        MS = MS{ ord };
        return;

      case {'shortest'}
        if numel( varargin )
          error('no further options are allowed after this (terminal) selector.');
        end
        ord = zeros(size(MS));
        for s = 1:numel(MS)
          ord(s) = sum( meshQuality( MS{s} ,'size') );
        end
        [~,ord] = min( ord );
        MS = MS{ ord };
        return;

      case {'remove','delete'}
        w = varargin{1}; varargin(1) = [];
        if ~isa( w , 'function_handle' ), error('a predicate was expected'); end

        for s = 1:numel(MS)
          if feval( w , MS{s} ), MS{s} = []; end
        end
        MS( cellfun('isempty',MS) ) = [];

      case {'sort','order'}
        w = varargin{1}; varargin(1) = [];
        if ischar( w )
          switch lower(w)
            case {'ntri'}
              w = @(M)size(M.tri,1);
            case {'nxyz'}
              w = @(M)size(M.xyz,1);
            otherwise,error('not implemented ordering');
          end
        end
        if ~isa( w , 'function_handle' ), error('a scalar function was expected'); end

        ord = zeros(size(MS));
        for s = 1:numel(MS)
          ord(s) = feval( w , MS{s} );
        end
        [~,ord] = sort( ord );
        MS = MS(ord);

      case 'select'
        w = varargin{1}; varargin(1) = [];
        if isnumeric( w )
          w( w <= 0 ) = numel( MS ) + w( w <= 0 );
        elseif isa( w , 'function_handle')
          ww = false( numel(MS) ,1);
          for s = 1:numel(MS)
            ww(s) = feval( w , MS{s} );
          end
          if ~islogical( ww )
            error('function for the selection should return logicals');
          end
          w = find(ww);

        else
          error('unknown selection type');
        end

        MS = MS( w );

      case {'combine','append'}

        MS = MeshAppend( MS{:} );
        return;

      otherwise
        error('unknown operation');
    end
  end


end


function MS = buildPieces( M , fc , KeepNodes )
%split M into one sub-mesh per label in fc (one label per face). Unless KeepNodes,
%each piece drops the nodes it does not use and renumbers its .tri accordingly --
%a single-pass equivalent of the old per-piece MeshRemoveFaces + MeshTidy(NaN,
%false) loop, but touching only each piece's own nodes (O(nFaces+nNodes) total).
%Field handling mirrors those two: .tri is subset by face AND renumbered; other
%tri*/celltype fields are subset by face; xyz* fields are subset to the used
%nodes; any other field is copied unchanged. Non-finite nodes are rare and are
%handled by the exact old loop as a fallback.
  us = unique( fc(:) ).';

  if ~KeepNodes && ~all( isfinite( M.xyz(:) ) )
    MS = cell( numel(us) , 1 );
    for k = 1:numel(us)
      MS{k} = MeshTidy( MeshRemoveFaces( M , fc ~= us(k) ) , NaN , false );
    end
    return;
  end

  Fs       = fieldnames( M ).';
  classTRI = class( M.tri );
  hasCT    = isfield( M , 'celltype' ) && ~isscalar( M.celltype );
  MS       = cell( numel(us) , 1 );

  for k = 1:numel(us)
    F = find( fc == us(k) );

    if ~KeepNodes
      sub  = M.tri( F ,:);
      nz   = sub ~= 0;
      used = unique( sub(nz) );                      %sorted used node ids (0 = pad)
      remap = zeros( max([ used(:) ; 1 ]) , 1 );  remap( used ) = 1:numel(used);
    end

    P = struct();
    for f = Fs, f = f{1};
      if strcmp( f , 'tri' )
        if KeepNodes
          P.tri = M.tri( F ,:);
        else
          t = zeros( size(sub) );  t(nz) = remap( sub(nz) );
          P.tri = feval( classTRI , t );
        end
      elseif strncmp( f , 'tri' , 3 )                %per-face field
        v = M.(f); sz = size(v); sz(1) = numel(F);
        P.(f) = reshape( v( F ,:,:,:,:,:,:) , sz );
      elseif strcmp( f , 'celltype' ) && hasCT
        P.celltype = M.celltype( F ,: );
      elseif strncmp( f , 'xyz' , 3 )                %per-node field
        if KeepNodes, P.(f) = M.(f);
        else,         P.(f) = M.(f)( used ,:,:,:,:,:,: );
        end
      else
        P.(f) = M.(f);                               %anything else: copy as-is
      end
    end
    MS{k} = P;
  end
end


function comp = cellConnectivity( T , K )
%connected components of the cells (rows of T) where two cells are neighbours
%iff they share at least K vertices (K=2 -> a shared EDGE, K=3 -> a shared
%triangular FACE). Bipartite cell<->(sorted K-subset) graph + one conncomp; two
%cells that only touch at a single vertex (fewer than K shared nodes) fall in
%different components.
  nC  = size( T ,1);
  cmb = nchoosek( 1:size(T,2) , K );        %every K-subset of a cell's vertices
  nk  = size( cmb ,1);
  S   = zeros( nC*nk , K );
  cid = zeros( nC*nk , 1 );
  for i = 1:nk
    r = (i-1)*nC + (1:nC);
    S(r,:) = sort( T(:,cmb(i,:)) , 2 );      %the K-subset (as a sorted key)
    cid(r) = (1:nC).';                        %...belongs to this cell
  end
  [ ~ , ~ , sk ] = unique( S , 'rows' );      %shared-subset node id
  g  = graph( cid , nC + sk , [] , nC + max(sk) );   %bipartite: cells <-> subsets
  cc = conncomp( g );
  [ ~ , ~ , comp ] = unique( cc(1:nC).' );    %component label per cell (1..nComp)
end
