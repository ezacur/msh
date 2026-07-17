function V = meshF2V( M , F , mode )
%MESHF2V  Transfer / aggregate per-FACE data onto the VERTICES.
%
%   V = meshF2V( M , F )         area-weighted (default) face->vertex transfer
%   V = meshF2V( M , F , mode )  with the weighting / accumulator MODE below
%
%   For every vertex, V gathers the values F of its incident faces, combined
%   according to MODE. F gives one value (row) per face and may be:
%     * a numeric array  nF x k (or nF x ...) of per-face values
%     * a char NAME      -> uses M.(['tri' NAME]) if present, else M.(NAME)
%     * a function handle-> F = fun(M)  (or fun(M.tri))
%   V has the same trailing size as F but nV rows. Mirror of meshV2F.
%
%   MODE (a string):
%     'sum'                 plain sum of the incident face values
%     'mean' / 'average'    simple average (over the number of incident faces)
%     'area'   (default)    area-weighted SUM      (triangles, celltype 5)
%     'narea'               area-weighted AVERAGE  (divided by the sum of areas)
%     'length' / 'nlength'  length-weighted sum / average (segments, celltype 3)
%     'volume' / 'nvolume'  volume-weighted sum / average (tetrahedra, celltype 10)
%     'angles'              incident-ANGLE-weighted sum (triangles): a face's
%                           weight at a vertex is its corner angle there -- the
%                           pseudonormal weighting used by meshNormals.
%   The 'n...' variants divide by the total weight (a weighted average); the rest
%   are weighted sums. Each weighting checks the matching celltype.
%
%   MODE may instead be a FUNCTION HANDLE, applied per vertex over its incident
%   face values via accumarray (e.g. @max, @min, @median, @(x)...).
%
% See also meshV2F, meshQuality, meshNormals.

  if nargin < 3, mode = 'area'; end
  
  nF = size( M.tri ,1);
  nS = size( M.tri ,2);
  nV = size( M.xyz ,1);
 

  if isa( F , 'function_handle' )
    try, F = feval( F , M ); catch
    try, F = feval( F , M.tri ); catch
      error('invalid function to evaluate on mesh');
    end; end
  elseif ischar( F )
    try, F = M.(['tri',F]); catch
    try, F = M.(F); catch
      error('invalid attribute name.');
    end; end
  end
  if size( F ,1) ~= nF
    error('invalid per-faces-field');
  end

  sz = size( F );
  sz(1) = nV;
  F = F(:,:);

  if nF == 0
    %no faces: every vertex gathers nothing -> zeros (callers like meshNormals
    %normalize that to NaN). The weighted paths below would die on the EMPTY
    %meshQuality weights otherwise (the isempty(W) branch never assigned V).
    V = zeros( sz );
    return;
  end

  if ischar( mode )

    N = []; W = [];
    switch lower(mode)
      case {'s','sum'}
        W = 1;
        N = 1;

      case {'m','mean','average'}
        W = 1;
        N = accumarray(  double( M.tri(:) ) , 1 );

      case {'l','length'}
        if meshCelltype(M) ~= 3, error('invalid weighting for the celltype'); end
        W = meshQuality( M , 'length' );
        N = 1;

      case {'nl','nlength','normalizedl','normalizedlength'}
        if meshCelltype(M) ~= 3, error('invalid weighting for the celltype'); end
        W = meshQuality( M , 'length' );

      case {'a','area'}
        if meshCelltype(M) ~= 5, error('invalid weighting for the celltype'); end
        W = meshQuality( M , 'area' );
        N = 1;

      case {'na','narea','normalizeda','normalizedarea'}
        if meshCelltype(M) ~= 5, error('invalid weighting for the celltype'); end
        W = meshQuality( M , 'area' );

      case {'v','vol','volume'}
        if meshCelltype(M) ~= 10, error('invalid weighting for the celltype'); end
        W = meshQuality( M , 'volume' );
        N = 1;

      case {'nv','nvol','nvolume','normalizedv','normalizedvol','normalizedvolume'}
        if meshCelltype(M) ~= 10, error('invalid weighting for the celltype'); end
        W = meshQuality( M , 'volume' );

      case {'g','angles'}
        if meshCelltype(M) ~= 5, error('invalid weighting for the celltype'); end
        W = meshQuality( M , 'angles' );
        N = 1;

      otherwise, error('invalid weighting option');
    end

    if     isempty( W )
    elseif issparse( W ) && size( W ,2) == nF && size( W ,1) == nV
      V = W * F;
      if isempty(N), N = sum( W ,2); end
    elseif isscalar( W )
      Tid = ( 1:nF ).';
      W = sparse( double( M.tri ) , repmat( Tid ,nS,1) , W ,nV,nF);
      V = W * F;
      if isempty(N), N = sum( W ,2); end
    elseif size( W ,1) == nF
      Tid = ( 1:nF ).';
      W = repmat( double(W) ,1,nS/size(W,2));
      W = sparse( double(M.tri) , repmat( double(Tid) ,nS,1) , W ,nV,nF);
      V = W * double( F );
      if isempty(N), N = sum( W ,2); end
    else
      error('incorrect W matrix');
    end
    
    %if required, normalized by N
    if 0
    elseif isscalar( N ) && N == 1
    elseif size( N ,1) == nV && size( N ,2) == 2 && all( N == 1 )
    elseif size( N ,1) == nV && size( N ,2) == 1
      V = bsxfun( @rdivide , V , N );
    elseif isscalar( N )
      V = V * (1/N);
    else
      error('incorrect normalization step');
    end
    
  elseif isa( mode , 'function_handle' )
    
    nC = size( F ,2);
    V = NaN( nV , nC );
    for c = 1:nC
      V(:,c) = accumarray( double( M.tri(:) ) , repmat( F(:,c) ,nS,1) , [nV,1] , mode );
    end
    
  end
  
  V = reshape( V , sz );
end
