function F = meshV2F( M , V , mode )
%MESHV2F  Transfer / aggregate per-VERTEX data onto the FACES (cells).
%
%   F = meshV2F( M , V )         area-weighted (default) vertex->face transfer
%   F = meshV2F( M , V , mode )  with the weighting / accumulator MODE below
%
%   For every face, F combines the values V at its vertices according to MODE.
%   V gives one value (row) per vertex and may be:
%     * a numeric array  nV x k (or nV x ...) of per-vertex values
%     * a char NAME      -> uses M.(['xyz' NAME]) if present, else M.(NAME)
%     * a function handle-> V = fun(M)  (or fun(M.xyz))
%   F has the same trailing size as V but nF rows. Mirror of meshF2V.
%
%   MODE (a string):
%     'sum'                 sum of the face's vertex values
%     'mean' / 'average'    simple average over the face's vertices
%     'area'   (default)    area-weighted (triangles, celltype 5): area(f)*sum_v
%     'narea'               normalized -> the area cancels, so it equals 'mean'
%     'length' / 'nlength'  segments   (celltype 3)
%     'volume' / 'nvolume'  tetrahedra (celltype 10)
%     'angles'              each vertex weighted by its corner angle in the face
%   The 'n...' variants divide by the total weight. Note the per-face SCALAR
%   weights ('area'/'length'/'volume') are constant across a face, so their
%   normalized forms collapse to 'mean'. Each weighting checks the celltype.
%
%   MODE may instead be a FUNCTION HANDLE, applied per face over its vertex
%   values via accumarray (e.g. @max, @min, @median, @(x)...).
%
% See also meshF2V, meshQuality, meshNormals.

  if nargin < 3, mode = 'area'; end
  
  nF = size( M.tri ,1);
  nS = size( M.tri ,2);
  nV = size( M.xyz ,1);
 
  if isa( V , 'function_handle' )
    try, V = feval( V , M ); catch
    try, V = feval( V , M.xyz ); catch
      error('invalid function to evaluate on mesh');
    end; end
  elseif ischar( V )
    try, V = M.(['xyz',V]); catch
    try, V = M.(V); catch
      error('invalid attribute name.');
    end; end
  end
  if size( V ,1) ~= nV
    error('invalid per-vertices-field');
  end

  sz = size( V );
  sz(1) = nF;
  V = V(:,:);
  
  if ischar( mode )

    N = []; W = [];
    switch lower(mode)
      case {'s','sum'}
        W = 1;
        N = 1;

      case {'m','mean','average'}
        W = 1/size( M.tri ,2);
        N = 1;
        %N = accumarray(  M.tri(:) , 1 );
        

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
    elseif issparse( W ) && size( W ,2) == nV && size( W ,1) == nF
      F = W * V;
      if isempty(N), N = sum( W ,2); end
    elseif isscalar( W )
      Tid = ( 1:nF ).';
      W = sparse( repmat( Tid ,nS,1) , double( M.tri ) , W ,nF,nV);
      F = W * double( V );
      if isempty(N), N = sum( W ,2); end
    elseif size( W ,1) == nF
      Tid = ( 1:nF ).';
      W = repmat( double(W) ,1,nS/size(W,2));
      W = sparse( repmat( double(Tid) ,nS,1) , double(M.tri) , W ,nF,nV);
      F = W * double( V );
      if isempty(N), N = sum( W ,2); end
    else
      error('incorrect W matrix');
    end
    
    %if required, normalized by N
    if 0
    elseif isscalar( N ) && N == 1
    elseif size( N ,1) == nF && size( N ,2) == 2 && all( N == 1 )
    elseif size( N ,1) == nF && size( N ,2) == 1
      F = bsxfun( @rdivide , F , N );
    elseif isscalar( N )
      F = F * (1/N);
    else
      error('incorrect normalization step');
    end
    
  elseif isa( mode , 'function_handle' )

    nC = size( V ,2);
    F = NaN( nF , nC );
    
    Fid = ( 1:nF ).';
    Fid = repmat( Fid , nS , 1 );
    
    for c = 1:nC
      F(:,c) = accumarray( Fid , V( M.tri(:) ,c) , [nF,1] , mode );
    end
    
  end
  
  F = reshape( F , sz );
end
