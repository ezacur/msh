function [xyz, cell_id, t, ray_id] = bvhIntersectRay( M , ray , B , MODE )
%bvhIntersectRay  Ray (line) vs mesh-triangle intersections on the BVH blob.
%
%   [xyz,cell,t,rid] = bvhIntersectRay( M , ray )            builds/uses a BVH
%   [xyz,cell,t,rid] = bvhIntersectRay( M , ray , B )        reuses a blob
%   [xyz,cell,t,rid] = bvhIntersectRay( {M,B} , ray )        mesh+blob bundled
%   [ ... ]          = bvhIntersectRay( ...   , MODE )
%
%   Same ray conventions as tools/IntersectSurfaceRay: ray is 2x3 [p0;p1],
%   2x3xN pages, or N x 6 rows [p0,p1]; the LINE hit = p0 + t*(p1-p0) with t
%   UNBOUNDED (negatives = behind p0). MODE: 'first' (default, smallest t) |
%   'last' | 'all' (every hit, sorted by ray,t) | 'any' (occlusion query,
%   1e-9 < t < 1-1e-5 with early exit).
%
%   Uses the SAME cached/transformable acceleration blob as bvhClosestElement
%   (one structure per mesh serving both query types): similarity transforms
%   folded in the frame cost O(1), deformations refit, and the blob is a value
%   (msh-cacheable, save/load). Only TRIANGLE cells are tested; other celltypes
%   in mixed meshes are skipped. Requires the compiled bvhIntersectRay_mx and
%   a v2 packed blob (BVH default when BVH_mx is built); for machines
%   without the MEXes use tools/IntersectSurfaceRay (pure-MATLAB fallback).
%
% See also BVH, bvhClosestElement, IntersectSurfaceRay.

  %bundled form: bvhIntersectRay( {M,B} , ray , MODE )
  if iscell( M )
    if nargin > 2, MODE = B; end          %shift: 3rd arg is MODE in this form
    B = M{2};  M = M{1};
  elseif nargin < 3
    B = [];
  end
  if ~exist('MODE','var') || isempty( MODE ), MODE = 'first'; end
  switch lower( MODE )
    case 'first', modeI = 1;
    case 'last',  modeI = 2;
    case 'all',   modeI = 3;
    case 'any',   modeI = 4;
    otherwise, error('bvhIntersectRay:mode','MODE must be first|last|all|any.');
  end
  if exist( 'bvhIntersectRay_mx' ,'file') ~= 3
    error('bvhIntersectRay:mex','bvhIntersectRay_mx is not compiled (mex COMPFLAGS="$COMPFLAGS /openmp" -lut bvhIntersectRay_mx.cpp).');
  end

  %normalize the ray shapes to N x 6
  if size( ray ,1) == 2 && size( ray ,2) == 3
    ray = reshape( permute( ray , [2 1 3] ) , 6 , [] ).';
  elseif size( ray ,2) ~= 6
    error('bvhIntersectRay:ray','ray must be 2x3, 2x3xN or Nx6.');
  end
  ray = double( ray );

  if isempty( B ), B = BVH( M ); end
  if ~isfield( B , 'pkV' ), B = BVH( M ); end        %needs the v2 packed blob

  %staleness spot-check (same pattern as bvhClosestElement)
  Xw = double( M.xyz ); Xw(:,end+1:3) = 0;
  ok = size( Xw ,1) == size( B.X ,1) && isequal( size( double(M.tri) ) , size( B.Tri ) );
  if ok && size( Xw ,1) > 0
    ii = unique( round( linspace( 1 , size(Xw,1) , 4 ) ) );
    Yw = B.X(ii,:) * B.frame(1:3,1:3).' + B.frame(1:3,4).';
    ok = max(max( abs( Yw - Xw(ii,:) ) )) <= 1e-6 * max( 1 , max(max( abs( Xw(ii,:) ) )) );
  end
  if ~ok
    %ERROR on purpose (see bvhClosestElement: an in-call rebuild would be
    %discarded and silently repeated on every later call).
    error('bvhIntersectRay:staleBVH', ...
          'B does not match M (stale or foreign blob). Refit it (B = BVH(B,M)) or rebuild it (B = BVH(M)).');
  end

  %global frame: rays into build space (t is INVARIANT under the similarity)
  F    = B.frame;
  hasF = ~isequal( F , eye(4) );
  if hasF
    Af = F(1:3,1:3);  tf = F(1:3,4).';
    s2 = trace( Af.'*Af )/3;
    Ai = Af / s2;                          %(x-t)*inv(A).' = (x-t)*(A/s2)
    ray = [ ( ray(:,1:3) - tf ) * Ai , ( ray(:,4:6) - tf ) * Ai ];
  end

  [ xyz , cell_id , t , ray_id ] = bvhIntersectRay_mx( ray , B , modeI , maxNumCompThreads );

  if hasF                                  %hit points back to world space
    xyz = xyz * Af.' + tf;
  end

end
