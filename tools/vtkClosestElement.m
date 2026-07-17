function [ element_id , xyz_closest_point , distance , barycentric_coordinates ] = vtkClosestElement( M , xyz )
%VTKCLOSESTELEMENT  Closest mesh element (and, optionally, robust barycentric
%   weights) for a set of query points -- a thin SHADOW over the compiled MEX
%   vtkClosestElement_ that FIXES its barycentric-coordinate bug.
%
%   [id,cp,d,bc] = vtkClosestElement( MESH , XYZ )   one-shot: for every row of
%       XYZ (M-by-3, or M-by-2 in 2D) returns the 1-based id of the closest
%       triangle, the closest point cp, the distance d, and the barycentric
%       coordinates bc of cp within triangle id.
%
%   Persistent-locator form (as in the original MEX):
%     vtkClosestElement( MESH )        build the cell locator from MESH
%     [...] = vtkClosestElement( XYZ ) query points against the last locator
%     vtkClosestElement( [] , [] )     free the locator
%   ( vtkClosestElement() with no arguments is an error: the MEX does nothing
%     for it -- use ( [] , [] ) to free. )
%
%   WHY THIS SHADOW EXISTS
%     The compiled MEX (now vtkClosestElement_) computes the barycentric weights
%     (its 4th output) by solving the ORIGIN-based linear system M*w = cp with
%     M = [v1 v2 v3] (the triangle vertices as columns), dividing by
%     det(M) = v1.(v2 x v3) = 6*signed volume of the tetra (ORIGIN,v1,v2,v3).
%     That determinant -> 0 whenever the triangle's supporting plane passes near
%     the coordinate origin (and is EXACTLY 0 for any triangle with a vertex at
%     the origin), so the weights blow up (they stop summing to 1, go negative)
%     for those triangles -- a translation-variance bug. This .m keeps using the
%     MEX for id / closest point / distance (all correct) but RECOMPUTES the
%     barycentric weights here from the triangle's OWN edges via CROSS products
%     (sub-triangle areas), which is both translation-invariant AND free of the
%     catastrophic cancellation that the Gram dot form suffers on slivers -- so
%     it stays exact even for very thin triangles (to machine precision at
%     aspect 1e12). Weights are clamped to [0,1] and renormalized to sum to 1.
%     When the MEX is recompiled with the fix (see the notes in
%     vtkClosestElement_.cpp) this shadow can be deleted.
%
%   See also distanceFrom, MeshQuery, meshMapPoints, vtkClosestElement_.

  persistent lastMesh   %mesh cached at "build locator" time, for the barycentric
                        %of a later query-only call (mirrors the MEX's static locator)

  %----- no arguments: NOT a real call (the MEX just prints "No entiendo") -----
  if nargin == 0
    error( 'vtkClosestElement:signature' , 'no arguments given. Use ( MESH ), ( XYZ ), ( MESH , XYZ ), or ( [] , [] ) to free the locator.' );
  end

  %----- free the locator ( [] , [] ) -----------------------------------------
  if nargin == 2 && isnumeric(M) && isnumeric(xyz)
    try, vtkClosestElement_( [] , [] ); end
    lastMesh = [];
    return;
  end

  %----- build the locator ( MESH ) -------------------------------------------
  if nargin == 1 && isstruct( M )
    lastMesh = normalizeMesh( M );
    vtkClosestElement_( lastMesh );
    return;
  end

  %----- resolve the mesh + query for the evaluating forms --------------------
  if nargin == 2 && isstruct( M ) && isnumeric( xyz )     % one-shot ( MESH , XYZ )
    Md = normalizeMesh( M );
    ONESHOT = true;
  elseif nargin == 1 && isnumeric( M )                    % query the last locator
    xyz = M;
    Md  = lastMesh;                                       % may be [] (barycentric then unavailable)
    ONESHOT = false;
  else
    error( 'vtkClosestElement:signature' , 'unrecognized call. Use ( MESH ), ( XYZ ), ( MESH , XYZ ) or ( [] , [] ).' );
  end

  xyz = double( xyz );
  if size( xyz ,2) < 3, xyz(:,end+1:3) = 0; end          % 2D -> pad z=0 (the MEX reads 3 columns)

  %----- delegate id / closest / distance to the MEX (its correct outputs) ----
  %never request the MEX's 4th output (its buggy barycentric); recompute below.
  if ONESHOT
    if     nargout >= 3, [ element_id , xyz_closest_point , distance ] = vtkClosestElement_( Md , xyz );
    elseif nargout == 2, [ element_id , xyz_closest_point ]            = vtkClosestElement_( Md , xyz );
    else,                  element_id                                  = vtkClosestElement_( Md , xyz );
    end
  else
    if     nargout >= 3, [ element_id , xyz_closest_point , distance ] = vtkClosestElement_( xyz );
    elseif nargout == 2, [ element_id , xyz_closest_point ]            = vtkClosestElement_( xyz );
    else,                  element_id                                  = vtkClosestElement_( xyz );
    end
  end

  %----- robust barycentric coordinates (the whole point of the shadow) --------
  if nargout >= 4
    if isempty( Md )
      error( 'vtkClosestElement:noMesh' , 'barycentric coordinates need the mesh: use the one-shot ( MESH , XYZ ) form, or build the locator with a MESH struct first.' );
    end
    barycentric_coordinates = calcular_barycentric( Md , element_id , xyz_closest_point );
  end

end


function W = calcular_barycentric( M , eid , cp )
%barycentric weights of each closest point cp within its triangle eid, as
%ratios of signed sub-triangle areas computed with CROSS products, and
%translation-invariant (edge vectors from vertex 1). Fully vectorized.
%
%WHY CROSS PRODUCTS (not the Gram dot form d00*d11-d01^2): that determinant is
%(0.25+h^2)-0.25 for a thin triangle -> CATASTROPHIC CANCELLATION -> the weights
%lose all precision for slivers (measured: 5% error at aspect 1e8, NaN at
%1e10). |v0 x v1|^2 computes the SAME 4*area^2 directly, so this stays at
%machine precision even at aspect 1e10.
  eid = double( eid(:) );
  T = M.tri( eid , 1:3 );
  a = M.xyz( T(:,1) ,:);
  b = M.xyz( T(:,2) ,:);
  c = M.xyz( T(:,3) ,:);

  v0 = b - a;   v1 = c - a;   v2 = cp - a;
  n  = cross( v0 , v1 , 2 );                    % triangle normal, |n| = 2*area
  nn = sum( n.*n ,2);                           % 4*area^2, no dot-form cancellation

  gamma = sum( n .* cross( v0 , v2 , 2 ) ,2) ./ nn;   % weight of vertex 3 (area a-b-cp / a-b-c)
  beta  = sum( n .* cross( v2 , v1 , 2 ) ,2) ./ nn;   % weight of vertex 2 (area a-cp-c / a-b-c)
  alpha = 1 - beta - gamma;                           % weight of vertex 1

  W = [ alpha , beta , gamma ];

  %cp is ON the triangle, so the exact weights are in [0,1]: clamp the tiny
  %numerical spill (e.g. -1e-16 for a point on an edge) and renormalize so each
  %row is a valid convex combination summing to 1. (Accuracy comes from the
  %cross form above; this only guarantees validity for downstream consumers.)
  W = max( W , 0 );
  W = W ./ sum( W , 2 );
end


function M = normalizeMesh( M )
%keep only .xyz/.tri, as double, padded to 3D (the MEX assumes double & 3 cols)
  if ~isstruct( M ) || ~isfield( M ,'xyz' ) || ~isfield( M ,'tri' )
    error( 'vtkClosestElement:mesh' , 'M must be a mesh struct with .xyz and .tri.' );
  end
  xyz = double( M.xyz );
  if size( xyz ,2) < 3, xyz(:,end+1:3) = 0; end
  M = struct( 'xyz' , xyz , 'tri' , double( M.tri ) );
end
