function [ xyz , P_id , cell_id , t , ray_id ]  = IntersectSurfaceRay( P , ray , MODE )
%INTERSECTSURFACERAY  Intersections of rays (lines) with triangulated surfaces.
%
%   [xyz,P_id,cell_id] = IntersectSurfaceRay( P )               pick under the mouse
%   [xyz,P_id,cell_id] = IntersectSurfaceRay( P , ray )         one ray (2 x 3)
%   [xyz,P_id,cell_id] = IntersectSurfaceRay( P , rays )        MANY rays
%   [ ... ]            = IntersectSurfaceRay( P , ray , MODE )
%   [xyz,P_id,cell_id,t,ray_id] = IntersectSurfaceRay( ... )
%
%   P    surface(s): a mesh struct (.vertices/.faces or .xyz/.tri), a PATCH or
%        SURFACE handle (Visible off is skipped; hgtransform parents applied),
%        or an array / cell array of them. P_id says which one was hit.
%   ray  2 x 3      [p0;p1]  the LINE through p0 -> p1, parametrized
%                            hit = p0 + t*(p1-p0), t UNBOUNDED (negatives =
%                            behind p0). Empty/omitted: the axes CurrentPoint
%                            under the mouse (interactive picking).
%        2 x 3 x N           N rays at once, one [p0;p1] page each.
%        N x 6     [p0,p1]   N rays at once, one row each.
%                            (One vectorized call is MUCH faster than N calls.)
%   MODE 'first' (default) the hit with the smallest t, per ray
%        'last'            the hit with the largest t, per ray
%        'all'             every hit, sorted by ( ray , t )
%        'any'             OCCLUSION (shadow-segment) query: SOME hit with
%                          1e-9 < t < 1-1e-5, early exit -- the cheapest way
%                          to ask "is the open segment p0->p1 blocked?"; the
%                          guard bands let p1 lie exactly ON the surface
%                          without self-reporting. One row per ray, as 'first'.
%
%   OUTPUTS ('first'/'last': one row per ray; 'all': one row per hit)
%     xyz      x 3    intersection points          (NaN row  = ray misses)
%     P_id     x 1    index into P of the hit surface  (0    = ray misses)
%     cell_id  x 1    face row within that surface     (0    = ray misses)
%     t        x 1    ray parameter of the hit         (NaN  = ray misses)
%     ray_id   x 1    which ray produced each row (needed to read 'all' with
%                     several rays)
%   With no surfaces / no hits, 'first'/'last' return the NaN/0 rows above and
%   'all' returns empties. (The old version ERRORED on a complete miss.)
%
%   The heavy lifting runs in the IntersectSurfaceRay_mx MEX (median-split BVH
%   + Moller-Trumbore, single-threaded; compile once with
%       mex IntersectSurfaceRay_mx.cpp
%   ). The MEX is a DROP-IN with this very syntax for DATA inputs -- mesh
%   structs / cell arrays of them, ray, MODE, same outputs -- and it CACHES the
%   acceleration tree internally (fingerprint of the mesh bytes), so repeated
%   queries against the same mesh only pay the traversal; only graphics handles
%   and the interactive empty-ray form need this .m file. Without the MEX a
%   pure-MATLAB fallback (the original per-ray projection algorithm) gives the
%   same results, just slower.
%
%   Run IntersectSurfaceRay with no arguments for an interactive demo.
%
%   See also silhouette, clickOnMesh, meshIsInterior.

if nargin == 0
%%
  close all;

  [x,y,z] = sphere(50); M1 = surf2patch(x,y,z,'triangle');
  M1 = patch(M1,'facecolor','r','edgecolor',[0 0 0],'FaceVertexAlphaData',ones([size(M1.faces,1),1]),'facealpha','flat');

  [x,y,z] = peaks(20);  M2 = surf2patch(x,y,z/7,'triangle');
  M2 = patch(M2,'facecolor','b','edgecolor',[0 0 0],'FaceVertexAlphaData',ones([size(M2.faces,1),1]),'facealpha','flat');

  marker = line(NaN,NaN,'marker','o','markersize',10,'markerfacecolor',[0 1 0],'linestyle','none','color',[0 1 0]);
  ray    = line(NaN,NaN,'linestyle','-','color',[1 0 0]);
  set(gcf,'WindowButtonDownFcn',@(h,e)demo([M1;M2]));

  view(3)
  axis equal

%%
  return;
end

  if nargin < 1, return; end
  if nargin < 2, ray = []; end
  if nargin < 3 || isempty( MODE ), MODE = 'first'; end

  if isempty( ray )
    ray = hittest;
    if isempty( ray ) || ~ishandle(ray) || any(strcmp(get(ray,'type'),{'panel','figure'}))
      xyz = [NaN,NaN,NaN]; P_id = 0; cell_id = 0; t = NaN; ray_id = 1;
      return;
    end
    if ~strcmp( get(ray,'type') , 'axes' )
      ray = ancestor(ray,'axes');
    end
    ray = get( ray , 'currentpoint' );
  end

  % normalize the ray forms to N x 6 rows of [p0,p1]
  if ndims( ray ) == 3 && size(ray,1) == 2 && size(ray,2) == 3
    RAYS = double( [ reshape(ray(1,:,:),3,[]) ; reshape(ray(2,:,:),3,[]) ].' );
  elseif ismatrix( ray ) && size( ray ,2) == 6
    RAYS = double( ray );
  elseif isequal( size( ray ) , [2,3] )
    RAYS = double( [ ray(1,:) , ray(2,:) ] );
  else
    error('IntersectSurfaceRay:ray','ray must be 2x3 ([p0;p1]), 2x3xN, or N x 6 ([p0,p1] per row).');
  end
  nR = size( RAYS ,1);

  switch lower( MODE )
    case 'first', mode = 0;
    case 'last' , mode = 1;
    case 'all'  , mode = 2;
    case 'any'  , mode = 3;
    otherwise,    error('unknown mode');
  end

  % ---- collect the surfaces as a CELL of plain mesh structs -------------------
  % struct inputs pass through UNTOUCHED (no copy -- MATLAB shares the arrays
  % and the MEX reads them in place); graphics handles are converted (Visible
  % off skipped, hgtransform applied). SIDX maps every cell entry back to its
  % position in P, so P_id refers to the caller's ordering even when entries
  % were skipped.
  SURF = {};  SIDX = zeros(0,1);
  for p_i = 1:numel(P)
    M = P(p_i);
    if iscell( M ), M = M{1}; end
    if isstruct(M) && ( ( isfield(M,'vertices') && isfield(M,'faces') ) || ...
                        ( isfield(M,'xyz')      && isfield(M,'tri')   ) )
      SURF{end+1,1} = M;                                       %#ok<AGROW>
      SIDX(end+1,1) = p_i;                                     %#ok<AGROW>
    elseif ishandle(M) && strcmp( get(M,'type') , 'patch' )
      if ~onoff( M , 'Visible' ), continue; end

      XYZ_ = get( M ,'vertices' );
      TRI_ = get( M ,'faces' );

      parent = get(M,'Parent');
      while ~strcmp( get(parent,'type') , 'axes' )
        if strcmp( get(parent,'type') , 'hgtransform' )
          XYZ_ = transform( XYZ_ , get(parent,'Matrix') , 'rows' );
        end
        parent = get(parent,'Parent');
      end

      SURF{end+1,1} = struct( 'vertices',double(XYZ_) , 'faces',double(TRI_) );  %#ok<AGROW>
      SIDX(end+1,1) = p_i;                                     %#ok<AGROW>

    elseif ishandle(M) && strcmp( get(M,'type') , 'surface' )
      if ~onoff( M , 'Visible' ), continue; end

      [ TRI_ , XYZ_ ] = surf2patch( M , 'triangles' );

      parent = get( M , 'Parent' );
      while ~strcmp( get(parent,'type') , 'axes' )
        if strcmp( get(parent,'type') , 'hgtransform' )
          XYZ_ = transform( XYZ_ , get(parent,'Matrix') , 'rows' );
        end
        parent = get(parent,'Parent');
      end

      SURF{end+1,1} = struct( 'vertices',double(XYZ_) , 'faces',double(TRI_) );  %#ok<AGROW>
      SIDX(end+1,1) = p_i;                                     %#ok<AGROW>
    end
  end

  % ---- intersect ---------------------------------------------------------------
  if exist( 'IntersectSurfaceRay_mx' , 'file' ) == 3
    % the MEX takes the CELL directly (it concatenates the surfaces itself --
    % or not at all on an acceleration-cache hit) and already returns the final
    % outputs; only the surface ids need remapping to the caller's P positions.
    [ xyz , pid , cell_id , t , ray_id ] = IntersectSurfaceRay_mx( SURF , RAYS , mode );
    P_id = zeros( size(pid) );
    w = pid > 0;
    P_id(w) = SIDX( pid(w) );
    return;
  end

  % ---- pure-MATLAB fallback: flatten into one soup and use the original core --
  XYZ = zeros(0,3);
  TRI = zeros(0,5);
  for s = 1:numel(SURF)
    S = SURF{s};
    if isfield( S ,'vertices'), XYZ_ = S.vertices;  TRI_ = S.faces;
    else,                       XYZ_ = S.xyz;       TRI_ = S.tri;
    end
    if size( TRI_ , 2 ) ~= 3        , error('non triangular mesh, size(TRI,2) ~= 3'); end
    if size( XYZ_ , 2 ) ~= 3        , error('non triangular mesh, size(XYZ,2) ~= 3'); end
    if any( TRI_(:) < 1 )           , error('non triangular mesh, TRI < 1');           end
    if any( TRI_(:) > size(XYZ_,1) ), error('non triangular mesh, TRI > size(XYZ,1)'); end

    TRI_ = double( TRI_ ) + size( XYZ , 1 );
    TRI_(:,4) = SIDX(s);
    TRI_(:,5) = ( 1:size(TRI_,1) ).';

    XYZ = [ XYZ ; double(XYZ_) ];                              %#ok<AGROW>
    TRI = [ TRI ; TRI_ ];                                      %#ok<AGROW>
  end

  [ th , hh , rh ] = localIntersect( XYZ , TRI(:,1:3) , RAYS , mode );

  if mode == 2                                   % 'all': one row per hit
    w       = hh > 0;
    th      = th(w); hh = hh(w); rh = rh(w);
    if isempty( hh )                             % no hit at all -> empty outputs
      xyz = zeros(0,3); P_id = zeros(0,1); cell_id = zeros(0,1);
      t   = zeros(0,1); ray_id = zeros(0,1);
      return;
    end
    xyz     = RAYS(rh,1:3) + th .* ( RAYS(rh,4:6) - RAYS(rh,1:3) );
    P_id    = TRI( hh , 4 );
    cell_id = TRI( hh , 5 );
    t       = th;
    ray_id  = rh;
  else                                           % 'first'/'last': one row per ray
    xyz     = NaN( nR , 3 );
    P_id    = zeros( nR , 1 );
    cell_id = zeros( nR , 1 );
    t       = th;
    ray_id  = ( 1:nR ).';
    w = hh > 0;
    if any( w )
      xyz(w,:)   = RAYS(w,1:3) + th(w) .* ( RAYS(w,4:6) - RAYS(w,1:3) );
      P_id(w)    = TRI( hh(w) , 4 );
      cell_id(w) = TRI( hh(w) , 5 );
    end
  end

  function demo(hs)
    try,   ry = get(gca,'CurrentPoint');
    catch, ry = NaN(2,3);
    end
    set( ray , 'XData',ry(:,1), 'YData',ry(:,2), 'ZData',ry(:,3) );
    [ xyz , P_id , cell_id ]  = IntersectSurfaceRay( hs , [] , 'all' );
    set(marker,'XData',xyz(:,1),'YData',xyz(:,2),'ZData',xyz(:,3));
    arrayfun(@(h)set(h,'EdgeColor',[1 1 1]*0.6),hs);
    arrayfun(@(h)set(h,'FaceVertexAlphaData',zeros(size(get(h,'FaceVertexAlphaData')))),hs);
    if isempty( P_id ) || ~P_id(1), return; end
    set( hs(P_id(1)) , 'EdgeColor', [1 1 1]*0 );
    for i = 1:numel(P_id)
      fc = get( hs(P_id(i)) , 'FaceVertexAlphaData' );
      fc( cell_id(i) ) = 0.5;
      set( hs(P_id(i)) , 'FaceVertexAlphaData' , fc );
    end
  end

end


function [ tq , hq , rq ] = localIntersect( XYZ , TRI , RAYS , mode )
%pure-MATLAB fallback with the SAME contract as the MEX: per ray, the original
%projection algorithm (project the vertices onto the plane transverse to the
%ray, sign-prune the triangles straddling the two transverse planes, keep the
%ones whose 2D projection contains the origin, and intersect their planes).
  nR = size( RAYS ,1);
  tq = NaN( nR ,1); hq = zeros( nR ,1); rq = ( 1:nR ).';
  if mode == 2, tq = []; hq = []; rq = []; end

  for r = 1:nR
    p0 = RAYS(r,1:3);  p1 = RAYS(r,4:6);
    R  = null( p0 - p1 );
    XY = XYZ * R - p0 * R;

    F = ( 1:size(TRI,1) ).';
    F = F( mod( sum( sign( [ XY(TRI(F,1),1) , XY(TRI(F,2),1) , XY(TRI(F,3),1) ] ) , 2 ) , 3 ) > 0 );
    F = F( mod( sum( sign( [ XY(TRI(F,1),2) , XY(TRI(F,2),2) , XY(TRI(F,3),2) ] ) , 2 ) , 3 ) > 0 );
    F = F( mod( ...
           sign( ( XY(TRI(F,2),1)-XY(TRI(F,1),1) ) .* XY(TRI(F,1),2) - ( XY(TRI(F,2),2)-XY(TRI(F,1),2) ) .* XY(TRI(F,1),1) ) + ...
           sign( ( XY(TRI(F,3),1)-XY(TRI(F,2),1) ) .* XY(TRI(F,2),2) - ( XY(TRI(F,3),2)-XY(TRI(F,2),2) ) .* XY(TRI(F,2),1) ) + ...
           sign( ( XY(TRI(F,1),1)-XY(TRI(F,3),1) ) .* XY(TRI(F,3),2) - ( XY(TRI(F,1),2)-XY(TRI(F,3),2) ) .* XY(TRI(F,3),1) ) , ...
         3 ) == 0 );
    if isempty( F ), continue; end

    N  = cross( XYZ(TRI(F,2),:) - XYZ(TRI(F,1),:) , XYZ(TRI(F,3),:) - XYZ(TRI(F,1),:) , 2 );
    D  = sum( XYZ(TRI(F,1),:).*N , 2 );
    tt = ( D - N * p0.' )./( N * (p1-p0).' );

    w  = isfinite( tt );
    F  = F(w); tt = tt(w);
    if isempty( F ), continue; end

    switch mode
      case 0, [ tq(r) , i ] = min( tt );  hq(r) = F(i);
      case 1, [ tq(r) , i ] = max( tt );  hq(r) = F(i);
      case 2
        [ tt , o ] = sort( tt );
        tq = [ tq ; tt ]; hq = [ hq ; F(o) ]; rq = [ rq ; r*ones(numel(o),1) ]; %#ok
      case 3                                   % 'any': some hit inside the segment
        i = find( tt > 1e-9 & tt < 1-1e-5 , 1 );
        if ~isempty( i ), tq(r) = tt(i);  hq(r) = F(i); end
    end
  end
end
