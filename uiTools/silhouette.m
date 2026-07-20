function hS_ = silhouette( hP , varargin )
%SILHOUETTE  View-dependent outline of a triangular patch, auto-updated with the
%            camera.
%
%   hS = silhouette( hP )            attach a silhouette to the patch hP
%   hS = silhouette( M )             mesh struct (.xyz,.tri): an invisible
%                                    carrier patch is created for it
%   hS = silhouette( [h1,h2,...] )   one silhouette per patch (vector output)
%   hS = silhouette( hP , 'Prop',v , ... )   extra properties passed to the
%                                    silhouette patch (EdgeColor, LineWidth,
%                                    LineStyle, ...). Defaults: black, 1.2pt.
%   silhouette( [] , hS )            re-attach an EXISTING silhouette (used by
%                                    the CreateFcn on file load).
%
%   The silhouette is drawn as a line patch (2-vertex faces) holding, for the
%   CURRENT camera, the edges that separate camera-facing from back-facing
%   triangles, plus the open boundary of the mesh (always shown). Listeners on
%   the axes camera (CameraPosition/Target/ViewAngle, View, Projection) keep it
%   in sync; both orthographic and perspective projections are supported, and
%   hgtransform parents are accounted for (normals are re-rotated).
%
%   OCCLUSION ALPHAS -- set( hS , 'UserData' , [aHIDDEN aVISIBLE] ) draws the
%   silhouette vertices that are OCCLUDED by the mesh itself with alpha aHIDDEN
%   and the directly visible ones with aVISIBLE (EdgeAlpha 'interp'). A scalar
%   means [a 1]; empty (default) = uniform EdgeAlpha, no occlusion test.
%   The occlusion runs LIVE at every camera event: one segment query per
%   silhouette vertex against a per-mesh BVH blob (bvhIntersectRay_mx, engine
%   in BVH\, sub-ms at any mesh size; independent of the msh class and of any
%   Mesh*/mesh* function). Without the BVH engine it uses the legacy
%   IntersectSurfaceRay_mx (internally cached tree), and without any MEX it
%   falls back to a transverse-grid test (orthographic) / per-vertex rays
%   (perspective), also live. An ADAPTIVE SAFETY THROTTLE caps recomputes at
%   max( 33 ms , 3x the measured pass cost ) during event bursts (isolated
%   events are always synchronous; the duty cycle stays ~1/3 at any mesh size)
%   and a trailing pass lands the exact state ~50 ms after the stream pauses.
%
%   LIFETIME -- deleting the patch deletes its silhouette (and only then).
%   Visible follows the patch: hiding the patch hides the silhouette, and
%   re-showing it restores the silhouette to the USER's last explicit state --
%   a silhouette you turned off with set(hS,'Visible','off') STAYS off when the
%   patch is shown again (turn it back on explicitly). Editing the patch
%   Vertices / Faces re-caches and redraws. A
%   savefig/openfig round-trip re-attaches automatically (CreateFcn).
%   COPYOBJ does NOT fire CreateFcn nor copy appdata: a copied silhouette is a
%   dead patch -- re-attach with silhouette( copiedPatch ).
%
%   HOW IT WORKS / PERFORMANCE
%     updateMesh (per mesh edit)   caches face normals/centers and the STATIC
%       edge structure: interior edges + their 2 incident faces + boundary.
%     updateCamera (per camera event)  first a guard on the EFFECTIVE inputs:
%       in orthographic everything depends only on the view DIRECTION
%       (quantized 1e-12), so zoom / dolly / roll recompute NOTHING; in
%       perspective only on CameraPosition (zoom and CameraTarget are free).
%       Then faces are classified by sign(N.viewdir) and the silhouette is the
%       cached interior edges whose two faces fall on OPPOSITE sides -- pure
%       logical indexing (~6 ms @ 350k faces), with the Faces re-upload
%       skipped when the edge set did not change. Non-manifold meshes (edges
%       with >2 faces) use an exact per-event fallback. The occlusion test is
%       one viewer->vertex segment query per silhouette vertex against a BVH
%       blob built once per mesh in updateMesh (sub-ms per event even at 80k+
%       faces); without the MEXes, a transverse 2D grid (orthographic) or
%       per-vertex rays (perspective) do the same, slower. No forced drawnow inside the
%       camera callbacks: MATLAB repaints naturally between events.
%
%   NOTE: this shadows the Statistics Toolbox SILHOUETTE (cluster plot); the
%   CreateFcn guards against resolving to that one on load.
%
%   See also headlight, patchBorder, plotMESH, camlight.

  if numel( hP ) > 1
    try,    hS = gobjects( size( hP ) );
    catch,  hS = zeros( size( hP ) );
    end
    for i = 1:numel( hS )
      hS(i) = silhouette( hP(i) , varargin{:} );
    end
    return;
  end

  if isstruct( hP ) && isfield( hP ,'tri') && isfield( hP ,'xyz')
    hP = patch( 'Vertices' , hP.xyz , 'Faces' , hP.tri , 'FaceColor','none','EdgeColor','none' , 'HandleVisibility' , 'off' );
  end


  if isempty( hP ) && numel( varargin ) == 1 && ishandle( varargin{1} )
    hS = varargin{1};
  else

    if size( get( hP , 'Faces' ) ,2) ~= 3
      error('Only triangular meshes are allowed.');
    end

    hS = patch( 'Parent'            , get(hP,'Parent') , ...
                'Vertices'          , [] , ...
                'Faces'             , [] , ...
                'EdgeColor'         , [0,0,0] , ...
                'EdgeLighting'      , 'none' ,...
                'LineWidth'         , 1.2 , ...
                'Visible'           , get( hP , 'Visible') , ...
                'Clipping'          , get( hP , 'Clipping') , ...
                'HandleVisibility'  , 'off' , ...
                'Hittest'           , 'off' , ...
                'PickableParts'     , 'none' , ...
                'XLimInclude'       , get( hP , 'XLimInclude' ) ,...
                'YLimInclude'       , get( hP , 'YLimInclude' ) ,...
                'ZLimInclude'       , get( hP , 'ZLimInclude' ) ,...
                'Tag'               , [ 'silhouette_' , get( hP ,'Tag' ) ] ,...
                varargin{:} );
    try, set( hS ,'LineJoin','round'); end
    try, set( hS ,'LineSmoothing','on'); end
    try
      DisplayName_hP = get( hP , 'DisplayName' );
      DisplayName_hS = get( hS , 'DisplayName' );
      if isempty( DisplayName_hS ) && ~isempty( DisplayName_hP )
        set( hS , 'DisplayName' , [ 'silhouette_' , DisplayName_hP ] );
      end
    end
  
    setappdata( hS , 'silhouette_parentPatch' , hP );
    set( hS , 'CreateFcn' , @(h,~)( any( ~isempty(which('silhouette'))&&isempty(regexp(fileparts(which('silhouette.m')),'stats$','once'))&&double(silhouette( [] , h )) ) || any( set(h,'Visible','off') ) ) );
  end

  hAx = ancestor( hS , 'axes' ); hFig = ancestor( hAx ,'figure');
  hP  = getappdata( hS , 'silhouette_parentPatch' );

  if isnumeric( hAx ), hAx = handle( hAx ); end
  if isnumeric( hP  ), hP  = handle( hP  ); end
  if isnumeric( hS  ), hS  = handle( hS  ); end
  
  set( hP , 'DeleteFcn' , @(~,~)delete( hS ) );
  set( hS , 'DeleteFcn' , @(~,~)cleanupSilhouette( hP , hS ) );

  % SAFETY THROTTLE timer: isolated camera events update synchronously (real
  % time), but a stream firing faster than ~30 Hz skips the in-between
  % recomputes; this single-shot timer then runs one TRAILING pass ~50 ms after
  % the last skipped event, so the resting state is always exact. It lives in
  % the silhouette's appdata and cleanupSilhouette (the DeleteFcn above)
  % stops/deletes it -- no timer outlives its silhouette.
  tOld = getappdata( hS , 'silhouette_throttle' );
  if ~isempty( tOld ), try, stop( tOld ); delete( tOld ); end; end
  tOld = getappdata( hS , 'silhouette_timer' );          %legacy (pre-throttle) name,
  if ~isempty( tOld ), try, stop( tOld ); delete( tOld ); end; end   %e.g. old .fig files
  tD = timer( 'ExecutionMode','singleShot' , 'StartDelay',0.05 , ...
              'Name','silhouette_throttle' , 'ObjectVisibility','off' , ...
              'TimerFcn', @(~,~)trailingUpdate( hS , hAx ) );
  setappdata( hS , 'silhouette_throttle' , tD );
  setappdata( hS , 'silhouette_passCost' , 0 );   %measured cost of the last pass
                                                  %(feeds the ADAPTIVE throttle)

  % the USER's desired visibility of the silhouette. The parent coupling means:
  % patch off -> silhouette off (forced); patch back on -> silhouette returns to
  % THIS state, so a silhouette the user explicitly hid STAYS hidden. It is
  % updated by the hS Visible listener only when the change did not come from
  % the parent sync (see syncVisibleFromParent / visibleChanged).
  if isempty( getappdata( hS , 'silhouette_userVisible' ) )
    setappdata( hS , 'silhouette_userVisible' , get( hS , 'Visible' ) );
  end

  try, updateMesh( hS ); end
  setappdata( hS , 'silhouette_listeners' , {   newListener( hAx , 'CameraPosition'       , @(~,~)updateCamera( hS , hAx ) ) ;...
                                                newListener( hAx , 'CameraTarget'         , @(~,~)updateCamera( hS , hAx ) ) ;...
                                                newListener( hAx , 'CameraViewAngle'      , @(~,~)updateCamera( hS , hAx ) ) ;...
                                                newListener( hAx , 'Projection'           , @(~,~)updateCamera( hS , hAx ) ) ;...
                                                newListener( hAx , 'View'                 , @(~,~)updateCamera( hS , hAx ) ) ;...
                                                newListener( hS  , 'Visible'              , @(~,~)visibleChanged( hS , hAx ) ) ;...
                                                newListener( hS  , 'UserData'             , @(~,~)updateCamera( hS , hAx ) ) ;...
                                                newListener( hP  , 'Visible'              , @(~,~)syncVisibleFromParent( hP , hS ) ) ;...
                                                newListener( hP  , 'Vertices'             , @(~,~)updateMesh( hS ) ) ;...
                                                newListener( hP  , 'Faces'                , @(~,~)updateMesh( hS ) ) ;...
                                            } );
%                                                 newListener( hP  , 'FaceVertexAlphaData'  , @(~,~)set( hS ,'FaceVertexAlphaData' , get( hP , 'FaceVertexAlphaData' ) ) ) ;...
  try, updateCamera( hS , hAx ); end
  
  if nargout > 0, hS_ = hS; end
end

function updateMesh( hS )
%re-cache everything that depends ONLY on the mesh (not on the camera). Fired by
%the Vertices / Faces listeners of the parent patch, and once at attach time.
%Stores in the silhouette's appdata: the faces, the (unnormalized) face normals
%(front/back classification), the face centers (perspective classification), and
%the STATIC edge structure used by the per-event fast path. Finishes by clearing
%the dedup state and forcing one updateCamera pass.
  hP = getappdata( hS , 'silhouette_parentPatch' );
  syncVisibleFromParent( hP , hS );

  F = get( hP , 'Faces' );
  if size( F ,2) ~= 3
    warning( 'Only triangular meshes are allowed for showing silhouettes. Turning silhouette off.');
    set( hS ,'Visible', 'off' );
  end
%   if ~strcmp( get( get( hP , 'Parent' ) ,'Type' ) ,'axes' )
%     warning( 'Silhouettes can be wrong for hgtransform children.');
%   end
  
  V = get( hP , 'Vertices' );
  set( hS , 'Vertices' , V );
  
  setappdata( hS , 'silhouette_faces'       , F );
  setappdata( hS , 'silhouette_faceCenters' , meshFacesCenter( V , F ) );
  setappdata( hS , 'silhouette_normals'     , meshNormals( V , F ) );
  setappdata( hS , 'silhouette_lastBlocker' , zeros( size(V,1) ,1) );   %shadow cache (occlusion)

  % OCCLUSION ACCELERATOR: explicit BVH blob for this mesh (BVH engine in
  % BVH\ -- independent of the msh class and of any Mesh*/mesh* function).
  % Built HERE (attach/edit time) so the first occlusion camera event doesn't
  % hitch, and stored in the silhouette's appdata: it dies with the silhouette
  % and ANY number of silhouettes coexist (no shared cache to thrash --
  % IntersectSurfaceRay_mx's internal LRU has 4 slots: 5+ occlusion
  % silhouettes there mean a full tree REBUILD per camera event). 'noframe'
  % keeps the blob in WORLD coordinates, so updateCamera feeds the rays to
  % bvhIntersectRay_mx DIRECTLY (no wrapper, no frame folding).
  Bb = [];
  if exist( 'bvhIntersectRay_mx' ,'file') == 3 && exist( 'BVH_mx' ,'file') == 3 ...
                                               && exist( 'BVH' ,'file') == 2
    try, Bb = BVH( struct( 'xyz',double(V) , 'tri',double(F) ) , [] , 'noframe' ); end
  end
  setappdata( hS , 'silhouette_bvhBlob' , Bb );

  % FALLBACK (no BVH engine on the path): WARM the internal tree cache of
  % IntersectSurfaceRay_mx for this mesh -- the mex builds (and caches, keyed
  % by the mesh bytes) its BVH on the first heavy call; forcing that build
  % HERE moves the one-time cost to attach/edit time. The dummy rays start far
  % outside any bbox: each costs O(1) (root slab fails immediately). Meshes
  % under 1024 faces never get a tree (mex policy: always brute) -> no warm.
  if isempty( Bb ) && exist( 'IntersectSurfaceRay_mx' , 'file' ) == 3 && size( F ,1) >= 1024
    kw = ceil( 8.1e6 / max( size(F,1) , 1 ) ) + 1;      %enough rays to cross the
    try                                                 %build threshold
      IntersectSurfaceRay_mx( struct( 'vertices',double(V) , 'faces',double(F) ) , ...
                              repmat( [1e300 1e300 1e300 1e300 1e300 1e299] , kw , 1 ) , 'first' );
    end
  end

  % STATIC edge structure (built once per mesh, reused at every camera event):
  % unique undirected edges Eu, their (up to) two incident faces EF, and the open
  % boundary (1-face edges). With it the per-event silhouette is a pure logical
  % test front(EF1) ~= front(EF2) -- no unique/sort per event. Non-manifold
  % meshes (an edge with >2 faces) fall back to the exact per-event path.
  F3 = double( F );
  [ Eu , ~ , ec ] = unique( sort( [ F3(:,[1,2]) ; F3(:,[2,3]) ; F3(:,[3,1]) ] ,2) , 'rows' );
  cnt = accumarray( ec , 1 );
  fid = repmat( ( 1:size(F3,1) ).' , 3 , 1 );
  [ ecs , o ] = sort( ec );  fids = fid( o );
  isf = [ true ; diff( ecs ) ~= 0 ];                  %first occurrence of each edge
  iss = [ false ; isf(1:end-1) ] & ~isf;              %second occurrence
  EF  = zeros( size( Eu ,1) , 2 );
  EF( ecs(isf) , 1 ) = fids( isf );
  EF( ecs(iss) , 2 ) = fids( iss );
  wint = EF(:,2) > 0;                                 %interior (2-face) edges,
  setappdata( hS , 'silhouette_edgesInt'     , Eu( wint ,:) );   %pre-extracted so the
  setappdata( hS , 'silhouette_edgeFacesInt' , EF( wint ,:) );   %per-event path skips
  setappdata( hS , 'silhouette_manifold'  , all( cnt <= 2 ) );   %the interior masking
  setappdata( hS , 'silhouette_boundary'  , Eu( cnt == 1 ,:) );   % == meshBoundary( F )

  setappdata( hS , 'silhouette_lastState'   , [] );   %mesh changed: force the update
  setappdata( hS , 'silhouette_lastFaces'   , [] );   %...and never skip the re-upload
  updateCamera( hS , ancestor( hS ,'axes') );
end
function updateCamera( hS , hAx )
%recompute the silhouette for the CURRENT camera. Fired by the camera listeners
%(CameraPosition/Target/ViewAngle, View, Projection) and by the silhouette's own
%Visible / UserData listeners. The pipeline is:
%   1. dedup guard        skip if nothing relevant changed since the last pass
%   2. hgtransform walk   rotate the cached normals (and, in perspective, move
%                         the face centers) by the accumulated parent transforms
%   3. classification     angles = N . viewdir  (orthographic: one direction;
%                         perspective: per-face center-to-camera direction)
%   4. edge set           FAST path: cached interior edges whose two incident
%                         faces fall on opposite sides of the front/back
%                         partition, plus the fixed open boundary.
%                         FALLBACK (non-manifold): boundary of the front subset.
%   5. occlusion (opt.)   per-silhouette-vertex hidden test -> interp alphas,
%                         computed LIVE at every event: one segment query per
%                         vertex against the per-mesh BVH blob (sub-ms; legacy
%                         IntersectSurfaceRay_mx without the BVH engine; the
%                         no-MEX fallbacks -- transverse 2D grid in ortho,
%                         per-vertex rays in perspective -- are also live).
  if ~ishandle( hS  ), return; end
  if ~ishandle( hAx ), return; end
  if strcmp( get( hS ,'Visible') ,'off' ), return; end

  % dedup guard on the EFFECTIVE inputs (same idea as silhouetteByShading): in
  % ORTHOGRAPHIC both the edge set and the occlusion depend only on the view
  % DIRECTION (normalized, quantized at 1e-12 to absorb float noise), so zoom,
  % dolly and roll recompute NOTHING and a pan collapses to one throttled pass
  % (its two property events go through an intermediate direction); in
  % PERSPECTIVE they depend only on CameraPosition (zoom and CameraTarget
  % changes are free). UserData (the occlusion alphas) stays in the state so
  % toggling it always passes. updateMesh clears the state to force a pass.
  % (The very first guard here compared the figure's UpdateToken, which skipped
  % exactly the FIRST camera event after each render -> stale after a
  % programmatic view(); never compare render tokens.)
  if strcmp( get( hAx ,'Projection') ,'orthographic' )
    dirQ = get( hAx ,'CameraTarget') - get( hAx ,'CameraPosition');
    dirQ = dirQ ./ sqrt( sum( dirQ.^2 ) );
    dirQ = round( dirQ * 1e12 ) * 1e-12;
    state = { 'o' , dirQ , get( hS ,'UserData') };
  else
    state = { 'p' , get( hAx ,'CameraPosition') , get( hS ,'UserData') };
  end
  if isequal( state , getappdata( hS ,'silhouette_lastState') ), return; end

  % ADAPTIVE SAFETY THROTTLE: a real change, but the last full pass was less
  % than max( 33 ms , 3x the MEASURED pass cost ) ago (event burst) -> skip
  % this recompute and let the trailing timer land the exact state ~50 ms
  % after the stream pauses. The adaptive term bounds the silhouette duty
  % cycle at ~1/3 regardless of mesh size; isolated events always pass through
  % synchronously, so interaction still feels real-time.
  lastT = getappdata( hS , 'silhouette_lastUpdate' );
  pc    = getappdata( hS , 'silhouette_passCost' );  if isempty( pc ), pc = 0; end
  if ~isempty( lastT ) && toc( lastT ) < max( 0.033 , 3*pc )
    tD = getappdata( hS , 'silhouette_throttle' );
    if ~isempty( tD ) && isvalid( tD ), stop( tD ); start( tD ); end
    return;
  end
  tPass = tic;

  %SIL = get( hS , 'ApplicationData' );
  %try, if SIL.silhouette_stopped, return; end; end

  hP = getappdata( hS , 'silhouette_parentPatch' );
  if ~isequal( get( hS , 'Parent') , get( hP ,'Parent') )
    hPP = get( hP ,'Parent');
    set( hS , 'Parent' , hPP );
%     silhouette_listeners = getappdata( hS , 'silhouette_listeners' );
%     silhouette_listeners{end+1} = newListener( hPP , 'Matrix' , @(~,~)updateCamera( hS , hAx ) );
%     setappdata( hS , 'silhouette_listeners' , silhouette_listeners );
  end

  SILhouette_normals      = getappdata( hS , 'silhouette_normals' );
  SILhouette_faceCenters  = getappdata( hS , 'silhouette_faceCenters' );
  while 1
    hP = get( hP ,'Parent' ); if strcmp( get(hP,'Type') ,'axes' ), break; end
    if strcmp( get(hP,'Type') ,'hgtransform')
      TR = get( hP , 'Matrix' );
      SILhouette_normals = SILhouette_normals * TR(1:3,1:3).';
      if strcmp( get( hAx , 'Projection' ) , 'perspective' )
        SILhouette_faceCenters = bsxfun( @plus , SILhouette_faceCenters * TR(1:3,1:3).' , TR(1:3,4).' );
      end
    end
  end

  CameraPosition = get( hAx , 'CameraPosition' );
  switch get( hAx , 'Projection' )
    case 'orthographic'
      viewDIR = get( hAx , 'CameraTarget' ) - CameraPosition;
      angles  = SILhouette_normals * viewDIR(:);
    case 'perspective'
      viewDIR = bsxfun( @minus , SILhouette_faceCenters , CameraPosition );
      angles  = dot( SILhouette_normals , viewDIR ,2);
  end
  
  SILhouette_faces    = getappdata( hS , 'silhouette_faces' );
  SILhouette_boundary = getappdata( hS , 'silhouette_boundary' );
  if getappdata( hS , 'silhouette_manifold' )
    % fast path: silhouette edges = cached INTERIOR edges (pre-extracted in
    % updateMesh) whose two faces fall on OPPOSITE sides of the front/back
    % partition. Pure logical indexing, no per-event interior masking.
    Eui = getappdata( hS , 'silhouette_edgesInt' );
    EFi = getappdata( hS , 'silhouette_edgeFacesInt' );
    front = angles > 0;
    wi = front( EFi(:,1) ) ~= front( EFi(:,2) );
    FacesNew = [ SILhouette_boundary ; Eui( wi ,:) ];
  else
    % exact fallback (edges with >2 faces): boundary of the front subset per event
    S = meshBoundary( SILhouette_faces( angles > 0 ,:) );
    FacesNew = unique( [ SILhouette_boundary ; S ] ,'rows');
  end
  % upload only when the edge set actually CHANGED (small rotations often keep
  % it identical; the Faces set is the expensive graphics side of the pass)
  if ~isequal( FacesNew , getappdata( hS ,'silhouette_lastFaces') )
    set( hS , 'Faces' , FacesNew );
    setappdata( hS ,'silhouette_lastFaces' , FacesNew );
  end

  occludedAlpha = get( hS , 'UserData' ); %disp( occludedAlpha )
  if isempty( occludedAlpha ), occludedAlpha = [1,1]*get( hS ,'EdgeAlpha'); end
  if numel( occludedAlpha ) == 1, occludedAlpha(2) = 1; end
  if all( occludedAlpha == occludedAlpha(1) )
    set( hS , 'EdgeAlpha' , occludedAlpha(1) );
  else
    X = get(hS,'Vertices') ;  %mesh used for deep oclussion
    persp = strcmp( get( hAx , 'Projection' ) , 'perspective' );

    % CLOSED mesh (manifold, no boundary): any occluding crossing includes a
    % camera-FACING face, so testing only the front faces (angles<0) suffices --
    % halves the binning/testing work. Open meshes must test every face (a
    % back-facing sheet can be the only occluder).
    ON = getappdata( hS ,'silhouette_normals');
    if getappdata( hS ,'silhouette_manifold') && isempty( SILhouette_boundary )
      wf = angles < 0;
      OF = SILhouette_faces( wf ,:);  ON = ON( wf ,:);
    else
      OF = SILhouette_faces;
    end

    verts = unique( FacesNew );  verts = verts( verts > 0 );
    BLOB = getappdata( hS , 'silhouette_bvhBlob' );
    if ~isempty( BLOB ) || exist( 'IntersectSurfaceRay_mx' , 'file' ) == 3
      % MEX path: ONE occlusion ('any') query per silhouette vertex -- against
      % the explicit per-silhouette BVH blob (bvhIntersectRay_mx; 'noframe' at
      % build time = world-space rays go STRAIGHT to the mex) or, without the
      % BVH engine on the path, against IntersectSurfaceRay_mx's internally
      % cached tree (warmed per mesh in updateMesh). Both engines share the
      % exact same 'any' semantics: some hit strictly INSIDE the segment
      % viewer->vertex (1e-9 < t < 1-1e-5, so the vertex itself never
      % self-reports) = occluded; unordered traversal with early exit, cheaper
      % than 'first' and much cheaper than the old 'all'.
      %   orthographic: P0 = vertex - viewDIR (UNnormalized -> P0 falls outside
      %   the mesh, so every occluder lies at t > 0 and 'any' sees it);
      %   perspective : P0 = camera (hits behind it, t<=0, must not occlude --
      %   exactly what the 'any' window discards).
      Xv = X( verts ,:);
      VA = zeros( size(X,1) ,1) + occludedAlpha(2);
      if ~persp
        P0 = Xv - viewDIR;
      else
        P0 = repmat( CameraPosition , numel(verts) , 1 );
      end

      % SHADOW CACHE (temporal coherence): an occluded vertex tends to stay
      % occluded by the SAME triangle across camera events. Re-test each
      % vertex's last known blocker first (one vectorized Moller-Trumbore,
      % same guards as 'any'); only the unresolved rays go to the mex. A stale
      % blocker is harmless: its test fails and the mex re-decides.
      B = getappdata( hS , 'silhouette_lastBlocker' );
      if numel( B ) ~= size( X ,1), B = zeros( size(X,1) ,1); end
      occ = false( numel(verts) ,1);
      w = find( B( verts ) > 0 );
      if ~isempty( w )
        occ(w) = segmentHitsTri( X , SILhouette_faces( B(verts(w)) ,:) , P0(w,:) , Xv(w,:) );
      end
      q = find( ~occ );
      if ~isempty( q )
        if ~isempty( BLOB )
          [ ~ , cid ] = bvhIntersectRay_mx( [ P0(q,:) , Xv(q,:) ] , BLOB , 4 , 1 );
          hit = cid > 0;
        else
          SS = struct( 'vertices',X , 'faces',double(SILhouette_faces) );
          [ ~ , pid , cid ] = IntersectSurfaceRay_mx( SS , [ P0(q,:) , Xv(q,:) ] , 'any' );
          hit = pid > 0;
        end
        occ( q(hit) ) = true;
        B( verts( q( hit) ) ) = cid( hit );
        B( verts( q(~hit) ) ) = 0;
      end
      setappdata( hS , 'silhouette_lastBlocker' , B );
      VA( verts( occ ) ) = occludedAlpha(1);
    elseif ~persp
      % ORTHOGRAPHIC fallback (no MEX): project everything ONCE onto the
      % transverse plane and bin the faces on a 2D grid; each vertex is tested
      % only against the faces of its cell (exact superset). O(nF + nVerts*cand)
      % instead of the old per-vertex O(nF) ray intersection.
      VA = occlusionOrtho( X , OF , verts , viewDIR , CameraPosition , ON , occludedAlpha );
    else
      % PERSPECTIVE fallback (no MEX): rays differ per vertex; shoot each one
      % from the camera.
      VA = zeros( size(X,1) ,1) + occludedAlpha(2);
      for v = verts(:).'
        x = X(v,:);
        dv = x - CameraPosition;
        try
          y  = intersectSurfaceRay_( X , OF , [ x - dv ; x ] );
          d2x = sqrt( sum( ( CameraPosition - x ).^2 ) );
          d2y = sqrt( sum( ( CameraPosition - y ).^2 ) );
          if ( d2x - d2y )/d2x > 1e-5
            VA(v,:) = occludedAlpha(1);
          end
        end
      end
    end

    set( hS , 'EdgeAlpha' , 'interp' , 'AlphaDataMapping' , 'none' , 'FaceVertexAlphaData' , VA );
  end

  %NO drawnow here: at large meshes a forced render INSIDE the camera callback
  %adds its full cost to the interaction latency, and MATLAB repaints on its
  %own between events (the new silhouette rides that natural render). Scripted
  %camera loops must drawnow themselves.
  setappdata( hS , 'silhouette_lastState'  , state );
  setappdata( hS , 'silhouette_lastUpdate' , tic );      %throttle reference clock
  setappdata( hS , 'silhouette_passCost'   , toc( tPass ) );
end
function visibleChanged( hS , hAx )
%hS Visible listener: when the change is USER-initiated (not the parent sync),
%record it as the desired state -- so re-showing the patch will not resurrect a
%silhouette the user explicitly hid. Then refresh as usual.
  if ~ishandle( hS ), return; end
  if ~isequal( getappdata( hS , 'silhouette_syncingVisible' ) , true )
    setappdata( hS , 'silhouette_userVisible' , get( hS , 'Visible' ) );
  end
  updateCamera( hS , hAx );
end
function syncVisibleFromParent( hP , hS )
%propagate the PARENT patch visibility onto the silhouette: patch off forces the
%silhouette off; patch on restores the USER's desired state (appdata
%'silhouette_userVisible'). The syncing flag tells visibleChanged that this
%change is a propagation, NOT a user decision to be recorded.
  if ~ishandle( hP ) || ~ishandle( hS ), return; end
  if strcmp( get( hP , 'Visible' ) , 'on' )
    nv = getappdata( hS , 'silhouette_userVisible' );
    if isempty( nv ), nv = 'on'; end
  else
    nv = 'off';
  end
  if strcmp( get( hS , 'Visible' ) , nv ), return; end
  setappdata( hS , 'silhouette_syncingVisible' , true );
  try,   set( hS , 'Visible' , nv );  end
  setappdata( hS , 'silhouette_syncingVisible' , false );
end
function trailingUpdate( hS , hAx )
%throttle-timer callback: the event stream paused, run the pass that the
%throttle skipped so the RESTING state is exact. Goes through updateCamera's
%normal guards (dedup + throttle): if the last leading-edge pass was very
%recent it simply re-arms and converges one cycle later.
  if ~ishandle( hS ) || ~ishandle( hAx ), return; end
  updateCamera( hS , hAx );
end
function cleanupSilhouette( hP , hS )
%DeleteFcn of the silhouette: unhook the parent patch and stop/delete the
%throttle timer, so no timer ever outlives its silhouette.
  try, set( hP , 'DeleteFcn' , '' ); end
  try
    tD = getappdata( hS , 'silhouette_throttle' );
    if ~isempty( tD ) && isvalid( tD ), stop( tD ); delete( tD ); end
  end
end
function F = meshBoundary( F )
%edges appearing exactly ONCE among the faces F: the open boundary of a mesh (or,
%applied to the front-facing subset, its boundary = silhouette + open boundary).
%Only used by the non-manifold fallback; the fast path uses the cached edges.
  [u,~,c] = unique( sort( [ F(:,[1,2]) ; F(:,[2,3]) ; F(:,[3,1]) ] ,2) , 'rows' );
  F = u( accumarray( c(:) , 1) == 1 ,:);
end
function N = meshNormals( V , F )
%UNNORMALIZED face normals: only the SIGN of N.viewdir is used, so skipping the
%normalization is free accuracy (no division) and cheaper.
  x = V( F(:,1) ,:);
  N = cross( V( F(:,2) ,:) - x , V( F(:,3) ,:) - x , 2 );
end
function C = meshFacesCenter( V , F )
%face centroids: the per-face ray origins of the PERSPECTIVE classification.
  C = ( V( F(:,1) ,:) +  V( F(:,2) ,:) + V( F(:,3) ,:) )/3;
end

function occ = segmentHitsTri( X , TRI , P0 , P1 )
%one Moller-Trumbore per row: does the open segment P0(i,:)->P1(i,:) hit the
%triangle TRI(i,:)? Same guards as the mex 'any' mode (t in (1e-9,1-1e-5),
%inclusive barycentric tolerance) -- used by the occlusion shadow cache.
  d  = P1 - P0;
  v0 = X(TRI(:,1),:);  e1 = X(TRI(:,2),:)-v0;  e2 = X(TRI(:,3),:)-v0;
  p  = cross( d , e2 ,2);
  dt = dot( e1 , p ,2);
  s  = P0 - v0;
  iv = 1 ./ dt;
  u  = dot( s , p ,2).*iv;
  q  = cross( s , e1 ,2);
  v  = dot( d , q ,2).*iv;
  t  = dot( e2 , q ,2).*iv;
  ee = 1e-9;
  occ = abs(dt) >= 1e-300 & u >= -ee & u <= 1+ee & v >= -ee & u+v <= 1+ee ...
      & t > 1e-9 & t < 1-1e-5;
end
function VA = occlusionOrtho( X , F , verts , viewDIR , CameraPosition , N0 , occludedAlpha )
% per-vertex deep-occlusion alphas for PARALLEL (orthographic) rays. Project all
% vertices once onto the plane transverse to viewDIR, bin the faces on a uniform
% 2D grid with cell size >= the largest face extent (so a face lands in at most
% 4 cells) and test each silhouette vertex only against the faces of its cell:
% 2D point-in-triangle (same sign test as intersectSurfaceRay_) + depth of the
% plane hit along the ray (s < 0 means the hit is towards the camera = occluder).
  VA = zeros( size(X,1) ,1) + occludedAlpha(2);

  [ ~ , k ] = max( abs( viewDIR ) );                  %dominant component first
  pp = [ k , 1:k-1 , k+1:3 ];
  Z  = viewDIR( pp );  Xp = X(:,pp);
  R  = [ -Z(2)/Z(1) , -Z(3)/Z(1) ; 1,0 ; 0,1 ];
  PX = Xp * R(:,1);  PY = Xp * R(:,2);                %transverse coordinates

  F1 = F(:,1); F2 = F(:,2); F3 = F(:,3);
  bx0 = min( min( PX(F1) , PX(F2) ) , PX(F3) );  bx1 = max( max( PX(F1) , PX(F2) ) , PX(F3) );
  by0 = min( min( PY(F1) , PY(F2) ) , PY(F3) );  by1 = max( max( PY(F1) , PY(F2) ) , PY(F3) );

  x0 = min( PX );  y0 = min( PY );
  h  = max( max( bx1 - bx0 ) , max( by1 - by0 ) );    %cell >= largest face extent
  h  = max( h , eps( max( abs( [PX;PY] ) ) + 1 ) * 8 );
  ncx = floor( ( max(PX) - x0 ) / h ) + 1;
  ncy = floor( ( max(PY) - y0 ) / h ) + 1;
  cellx = @(p) min( floor( ( p - x0 ) / h ) , ncx-1 );
  celly = @(p) min( floor( ( p - y0 ) / h ) , ncy-1 );

  cx0 = cellx( bx0 );  cx1 = cellx( bx1 );            %<=2 cells per axis per face
  cy0 = celly( by0 );  cy1 = celly( by1 );
  lin = [ cx0*ncy+cy0 ; cx0*ncy+cy1 ; cx1*ncy+cy0 ; cx1*ncy+cy1 ] + 1;
  fid = repmat( ( 1:size(F,1) ).' , 4 , 1 );
  [ lin , o ] = sort( lin );  fid = fid( o );
  kp  = [ true ; diff(lin) ~= 0 | diff(fid) ~= 0 ];   %dedup (cell,face) pairs
  lin = lin(kp);  fid = fid(kp);
  last   = [ find( diff( lin ) ~= 0 ) ; numel(lin) ];
  cnt    = diff( [ 0 ; last ] );
  cstart = [ 1 ; last(1:end-1) + 1 ];
  map    = sparse( lin(last) , 1 , 1:numel(last) , ncx*ncy , 1 );

  d  = sum( N0 .* X( F1 ,:) , 2 );                    %face plane: N0.p = d
  nu = N0 * viewDIR(:);

  for v = verts(:).'
    px = PX(v);  py = PY(v);
    slot = full( map( ( cellx(px) )*ncy + celly(py) + 1 ) );
    if ~slot, continue; end
    c = fid( cstart(slot) : cstart(slot) + cnt(slot) - 1 );

    ax_ = PX(F1(c)) - px;  ay_ = PY(F1(c)) - py;
    bx_ = PX(F2(c)) - px;  by_ = PY(F2(c)) - py;
    cx_ = PX(F3(c)) - px;  cy_ = PY(F3(c)) - py;
    ins = mod( sign( ( bx_-ax_ ).*ay_ - ( by_-ay_ ).*ax_ ) + ...
               sign( ( cx_-bx_ ).*by_ - ( cy_-by_ ).*bx_ ) + ...
               sign( ( ax_-cx_ ).*cy_ - ( ay_-cy_ ).*cx_ ) , 3 ) == 0;
    if ~any( ins ), continue; end
    c = c( ins );

    x = X(v,:);
    s = ( d(c) - ( N0(c,:) * x(:) ) ) ./ nu(c);       %hit = x + s*viewDIR
    tol = 1e-5 * sqrt( sum( ( CameraPosition - x ).^2 ) ) / sqrt( sum( viewDIR.^2 ) );
    if any( s < -tol ), VA(v) = occludedAlpha(1); end %a hit towards the camera
  end
end

function nL = newListener( hh , prop , fcn )
%HG1 / HG2 property-PostSet listener factory (same idiom as headlight.m). The
%returned listener object must be KEPT REFERENCED (they live in the silhouette's
%appdata) or it stops listening. NB the version scalar: 8.4.0.build maps to
%804.000x > 804, so R2014b+ (HG2) correctly picks event.proplistener.
  persistent matlabV
  if isempty( matlabV )
    matlabV = sscanf(version,'%d.%d.%d.%d.%d',5); matlabV=[100,1,1e-2,1e-9,1e-13]*[ matlabV(1:min(5,end)) ; zeros(5-numel(matlabV),1) ];
  end

  if matlabV <= 804, nL = handle.listener(    handle( hh ) , findprop( handle(hh) , prop ) , 'PropertyPostSet' , fcn );    
  else,              nL = event.proplistener(         hh   , findprop(        hh  , prop ) , 'PostSet'         , fcn );
  end
end

function [ xyz ]  = intersectSurfaceRay_( V , F , ray )
%first intersection of the segment ray(1,:) -> ray(2,:) with the surface (V,F);
%used by the PERSPECTIVE occlusion test (one call per silhouette vertex; the
%orthographic path uses the gridded occlusionOrtho instead). Three culling
%passes in the plane transverse to the ray (X-slab, Y-slab, 2D point-in-triangle
%by orientation signs), then the earliest plane hit among the survivors. ERRORS
%when nothing is hit -- callers wrap it in try (no hit = nothing occludes).
  Z0 = ray(1,:);
  Z  = ray(2,:) - Z0;

  % put the DOMINANT ray component first: R below divides by Z(1), so an
  % axis-aligned ray (e.g. view(0,90) -> Z = [0 0 -d]) used to produce Inf/NaN
  % and silently disable the occlusion test. The permutation is undone at the end.
  [ ~ , k ] = max( abs( Z ) );
  pp = [ k , 1:k-1 , k+1:3 ];
  V  = V(:,pp);  Z0 = Z0(pp);  Z = Z(pp);

  R = [ -Z(2)/Z(1) , -Z(3)/Z(1) ; 1,0 ; 0,1 ];
  Z0R = Z0 * R;

  
  X = ( V * R(:,1) ) - Z0R(1);
  s = [ X(F(:,1)) , X(F(:,2)) , X(F(:,3)) ];
  w = ~( all( s > 0 ,2) | all( s < 0 ,2) );
  F = F( w , : );
  if isempty( F ), error('no Fs 1.'); end

  Y = ( V * R(:,2) ) - Z0R(2);
  s = [ Y(F(:,1)) , Y(F(:,2)) , Y(F(:,3)) ];
  w = ~( all( s > 0 ,2) | all( s < 0 ,2) );
  F = F( w , : );
  if isempty( F ), error('no Fs 1.'); end


  Ax = X(F(:,1));   Ay = Y(F(:,1));
  Bx = X(F(:,2));   By = Y(F(:,2));
  Cx = X(F(:,3));   Cy = Y(F(:,3));
  F = F( mod( ...
           sign( ( Bx-Ax ) .* Ay - ( By-Ay ) .* Ax ) + ...
           sign( ( Cx-Bx ) .* By - ( Cy-By ) .* Bx ) + ... 
           sign( ( Ax-Cx ) .* Cy - ( Ay-Cy ) .* Cx ) , ...
         3 ) == 0 , : );
  if isempty( F ), error('no Fs 2.'); end

  V1 = V( F(:,1) ,:);
  N = cross( V(F(:,2),:) - V1 , V(F(:,3),:) - V1 , 2 );

  D = dot( V1 , N , 2 );
  t = ( D - N * Z0(:) )./( N * Z(:) );
  t = min( t );
  xyz = Z0 + ( t .* Z );
  xyz( pp ) = xyz;                %undo the dominant-component permutation

end

