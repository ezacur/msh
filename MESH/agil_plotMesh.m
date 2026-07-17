function [hP_,hX_] = agil_plotMesh( M , varargin )
%AGIL_PLOTMESH  Plot (or agilize) a huge mesh with a decimated INTERACTION PROXY.
%
%   hP = agil_plotMesh( M )                 plot M (via plotMESH) and attach
%   hP = agil_plotMesh( M , NFACES , ... )  proxy target size (default 75000)
%   hP = agil_plotMesh( hP , ... )          agilize an ALREADY plotted patch
%   [hP,hX] = agil_plotMesh( ... )          also return the proxy patch hX
%
%   While the camera moves, the FULL patch is hidden and a DECIMATED proxy
%   (~NFACES triangles, built once with reducepatch) is shown instead, so the
%   renderer never has to draw millions of triangles inside an interaction
%   gesture; ~0.2 s after the last camera event a single-shot timer restores
%   the full mesh. Isolated camera events behave the same (brief proxy flash,
%   then the exact full render lands at rest).
%
%   M        a mesh struct (plotted with plotMESH; extra arguments after
%            NFACES are forwarded to it) or an existing PATCH handle (extra
%            arguments are set() on it).
%   NFACES   optional numeric scalar: target triangle count of the proxy
%            (default 75000). If the mesh has <= 2*NFACES triangles NO proxy
%            machinery is attached (it would not pay) and the plain patch is
%            returned.
%
%   WHAT THE PROXY MIRRORS -- at every gesture START the proxy re-copies the
%   full patch's scalar appearance (FaceColor/EdgeColor/alphas/lighting/
%   material/line properties) and, when the full patch carries PER-VERTEX
%   FaceVertexCData / FaceVertexAlphaData, gathers them onto the proxy
%   vertices through a nearest-vertex map (knnsearch, cached at build time --
%   the per-gesture gather is instantaneous and reflects the CURRENT data,
%   e.g. a live EP map). PER-FACE CData is NOT mapped (no face correspondence
%   survives the decimation): the proxy then falls back to a neutral flat
%   color.
%
%   LIFETIME
%     - deleting hX detaches everything and restores the full patch;
%     - deleting hP deletes hX (and the timer / listeners);
%     - re-attaching (calling agil_plotMesh(hP) again) replaces cleanly;
%     - editing hP's Vertices/Faces rebuilds the proxy (eagerly: reducepatch
%       runs at edit time, never inside a camera callback);
%     - hiding hP yourself disables the proxy until you show it again.
%
%   NOTES
%     - the proxy lives in the same Parent (hgtransform-safe) with
%       PickableParts 'none' and hidden handle;
%     - scripted figure export: view(...); print(...) back to back may catch
%       the proxy frame -- do  pause(0.25); drawnow;  before printing (same
%       caveat as the silhouette/backFaceCulling trackers);
%     - camera trackers attached to the full patch (silhouette, shading, ...)
%       follow its visibility automatically and pause during the gesture.
%
%   Example:
%     M = sphereMesh( 8 );                    % ~1.3M triangles
%     agil_plotMesh( M , 'FaceColor',[.8 .8 .9] , 'EdgeColor','none' );
%     % rotate: the gesture moves a 75k proxy; the 1.3M mesh lands at rest
%
%   See also plotMESH, reducepatch, silhouette, backFaceCullingSplit.

if 0

  M = sphereMesh( 8 );                                   %#ok<UNRCH> ~1.3M tris
  clf
  hP = agil_plotMesh( M , 'FaceColor',[.8 .8 .9] ,'EdgeColor','none','nice');
  headlight

%%
end

  TRAIL = 0.20;                 %seconds after the last camera event -> full mesh

  %-- optional leading numeric scalar = proxy target faces
  NFACES = 75000;
  if numel( varargin ) && isnumeric( varargin{1} ) && isscalar( varargin{1} )
    NFACES = double( varargin{1} );  varargin(1) = [];
    if ~isfinite( NFACES ) || NFACES < 100
      error( 'agil_plotMesh:nfaces' , 'NFACES must be a finite scalar >= 100.' );
    end
  end

  %-- the FULL patch: plot the struct, or take the given patch
  if isstruct( M )
    hP = plotMESH( M , varargin{:} );
    if ~strcmp( get( hP ,'Type') ,'patch' )
      hP = findobj( hP , 'Type','patch' );
      if isempty( hP ), error( 'agil_plotMesh:plot' , 'plotMESH did not return a patch.' ); end
      hP = hP(1);
    end
  elseif isgraphics( M ) && strcmp( get( M ,'Type') ,'patch' )
    hP = M;
    if numel( varargin ), set( hP , varargin{:} ); end
  else
    error( 'agil_plotMesh:input' , 'A mesh struct or a patch handle was expected.' );
  end
  if size( get( hP ,'Faces') ,2) ~= 3
    error( 'agil_plotMesh:tri' , 'Only triangular meshes are allowed.' );
  end

  %-- re-attach: deleting the previous proxy detaches its whole machinery
  try, delete( getappdata( hP , 'agil_proxy' ) ); end

  %-- small meshes: not worth a proxy
  if size( get( hP ,'Faces') ,1) <= 2*NFACES
    if nargout,     hP_ = hP;          end
    if nargout > 1, hX_ = gobjects(0); end
    return;
  end

  hAx = ancestor( hP ,'axes' );

  %-- the proxy patch: same parent (hgtransform-safe), invisible, unpickable
  hX = patch( 'Parent'          , get( hP ,'Parent') , ...
              'Vertices'        , [] , 'Faces' , [] , ...
              'Visible'         , 'off' , ...
              'HandleVisibility', 'off' , ...
              'PickableParts'   , 'none' , ...
              'Tag'             , 'agil_proxy' );

  setappdata( hP , 'agil_proxy'   , hX );
  setappdata( hP , 'agil_nfaces'  , NFACES );
  setappdata( hP , 'agil_inProxy' , false );
  setappdata( hX , 'agil_full'    , hP );

  rebuildProxy( hP );                                  %reducepatch + knnsearch map

  %-- trailing timer: lands the full mesh when the gesture pauses
  tD = timer( 'ExecutionMode','singleShot' , 'StartDelay', TRAIL , ...
              'Name','agil_trail' , 'ObjectVisibility','off' , ...
              'TimerFcn', @(~,~)exitProxy( hP ) );
  setappdata( hP , 'agil_timer' , tD );

  setappdata( hP , 'agil_listeners' , { ...
      newListener( hAx , 'CameraPosition'  , @(~,~)onCamera( hP ) ) ;...
      newListener( hAx , 'CameraTarget'    , @(~,~)onCamera( hP ) ) ;...
      newListener( hAx , 'CameraViewAngle' , @(~,~)onCamera( hP ) ) ;...
      newListener( hAx , 'Projection'      , @(~,~)onCamera( hP ) ) ;...
      newListener( hAx , 'View'            , @(~,~)onCamera( hP ) ) ;...
      newListener( hP  , 'Vertices'        , @(~,~)rebuildProxy( hP ) ) ;...
      newListener( hP  , 'Faces'           , @(~,~)rebuildProxy( hP ) ) ;...
      } );

  set( hX , 'DeleteFcn' , @(~,~)agil_detach( hX ) );   %delete(hX) = undo
  try, addlistener( hP , 'ObjectBeingDestroyed' , @(~,~)killAll( hP ) ); end

  if nargout,     hP_ = hP; end
  if nargout > 1, hX_ = hX; end
end


function rebuildProxy( hP )
%(re)build the decimated proxy: reducepatch to the target count + the nearest-
%vertex map used to gather per-vertex data. Runs at attach and EAGERLY on
%Vertices/Faces edits -- never inside a camera callback (reducepatch on
%millions of triangles takes seconds; camera passes must stay instant).
  if ~ishandle( hP ), return; end
  hX = getappdata( hP , 'agil_proxy' );
  if isempty( hX ) || ~ishandle( hX ), return; end

  V = get( hP ,'Vertices' );
  F = get( hP ,'Faces' );
  [ Fp , Vp ] = reducepatch( F , V , getappdata( hP ,'agil_nfaces' ) );
  set( hX , 'Vertices' , Vp , 'Faces' , Fp );
  setappdata( hP , 'agil_idx' , knnsearch( V , Vp ) );  %proxy vtx -> nearest full vtx
end

function onCamera( hP )
%camera event: first one of a gesture swaps full -> proxy (after mirroring the
%CURRENT appearance); every one re-arms the trailing timer.
  if ~ishandle( hP ), return; end
  hX = getappdata( hP , 'agil_proxy' );
  if isempty( hX ) || ~ishandle( hX ), return; end

  if ~isequal( getappdata( hP , 'agil_inProxy' ) , true )
    if ~strcmp( get( hP ,'Visible') ,'on' ), return; end    %user hid the mesh
    mirrorAppearance( hP , hX );
    set( hX , 'Visible' , 'on'  );
    set( hP , 'Visible' , 'off' );
    setappdata( hP , 'agil_inProxy' , true );
  end

  % re-arm the trailing timer, but at most every 50 ms: a stop/start pair
  % costs ~1.5 ms and a drag fires ~60 events/s -- skipping re-arms inside
  % that window only stretches the effective trail to <= TRAIL+50 ms.
  la = getappdata( hP , 'agil_lastArm' );
  if ~isempty( la ) && toc( la ) < 0.05, return; end
  tD = getappdata( hP , 'agil_timer' );
  if ~isempty( tD ) && isvalid( tD ), stop( tD ); start( tD ); end
  setappdata( hP , 'agil_lastArm' , tic );
end

function exitProxy( hP )
%trailing-timer callback: the gesture paused -> restore the full mesh.
  if ~ishandle( hP ), return; end
  hX = getappdata( hP , 'agil_proxy' );
  if ~isequal( getappdata( hP , 'agil_inProxy' ) , true ), return; end
  setappdata( hP , 'agil_inProxy' , false );
  try, set( hP , 'Visible' , 'on'  ); end
  if ~isempty( hX ) && ishandle( hX ), set( hX , 'Visible' , 'off' ); end
end

function mirrorAppearance( hP , hX )
%copy the full patch's CURRENT look onto the proxy: scalar appearance props,
%plus per-vertex color/alpha data gathered through the cached nearest map
%(instantaneous: one 75k-row indexing). Per-FACE data cannot be mapped (no
%face correspondence survives reducepatch) -> neutral flat color instead.
  for p = { 'EdgeColor','LineStyle','LineWidth','Marker','MarkerSize' , ...
            'FaceLighting','EdgeLighting','BackFaceLighting' , ...
            'AmbientStrength','DiffuseStrength','SpecularStrength' , ...
            'SpecularExponent','SpecularColorReflectance' , ...
            'AlphaDataMapping','CDataMapping' }
    try, set( hX , p{1} , get( hP , p{1} ) ); end
  end

  idx = getappdata( hP , 'agil_idx' );
  nV  = size( get( hP ,'Vertices') ,1);

  FC = get( hP ,'FaceColor' );
  CD = get( hP ,'FaceVertexCData' );
  if ischar( FC ) && any( strcmp( FC , {'interp','flat'} ) )
    if ~isempty( CD ) && size( CD ,1) == nV
      set( hX , 'FaceVertexCData' , CD( idx ,:) , 'FaceColor' , FC );
    else                                         %per-face data: not mappable
      set( hX , 'FaceColor' , [0.7 0.7 0.7] );
    end
  else
    try, set( hX , 'FaceColor' , FC ); end
  end

  FA = get( hP ,'FaceAlpha' );
  AD = get( hP ,'FaceVertexAlphaData' );
  if ischar( FA ) && any( strcmp( FA , {'interp','flat'} ) )
    if ~isempty( AD ) && size( AD ,1) == nV
      set( hX , 'FaceVertexAlphaData' , AD( idx ,:) , 'FaceAlpha' , FA );
    else
      set( hX , 'FaceAlpha' , 1 );
    end
  else
    try, set( hX , 'FaceAlpha' , FA ); end
  end
end

function agil_detach( hX )
%DeleteFcn of the proxy = UNDO: restore the full patch, drop the machinery.
  hP = getappdata( hX , 'agil_full' );
  if ~ishandle( hP ), return; end
  try
    tD = getappdata( hP , 'agil_timer' );
    if ~isempty( tD ) && isvalid( tD ), stop( tD ); delete( tD ); end
  end
  try, set( hP , 'Visible' , 'on' ); end
  try, rmappdata( hP , 'agil_proxy'     ); end
  try, rmappdata( hP , 'agil_idx'       ); end
  try, rmappdata( hP , 'agil_nfaces'    ); end
  try, rmappdata( hP , 'agil_inProxy'   ); end
  try, rmappdata( hP , 'agil_timer'     ); end
  try, rmappdata( hP , 'agil_lastArm'   ); end
  try, rmappdata( hP , 'agil_listeners' ); end
end

function killAll( hP )
%full patch dying: take the proxy (whose DeleteFcn cleans the rest) with it.
  try
    hX = getappdata( hP , 'agil_proxy' );
    if ~isempty( hX ) && ishandle( hX ), delete( hX ); end
  end
  try
    tD = getappdata( hP , 'agil_timer' );
    if ~isempty( tD ) && isvalid( tD ), stop( tD ); delete( tD ); end
  end
end

function nL = newListener( hh , prop , fcn )
%HG1 / HG2 property-PostSet listener factory (same idiom as silhouette.m). The
%returned listeners must be KEPT REFERENCED (they live in the patch appdata).
  persistent matlabV
  if isempty( matlabV )
    matlabV = sscanf(version,'%d.%d.%d.%d.%d',5); matlabV=[100,1,1e-2,1e-9,1e-13]*[ matlabV(1:min(5,end)) ; zeros(5-numel(matlabV),1) ];
  end
  if matlabV <= 804, nL = handle.listener(    handle( hh ) , findprop( handle(hh) , prop ) , 'PropertyPostSet' , fcn );
  else,              nL = event.proplistener(         hh   , findprop(        hh  , prop ) , 'PostSet'         , fcn );
  end
end
