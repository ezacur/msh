function hP_ = silhouetteByShading( hP , shadingFCN , MODE )
%SILHOUETTEBYSHADING  View-dependent alpha shading of a patch (sketch look),
%                     auto-updated with the camera.
%
%   silhouetteByShading( hP )               default 'sketch' shading
%   silhouetteByShading( hP , style )       style: 'sketch' | 'fresnel' |
%                                           'glass0' | 'glass' | 'dirtyglass'
%   silhouetteByShading( hP , @(n)... )     custom shading: receives
%                                           n = abs( cos( normal , viewdir ) )
%                                           per vertex (in [0,1]) and returns
%                                           the per-vertex FaceAlpha.
%   silhouetteByShading( hP , style ,'opaque' )   OPAQUE mode: the
%       alpha-over-background blend is baked into per-vertex COLORS
%       ( FaceVertexCData = a*FaceColor + (1-a)*background , FaceAlpha 1 ), so
%       the renderer skips the depth-sorted transparency of FaceAlpha 'interp'
%       -- MUCH cheaper per frame on big meshes. Identical look wherever the
%       mesh does not overlap itself on screen; where it does, you see the
%       front surface only (the see-through 'glass' styles lose their point:
%       use it with 'sketch' / 'fresnel'). The patch FaceColor and the axes
%       (or figure) background are captured AT ATTACH time. Default mode:
%       'alpha' (the classic translucent behavior).
%
%   The patch is set to FaceAlpha 'interp' and its FaceVertexAlphaData is
%   recomputed from the CURRENT camera: n = VertexNormals . viewdir (one global
%   direction in orthographic, per-vertex camera directions in perspective),
%   mapped through the shading function. With BackFaceLighting 'unlit' the
%   back-facing vertices get alpha 0. Listeners on the camera keep it live.
%
%   AGILITY (same ideas as silhouette / backFaceCullingSplit):
%     - dedup guard on the EFFECTIVE inputs: orthographic shading depends ONLY
%       on the normalized view direction (quantized at 1e-12 to absorb float
%       noise), so zoom / dolly / roll recompute NOTHING and pan collapses to
%       at most one throttled recompute (its two property events go through an
%       intermediate direction); perspective depends only on CameraPosition
%       (zoom and CameraTarget changes are free). The guard also absorbs the
%       double View+CameraPosition firing per gesture. (It replaces the old
%       figure-UpdateToken guard, which skipped exactly the FIRST event after
%       each render: a single programmatic view() change left the shading
%       STALE until the next event.)
%     - ADAPTIVE SAFETY THROTTLE: event bursts recompute at most every
%       max( 33 ms , 3x the measured pass cost ) -- the shading duty cycle is
%       bounded at ~1/3 regardless of mesh size; a single-shot trailing timer
%       (~50 ms) lands the exact resting state. Isolated events stay
%       synchronous (real time).
%     - the numeric pipeline runs in SINGLE (normals cached as single, alphas
%       stored as single by the patch): ~30% less memory traffic, alpha
%       precision 1e-7 (far beyond visible). Custom shadingFCNs receive a
%       single n and their arithmetic propagates it -- fine for any usual math.
%     - no forced drawnow inside the camera callbacks: at large meshes it
%       serialized a full render into the interaction latency; MATLAB repaints
%       naturally between events.
%
%   Editing Vertices / Faces / VertexNormals / BackFaceLighting forces the next
%   pass. Re-attaching (calling this again on the same patch) replaces the
%   previous shading cleanly. NOTE: hgtransform parents are NOT accounted for
%   (unlike silhouette): normals are used as stored.
%
%   See also silhouette, backFaceCullingSplit, plotMESH.

if 0

[n,y] = meshgrid( linspace(1,50,100) , linspace(1,50,100) );

M = Mesh( [n(:),y(:), vec( peaks( n(:) , y(:) ) ) ] , delaunay(n,y) );

silhouetteByShading( plotMESH( M ,'k') ); axis('off');set(gcf,'Color','w')

%%
end

  if ~isgraphics( hP ) || ~strcmp( get( hP ,'Type') ,'patch')
    error('Patch graphic objects were expected.');
  end
  if size( get( hP , 'Faces' ) ,2) ~= 3
    error('Only triangular meshes are allowed.');
  end

  if nargin < 2
    shadingFCN = 'sketch';
  end
  if ischar( shadingFCN )
    switch lower( shadingFCN )
      case 'sketch'
        shadingFCN = @(n)power( round( cos( n * pi/2 ) * 5 )/5 ,8);
        %C = @(x)( ( round( ( abs( cos( abs(x)*pi/2 ) ).^1.0 .* ( x > 0 ) ) * 5 )/5 ).^8 ) * 1;
        %n = C( n );

      case 'fresnel'
        shadingFCN = @(n) power( 1 - n ,4);

      case 'glass0'
        below = @(x,t) x .* ( x > t );
        shadingFCN = @(n) min( 1 - below( 1-below( power( 1 - n ,4) ,0.2) ,0.4) , 0.9 );

      case 'glass'
        shadingFCN = @(n) ( 1/2 - sin( ( n - 0.5 )*pi )/2 ).^8;

      case 'dirtyglass'
        shadingFCN = @(n) max( 1/2 - sin( ( n - 0.5 )*pi )/2 ,0.2);

    end
  end
  if ~isa( shadingFCN ,'function_handle')
    error('A shadingFCN was expected.');
  end

  if nargin < 3, MODE = 'alpha'; end
  MODE = lower( MODE );
  if ~any( strcmp( MODE , {'alpha','opaque'} ) )
    error( 'silhouetteByShading:badMode' , 'MODE must be ''alpha'' (default) or ''opaque''.' );
  end

  hAx = ancestor( hP ,'axes'); hFig = ancestor( hAx ,'figure');
  set( hAx ,'ALim',[0,1]);

  %normals FIRST: Mesh(hP) must read the patch BEFORE the opaque mode switches
  %FaceColor to 'interp' with a still-empty FaceVertexCData (it would warn)
  if isempty( get( hP ,'VertexNormals') )
    set( hP ,'VertexNormals' , meshNormals( Mesh( hP ) ,'a') );
  end

  %re-attach in the OTHER mode: undo what the previous mode took over
  prevMODE = getappdata( hP ,'shading_MODE' );
  if ~isempty( prevMODE ) && ~strcmp( prevMODE , MODE )
    if strcmp( prevMODE ,'opaque' )
      try, set( hP ,'FaceColor', getappdata( hP ,'shading_FaceColor0' ) ); end
    else
      set( hP ,'FaceAlpha',1 );
    end
  end
  setappdata( hP ,'shading_MODE' , MODE );

  set( hP ,'FaceLighting','none');
  if strcmp( MODE ,'alpha' )
    set( hP ,'AlphaDataMapping','none');
    set( hP ,'FaceAlpha','interp');
  else
    %OPAQUE: capture the blend endpoints. FaceColor may already be 'interp'
    %when re-attaching in opaque mode: keep the stored C1 from the previous
    %attach in that case.
    FC = get( hP ,'FaceColor' );
    if isnumeric( FC ) && numel( FC ) == 3
      setappdata( hP ,'shading_FaceColor0' , FC );
      setappdata( hP ,'shading_C1' , single( FC(:).' ) );
    elseif isempty( getappdata( hP ,'shading_C1' ) )
      setappdata( hP ,'shading_FaceColor0' , [0 0 0] );
      setappdata( hP ,'shading_C1' , single( [0 0 0] ) );
    end
    BG = get( hAx ,'Color' );
    if ~isnumeric( BG ) || numel( BG ) ~= 3, BG = get( hFig ,'Color' ); end
    if ~isnumeric( BG ) || numel( BG ) ~= 3, BG = [1 1 1]; end
    setappdata( hP ,'shading_C0' , single( BG(:).' ) );
    set( hP ,'FaceAlpha',1 ,'FaceColor','interp' );
  end

  % SAFETY THROTTLE timer (same pattern as silhouette / backFaceCullingSplit):
  % isolated camera events update synchronously; a stream faster than ~30 Hz
  % skips the in-between recomputes and this single-shot timer lands the exact
  % shading ~50 ms after the stream pauses. Re-attaching kills the previous one;
  % the ObjectBeingDestroyed listener below kills it with the patch.
  tOld = getappdata( hP , 'shading_throttle' );
  if ~isempty( tOld ), try, stop( tOld ); delete( tOld ); end; end
  tD = timer( 'ExecutionMode','singleShot' , 'StartDelay',0.05 , ...
              'Name','shading_throttle' , 'ObjectVisibility','off' , ...
              'TimerFcn', @(~,~)updateShading( ) );
  setappdata( hP , 'shading_throttle' , tD );

  setappdata( hP , 'shading_FCN'        , shadingFCN );
  setappdata( hP , 'shading_lastState'  , [] );
  setappdata( hP , 'shading_lastUpdate' , [] );
  setappdata( hP , 'shading_lastN'      , [] );   %last UPLOADED mapped values
  setappdata( hP , 'shading_passCost'   , 0  );   %measured cost of the last pass
  setappdata( hP , 'shading_dirty'      , 2  );   %passes left re-reading normals

  setappdata( hP , 'shading_listeners' , {  newListener( hAx , 'CameraPosition'       , @(~,~)updateShading( ) ) ;...
                                            newListener( hAx , 'CameraTarget'         , @(~,~)updateShading( ) ) ;...
                                            newListener( hAx , 'CameraViewAngle'      , @(~,~)updateShading( ) ) ;...
                                            newListener( hAx , 'Projection'           , @(~,~)updateShading( ) ) ;...
                                            newListener( hAx , 'View'                 , @(~,~)updateShading( ) ) ;...
                                            newListener( hP  , 'Vertices'             , @(~,~)meshChanged( ) )   ;...
                                            newListener( hP  , 'Faces'                , @(~,~)meshChanged( ) )   ;...
                                            newListener( hP  , 'VertexNormals'        , @(~,~)meshChanged( ) )   ;...
                                            newListener( hP  , 'BackFaceLighting'     , @(~,~)meshChanged( ) )   ;...
                                            newListener( hP  , 'Visible'              , @(~,~)updateShading( ) ) ;...
                                         } );
  try, addlistener( hP , 'ObjectBeingDestroyed' , @(~,~)killThrottle( ) ); end
  try, updateShading( ); end

  function killThrottle()
    try
      t = getappdata( hP , 'shading_throttle' );
      if ~isempty( t ) && isvalid( t ), stop( t ); delete( t ); end
    end
  end

  function meshChanged()
  %mesh-ish change (Vertices/Faces/VertexNormals/BackFaceLighting): the cached
  %state no longer proves anything -- clear it so the next pass cannot be
  %dedup-skipped, and mark the normals cache dirty for TWO passes (auto-mode
  %normals may be recomputed by MATLAB only at the next render, i.e. AFTER our
  %first re-read: the second pass heals that, same self-healing the un-cached
  %per-event get() used to give), then refresh.
    setappdata( hP , 'shading_dirty'     , 2  );
    setappdata( hP , 'shading_lastState' , [] );
    updateShading( );
  end

  function updateShading()
    if ~ishandle( hP ), return; end                   %the trailing timer may
    if ~ishandle( hAx ), return; end                  %outlive a deleted patch
    if strcmp( get( hP ,'Visible') ,'off'), return; end

    % dedup guard on the EFFECTIVE inputs of the shading: in ORTHOGRAPHIC the
    % alphas depend only on the normalized view direction -- zoom, pan, dolly
    % and roll are skipped entirely; in PERSPECTIVE only on CameraPosition --
    % zoom is skipped. This also absorbs the double View+CameraPosition firing
    % per gesture. (The old guard compared the figure's UpdateToken, which
    % skipped exactly the FIRST event after each render: one programmatic
    % view() change left the shading STALE until the next event.)
    CameraPosition = get( hAx , 'CameraPosition' );
    if strcmp( get( hAx ,'Projection') ,'orthographic')
      viewDirection = CameraPosition - get( hAx , 'CameraTarget' );
      viewDirection = viewDirection ./ sqrt( sum( viewDirection.^2 ) );
      %quantize: a dolly recomputes the SAME direction with last-ulp noise
      %(x/|x| vs 2.5x/|2.5x|) and would beat the guard; a cosine change of
      %1e-12 is invisible, so snapping there makes the skip robust.
      viewDirection = round( viewDirection * 1e12 ) * 1e-12;
      state = { 'o' , viewDirection };
    else
      state = { 'p' , CameraPosition };
    end
    if isequal( state , getappdata( hP ,'shading_lastState') ), return; end

    % ADAPTIVE SAFETY THROTTLE: a real change but the last full pass was less
    % than max( 33 ms , 3x the MEASURED pass cost ) ago (event burst) -> skip
    % and let the trailing timer land the resting state. The adaptive term
    % bounds the shading duty cycle at ~1/3 by construction: a mesh whose pass
    % costs 30 ms recomputes at most every 90 ms instead of eating the frame.
    lastT = getappdata( hP , 'shading_lastUpdate' );
    if ~isempty( lastT ) && toc( lastT ) < max( 0.033 , 3 * getappdata( hP ,'shading_passCost' ) )
      t = getappdata( hP , 'shading_throttle' );
      if ~isempty( t ) && isvalid( t ), stop( t ); start( t ); end
      return;
    end
    tPass = tic;

    % vertex normals, cached as SINGLE (~30% less memory traffic through the
    % whole GEMV/FCN/upload chain; the patch stores single alphas natively and
    % an alpha precision of 1e-7 is far beyond visible). Re-read while 'dirty'
    % (two passes after every mesh-ish event, see meshChanged).
    dirty = getappdata( hP , 'shading_dirty' );
    if dirty > 0
      NORMALS  = get( hP , 'VertexNormals' );
      vertices = get( hP , 'Vertices' );
      if size( vertices ,1) ~= size( NORMALS ,1), return; end
      NORMALS = single( NORMALS );
      setappdata( hP , 'shading_NORMALS'  , NORMALS );
      setappdata( hP , 'shading_VERTICES' , single( vertices ) );  %cached too: the
      setappdata( hP , 'shading_dirty'    , dirty - 1 );           %perspective pass
    else                                                           %used to re-read +
      NORMALS = getappdata( hP , 'shading_NORMALS' );              %re-convert nV x 3
    end                                                            %on EVERY event

    if state{1} == 'o'
      n = NORMALS * single( state{2}(:) );
    else
      D = single( CameraPosition ) - getappdata( hP , 'shading_VERTICES' );
      %normalization FUSED into the scalar dot: dividing the nV x 1 dot by
      %|D| costs one pass over nV instead of normalizing all 3 components
      %of D (one nV x 3 divide + temporary saved per event)
      n = sum( NORMALS .* D ,2) ./ sqrt( sum( D.^2 ,2) );
    end
    %n0 = n;
    w = n < 0;

    n = abs(n);
    n = feval( getappdata( hP , 'shading_FCN' ) , n );

    if isequal( get( hP ,'BackFaceLighting') ,'unlit')
      n( w ) = 0;
    end
    % upload only when the mapped values actually CHANGED: the quantized
    % 'sketch' preset lands on 6 discrete levels, so slow rotations often
    % produce an IDENTICAL vector and skip the upload + attribute rebuild
    % (smooth FCNs almost never hit this; the compare costs ~0.2 ms).
    if ~isequal( n , getappdata( hP ,'shading_lastN' ) )
      setappdata( hP ,'shading_lastN' , n );
      if isequal( getappdata( hP ,'shading_MODE' ) ,'opaque' )
        C0 = getappdata( hP ,'shading_C0' );
        C1 = getappdata( hP ,'shading_C1' );
        set( hP , 'FaceVertexCData' , n*(C1-C0) + C0 );   %a*C1 + (1-a)*C0
      else
        set( hP , 'FaceVertexAlphaData' , n );
      end
    end
    %NO drawnow here: at 350k translucent tris a forced render INSIDE the
    %camera callback adds its full cost to the interaction latency, and MATLAB
    %repaints on its own between events anyway (the new alphas ride that
    %natural render). Scripted camera loops must drawnow themselves.

    setappdata( hP , 'shading_lastState'  , state );
    setappdata( hP , 'shading_lastUpdate' , tic );
    setappdata( hP , 'shading_passCost'   , toc( tPass ) );
  end

  if nargout, hP_ = hP; end

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
