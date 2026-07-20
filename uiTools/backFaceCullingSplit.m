function hB_ = backFaceCullingSplit( hF , varargin )
%BACKFACECULLINGSPLIT  Split a patch into camera-facing and back-facing halves.
%
%   hB = backFaceCullingSplit( hF )
%   hB = backFaceCullingSplit( hF , 'Prop',val , ... )
%
%   hF keeps only its camera-FACING faces and a NEW patch hB (a clone of hF's
%   properties, then the given overrides) receives the BACK-facing ones; the
%   partition follows the camera (listeners on the axes, both projections).
%   E.g.  backFaceCullingSplit( hF ,'FaceColor','y','EdgeColor','none' )
%   paints the inside of an open shell yellow. While split, get(hF,'Faces')
%   returns only the front subset (the full set is restored on un-split).
%
%   LIFETIME -- deleting hB UNDOES the split (hF gets its full Faces back);
%   deleting hF deletes hB. Editing hF's Vertices re-caches and repartitions.
%
%   PERFORMANCE
%     - dedup guard on the EFFECTIVE inputs: in orthographic the partition
%       depends only on the (quantized) view direction, so zoom / dolly / roll
%       skip even the classification GEMV and a pan collapses to one throttled
%       pass; in perspective it depends only on CameraPosition (zoom and
%       CameraTarget changes are free). The guard also absorbs the double
%       View+CameraPosition firing per gesture;
%     - per (surviving) camera event one GEMV classifies the faces; the
%       (heavy) two Faces re-uploads happen only when the partition actually
%       CHANGED;
%     - while hB is INVISIBLE (the backFaceCulling case) its Faces are not
%       updated at all; they refresh automatically when it is shown;
%     - an ADAPTIVE SAFETY THROTTLE caps recomputes at max( 33 ms , 3x the
%       measured pass cost ) during event bursts (duty cycle ~1/3 at any mesh
%       size), with a trailing pass (single-shot timer owned by the split,
%       deleted with it) landing the exact state when the stream pauses;
%     - no forced drawnow inside the camera callbacks: MATLAB repaints
%       naturally between events (scripted camera loops must drawnow).
%
%   See also backFaceCulling, silhouette, plotMESH.

if 0

  M = MeshClip( sphereMesh(5) , eye(4) );
  clf
  hF = plotMESH( M ,'b','EdgeColor','k','nice','gouraud');
  hB = backFaceCullingSplit( hF ,'FaceColor','y','EdgeColor','none');

%%

end

  if size( get( hF , 'Faces' ) ,2) ~= 3
    error('Only triangular meshes are allowed.');
  end

  try
    delete( getappdata( hF , 'bfc_backPatch' ) );   %re-split: undo a previous one
  end

  hB = patch( 'Parent'   , get(hF,'Parent') , ...
              'Vertices' , []  ,...
              'Faces'    , []  );
  props = get( hF );
  for p = fieldnames(props).', p = p{1};
    if any( strcmp( p , {'Faces','Vertices','XData','YData','ZData','Children','Parent','HandleVisibility','CreateFcn','DeleteFcn','Type'} ) ), continue; end
    try, set( hB , p , props.(p) ); end
  end
  set( hB , 'HandleVisibility' ,'off' );
  if numel(varargin), set( hB , varargin{:} ); end

  setappdata( hB , 'bfc_frontalPatch' , hF );
  set( hB , 'DeleteFcn' , @(~,~)bfc_delete( hB ) );

  setappdata( hF , 'bfc_backPatch' , hB );
  set( hF , 'DeleteFcn' , @(~,~)delete( hB ) );

  matlabV = sscanf(version,'%d.%d.%d.%d.%d',5); matlabV=[100,1,1e-2,1e-9,1e-13]*[ matlabV(1:min(5,end)) ; zeros(5-numel(matlabV),1) ];
  if matlabV <= 804, newListener = @(hh,prop,fcn)handle.listener(    handle( hh ) , findprop( handle(hh) , prop ) , 'PropertyPostSet' , fcn );
  else,              newListener = @(hh,prop,fcn)event.proplistener(         hh   , findprop(        hh  , prop ) , 'PostSet'         , fcn );
  end

  hAx = ancestor( hB , 'axes' );

  if isnumeric( hAx ), hAx = handle( hAx ); end
  if isnumeric( hF  ), hF  = handle( hF  ); end
  if isnumeric( hB  ), hB  = handle( hB  ); end

  % SAFETY THROTTLE timer (same pattern as silhouette): isolated camera events
  % repartition synchronously; a stream faster than ~30 Hz skips the in-between
  % recomputes and this single-shot timer lands the exact partition ~50 ms
  % after the stream pauses. bfc_delete stops/deletes it.
  tOld = getappdata( hF , 'bfc_throttle' );
  if ~isempty( tOld ), try, stop( tOld ); delete( tOld ); end; end
  tD = timer( 'ExecutionMode','singleShot' , 'StartDelay',0.05 , ...
              'Name','bfc_throttle' , 'ObjectVisibility','off' , ...
              'TimerFcn', @(~,~)trailingUpdate( hF , hAx ) );
  setappdata( hF , 'bfc_throttle' , tD );
  setappdata( hF , 'bfc_passCost' , 0 );    %measured cost of the last pass
                                            %(feeds the ADAPTIVE throttle)

  try, updateMesh( hF ); end
  setappdata( hF , 'bfc_listeners' , {  newListener( hAx , 'CameraPosition'  , @(~,~)updateCamera( hF , hAx ) ) ;...
                                        newListener( hAx , 'CameraTarget'    , @(~,~)updateCamera( hF , hAx ) ) ;...
                                        newListener( hAx , 'CameraViewAngle' , @(~,~)updateCamera( hF , hAx ) ) ;...
                                        newListener( hAx , 'Projection'      , @(~,~)updateCamera( hF , hAx ) ) ;...
                                        newListener( hAx , 'View'            , @(~,~)updateCamera( hF , hAx ) ) ;...
                                        newListener( hF  , 'Vertices'        , @(~,~)updateMesh( hF ) )         ;...
                                        newListener( hB  , 'Visible'         , @(~,~)backShown( hF ) )          ;...
                                    } );

  if nargout > 0, hB_ = hB; end
end


function bfc_delete( hB )
%DeleteFcn of the back patch = UNDO the split: restore the frontal patch's full
%Faces, drop every cache/listener/timer. Also reached when hF is deleted (its
%DeleteFcn deletes hB) -- then the restores on the dying hF just no-op in try.
  hF = getappdata( hB , 'bfc_frontalPatch' );
  delete( hB );

  try
    tD = getappdata( hF , 'bfc_throttle' );
    if ~isempty( tD ) && isvalid( tD ), stop( tD ); delete( tD ); end
  end

  F = getappdata( hF , 'bfc_faces' );

  try, rmappdata( hF , 'bfc_backPatch' );     end
  try, rmappdata( hF , 'bfc_listeners' );     end
  try, rmappdata( hF , 'bfc_faces' );         end
  try, rmappdata( hF , 'bfc_center' );        end
  try, rmappdata( hF , 'bfc_normalsS' );      end
  try, rmappdata( hF , 'bfc_NdotCrel' );      end
  try, rmappdata( hF , 'bfc_partition' );     end
  try, rmappdata( hF , 'bfc_backStale' );     end
  try, rmappdata( hF , 'bfc_lastState' );     end
  try, rmappdata( hF , 'bfc_lastUpdate' );    end
  try, rmappdata( hF , 'bfc_passCost' );      end
  try, rmappdata( hF , 'bfc_throttle' );      end

  try, set( hF , 'Faces' , F ); end
  try, set( hF , 'DeleteFcn' , '' ); end
end

function trailingUpdate( hF , hAx )
%throttle-timer callback: run the pass that the throttle skipped, so the
%RESTING partition is exact. Re-enters updateCamera through its normal guards.
  if ~ishandle( hF ) || ~ishandle( hAx ), return; end
  updateCamera( hF , hAx );
end

function backShown( hF )
%hB Visible listener: while hidden the back patch's Faces are NOT updated (the
%backFaceCulling case pays nothing for it); on showing it, refresh if stale.
  if ~ishandle( hF ), return; end
  hB = getappdata( hF , 'bfc_backPatch' );
  if isempty( hB ) || ~ishandle( hB ) || ~strcmp( get( hB ,'Visible') ,'on' ), return; end
  if isequal( getappdata( hF , 'bfc_backStale' ) , true )
    F = getappdata( hF , 'bfc_faces' );
    w = getappdata( hF , 'bfc_partition' );
    if ~isempty( w ), set( hB , 'Faces' , F( w ,:) ); end
    setappdata( hF , 'bfc_backStale' , false );
  end
end

function updateMesh( hF )
%re-cache everything that depends only on the mesh (fired by the Vertices
%listener and once at split time), then force a repartition pass.
  if ~strcmp( get( get( hF , 'Parent' ) ,'Type' ) ,'axes' )
    warning( 'Culling of backFace can be wrong for hgtransform children.');
  end

  hB = getappdata( hF , 'bfc_backPatch' );

  F = get( hF , 'Faces' );
  if size( F ,2) ~= 3, warning( 'Only triangular meshes are allowed for backFaceCullingSplit.'); end

  V = get( hF , 'Vertices' );
  set( hB , 'Vertices' , V );

  % classification caches: SINGLE precision (only the SIGN of the angles
  % matters, and single halves the GEMV memory traffic). For the perspective
  % path,  dot(N, C-P) = dot(N, C-g) - dot(N, P-g)  with g = mesh centroid:
  % the per-face term is cached (bfc_NdotCrel) so each event needs ONE GEMV
  % and a subtraction -- no nF x 3 temp; centering at g keeps the single-
  % precision cancellation mesh-local (safe for meshes far from the origin).
  N = meshNormals( V , F );
  g = mean( V , 1 );
  C = meshFacesCenter( V , F );
  setappdata( hF , 'bfc_faces'    , F );
  setappdata( hF , 'bfc_center'   , g );
  setappdata( hF , 'bfc_normalsS' , single( N ) );
  setappdata( hF , 'bfc_NdotCrel' , single( sum( N .* bsxfun( @minus , C , g ) ,2) ) );

  setappdata( hF , 'bfc_partition' , [] );            %mesh changed: force the pass
  setappdata( hF , 'bfc_lastState' , [] );
  updateCamera( hF , ancestor( hF ,'axes') );
end

function updateCamera( hF , hAx )
%repartition for the CURRENT camera:
%   1. dedup guard      skip if the camera state is unchanged (absorbs the
%                       double View+CameraPosition firing per gesture)
%   2. safety throttle  a real change < 33 ms after the last full pass is
%                       skipped; the trailing timer lands it when the burst ends
%   3. classification   angles = N . viewdir (ortho: one direction;
%                       perspective: per-face center-to-camera)
%   4. re-upload        ONLY if the front/back partition actually changed --
%                       and only onto the back patch when it is visible.
  if ~ishandle( hF  ), return; end
  if ~ishandle( hAx ), return; end

  % dedup guard on the EFFECTIVE inputs (same idea as silhouette and
  % silhouetteByShading): in ORTHOGRAPHIC the partition depends only on the
  % view DIRECTION (normalized, quantized 1e-12 to absorb float noise), so
  % zoom / dolly / roll skip even the classification GEMV, and a pan collapses
  % to one throttled pass; in PERSPECTIVE it depends only on CameraPosition
  % (zoom and CameraTarget changes are free). updateMesh clears the state.
  if strcmp( get( hAx ,'Projection') ,'orthographic' )
    dirQ = get( hAx ,'CameraTarget') - get( hAx ,'CameraPosition');
    dirQ = dirQ ./ sqrt( sum( dirQ.^2 ) );
    dirQ = round( dirQ * 1e12 ) * 1e-12;
    state = { 'o' , dirQ };
  else
    state = { 'p' , get( hAx ,'CameraPosition') };
  end
  if isequal( state , getappdata( hF ,'bfc_lastState') ), return; end

  % ADAPTIVE SAFETY THROTTLE: recompute at most every max( 33 ms , 3x the
  % MEASURED pass cost ) during event bursts -- the duty cycle stays ~1/3 at
  % any mesh size; the trailing timer lands the exact resting partition.
  lastT = getappdata( hF , 'bfc_lastUpdate' );
  pc    = getappdata( hF , 'bfc_passCost' );  if isempty( pc ), pc = 0; end
  if ~isempty( lastT ) && toc( lastT ) < max( 0.033 , 3*pc )
    tD = getappdata( hF , 'bfc_throttle' );
    if ~isempty( tD ) && isvalid( tD ), stop( tD ); start( tD ); end
    return;
  end
  tPass = tic;

  hB = getappdata( hF , 'bfc_backPatch' );

  CameraPosition = get( hAx , 'CameraPosition' );
  Ns = getappdata( hF ,'bfc_normalsS');               %single: sign-only use
  switch get( hAx , 'Projection' )
    case 'orthographic'
      viewDIR = get( hAx , 'CameraTarget' ) - CameraPosition;
      angles  = Ns * single( viewDIR(:) );
    case 'perspective'
      % dot(N, C-P) = (cached) dot(N, C-g) - dot(N, P-g): one GEMV + one
      % subtraction, no nF x 3 temporary (see the cache note in updateMesh)
      p = CameraPosition - getappdata( hF ,'bfc_center');
      angles = getappdata( hF ,'bfc_NdotCrel') - Ns * single( p(:) );
  end
  w = angles > 0;                                     %back-facing faces

  if ~isequal( w , getappdata( hF , 'bfc_partition' ) )
    F = getappdata( hF , 'bfc_faces' );
    set( hF , 'Faces' , F(~w,:) );
    if strcmp( get( hB ,'Visible') ,'on' )
      set( hB , 'Faces' , F( w,:) );
      setappdata( hF , 'bfc_backStale' , false );
    else
      setappdata( hF , 'bfc_backStale' , true );      %refresh on show (backShown)
    end
    setappdata( hF , 'bfc_partition' , w );
    %NO drawnow here: at large meshes a forced render INSIDE the camera
    %callback adds its full cost to the interaction latency, and MATLAB
    %repaints on its own between events (the new partition rides that natural
    %render). Scripted camera loops must drawnow themselves.
  end

  setappdata( hF , 'bfc_lastState'  , state );
  setappdata( hF , 'bfc_lastUpdate' , tic );
  setappdata( hF , 'bfc_passCost'   , toc( tPass ) );
end

function N = meshNormals( V , F )
%UNNORMALIZED face normals: only the sign of N.viewdir matters.
  x = V( F(:,1) ,:);
  N = cross( V( F(:,2) ,:) - x , V( F(:,3) ,:) - x , 2 );
end
function C = meshFacesCenter( V , F )
%face centroids: the perspective classification ray origins.
  C = ( V( F(:,1) ,:) +  V( F(:,2) ,:) + V( F(:,3) ,:) )/3;
end
