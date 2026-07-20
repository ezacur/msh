function hL_ = headlight( varargin )
%HEADLIGHT  A camlight that keeps following the camera.
%
%   headlight                creates a light AT the camera ('headlight') that
%                            tracks the camera as the view changes.
%   headlight('left')        light left-and-up of the camera, tracking it.
%   headlight('right')       light right-and-up of the camera, tracking it.
%   headlight('headlight')   same as no argument (light at the camera).
%   headlight(az,el)         light at azimuth AZ, elevation EL (degrees) from the
%                            camera, tracking it.
%   headlight(...,'infinite')  directional light (parallel rays).
%   headlight(...,'local')     point light (the camlight default).
%   headlight(hL,...)        (re)position an EXISTING light hL and make it track.
%   headlight(ax,...)        create the tracking light in axes AX.
%   hL = headlight(...)      returns the light handle.
%
%   The position / style arguments are passed straight to CAMLIGHT (headlight's
%   default is 'headlight', i.e. AZ = EL = 0 -- NOT camlight's own default of
%   'right'; a spec with ONLY a style, e.g. headlight('local'), also keeps the
%   'headlight' position).  headlight only adds:  (1) listeners on the axes
%   camera so the light is re-placed with camlight(hL,...) whenever View or
%   CameraPosition changes, and (2) a CreateFcn so the behaviour survives a
%   FILE round-trip (savefig / hgsave -> openfig / hgload fire it).  COPYOBJ
%   does NOT fire CreateFcn nor copy appdata (verified R2022a): a copied light
%   is a plain static camlight -- call headlight( copiedLight , spec... ) to
%   re-attach.
%
%   AZ,EL are RELATIVE TO THE CURRENT VIEW (never a global position): camlight
%   rotates the camera position about the CameraTarget by AZ degrees
%   horizontally (+AZ = to the RIGHT of the camera) and EL vertically
%   (+EL = above), in visually-corrected (DataAspectRatio-aware) space. The
%   named presets are plain aliases:  'headlight' == (0,0) ,
%   'right' == (30,30) , 'left' == (-30,30).  On every camera change the SAME
%   offset is re-applied relative to the NEW camera -- the light behaves like
%   a lamp mounted on the camera rig.
%
%   BOTH styles track the camera and share the same angular placement:
%     'local' (default)  point source located ON the camera's orbit sphere
%                        (same distance to the target as the camera; for
%                        'headlight' that is exactly CameraPosition): rays
%                        diverge from it, so nearby geometry is lit unevenly.
%     'infinite'         directional source: Position becomes a UNIT vector
%                        (the same direction the 'local' point would have,
%                        seen from the target) and the rays are PARALLEL --
%                        illumination depends only on the surface normals,
%                        which gives the more uniform, view-stable look for a
%                        tracking key light.
%
%   See also CAMLIGHT, LIGHT, LIGHTANGLE.

  % ---- separate an optional leading light/axes handle from the camlight spec ---
  % everything after it is passed VERBATIM to camlight ( 'headlight'|'left'|'right'
  % | az,el , plus an optional 'local'|'infinite' ).
  hL = [];  ax = [];  spec = varargin;
  if ~isempty( spec ) && isscalar( spec{1} ) && ishghandle( spec{1} )
    switch get( spec{1} , 'Type' )
      case 'light', hL = spec{1};  spec(1) = [];
      case 'axes',  ax = spec{1};  spec(1) = [];
    end
  end

  % re-attach on an existing light with no new spec (e.g. from the CreateFcn):
  % recover the offset stored on the light so a copy keeps its 'left'/'right'/...
  if isempty( spec ) && ~isempty( hL ) && isappdata( hL , 'headlight_spec' )
    spec = getappdata( hL , 'headlight_spec' );
  end
  if isempty( spec ), spec = { 'headlight' }; end        % default: AT the camera
  % a spec with ONLY style token(s) (e.g. headlight('local')) must keep
  % headlight's documented default POSITION: bare camlight(hL,'local') would
  % silently fall back to camlight's own default 'right' (probe-verified)
  if all( cellfun( @(s)ischar(s) && any(strcmpi(s,{'local','infinite'})) , spec ) )
    spec = [ { 'headlight' } , spec ];
  end

  % ---- create the light, or (re)position an existing one, per the spec ---------
  if     ~isempty( hL ), camlight( hL , spec{:} );
  elseif ~isempty( ax ), hL = camlight( ax , spec{:} );
  else,                  hL = camlight( spec{:} );
  end
  if nargout > 0, hL_ = hL; end
  setappdata( hL , 'headlight_spec' , spec );

  % ---- HG1 / HG2 property-listener factory -------------------------------------
  matlabV = sscanf(version,'%d.%d.%d.%d.%d',5); matlabV=[100,1,1e-2,1e-9,1e-13]*[ matlabV(1:min(5,end)) ; zeros(5-numel(matlabV),1) ];
  if matlabV <= 804, newListener = @(hh,prop,fcn)handle.listener(    handle( hh ) , findprop( handle(hh) , prop ) , 'PropertyPostSet' , fcn );
  else,              newListener = @(hh,prop,fcn)event.proplistener(         hh   , findprop(        hh  , prop ) , 'PostSet'         , fcn );
  end

  hAx = ancestor( hL , 'axes' );

  % re-place the light relative to the CURRENT camera on any camera change. camlight
  % only READS the camera and SETS the light, so this never re-triggers the listeners.
  follow = @() camlight( hL , spec{:} );
  setappdata( hL , 'headlight_listeners' , { newListener( hAx , 'View'           , @(~,~)follow() ) ; ...
                                             newListener( hAx , 'CameraPosition' , @(~,~)follow() ) } );

  % survive copyobj / hgload / print: re-attach WITH the same spec (captured here)
  set( hL , 'CreateFcn' , @(h,~)~isempty(which('headlight'))&&~~double(headlight(h,spec{:})) );

end
