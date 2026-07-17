function hP = exploreMESH( M , varargin )

  if iscell( M )
    M = MeshAppend( M );
  end

  if isstruct( M ) && ~isfield( M ,'xyz' )
    A = [];
    for f = fieldnames(M).', f = f{1};
      A = MeshAppend( A , MeshAddField( M.(f) ,'triPIECE' , f ) );
    end
    exploreMESH( A , varargin{:} );
    return;
  end

  if size( M.tri ,2) == 4 && ~any( M.tri(:,4) ),        M.tri = M.tri(:,1:3); end
  M.meshCelltype = meshCelltype( M );

  try, M.triLABELS; f = 'triLABELS'; end
  try, M.triLABEL;  f = 'triLABEL';  end
  try, M.triPIECE;  f = 'triPIECE';  end
  F = M.(f);

  W = M;
  if numel( varargin ) && isstruct( varargin{1} )
    W = varargin{1}; varargin(1) = [];
  end

  hFig = figure('Color',[1,1,1]);
  hW = plotMESH( W , 'ne','gouraud','FaceAlpha',0.125,'FaceColor','w','HandleVisibility','off','nice','Hittest','off','PickableParts','none','Clip','off','Tag','skeleton');
  set( hW ,'FaceAlpha',0.05,'AmbientStrength',0,'DiffuseStrength',1,'SpecularStrength',1,'SpecularExponent',10^-1);
  set( gca ,'Visible','off');
  axis(gca ,'vis3d');

  hB = patch('Vertices',M.xyz,'Faces',zeros(0,2),'LineWidth',3,'EdgeColor','k','Hittest','off','PickableParts','none','Clip','off','Tag','boundaries');
  if M.meshCelltype == 3
    set( hB ,'Marker','o','LineWidth',1,'MarkerFaceColor','k');
  end

  P = unique( F ,'stable');

  cols = jet( numel(P) )/2 + 0.5;
  rng_state = rng; rng( numel(P) ); rp = randperm( size( cols ,1) ); rng( rng_state );
  cols = cols( rp ,:);

  hP = gobjects( numel(P) ,1);
  for p = 1:numel(P)
    if isempty( P{p} ), P{p} = 'no-PIECEname'; end
    if M.meshCelltype == 5
      hP(p) = hplotMESH( MeshRemoveFaces( M , ~strcmp( F , P{p} ) , true ) ,'FaceColor',cols(p,:),'DisplayName',P{p},'ne','nice',varargin{:},'Clip','off','Tag',P{p});
    elseif M.meshCelltype == 3
      hP(p) = hplotMESH( MeshRemoveFaces( M , ~strcmp( F , P{p} ) , true ) ,'EdgeColor',cols(p,:),'DisplayName',P{p},varargin{:},'Clip','off','Tag',P{p});
    end
  end
  updateExplorer( );

  debouncerVisibility_timer = timer('ExecutionMode','singleShot','StartDelay',0.05,'TimerFcn', @(~,~) updateExplorer( ) );
  setappdata( hFig , 'debouncerVisibility_timer_cleanup' , onCleanup( @()deleteTimer(debouncerVisibility_timer) ) );

  hLeg = legend_();
  hLeg.UserData = hLeg.String(:);
  hLeg.Interpreter = 'tex';
  hLeg.String = strrep( hLeg.UserData , '_' , '\_' );
  hLeg.FontSize = 10;
  hLeg.FontName = 'Segoe UI';
  hLeg.NumColumns = ceil( numel(P) / 43 );
  
  Ls = cell( size(hP) );
  for p = 1:numel(P)
    set( hP(p) ,'ButtonDownFcn' , @(h,e)clickOnPatch( h ) );
    Ls{p} = newListener( hP(p) , 'Visible' , @(~,~)debounceUpdate( ) );
  end
  setappdata( hFig ,'explorer_listeners' , Ls );
  
  
  hFilter = uicontrol('Style','edit','Position',[40 40 400 60],'FontSize',25,'String','','HorizontalAlignment','left');
  highlightItems( get( hFilter , 'String' ) );
  %set( hFilter ,'Callback', @(h,e)highlightItems( hLeg , get( h ,'String' ) ) );

  oldWarn = warning( 'off' , 'MATLAB:ui:javaframe:PropertyToBeRemoved' ); CLEAN = onCleanup( @()warning( oldWarn ) );
  jFilter = findjobj( hFilter );
  set( jFilter , 'KeyPressedCallback' , @(j,e)highlightItems( jFilter.Text ) );

  uicontrol( 'Style' , 'pushbutton' , 'Position' , [ 450 50 40 25 ] ,'String','Show'  , 'Callback' , @(h,e)fShow()   );
  uicontrol( 'Style' , 'pushbutton' , 'Position' , [ 500 50 40 25 ] ,'String','Hide'  , 'Callback' , @(h,e)fHide()   );
  uicontrol( 'Style' , 'pushbutton' , 'Position' , [ 550 50 40 25 ] ,'String','TOP'   , 'Callback' , @(h,e)fTop()    );
  uicontrol( 'Style' , 'pushbutton' , 'Position' , [ 600 50 40 25 ] ,'String','BOTTOM', 'Callback' , @(h,e)fBottom() );

  function fShow()
    set( hP( ~cellfun( 'isempty' , regexp( hLeg.String , '\\color{red}' , 'once' ) ) ) , 'Visible','on');
  end
  function fHide()
    set( hP( ~cellfun( 'isempty' , regexp( hLeg.String , '\\color{red}' , 'once' ) ) ) , 'Visible','off');
  end
  function fTop()
    reorderLegend( ~cellfun( 'isempty' , regexp( hLeg.String , '\\color{red}' , 'once' ) ) ,  1 );
  end
  function fBottom()
    reorderLegend( ~cellfun( 'isempty' , regexp( hLeg.String , '\\color{red}' , 'once' ) ) , -1 );
  end

%   set( hFig , 'WindowButtonMotionFcn' ,@(h,e)MovingMouse( e ) );


  hAlpha = uicontrol( 'Style','slider','Position',[10,10,150,16],'Value',get( hW ,'FaceAlpha'),'Callback',@(h,e)set(hW,'FaceAlpha',h.Value) );
  addlistener( hAlpha ,'ContinuousValueChange', @(varargin)feval(get(hAlpha,'Callback'),hAlpha,[]) );


  hLeg.Units = 'pixels';
  hLeg_pos = hLeg.Position;
  hLeg.Position = [10,130,hLeg_pos(3:4)];
  SetPosition( gca , [ hLeg_pos(3)+10 , 0 , -hLeg_pos(3)-1 , -1 ] ,true);
  


  
  hN = uicontrol( 'Style','popupmenu' ,'Position',[ 450 , 80 , 120 , 23 ] );
  try, hN.String = [ {''} ; getH4() ]; end
  hN.Callback = @(h,e)getH4_seleccion(hN,e);
  function getH4_seleccion(hN,e)
    name = hN.String{ hN.Value };
    if ~isempty( name )
      name = getH4( string( name ) );
      name = strjoin( name , '|' );
    end
    hFilter.String = name;
    highlightItems( hFilter.String );
  end








  function highlightItems( idx )

    if ischar( idx ) && isempty( idx )
      idx = [];
    end

    if ischar( idx )
      idx = ~cellfun( 'isempty' , regexp( hLeg.UserData , [ '^' , idx ] ,'once') );
    end

    if islogical( idx )
      idx = find(idx);
    end

    S = strrep( hLeg.UserData ,'_','\_');
    for i = idx(:).'
      S{i} = [ '\color{red}' , S{i} ];
    end
    hLeg.String = S;

  end

  function clickOnPatch( h )
    pk = pressedkeys(); if ~isempty( pk ), return; end
    pb = pressedkeys(3); if strcmp( get( hFig ,'SelectionType' ) ,'open'), pb = pb*10; end
    if pb == 10
      h.Visible = false;
    elseif pb == 2
      reorderLegend( h , 1);
    elseif pb == 4
      reorderLegend( h ,-1);
    elseif pb == 20
      h.Visible = false;
      reorderLegend( h , 1);
    elseif pb == 40
      h.Visible = false;
      reorderLegend( h ,-1);
    end
  end
  
  function reorderLegend( w , pos )
    if nargin < 2, pos = 1; end
    if ishghandle( w )
      w = hP == w;
    elseif ischar( w )
      w = ~cellfun( 'isempty' , regexp( hLeg.UserData , [ '^' , w ] ,'once') );
    elseif iscell( w )
      w = ismember( hLeg.UserData , w );
    end
    if pos < 0, w = ~w; end
    if numel( w ) ~= numel( hP )
      keyboard;
    end

    hP            = [            hP( w ,:) ;            hP( ~w ,:) ];
    hLeg.UserData = [ hLeg.UserData( w ,:) ; hLeg.UserData( ~w ,:) ];
    legend( hP );
  end

  function MovingMouse( e )
    ho = e.HitObject;
    hLeg.String = regexprep( hLeg.String , '^(?:\\bf)*' , '' );
    if strcmp( ho.Type ,'patch')
      ho.DisplayName = [ '\bf' , ho.DisplayName ];
    end
  end

  function debounceUpdate( )
    if strcmp( debouncerVisibility_timer.Running , 'on' ), stop( debouncerVisibility_timer ); end
    start( debouncerVisibility_timer );
  end

  function updateExplorer( )
    ONs = findall( hP ,'Visible','on');
    if ~isempty( ONs )
      ONs = get( ONs ,'DisplayName');
      ONs = regexprep( ONs , '.*\color\{.*\}' , '' );  % elimina prefijo TeX
      ONs = strrep( ONs , '\_' , '_' );         % des-escapea underscores
    end
  
    if isfield( W ,'triPIECE')
      if isempty( ONs )
         w = W;
      else
        w = MeshRemoveFaces( W , ismember( W.triPIECE , ONs ) );
      end
      set( hW ,'Faces' , w.tri );
    end
  
    if ~isempty( ONs )
      B = MeshBoundary( MeshRemoveFaces( M , ~ismember( F , ONs ) ) );
    else
      B.tri = zeros(0,2);
    end
    set( hB ,'Faces' , B.tri );
  end

end



function nL = newListener( hh , prop , fcn )
  persistent matlabV
  if isempty( matlabV )
    matlabV = sscanf(version,'%d.%d.%d.%d.%d',5); matlabV=[100,1,1e-2,1e-9,1e-13]*[ matlabV(1:min(5,end)) ; zeros(5-numel(matlabV),1) ];
  end

  if matlabV <= 804, nL = handle.listener(    handle( hh ) , findprop( handle(hh) , prop ) , 'PropertyPostSet' , fcn );    
  else,              nL = event.proplistener(         hh   , findprop(        hh  , prop ) , 'PostSet'         , fcn );
  end
end

function deleteTimer( t )
  if isvalid(t), stop(t); delete(t); end
end
