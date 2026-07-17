function hTr = interactivePoseMesh( hM , undoFcn )
if 0
  M = loadv( 'C:\Dropbox\mTools\Corify_tools\RESOURCES\MARKERS3D.mat' , 'MARKERS3D' );
  close all
  hM = plotMESH( M ,'ne','nice','sil'); ze(3)
  hTr = interactivePoseMesh( hM );

  %%
end

  if nargin < 1, undoFcn = []; end
  if ~strcmp( get(hM,'Type') , 'patch' ), error('A patch handle was expected.'); end

  if ~strcmp( get(get(hM,'Parent'),'Type') , 'hgtransform' )
    set( hM , 'Parent' , hgtransform( 'Parent' , ancestor( hM ,'axes' ) ,'Hittest','off') );
  end
  hTr = get( hM ,'Parent' );

  M = Mesh( hM ,0);
  M = MeshTidy(M,0,true);
  M = meshSeparate( M ,'largest');
  set( hTr ,'UserData' , miniball( M.xyz ) );

  %set( hM , 'PickableParts' , 'all' );
  set( hM , 'ButtonDownFcn' , @(~,~)ButtonDownFcn_on_hM( hM , undoFcn ) );
end


function out = ButtonDownFcn_on_hM( hM , undoFcn )
  out = false;
  pk = pressedkeys(1);
  w = strcmp( pk , 'SPACE'    ); if  any(w), return; end
  w = strcmp( pk , 'LCONTROL' ); if ~any(w), return; end; pk(w) = [];
  hFig = ancestor( hM , 'figure' );
  if strcmp( get(hFig,'SelectionType') , 'open' ), return; end
  %if numel( pk ) ~= 1, return; end
  out = true;
  
  %SCALE action
  if false && numel( pk ) == 1 && strcmp( pk{1} , 'BUTTON3' )
    STOP_figure_actions( hFig );
  
    hTr = get( hM , 'Parent' );
    TR  = get( hTr , 'Matrix' ); try, undoFcn( hTr ); end
    C   = get( hTr , 'UserData' ); C = TR(1:3,1:3) * C(:) + TR(1:3,4);
    CC  = [ eye(3) , C ; 0 , 0 , 0 , 1 ];
    iCC = [ eye(3) , -CC(1:3,4) ; 0 , 0 , 0 , 1 ];
  
    FigurePoint = get( hFig , 'CurrentPoint' );
  
    set( hFig , 'WindowButtonMotionFcn' , @(~,~) START_SCALE() );
    return;
  end
  function START_SCALE( )
    newFigurePoint = get( hFig , 'CurrentPoint' );
    F = exp( -( FigurePoint(2) - newFigurePoint(2) )/150 );

    set( hTr , 'Matrix' , CC * [F,0,0,0;0,F,0,0;0,0,F,0;0,0,0,1] * iCC * TR );
  end
  

  %TRANSLATION action
  if numel( pk ) == 1 && strcmp( pk{1} , 'BUTTON2' )
    STOP_figure_actions( hFig );

    hTr = get( hM , 'Parent' );
    TR  = get( hTr , 'Matrix' ); try, undoFcn( hTr ); end
    hAx = ancestor( hTr , 'axes');
    C   = get( hTr , 'UserData' ); C = TR(1:3,1:3) * C(:) + TR(1:3,4);

    CP = get( hAx , 'CurrentPoint' );  %CurrentPoint
    [PL,iPL] = getPlane( [ C(:).' ; diff( CP , 1 , 1 ) ] );
    CP   = mean( CP , 1 );
    CP2d = diag([1,1,0]) * ( iPL(1:3,1:3) * CP(:) + iPL(1:3,4) );

    set( hFig , 'WindowButtonMotionFcn' , @(~,~) START_TRANSLATE() );
    return;
  end
  function START_TRANSLATE( )
    nCP = get( hAx , 'CurrentPoint' );
    nCP = mean( nCP , 1 );
    nCP2d = [1,0,0;0,1,0;0,0,0] * ( iPL(1:3,1:3) * nCP(:) + iPL(1:3,4) );

    t = PL(1:3,1:3) * ( nCP2d - CP2d );

    set( hTr , 'Matrix' , [1,0,0,t(1);0,1,0,t(2);0,0,1,t(3);0,0,0,1] * TR );
  end


  %ROTATE action
  if numel( pk ) == 1 && strcmp( pk{1} , 'BUTTON1' )
    STOP_figure_actions( hFig );

    hTr = get( hM , 'Parent' );
    TR  = get( hTr , 'Matrix' ); try, undoFcn( hTr ); end
    hAx = ancestor( hTr , 'axes');
    C   = get( hTr , 'UserData' ); C = TR(1:3,1:3) * C(:) + TR(1:3,4);

    CP = get( hAx , 'CurrentPoint' );  %CurrentPoint
    [PL,iPL] = getPlane( [ C(:).' ; diff( CP , 1 , 1 ) ] );
    CP   = mean( CP , 1 );
    CP2d = diag([1,1,0]) * ( iPL(1:3,1:3) * CP(:) + iPL(1:3,4) );
    
    set( hFig , 'WindowButtonMotionFcn' , @(~,~) START_ROTATE() );
    return;
  end
  function START_ROTATE( )
    nCP = get( hAx , 'CurrentPoint' );
    nCP = mean( nCP , 1 );
    nCP2d = diag([1,1,0]) * ( iPL(1:3,1:3) * nCP(:) + iPL(1:3,4) );

    a = - atan2d( CP2d(2) , CP2d(1) ) + atan2d( nCP2d(2) , nCP2d(1) );
    R = [ cosd(a) , -sind(a) , 0 , 0 ; sind(a) , cosd(a)  , 0, 0 ; 0 , 0 , 1  , 0 ; 0 , 0 , 0 , 1 ];
    set( hTr , 'Matrix' , PL * R * iPL * TR );
  end
  
end
function out = STOP_figure_actions( hFig )
  out = true;
  oldFCNs = getappdata( hFig , 'oldFCNs' );
  if isempty( oldFCNs )
    oldFCNs = struct( 'WBD' , get( hFig , 'WindowButtonDownFcn' ) , 'WBM' , get( hFig , 'WindowButtonMotionFcn' ) , 'WBU' , get( hFig , 'WindowButtonUpFcn' ) , 'KP' , get( hFig , 'KeyPressFcn' ) );
    setappdata( hFig , 'oldFCNs' , oldFCNs );
    set( hFig , 'WindowButtonDownFcn'   , '' );
    set( hFig , 'WindowButtonMotionFcn' , '' );
    set( hFig , 'WindowButtonUpFcn'     , @(~,~)STOP_figure_actions( hFig ) );
    set( hFig , 'KeyPressFcn'           , '' );
  else
    rmappdata( hFig , 'oldFCNs' );
    set( hFig , 'WindowButtonDownFcn'   , oldFCNs.WBD );
    set( hFig , 'WindowButtonMotionFcn' , oldFCNs.WBM );
    set( hFig , 'WindowButtonUpFcn'     , oldFCNs.WBU );
    set( hFig , 'KeyPressFcn'           , oldFCNs.KP  );

    h2d = findall( hFig , 'Tag','_auxilar_graphics_' );
    if ~isempty( h2d ), delete( h2d ); end
  end
end
