function clickOnMesh( hM , vn )

  if nargin == 1, vn = 'onMesh'; end

  M = Mesh( hM );
  
  hAxe = ancestor( hM , 'axes' );
  hFig = ancestor( hAxe , 'figure' );
  
  hC = line(NaN,NaN,NaN,'Parent',hAxe,'Color','r','Marker','o','MarkerFaceColor','r','MarkerSize',5,'Hittest','off');
  hP = line(NaN,NaN,NaN,'Parent',hAxe,'Color','k','Marker','o','MarkerFaceColor','y','MarkerSize',12,'LineStyle','none','Hittest','off');
  hT = [];


  set( hM   , 'ButtonDownFcn'         , @(h,e)addPoint );
  set( hFig , 'WindowButtonMotionFcn' , @(h,e)setCursor );
  
  function setCursor
    ray = get(hAxe,'currentpoint');
    try
      xyz = IntersectSurfaceRay( M , ray , 'first' );
      set( hC ,'XData',xyz(1,1),'YData',xyz(1,2),'ZData',xyz(1,3));
    catch
      set( hC ,'XData',[],'YData',[],'ZData',[])
    end
  end
  
  function addPoint
    if ~pressedkeys(3) == 1 || ~strcmp( get( hFig , 'SelectionType' ) , 'open' )
      return;
    end
    
    ray = get( hAxe , 'currentpoint' );
    xyz = IntersectSurfaceRay( M , ray , 'first' );
    
    x = get( hP , 'XData' );
    y = get( hP , 'YData' );
    z = get( hP , 'ZData' );
    if numel(x) == 1 && isnan( x )
      x = []; y = []; z = [];
    end

    x = [ x(:) ; xyz(1,1) ];
    y = [ y(:) ; xyz(1,2) ];
    z = [ z(:) ; xyz(1,3) ];
    
    set( hP , 'XData' , x , 'YData' , y , 'ZData' , z );
    
    if ~isempty( vn )
      assignin( 'base' , vn , [ x , y , z ] );
    end
    
    hT = [ hT ; text(x(end),y(end),z(end), sprintf('%d',numel(hT)+1) ,'Parent',hAxe,'VerticalAlignment','middle','HorizontalAlignment','center','FontWeight','bold','FontSize',8) ];
    
  end
  
end
