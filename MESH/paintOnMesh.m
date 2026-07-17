function paintOnMesh( hM , hC , fcn )

  if nargin < 3, fcn = []; end
  if nargin < 2, hC  = []; end
  

  M = Mesh( hM );
  
  hAxe = ancestor( hM , 'axes' );
  hFig = ancestor( hAxe , 'figure' );

  if isempty( hC )
    hC = line( NaN , NaN , NaN ,'Parent',hAxe,'Color','r','Marker','o','MarkerSize',4,'Hittest','off');
    %set( ancestor( hAxe ,'figure') , 'CurrentObject' , hC );
  end

  set( hFig , 'KeyPressFcn' , @(h,e)deletePoint(e) );
  set( hM , 'ButtonDownFcn' , @(h,e)clickOnM );
  
  function deletePoint(e)
    if strcmp( e.Key , 'backspace' ) && isempty( e.Modifier )
      
      C = getXYZ( hC );
      if isempty( C ), return; end
      
      C(end,:) = [];
      setXYZ( hC , C );
      if ~isempty( fcn )
        fcn( C );
      end
    end
  end
  
  
  function clickOnM
    pk = pressedkeys(1);
    if numel( pk ) == 2
      if strcmp( pk{2} , 'SPACE' ) 
%         ObjectViewRotate( hM );
      end
      return;
    end

    
    if numel( pk ) ~= 1, return; end
    if ~strcmp( pk{1} , 'BUTTON1' ), return; end
    
    xyz = [];
    
    try, xyz = IntersectSurfaceRay( M , get(hAxe,'currentpoint') , 'first' ); end
    
    if isempty( xyz ), return; end
    oldMoving = [];
    oldMoving = get( hFig , 'WindowButtonMotionFcn' );
    set( hFig , 'WindowButtonUpFcn'     , @(h,e)STOP() );
    
    C = getXYZ( hC );
    addPoint( NaN(1,3) );
    addPoint( xyz );

    set( hFig , 'WindowButtonMotionFcn' , @(h,e)DRAW() );
    
    function DRAW()
      xyz = [];
      try, xyz = IntersectSurfaceRay( M , get(hAxe,'currentpoint') , 'first' ); end
      
      if isempty( xyz )
        STOP();
        return;
      end
      
      addPoint( xyz );
    end
    
    function STOP()
      addPoint( NaN(1,3) );
      set( hFig , 'WindowButtonMotionFcn' , oldMoving );
    end
    
    function addPoint( xyz )
      if any( isnan( xyz ) ) && any( isnan( C(end,:) ) ), return; end
      if isequal( C(end,:) , xyz ), return; end

      C = [ C ; xyz ];

      setXYZ( hC , C );
      if ~isempty( fcn )
        fcn( C );
      end
    end
  end


end
function setXYZ( h , xyz )

  if size( xyz ,2) == 2
    set( h , 'XData' , xyz(:,1) , 'YData' , xyz(:,2) , 'ZData' , zeros(size(xyz,1),1) );
    
  elseif size( xyz ,2) == 3
    
    set( h , 'XData' , xyz(:,1) , 'YData' , xyz(:,2) , 'ZData' , xyz(:,3) );

  end

end
function xyz = getXYZ( h )

  x = get( h , 'XData' );
  y = get( h , 'YData' );
  z = get( h , 'ZData' );
  
  xyz = [ x(:) , y(:) , z(:) ];

end


