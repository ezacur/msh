function paintMesh( hM )
  M = Mesh( Mesh( hM ) ,0);
  M.triID = ( 1:size(M.tri,1) ).';
  
  if isempty( get( hM , 'FaceVertexCData' ) )
    set( hM ,'FaceVertexCData', zeros( size( M.tri ,1) ,1) );
  end
  set( hM , 'FaceColor','flat');
  
  hAxe = ancestor( hM , 'axes' );
  set( hAxe , 'CLim' , [0 21] );
  
  hFig = ancestor( hAxe , 'figure' );
  set( hFig , 'Colormap' , [ 1 , 1 , 1 ; color( 1:20 ) ] );

  set( hFig , 'CurrentObject' , hM );

  set( hFig , 'KeyPressFcn'          , @(h,e)KeyPressFcn(         hFig , e , hAxe , hM ) );
  set( hM   , 'ButtonDownFcn'        , @(h,e)ButtonDownFcn_on_hM( hFig , e , hAxe , hM , M ) );

end

function ButtonDownFcn_on_hM( hFig , e , hAxe , hM , M )
  pk = pressedkeys_win( 1 );
  if any(strcmp(pk,'SPACE')), return; end
  pk = OPZ.pressedK( pk , hFig );

  if numel( pk ) == 1 && strcmp( pk{1} , 'BUTTON10' )
    xyz = []; try, xyz = IntersectSurfaceRay( M , get(hAxe,'CurrentPoint') , 'first' ); end; if isempty( xyz ), return; end
  
    LS = getappdata( hM , 'LANDMARKS' );
    if isempty( LS )
      LS = line( 'Parent' , hAxe ,...
        'XData' , [] , 'YData' , [] , 'ZData' , [] ,...
        'Marker','o','LineStyle','none','LineWidth',2,'MarkerSize',10,'MarkerFaceColor',[1,.5,0] );
      setappdata( hM , 'LANDMARKS' , LS );
    end
  
    XYZ = [ vec( get(LS,'XData') ) , vec( get(LS,'YData') ) , vec( get(LS,'ZData') ) ];
    XYZ = [ XYZ ; xyz ];
    set( LS , 'XData' , XYZ(:,1) , 'YData' , XYZ(:,2) , 'ZData' , XYZ(:,3) );
    return;
  end

  if numel( pk ) == 2 && strcmp( pk{1} , 'BUTTON1' ) && numel( pk{2} ) == 1 && any( pk{2} == '123456789' )
    OPZ.STOP( hFig ); OPZ.setPointer( hFig , 'point' );
    C = get( hM , 'FaceVertexCData' ); LP = NaN(1,2);
    LABEL = str2double( pk{2} ); DO_PAINT( LABEL );
    set( hFig , 'WindowButtonMotionFcn' , @(h,e)DO_PAINT( LABEL ) );
    return;
  end
  if ( numel( pk ) == 1 && strcmp( pk{1} , 'BUTTON3' ) ) ||...
     ( numel( pk ) == 2 && strcmp( pk{1} , 'BUTTON3' ) && numel( pk{2} ) == 1 && any( pk{2} == '123456789' ) )
    OPZ.STOP( hFig ); OPZ.setPointer( hFig , 'point' );
    C = get( hM , 'FaceVertexCData' ); LP = NaN(1,2);
    LABEL = 0; DO_PAINT( LABEL );
    set( hFig , 'WindowButtonMotionFcn' , @(h,e)DO_PAINT( LABEL ) );
    return;
  end
  function DO_PAINT( LABEL )
    lp = get( hFig ,'CurrentPoint'); %disp( max( abs( LP - lp ) ) )
    if max( abs( LP - lp ) ) < 5, return; end
    LP = lp;

    try
      [~,~,cid] = IntersectSurfaceRay( M , get(hAxe,'CurrentPoint') , 'first' ); if isempty(cid), return; end
      cLABEL = C(cid);
      if cLABEL ~= LABEL
        C( cid ) = LABEL;
        set( hM ,'FaceVertexCData', C ); drawnow('update');
      end
    end
  end

  
  if numel( pk ) == 3 && strcmp( pk{1} , 'BUTTON1' ) && numel( pk{2} ) == 1 && any( pk{2} == '123456789' ) && strcmp( pk{3} ,'LSHIFT')
    try, DO_FILL( str2double( pk{2} ) ); end
    return;
  end
  if ( numel( pk ) == 2 && strcmp( pk{1} , 'BUTTON3' ) && strcmp( pk{2} ,'LSHIFT') ) || ...
     ( numel( pk ) == 3 && strcmp( pk{1} , 'BUTTON3' ) && strcmp( pk{3} ,'LSHIFT') && numel( pk{2} ) == 1 && any( pk{2} == '123456789' ) )
    try, DO_FILL( 0 ); end
    return;
  end
  function DO_FILL( LABEL )
    [~,~,cid] = IntersectSurfaceRay( M , get(hAxe,'CurrentPoint') , 'first' ); if isempty(cid), return; end
    C = get( hM , 'FaceVertexCData' );
    LS = meshSeparate( M ,C);
    LS = fun(  @meshSeparate ,LS ,'un',0);
    LS = cat(1,LS{:});
    LS = LS{ fun( @(m)any( m.triID == cid ) ,LS) };

    C( LS.triID ) = LABEL;
    set( hM ,'FaceVertexCData', C );
  end

end
function KeyPressFcn( hFig , ev , hAx , hM )
  pk = OPZ.pressedK( 0 );
  if any( strcmp( pk , 'SPACE' ) ), return; end
  K = upper(ev.Key);
  K = [ pk{ ~strcmp(pk,K) } , K ];
  
  switch K
    case 'LSHIFTE',
      ev = get( hM ,'EdgeColor');
      if strcmp( ev ,'none'), set( hM ,'EdgeColor',[0,0,0]);
      else,                   set( hM ,'EdgeColor','none');
      end
  end
end
