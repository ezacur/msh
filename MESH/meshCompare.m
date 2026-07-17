function [areeq,a,b] = meshCompare( A , B , th )

  if nargin < 3, th = 'plot'; end

  if ischar( th ) && strcmp( th , 'plot' )
    try, A = toMesh( A ); end
    try, B = toMesh( B ); end
%     A = MeshTidy( Mesh(A,0) ,0,true);
%     B = MeshTidy( Mesh(B,0) ,0,true);
    [~,a,b] = meshCompare( A , B , 1e-12 );

    A.xyz = A.xyz + ( rand(size(A.xyz))*2 - 1 )*1e-3;
    B.xyz = B.xyz + ( rand(size(B.xyz))*2 - 1 )*1e-3;

    figure;
     plotMESH( A ,'FaceColor','r','FaceAlpha',0.05,'EdgeColor',[0.50,0,0],'EdgeAlpha',0.05,'DisplayName','A','gouraud','dull');
    hplotMESH( B ,'FaceColor','b','FaceAlpha',0.05,'EdgeColor',[0,0,0.50],'EdgeAlpha',0.05,'DisplayName','B','gouraud','dull');
    hplotMESH( a ,'FaceColor','r','FaceAlpha',0.80,'EdgeColor',[0.25,0,0],'EdgeAlpha',1.00,'DisplayName','A-B','LineWidth',2,'gouraud','dull');
    hplotMESH( b ,'FaceColor','b','FaceAlpha',0.80,'EdgeColor',[0,0,0.25],'EdgeAlpha',1.00,'DisplayName','B-A','LineWidth',2,'gouraud','dull');
    headlight();
    ze;
    legend_();
    return;

  end


  if isequal( th , 0 )

    A = MeshTidy( Mesh(A,0) ,0,true);
    B = MeshTidy( Mesh(B,0) ,0,true);
    areeq = isequal( A.xyz , B.xyz ) && isequal( A.tri , B.tri );
    if areeq
      areeq = 0;
      return;
    end


    FACTOR = 1/2;

    th = 1;
    while 1
      areeq = meshCompare( A , B , th );
      if th < 1e-12
        areeq = th;
        return
      end
      if ~areeq
        areeq = th / FACTOR;
        return;
      end
      th = th * FACTOR;
    end
  end


  areeq = false;
  A = MeshTidy( Mesh(A,0) ,0,true);
  B = MeshTidy( Mesh(B,0) ,0,true);

  a = MeshRemoveFaces( A , stickToPoints( A , B , th ) );
  if ~isempty( a.tri ) && nargout < 3, return; end

  b = MeshRemoveFaces( B , stickToPoints( B , A , th ) );
  if ~isempty( b.tri ), return; end
  
  areeq = isempty( a.tri ) && isempty( b.tri );

end