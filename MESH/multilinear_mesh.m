
addCell = @(M,xyz)struct('xyz',[M.xyz;bsxfun(@plus,ndmat([0 1],[0 1],0),xyz)],'tri',[M.tri;[1 2 4 3]+size(M.xyz,1)]);

M = struct('xyz',[],'tri',[]);

for x = 0:3, for y = 0:4,
  M = addCell(M,[x y 0]);
end; end

M.xyz = M.xyz + rand( size(M.xyz) )/10;
M.xyz(:,3) = 0;

M = TidyMesh( M , 0.2 );

%%

cla
% subplot(1,2,1); 
patch('vertices',M.xyz,'faces',M.tri,'facecolor','none','edgecolor',[1 0 0],'marker','o','markersize',20,'markerfacecolor',[1 0.5 0.5],'linewidth',3)
text( M.xyz(:,1) , M.xyz(:,2) , M.xyz(:,3) , arrayfun(@(n)sprintf('%d',n),1:size(M.xyz),'un',0) , 'horizontalAlignment','center','verticalAlignment','middle')
C = [ mean( reshape( M.xyz( M.tri , 1 ) , size(M.tri) ) , 2 ) , mean( reshape( M.xyz( M.tri , 2 ) , size(M.tri) ) , 2 ) , mean( reshape( M.xyz( M.tri , 3 ) , size(M.tri) ) , 2 ) ];
line( C(:,1) , C(:,2) , C(:,3) , 'linestyle','none','marker','s','markersize',20,'markerfacecolor',[1 1 0]);
text( C(:,1) , C(:,2) , C(:,3) , arrayfun(@(n)sprintf('%d',n),1:size(M.tri),'un',0) , 'horizontalAlignment','center','verticalAlignment','middle','fontweight','bold')
axis equal

%%

MP = mvpoly( zeros(1,3) , 0 , 2 );
for i = 1:size(M.tri,1)
  MP = [ MP ; mvpoly( [1 0 0 0;-1 1 0 0;-1 0 1 0;1 -1 -1 1] * M.xyz( M.tri(i,[1;2;4;3]) , : ) , [1 1] ) ];
end
MP = MP([2:end]);

hplot( affCoordinates( MP.c([1 2]) ,[0 1] , [0 1] ) , 'rc' ,'res',[5 10],'facealpha',0.2)

L = [];
for i = 1:size(M.tri,1)
  L = [ L ; kron(  eye(3) , [1 0 0 0;-1 1 0 0;-1 0 1 0;1 -1 -1 1] ) * kron( eye(3) , sparse( 1:4 , M.tri(i,[1;2;4;3]) , 1 , 4 , size(M.xyz,1) ) ) ];
end
maxnorm( L * M.xyz(:) , MP.C(:) )


%%

M.xyz = M.xyz(:,1:2);
[MP,L] = mesh2mvpoly( M , 'bilinear' );
maxnorm( L * M.xyz(:) , MP.C(:) )

%%


addCell = @(M,xyz)struct('xyz',[M.xyz;bsxfun(@plus,ndmat([0 1],[0 1],0),xyz)],'tri',[M.tri;[1 2 4 3]+size(M.xyz,1)]);

M = struct('xyz',[],'tri',[]);
for x = 0:5, for y = 0:5,
  M = addCell(M,[x y 0]);
end; end
%plotMESH( M )

M.r = M.xyz(:,2);
M.t = M.xyz(:,1)/max(M.xyz(:,1))*2*pi;

M.xyz(:,1) = M.r .* cos( M.t );
M.xyz(:,2) = M.r .* sin( M.t );
M.xyz(:,3) = max(M.r)-sqrt(max(M.r)^2-M.r.^2);


M = TidyMesh( M , 1e-5 , true );
plotMESH( M ); view(2)

[MP,L] = mesh2mvpoly( M , 'bilinear' );
maxnorm( L * M.xyz(:) , MP.C(:) )

C = MP.C(:) + randn(size(MP.C(:)))/1.020;
MP.C(:) = L*pinv(L)*C;
M.xyz(:) = pinv(L)*C;

plot( MP , 'randcolor' , 'edgecolor','none' );
hplotMESH( M ,'facecolor','none','edgecolor','r','linewidth',3);

