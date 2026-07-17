function MP = CubicMesh_2_mvpoly( CH )
%{

CH = CubicMeshClass( 'Mesh.exnode' , 'Mesh.exelem' );
MP = CubicMesh_2_mvpoly( CH );
plot( MP , 'randcolor' ,'bound','on','edgecolor','none','facealpha',0.5)



CH = CubicMeshClass( 'biVentricular.exnode' , 'biVentricular.exelem' );
MP = CubicMesh_2_mvpoly( CH );
plot( MP , 'randcolor' ,'bound','on','edgecolor','none','facealpha',0.5)

%}

  sz = [ 64 , 3 , CH.nElems ];
  C = zeros(sz);
  for e = 1:sz(3)
    %disp( CH.liNofElem(e,:) );
    id = CH.giNofElem( e , CH.liNofElem(e,:) );
    C(:,:,e) = cat(1, CH.coordinates(         id ,:,1) ,...
                      CH.derivatives.duds1(   id ,:,1) ,...
                      CH.derivatives.duds2(   id ,:,1) ,...
                      CH.derivatives.duds3(   id ,:,1) ,...
                      CH.derivatives.duds12(  id ,:,1) ,...
                      CH.derivatives.duds13(  id ,:,1) ,...
                      CH.derivatives.duds23(  id ,:,1) ,...
                      CH.derivatives.duds123( id ,:,1) );
  end
  
  o   = [1,2,9,10,3,4,11,12,17,18,33,34,19,20,35,36,5,6,13,14,7,8,15,16,21,22,37,38,23,24,39,40,25,26,41,42,27,28,43,44,49,50,57,58,51,52,59,60,29,30,45,46,31,32,47,48,53,54,61,62,55,56,63,64];
  H2P = [1,0,-0,-0;-0,0,1,-0;-3,3,-2,-1;2,-2,1,1];

  C = C(o,:,:);
  C = kron(H2P,H2P,H2P)*C(:,:);
  C = reshape( C , sz );
  
  MP = mvpoly( C , [3 3 3] , 3 );
end

% % CH = CubicMeshClass( 'Mesh.exnode' , 'Mesh.exelem' );
% CH = CHo;
% 
% iElem = 144;
% 
% M0 = [];for l=1:6, [ M.tri , M.xyz ] = CH.tesselateElementFace( iElem , l , 10 ); M0 = AppendMeshes( M0 , M );end
% clf; plotMESH( M0 );
% 
% CH.coordinates         = CH.coordinates         + randn( size( CH.coordinates         ) )*1;
% CH.derivatives.duds1   = CH.derivatives.duds1   + randn( size( CH.derivatives.duds1   ) )*1;
% CH.derivatives.duds2   = CH.derivatives.duds2   + randn( size( CH.derivatives.duds2   ) )*1;
% CH.derivatives.duds3   = CH.derivatives.duds3   + randn( size( CH.derivatives.duds3   ) )*1;
% CH.derivatives.duds12  = CH.derivatives.duds12  + randn( size( CH.derivatives.duds12  ) )*1;
% CH.derivatives.duds23  = CH.derivatives.duds23  + randn( size( CH.derivatives.duds23  ) )*1;
% CH.derivatives.duds13  = CH.derivatives.duds13  + randn( size( CH.derivatives.duds13  ) )*1;
% CH.derivatives.duds123 = CH.derivatives.duds123 + randn( size( CH.derivatives.duds123 ) )*1;
% 
% M1 = [];for l=1:6, [ M.tri , M.xyz ] = CH.tesselateElementFace( iElem , l , 10 ); M1 = AppendMeshes( M1 , M );end
% hplotMESH( M1 ,'facecolor','r');
% 
% 
% C = ndmat( linspace(0,1,5) , linspace(0,1,5) , linspace(0,1,5) );
% XYZ = zeros(size(C,1),3);
% for c=1:size(C,1)
%   XYZ(c,:) = CH.fastCoordinateEvaluation( iElem , C(c,:) ).';
% end
% MP = fit( mvpoly( zeros(64,3,1) , [3 3 3] ) , XYZ(:,:) , C );
% hplot( MP );
% 
% id = CH.giNofElem( iElem , CH.liNofElem(iElem,:) );
% u = cat(1, CH.coordinates(         id ,:,1) ,...
%            CH.derivatives.duds1(   id ,:,1) ,...
%            CH.derivatives.duds2(   id ,:,1) ,...
%            CH.derivatives.duds3(   id ,:,1) ,...
%            CH.derivatives.duds12(  id ,:,1) ,...
%            CH.derivatives.duds13(  id ,:,1) ,...
%            CH.derivatives.duds23(  id ,:,1) ,...
%            CH.derivatives.duds123( id ,:,1) );
% 
% H2C = [1,0,-0,-0;-0,0,1,-0;-3,3,-2,-1;2,-2,1,1];
% C2H = [1 0 0 0 ; 1 1 1 1 ; 0 1 0 0 ; 0 1 2 3 ];
% 
% kC = kron(C2H,C2H,C2H)*reshape( MP.C , [] , 3  );
% 
% o = []; for i=1:64, o = [ o , val2ind( u(:,1) , kC(i,1) )]; end
% io(o) = 1:64;
% maxnorm( kC(:,:) , u(o,:) )
% maxnorm( kC(io,:) , u(:,:) )
% 
% maxnorm( reshape(MP.C,[],3) , kron(H2C,H2C,H2C)*u(o,:) )
% 
% 
% ku = kron(H2C,H2C,H2C) * u;
% 
% [ sort(ku(:)) , sort( MP.C(:) ) ]
% 
% 
% o = [];for i=1:64, 
% 
% 
% ku(o,:) - reshape(MP.C,[],3)
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% % function MP = CH2mvpoly( CH )
% 
% %{
% 
% if isempty(which('CubicMeshClass')), cwd = pwd; cd('c:\Work\computationalcardiacanatomy\'); getCubicHermiteHomeDir(); cd(cwd); clearvars(cwd); end
% 
% CH = CubicMeshClass( 'Mesh.exnode' , 'Mesh.exelem' );
% 
% %}
% 
% %%
% 
% 
% 
% 
% 
% %%
% 
% 
% 
% 
% MP = mvpoly( zeros(64,3) , 3 );
% MP = MP([]);
% for iElem = 1:CH.nElems
%   C = ndmat( linspace(0,1,4) , linspace(0,1,4) , linspace(0,1,4) );
%   XYZ = zeros(size(C,1),3);
%   for c=1:size(C,1)
%     XYZ(c,:) = CH.fastCoordinateEvaluation( iElem , C(c,:) ).';
%   end
%   MP = [ MP ; fit( mvpoly( zeros(64,3,1) , [3 3 3] ) , XYZ(:,:) , C ) ];
%   
% %   C = rand(1,3); maxnorm( MP.v(C), CH.fastCoordinateEvaluation( iElem , C ).' )
% end
% 
% %%
% B   = [1,0,-3,2;0,0,3,-2;0,1,-2,1;0,0,-1,1];
% iB  = pinv(B);
% rp  = [1,17,5,21,2,18,6,22,33,49,37,53,34,50,38,54,9,25,13,29,10,26,14,30,3,19,7,23,4,20,8,24,41,57,45,61,42,58,46,62,35,51,39,55,36,52,40,56,11,27,15,31,12,28,16,32,43,59,47,63,44,60,48,64];
% MAT = comm(8) * sparse(1:64,rp,1) * kron(B,B,B);
% 
% for iElem = 1:12;
% %%
% clc
% 
%   id = CH.giNofElem( iElem , CH.liNofElem(iElem,:) );
%   u = cat(1, CH.coordinates(         id ,:,1) ,...
%              CH.derivatives.duds1(   id ,:,1) ,...
%              CH.derivatives.duds2(   id ,:,1) ,...
%              CH.derivatives.duds3(   id ,:,1) ,...
%              CH.derivatives.duds12(  id ,:,1) ,...
%              CH.derivatives.duds13(  id ,:,1) ,...
%              CH.derivatives.duds23(  id ,:,1) ,...
%              CH.derivatives.duds123( id ,:,1) );
%   
% B = kron(C2H,C2H,C2H) * reshape(MP.C,[],3);
%   
% 
% o=[];
% for i=1:64
%   o = [ o , val2ind( u(:,1) , B(i,1) ) ];
% end
% io(o) = 1:64;
% 
% uu = kron(H2C,H2C,H2C)*u;
% maxnorm( uu(o,:) , reshape(MP.C,[],3) )
% 
% 
%            
%   extent([ sort(MP.C(:)) .\ sort( vec( kron(  B  ,  B  ,  B  ) * reshape( permute(u,[1 3 2]) ,[],3) )) ])
%   extent([ sort(MP.C(:)) .\ sort( vec( kron(  B  ,  B  ,  B  ) * reshape( permute(u,[3 1 2]) ,[],3) )) ])
%   extent([ sort(MP.C(:)) .\ sort( vec( kron(  B' ,  B' ,  B' ) * reshape( permute(u,[1 3 2]) ,[],3) )) ])
%   extent([ sort(MP.C(:)) .\ sort( vec( kron(  B' ,  B' ,  B' ) * reshape( permute(u,[3 1 2]) ,[],3) )) ])
%   extent([ sort(MP.C(:)) .\ sort( vec( kron( iB  , iB  , iB  ) * reshape( permute(u,[1 3 2]) ,[],3) )) ])
%   extent([ sort(MP.C(:)) .\ sort( vec( kron( iB  , iB  , iB  ) * reshape( permute(u,[3 1 2]) ,[],3) )) ])
%   extent([ sort(MP.C(:)) .\ sort( vec( kron( iB' , iB' , iB' ) * reshape( permute(u,[1 3 2]) ,[],3) )) ])
%   extent([ sort(MP.C(:)) .\ sort( vec( kron( iB' , iB' , iB' ) * reshape( permute(u,[3 1 2]) ,[],3) )) ])
% %%
%             
% end
% 
% MP = mvpoly( permute(U,[2 1 3]) , 3 );
% 
% 
% 
% 
% % end
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% %%
% 
% MP = mvpoly( zeros([0,1]) , [3 3] , 2 );
% MP.C(:) = 1:numel(MP.C);
% 
% C2H = [1 0 0 0 ; 1 1 1 1 ; 0 1 0 0 ; 0 1 2 3 ];
% 
% V = ndmat([0 1],[0 1]);
% [ (1:16).' , ([...
%    evaluate( MP        , 0 , 0 ) ;...   %1
%    evaluate( MP        , 1 , 0 ) ;...   %2
%    evaluate( MP        , 0 , 1 ) ;...   %3
%    evaluate( MP        , 1 , 1 ) ;...   %4
%    evaluate( MP.d(1)   , 0 , 0 ) ;...   %5
%    evaluate( MP.d(1)   , 1 , 0 ) ;...   %6
%    evaluate( MP.d(1)   , 0 , 1 ) ;...   %7
%    evaluate( MP.d(1)   , 1 , 1 ) ;...   %8
%    evaluate( MP.d(2)   , 0 , 0 ) ;...   %9
%    evaluate( MP.d(2)   , 1 , 0 ) ;...   %10
%    evaluate( MP.d(2)   , 0 , 1 ) ;...   %11
%    evaluate( MP.d(2)   , 1 , 1 ) ;...   %12
%    evaluate( MP.d(1,2) , 0 , 0 ) ;...   %13
%    evaluate( MP.d(1,2) , 1 , 0 ) ;...   %14
%    evaluate( MP.d(1,2) , 0 , 1 ) ;...   %15
%    evaluate( MP.d(1,2) , 1 , 1 ) ;...   %16
%   ]) , ([...
%    evaluate( MP        , V ) ;...   %1
%    evaluate( MP.d(1)   , V ) ;...   %5
%    evaluate( MP.d(2)   , V ) ;...   %12
%    evaluate( MP.d(1,2) , V ) ;...   %13
%   ]) , ...
% getv( kron(C2H,C2H) * MP.C(:) , vec( permute( reshape( 1:16 , [2 2 2 2] ) , [1 3 2 4] ) )' ) ]
% 
% 
% 
