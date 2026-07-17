function M = MeshUniform( M , el , Nits )

  M = Mesh( M , 0 );

  M0 = M;
  [ M.xyz , M.tri ] = GEOMETRY_subdivideTri( M.xyz , M.tri , el , M.xyz , M.tri );
  
  for it = 1 : Nits
  
    M = GEOMETRY_mergeEdges0( M , el , M0 );
    [ M.xyz, M.tri ] = GEOMETRY_removeBadTri( M.xyz , M.tri ,      M0.xyz , M0.tri );
    [ M.xyz, M.tri ] = GEOMETRY_subdivideTri( M.xyz , M.tri ,  0 , M0.xyz , M0.tri );
    [ M.xyz, M.tri ] = GEOMETRY_removeBadTri( M.xyz , M.tri ,      M0.xyz , M0.tri );
    
    fprintf( 'Iteration:  %2d    Output mesh: ( %6d triangles , %6d vertices ).\n' , it , size( M.tri ,1) , size( M.xyz ,1) );
  end

end

function M = GEOMETRY_mergeEdges0( M , el , M0 )

  if nargin < 3, M0 = M; end

  cases = M.xyz;
  el = el/sqrt(3);

  while size(cases,1)>0
  
    fk1 = M.tri(:,1);
    fk2 = M.tri(:,2);
    fk3 = M.tri(:,3);
  
    numfaces = (1:size(M.tri,1))';
  
    e1=sqrt(sum((M.xyz(fk1,:)-M.xyz(fk2,:)).^2,2));
    e2=sqrt(sum((M.xyz(fk1,:)-M.xyz(fk3,:)).^2,2));
    e3=sqrt(sum((M.xyz(fk2,:)-M.xyz(fk3,:)).^2,2));
  
    temp(:,1)=mean(e1);
    temp(:,2)=std(e1);
  
    ed1=sort([fk1 fk2 ]')';
    ed2=sort([fk1 fk3 ]')';
    ed3=sort([fk2 fk3 ]')';
  
    e1=[e1 numfaces ed1 ];
    e2=[e2 numfaces ed2 ];
    e3=[e3 numfaces ed3 ];
  
    e=[e1 ; e2 ; e3];
    e=e(e(:,1)<el,:);
  
    e=sortrows(e,1);
  
  
    [etemp,ia,ic]=unique(e(:,3),'rows','stable');
    e=e(ia,:);
    [etemp,ia,ic]=unique(e(:,4),'rows','stable');
    e=e(ia,:);
    [test,ia,ic]=unique(e(:,2),'rows','stable');
    e=e(ia,:);
  
    ind1=(1:2:(2*size(ia,1)-1))';
    ind2=(2:2:(2*size(ia,1)))';
    ind3=(1:2*size(ia,1))';
  
    test1=ones(2*size(ia,1),1);
    test1(ind1)=e(:,3);
    test1(ind2)=e(:,4);
  
    test1(:,2)=ones;
    test1(ind1,2)=(1:size(ia))';
    test1(ind2,2)=(1:size(ia))';
  
    [etemp1,ia,ic]=unique(test1(:,1),'stable');
  
    test1=(test1(ia,:));
  
    test1(:,3)=ones;
    test1(2:end,3)=test1(1:end-1,2);
  
    test1(:,4)=test1(:,3)-test1(:,2);
  
    indicesseries= test1(test1(:,4)==0,2);
    indicesseries=unique(indicesseries,'stable');
  
  
    e=e(indicesseries,:);
  
    cases=e(:,3:4);
  
    averages=(M.xyz(cases(:,1),:)+M.xyz(cases(:,2),:)).*0.5;
    M.xyz(cases(:,1),:)=averages;
    M.xyz(cases(:,2),:)=averages;

    [ M.xyz , M.tri ] = GEOMETRY_cleanPatch( M.xyz , M.tri );

    M = GEOMETRY_project0( M , M0.xyz );
    %M.xyz = GEOMETRY_project( M.xyz , M.tri , M0.xyz , M0.tri );

  end

end

function M = GEOMETRY_project0( M , X )

  TRS = triangulation( M.tri , M.xyz );
  normalsS = vertexNormal(TRS);
  
  [IDXsource,~] = knnsearch( X , M.xyz );
  
  vector_s_to_t = X(IDXsource,:) - M.xyz;
  
  M.xyz = M.xyz + [  sum( vector_s_to_t.*normalsS ,2)./( norm(normalsS).^2 ).*normalsS(:,1) ,...
                     sum( vector_s_to_t.*normalsS ,2)./( norm(normalsS).^2 ).*normalsS(:,2) ,...
                     sum( vector_s_to_t.*normalsS ,2)./( norm(normalsS).^2 ).*normalsS(:,3) ];

end
