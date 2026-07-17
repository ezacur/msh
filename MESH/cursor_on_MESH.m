function txt = cursor_on_MESH( obj , ev )

if nargin == 1 && ishandle( obj )
  hFig = ancestor( obj , 'figure' );
  dcm_obj = datacursormode( hFig );
  set( dcm_obj ,'DisplayStyle','datatip',...
    'SnapToDataVertex','on',...
    'Enable','on',...
    'UpdateFcn',@(obj,ev)cursor_on_MESH(obj,ev) );
  if nargout
    txt = dcm_obj;
  end
  return;
end


txt = {};

T = get( ev , 'Target' );
if ~strcmp( get(T,'Type') , 'patch' )
  return;
end

M = Mesh( T );

xyz = get( ev , 'Position' );


txt{end+1} = sprintf('[%g,%g,%g]', xyz );
txt{end} = [ 'position: ' , txt{end} ];



node_ids = find( all( bsxfun( @eq , M.xyz , xyz ) ,2) );
txt{end+1} = sprintf('%d ,', node_ids );
txt{end} = [ 'nodes: ' , txt{end}(1:end-1) ];

face_ids = find( any( ismember( M.tri , node_ids ) ,2) );
txt{end+1} = sprintf('%d ,', face_ids );
txt{end} = [ 'faces: ' , txt{end}(1:end-1) ];

%%

delete( findall( get(T,'Parent') , 'Tag','cursor_on_MESH_face_ids') );

txyz = meshFacesCenter( M );
txyz = txyz( face_ids ,:);
text( txyz(:,1) , txyz(:,2) , txyz(:,3) , ...
  arrayfun( @num2str , face_ids , 'un' , 0 ) ,...
  'Parent' , get(T,'Parent') , ...
  'HorizontalAlignment','c','VerticalAlignment','m','FontWeight','bold',...
  'BackgroundColor',[0 1 1],'EdgeColor','k','Margin',3,'Color','r' ,...
  'Tag','cursor_on_MESH_face_ids');




delete( findall( get(T,'Parent') , 'Tag','cursor_on_MESH_node_ids') );

node_ids = setdiff( M.tri( face_ids ,:) , node_ids );
txyz = M.xyz( node_ids ,:);
text( txyz(:,1) , txyz(:,2) , txyz(:,3) , ...
  arrayfun( @num2str , node_ids , 'un' , 0 ) ,...
  'Parent' , get(T,'Parent') , ...
  'HorizontalAlignment','c','VerticalAlignment','m','FontWeight','demi',...
  'BackgroundColor',[1 1 0],'EdgeColor','k','Margin',3 ,...
  'Tag','cursor_on_MESH_node_ids');

end
