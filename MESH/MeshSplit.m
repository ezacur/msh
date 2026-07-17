function M = MeshSplit( M , SE )
%MESHSPLIT  Cut a mesh open along edges (or a polyline at vertices) by
%           DUPLICATING the nodes on the cut, so the pieces disconnect.
%
%   M = MeshSplit( M , SE )              triangle mesh: SE is n x 2 EDGES (node
%                                        id pairs) to cut along. Nodes on the
%                                        cut are duplicated (interior cut ->
%                                        the two sides separate; the geometry
%                                        is untouched, only connectivity).
%   M = MeshSplit( M , -ang )            triangle mesh, scalar <= 0: cut every
%                                        edge whose inter-face NORMAL angle is
%                                        >= ang degrees (creases; see
%                                        meshCellsContact).
%   M = MeshSplit( M , 'nonmanifold' )   cut every edge shared by MORE than 2
%                                        faces (each face keeps its own copy).
%   M = MeshSplit( M , V )               polyline (celltype 3): V is a list of
%                                        VERTEX ids; every extra segment
%                                        incident to each v gets its own copy
%                                        of v (a chain splits in two, a Y in
%                                        three).
%
%   Every xyz* field is duplicated along with the nodes. Original nodes keep
%   their ids; the copies are APPENDED at the end. Uses Alec Jacobson's
%   cut_edges construction (sparse corner graph + one conncomp).
%
% See also meshCellsContact, meshSeparate, MeshBoundary, MeshReOrderNodes.

  M.celltype = meshCelltype( M );

  if M.celltype == 5 && isscalar( SE ) && SE <= 0
    [E,~,A] = meshCellsContact( M );
    M = MeshSplit( M , E( A >= -SE ,:) );
    
    return;
    
  elseif ischar( SE ) && strcmpi( SE , 'nonmanifold' )
    [E,C] = meshCellsContact( M );
    M = MeshSplit( M , E( cellfun('prodofsize',C) > 2 ,:) );
    
    return;
  end


  if isempty( SE )
    return;
  end
  if ~isnumeric( SE ) || any( SE(:) < 1 ) || any( SE(:) > size( M.xyz ,1) ) || any( SE(:) ~= round( SE(:) ) )
    error('MeshSplit:SE','SE must contain valid node ids (integers in 1..%d).', size( M.xyz ,1) );
  end

  switch M.celltype
    case 3
      for v = SE(:).'
        e = find( any( M.tri == v ,2) ).';
        for ee = e(2:end)
          for f = fieldnames( M ,'^xyz').', M.(f{1}) = M.(f{1})([1:end,v],:,:,:,:,:,:); end
          M.tri( ee , M.tri( ee ,:) == v ) = size( M.xyz ,1);
        end
      end
      return;
      

      
    case 5
      if 0
        rand('seed',0);
        M = struct();
        M.xyz = [ rand(10,2) ];
        M.tri = delaunayn( M.xyz(1:10,:) );
        M.xyz = [ M.xyz ; bsxfun( @minus , rand(5,2) ,[1 0] ) ];
        
        M = MeshReOrderNodes( M , randperm(size(M.xyz,1)) );
        M.tri = [ M.tri ; 2 15 3 ; 4 2 13 ; 3 13 2 ];
        
        
        clf;set(gcf,'Position',[965,49,952,964])
        subplot(211);
        plot3d( M.xyz , 'o1kr7','eq')
        hplotMESH( M , 'textpoint','facealpha',0.2,'td',meshFacesConnectivity(M));
        
        
        SE = [8 4 2 15 14 10];
        % SE = [4 2 15];
        
        hplot3d( M.xyz(SE,:) , 'r2' )
        
        SE = [ SE(1:end-1).' , SE(2:end).' ];
        MM = MeshSplit( M , SE )
        
%         MM = Mesh( MeshRelax( MeshSmooth( MM ,10)  ))
        
        subplot(212);
        plot3d( MM.xyz , 'o1kr7','eq')
        hplotMESH( MM , 'textpoint','facealpha',0.2,'td',meshFacesConnectivity(MM));
        hplotMESH( MeshBoundary(MM) ,'edgecolor','r','linewidth',2)
        
        %%
      end
      
      if size( SE ,2) ~= 2
        error('splitting edges should be n x 2');
      end
      
      T = M.tri;

      nP = size( M.xyz ,1);
      Oid = ( 1:nP ).';

      a = unique( T(:) );
      Oid = Oid( a );
      a( a ) = 1:numel(a);
      T = a(T);
      SE = reshape( a(SE) , size(SE) );


      %%from Alec Jacobson's cut_edges function!!
      nT  = size(T,1);
      nT3 = 3*nT;
      F = reshape( 1:nT3 , nT , 3 );

      allE = sort( [ T(:,[2,3]); T(:,[3,1]); T(:,[1,2]) ] ,2);
      T = double( T(:) );

      [E,~,IC] = unique( allE ,'rows');
      nE = size( E ,1);

      [~,P] = setdiff( E , sort(SE,2) ,'rows');
      if size(P,1) == nE, return; end
      A = sparse( P , P , 1 , nE , nE );
      B = sparse( IC , F , 1 , nE , nT3 );
      C = sparse( F(:,[1,2,3,1,2,3]) , F(:,[2,3,1,3,1,2]) , 1 , nT3 , nT3 );
      D = sparse( F , T , 1 , nT3 , nP );

      G = ( C * ( B.' * A * B ) * C.' )  &  ( D * D.' );
      [~,J] = conncomp(G);

      F = J( F );
      Pid( J ,1) = T;
      %%thanks Alec.

      Pid = Oid( Pid );

      
    case 10
      error('not implemented for this celltype... and maybe it will never be implemented');

  end


  ID  = [ Pid ; setdiff( 1:nP , Pid ).' ];
  ord = zeros( 1 , numel(ID) );
  Z = false( 1 , numel(ID) );
  m = nP;
  for i = 1:numel(ID)
    if Z(ID(i)), m = m + 1; ord(i) = m;
    else,        ord(i) = ID(i);
    end
    Z(ord(i)) = true;
  end
  [~,ord] = sort(ord);
  F = iperm( ord , F );
  ID = ID( ord );
      
  M.tri = F;
  for f = fieldnames( M ).', f = f{1};
    if ~strncmp( f , 'xyz',3), continue; end
    M.(f) = M.(f)( ID ,:,:,:,:,:,:);   %extra colons: keep >2-dim fields' shape
  end

end

function [S,C] = conncomp(G)
  % CONNCOMP Drop in replacement for graphconncomp.m from the bioinformatics
  % toobox. G is an n by n adjacency matrix, then this identifies the S
  % connected components C. This is also an order of magnitude faster.
  %
  % [S,C] = conncomp(G)
  %
  % Inputs:
  %   G  n by n adjacency matrix
  % Outputs:
  %   S  scalar number of connected components
  %   C  

  % Transpose to match graphconncomp
  G = G';

  [p,q,r] = dmperm(G+speye(size(G)));
  S = numel(r)-1;
  C = cumsum(full(sparse(1,r(1:end-1),1,1,size(G,1))));
  C(p) = C;
end

%(a dead local inter_faces_angle used to live here: never called, and it
% computed N(E(:,3))-N(E(:,3)) = angle 0 always; meshCellsContact provides the
% correct inter-face angle used by the scalar-threshold form above)