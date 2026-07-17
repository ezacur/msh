function [ M , ids ] = MeshForceNode( M , X )
%MESHFORCENODE  Insert points into a triangular mesh, forcing them to be nodes.
%
%   [ M , ids ] = MeshForceNode( M , X )
%
%   Modifies the triangular surface mesh M so that every point of X becomes a
%   VERTEX of the mesh. Each point is first PROJECTED onto the surface of M
%   (via distanceFrom), and then inserted according to where its projection
%   lands:
%     * ON an existing vertex  -> reused, no new geometry.
%     * INSIDE a face          -> the face is split into 3 (a fan around the
%                                 new vertex): 1 triangle -> 3 triangles.
%     * ON an edge             -> the (up to 2) faces sharing that edge are
%                                 each split in 2.
%   So X are snapped to the surface, NOT kept at their input coordinates; the
%   inserted node sits at the projected location  M.xyz( abs(ids) , : ).
%
%   INPUTS
%     M : struct with fields .xyz (nV x 3) and .tri (nF x 3). Any extra field
%         named  xyz*  (e.g. xyzF, xyzRGB) is treated as a PER-VERTEX attribute
%         and interpolated barycentrically at inserted vertices; any field
%         named  tri*  (e.g. triF) is a PER-FACE attribute and is copied from
%         the parent face onto the child faces. (The plain 'xyz'/'tri' — 3-char
%         names — are the geometry; other 3-char fields are ignored.)
%     X : P x 3 list of points to force into the mesh.
%
%   OUTPUTS
%     M   : the densified mesh (more vertices / faces than the input).
%     ids : P x 1. For each input point i, WHERE it ended up. The node itself
%           is always  M.xyz( abs(ids(i)) , : ). The value encodes the kind:
%
%             ids(i) > 0  and <= nV0  : coincided with an EXISTING vertex
%                                       (nV0 = size of the ORIGINAL M.xyz;
%                                        record it before calling if you need
%                                        to tell this case apart from the next)
%             ids(i) >  nV0           : inserted INSIDE a face (face split 1->3)
%             ids(i) <  0             : inserted ON an edge; its vertex index
%                                       is  -ids(i)  (edge faces split in 2)
%
%           i.e. the SIGN flags edge-insertions (negative) vs vertex/interior
%           (positive), and the MAGNITUDE is always the vertex index in the
%           returned M. On return every ids(i) is nonzero and finite (see the
%           THRESHOLD note below).
%
%   THRESHOLD (hard-coded at the top of this file) is the tolerance used both
%   for vertex coincidence (a point closer than THRESHOLD to a vertex reuses
%   it) and for the interior-vs-edge decision (a barycentric coordinate below
%   THRESHOLD means the point lies on that edge). Points are inserted one per
%   face/edge per pass (so several points on the same face/edge are handled in
%   successive passes). If a whole pass inserts NOTHING (a numerically stuck /
%   degenerate case), THRESHOLD is relaxed x10 with a warning, which guarantees
%   termination: once large enough, remaining points simply snap to their
%   nearest vertex.
%
%   Requires VTK (run enableVTK once): uses vtkClosestPoint and distanceFrom.
%
% See also distanceFrom, vtkClosestPoint, MeshSubdivide, sphereMesh.

  THRESHOLD = 1e-10;   % tolerance: vertex-coincidence distance AND barycentric
                       % interior/edge split. Hard-coded on purpose; tune here.

if 0
  %%
  X = randn( 2000  ,3); X = bsxfun( @rdivide , X , fro( X ,2) )*1.15;
  M = sphereMesh( 3  );
  M.triF = rand( size(M.tri,1) ,1);
  M.xyzF = M.xyz(:,3); 
  [~,X] = distanceFrom( X , M ); 

  plot3d( X(:,1:3) ,'1okw6','eq');  hplotMESH( M ,'EdgeColor','r', 'LineWidth',2,'nf');
  [ MM , ids ] = MeshForceNode( M , X );
  hplotMESH( MM ,'td','F'); colormap jet
%   hplot3d( MM.xyz( abs(ids) ,:) , X , '.-' )
   
  w = ids > size(M. xyz,1);             hplot3d( X(w,:) , MM.xyz(  ids(w) ,:) , '1obb4-');
  w = ids > 0 & ids <= size(M.xyz,1);   hplot3d( X(w,:) , MM.xyz(  ids(w) ,:) , '1okg8-');
  w = ids < 0 & isfinite( ids );        hplot3d( X(w,:) , MM.xyz( -ids(w) ,:) , '1okw6-');
  ze
  %%

end

  xyzAtts = {}; triAtts = {};
  for f = fieldnames( M ).', f = f{1};
    if numel( f ) <= 3, continue; end
    if strncmp( f , 'xyz' ,3)
      xyzAtts{end+1,1} = f;
      if isnumeric( M.(f) ), M.(f) = double( M.(f) ); end
    end
    if strncmp( f , 'tri' ,3)
      triAtts{end+1,1} = f;
    end
  end

  [ ~ , X ] = distanceFrom( X , M );
  
%   [ ~ , ~ , ~ , B ] = distanceFrom( X , M );
%   X( all( B > 1e-6 ,2) ,:) = [];


  X(:,4) = 0;
  prevLeft = Inf;                 % #unassigned points at the start of last pass
  while ~all( X(:,4) )
    nLeft = sum( ~X(:,4) );
    if nLeft >= prevLeft
      % the whole previous pass inserted NOTHING: borderline/stuck points.
      % Relax THRESHOLD so they get classified (eventually they snap to their
      % nearest existing vertex, which guarantees the loop terminates).
      THRESHOLD = THRESHOLD * 10;
      warning( 'MeshForceNode:noProgress' , ...
        'no point could be inserted; relaxing THRESHOLD to %g (%d point(s) left).' , THRESHOLD , nLeft );
      if THRESHOLD > 1e40      % backstop: 1e-10*10^50 dwarfs any coordinate
        error( 'MeshForceNode:stuck' , 'cannot insert %d remaining point(s).' , nLeft );
      end
    end
    prevLeft = nLeft;

    p = vec( find( ~X(:,4) ) );

    %% on vertices
    [v,~,d] = vtkClosestPoint( struct('xyz',double(M.xyz)) , X(p,1:3) );
    w = d < THRESHOLD;
    if any(w)
      pp = p(w,:);  p(w,:) = [];

      X(pp,4) = v(w,:);
    end
    if isempty( p ), continue; end

    %% interior in faces
    [ ~ , ~ , F , B ] = distanceFrom( X(p,1:3) , M );
    if numel(p) > 1
      [~,o] = sort( fro( B-1/3 ,2) );
      p = p(o,:);
      F = F(o,:);
      B = B(o,:);
    end

    w = all( B > THRESHOLD ,2);
    if any(w)
      pp = p(w,:);  p(w,:) = [];
      FF = F(w,:);  F(w,:) = [];
      BB = B(w,:);  B(w,:) = [];
  
      [~,w] = unique( FF ,'first'); w = sort(w);
      pp = pp(w,:);     
      FF = FF(w,:);
      BB = BB(w,:);
  
      nV = size( M.xyz ,1);
      nP = numel(pp);
  
      X(pp,4) = nV+(1:nP);
  
      M.xyz = [ M.xyz ; X(pp,1:3) ];
      for a = xyzAtts(:).', a = a{1};
        M.(a)( nV+(1:nP) ,:,:,:,:) = M.(a)( M.tri(FF,1) ,:,:,:,:) .* BB(:,1) + ...
                                     M.(a)( M.tri(FF,2) ,:,:,:,:) .* BB(:,2) + ...
                                     M.(a)( M.tri(FF,3) ,:,:,:,:) .* BB(:,3);
      end 
      
      M.tri = [ M.tri ; M.tri(FF,[2,3]) , X(pp,4) ; M.tri(FF,[3,1]) , X(pp,4) ];
      M.tri(FF,3) = X(pp,4);
      for a = triAtts(:).', a = a{1};
        M.(a) = [ M.(a) ; M.(a)( FF ,:,:,:) ; M.(a)( FF ,:,:,:) ];
      end
      %figure; plotMESH( M ); hplot3d( X(:,1:3) ,'1okr4');
    end
    if isempty( p ), continue; end
    
    %% on edges
    nV = size( M.xyz ,1);
    [ ~ , ~ , F , B ] = distanceFrom( X(p,1:3) , M );
    if numel(p) > 1
      [~,o] = sort( fro( B-1/3 ,2) );
      p = p(o,:);
      F = F(o,:);
      B = B(o,:);
    end

    %%%%%%%M.tri = double( M.tri );
    M.tri = double( M.tri );
    E = sort( M.tri( F ,:) .* bsxfun( @lt , min(B,[],2) , B ) ,2); E = E(:,2:3);
    [~,w] = unique( E ,'rows','first'); w = sort(w);
    pp = p(w);
    EE = E(w,:); EE(:,3) = nV + ( 1:size(EE,1) ).';

    X(pp,4) = -(EE(:,3));
    M.xyz = [ M.xyz ; X(pp,1:3) ];
    B = fro( X(pp,1:3) - M.xyz( EE(:,1) ,:) ,2) ./ fro( M.xyz( EE(:,2) ,:) - M.xyz( EE(:,1) ,:) ,2);
    for a = xyzAtts(:).', a = a{1};
      M.(a)( EE(:,3) ,:,:,:,:) = M.(a)( EE(:,1) ,:,:,:,:) .* (1-B) + ...
                                 M.(a)( EE(:,2) ,:,:,:,:) .*    B;
    end 



    for e = 1:numel( pp )
      fs = sum( ismember( M.tri , EE(e,1:2) ) ,2) == 2;
      FS0 = M.tri( fs ,:); FS1 = FS0;
      FS0( FS0 == EE(e,1) ) = EE(e,3);
      M.tri( fs ,:) = FS0;

      FS1( FS1 == EE(e,2) ) = EE(e,3);
      M.tri = [ M.tri ; FS1 ];
      for a = triAtts(:).', a = a{1};
        M.(a) = [ M.(a) ; M.(a)( fs ,:,:,:) ];
      end
    end

%     b = arrayfun( @(e)find(any(b==e,2)).' , 1:size(EE) ,'un',0).';
%     try
%       b = cell2mat( b );
%     catch
%       n = cellfun('prodofsize',b); m = max(n);
%       b( n < m ) = cellfun( @(e)[e,zeros(1,m-numel(e))] , b( n < m ) ,'un',0);
%       b = cell2mat( b );
%     end
% 
%     
%     for c = 1:size( b ,2)
%       w = ~~b(:,c);
% 
%     end
% 
% 
%     Fid = find( ~~b );
%     Xid =    b( ~~b );
%     EE  = EE( b( ~~b ) ,1:2);
%     [~,o] = sort( Xid );
%     Xid = Xid(o);
%     Fid = Fid(o);
%     EE  = EE( o ,:);
% 
%     nV = size( M.xyz ,1);
%     Xid = Xid + nV;
%     M.xyz = [ M.xyz ; X(pp,1:3) ];
%     for a = xyzAtts(:).', a = a{1};
%       M.(a)( nV+(1:numel(pp)) ,:,:,:,:) = NaN;
%     end
%     X(pp,4) = -unique( Xid );
%     for x = 1:numel( Xid )
%       F = M.tri( Fid(x) ,:);
%       w = find( ismember( F , EE(x,:) ) );
%       F1 = F; F1( w(2) ) = Xid(x);
%       F2 = F; F2( w(1) ) = Xid(x);
%       M.tri = [ M.tri ; F1 ; F2 ];
%       for a = triAtts(:).', a = a{1};
%         M.(a) = [ M.(a) ; M.(a)( Fid(x) ,:,:,:) ; M.(a)( Fid(x) ,:,:,:) ];
%       end
%     end
%     M.tri( Fid ,:) = [];
%     for a = triAtts(:).', a = a{1};
%       M.(a)( Fid ,:) = [];
%     end


%     F = M.tri; F(:,end+1) = ( 1:size(M.tri,1) ).';
%     F = [ F( ~~b ,:) , b( ~~b ) , E( b( ~~b ) ,1:2) ];
% 
%     X( pp ,4) = -Inf;
  end

  ids = X(:,4);

end
