function C = MeshClip( M , V , varargin )
% clip out the positives

  if isempty( M.xyz ) || isempty( M.tri )
    C = M;
    return;
  end

  if numel(varargin) == 1 && isfloat(varargin{1}) && numel(varargin{1}) > 1
    M.xyz______V______ = V;
    z = unique( varargin{1} );

%     for z = z(:).'
%       M = MeshClip( M , M.xyz______V______ - z , 'both' );
%     end
    
    for z = unique( z(:).' )
      c1 = MeshClip( M , M.xyz______V______ - z , false); if isempty( c1.tri ), continue; end
      c2 = MeshClip( M , M.xyz______V______ - z , true); if isempty( c2.tri ), continue; end
      M = MeshAppend( c1 , c2 );
    end
    
    C = rmfields( M , 'xyz______V______' );
    return;
  end


  insideOut = false;
  KP = false;
  for v = 1:numel(varargin)
    if ischar( varargin{v} ) && ...
       ( strcmpi( varargin{v} ,'kp' ) || strcmpi( varargin{v} ,'keepparent' ) || strcmpi( varargin{v} ,'keepparentedge' ) )
      KP = true; continue;
    end
    insideOut = varargin{v};
  end

  if isa( V , 'function_handle' )
    try, V = feval( V , M ); catch
    try, V = feval( V , M.xyz ); catch
      error('invalid function to evaluate on mesh');
    end; end
  elseif ischar( V )
    try, V = M.(['xyz',V]); catch
    try, V = M.(V); catch
      error('invalid attribute name.');
    end; end
  elseif isnumeric( V ) && isequal( size( V ) , [4 4] )
    V = distance2Plane( M.xyz , V , true );
  end
  if numel( V ) ~= size( M.xyz ,1)
    error('invalid scalar for contouring');
  end
  
  BOTH = false; SPLIT = false;
  if isequal( insideOut , 2), insideOut = 'both'; end
  if isequal( insideOut , 3), insideOut = 'split'; end
  if ~ischar( insideOut ) && numel( insideOut ) == 1 && ismember( insideOut ,[0,1] )
    insideOut = ~~insideOut;
  elseif ischar( insideOut )
  else, error('Invalid insideOut.');
  end
  if 0
  elseif islogical( insideOut ) && numel( insideOut ) == 1 &&  insideOut
    V = -V;
  elseif islogical( insideOut ) && numel( insideOut ) == 1 && ~insideOut
  elseif ischar( insideOut ) && ( strcmp( insideOut , 'both' ) || strcmp( insideOut , 'b' ) )
    BOTH = true;
  elseif ischar( insideOut ) && ( strcmp( insideOut , 'split' ) || strcmp( insideOut , 's' ) )
    BOTH = true; SPLIT = true;
  else, error( 'invalid option' );
  end


  if ~~BOTH && ~~SPLIT
    C = MeshAppend( MeshClip( M , V , true ) , MeshClip( M , V , false ) );
    return;
  end



  try, M = MeshOrderCells( M , V ); end
  s = reshape( sign( V( M.tri ) ) ,size(M.tri) ); s( ~s ) = 1;

  C = struct(); if isfield( M , 'texture' ), C.texture = M.texture; end

  P   = [];
  T   = []; Tid = [];
  W   = []; Wid = [];
  
  M.celltype = meshCelltype( M );
  switch M.celltype
    case 3

      t = all( bsxfun( @eq , s , [ -1 ,  1 ] ) ,2); 
      addCELLS( t ,2); P = [ P ; E(1,1) ; E(1,2) ];
      if BOTH
      addCELLS( t ,2); P = [ P ; E(1,2) ; E(2,2) ];
      end
      
      t = all( bsxfun( @eq , s , [  1 , -1 ] ) ,2); 
      addCELLS( t ,2); P = [ P ; E(1,2) ; E(2,2) ];
      if BOTH
      addCELLS( t ,2); P = [ P ; E(1,1) ; E(1,2) ];
      end
      
    case 5

      t = all( bsxfun( @eq , s , [ -1 , -1 ,  1 ] ) ,2); 
      addCELLS( t ,3); P = [ P ; E(1,1) ; E(2,3) ; E(1,3) ];
      addCELLS( t ,3); P = [ P ; E(1,1) ; E(2,2) ; E(2,3) ];
      if BOTH
      addCELLS( t ,3); P = [ P ; E(1,3) ; E(2,3) ; E(3,3) ];
      end

        
      t = all( bsxfun( @eq , s , [ -1 ,  1 , -1 ] ) ,2); 
      addCELLS( t ,3); P = [ P ; E(1,1) ; E(1,2) ; E(2,3) ];
      addCELLS( t ,3); P = [ P ; E(1,1) ; E(2,3) ; E(3,3) ];
      if BOTH
      addCELLS( t ,3); P = [ P ; E(1,2) ; E(2,2) ; E(2,3) ];
      end

      t = all( bsxfun( @eq , s , [ -1 ,  1 ,  1 ] ) ,2); 
      addCELLS( t ,3); P = [ P ; E(1,1) ; E(1,2) ; E(1,3) ];
      if BOTH
      addCELLS( t ,3); P = [ P ; E(1,2) ; E(2,2) ; E(1,3) ];
      addCELLS( t ,3); P = [ P ; E(1,3) ; E(2,2) ; E(3,3) ];
      end
      
    case 10
      
      t = all( bsxfun( @eq , s , [ -1 , -1 , -1 ,  1 ] ) ,2);
      addWEDGES( t );  P = [ P ; E(2,2) ; E(1,1) ; E(3,3) ; E(2,4) ; E(1,4) ; E(3,4) ];
      if BOTH
      addCELLS( t ,4); P = [ P ; E(1,4) ; E(2,4) ; E(3,4) ; E(4,4) ];
      end
      
      t = all( bsxfun( @eq , s , [ -1 , -1 ,  1 , -1 ] ) ,2); 
      addWEDGES( t );  P = [ P ; E(1,1) ; E(2,2) ; E(4,4) ; E(1,3) ; E(2,3) ; E(3,4) ];
      if BOTH
      addCELLS( t ,4); P = [ P ; E(1,3) ; E(2,3) ; E(3,3) ; E(3,4) ];
      end
      
      t = all( bsxfun( @eq , s , [ -1 , -1 ,  1 ,  1 ] ) ,2);
      addWEDGES( t );  P = [ P ; E(1,4) ; E(1,1) ; E(1,3) ; E(2,4) ; E(2,2) ; E(2,3) ];
      if BOTH
      addWEDGES( t );  P = [ P ; E(1,4) ; E(2,4) ; E(4,4) ; E(1,3) ; E(2,3) ; E(3,3) ];
      end
      
      t = all( bsxfun( @eq , s , [ -1 ,  1 ,  1 ,  1 ] ) ,2); 
      addCELLS( t ,4); P = [ P ; E(1,1) ; E(1,2) ; E(1,3) ; E(1,4) ];
      if BOTH
      addWEDGES( t );  P = [ P ; E(1,3) ; E(1,2) ; E(1,4) ; E(3,3) ; E(2,2) ; E(4,4) ];
      end
      
    otherwise, error('not implemented for this celltype');
  end
  

  P = double( P );
  P = sort( P , 2 );
  P = [ P , V( P(:,1) ) ./ ( V( P(:,1) ) - V( P(:,2) ) ) ];
  P( ~isfinite( P(:,3) ) ,3) = 0;
  w   = P(:,3) == 1; P(w,1) = P(w,2); P(w,3) = 0;
  w   = P(:,3) == 0; P(w,2) = 1;

  [P,~,b] = unique( P , 'rows' , 'stable' ); T = reshape( b( T ) ,size(T) );
  
  if ~isempty( W )
    W = reshape( b( W ) ,size(W) );

    W = struct( 'tri' , W , 'triWid' , Wid , 'celltype' , 13 );
    % real coordinates of every (interpolated) node.  Used ONLY to orient the
    % resulting tets, never to choose the split diagonals -- those come purely
    % from node indices so that two wedges sharing a quad face always split it
    % along the same diagonal ( => conforming tet mesh, no cracks/T-junctions ).
    XYZ = bsxfun( @times , 1-P(:,3) , M.xyz( P(:,1) ,:) ) + ...
          bsxfun( @times ,   P(:,3) , M.xyz( P(:,2) ,:) );
    
    W = MeshClip_wedge2tets( W , XYZ );

    Tid = [ Tid ; W.triWid ];
    T   = [ T   ; W.tri    ];
  end
    
  for i = 1:size(T,2)-1
    for j = i+1:size(T,2)
      w = T(:,i) == T(:,j);
      T(w,:) = [];
      Tid(w) = [];
    end
  end

  %inclussion of the non-clipped original cells
  oP = ( 1:size(M.xyz,1) ).'; oP(:,2) = 1; oP(:,3) = 0;
  oT = all( bsxfun( @eq , s , -1 ) ,2);
  if BOTH, oT = oT | all( bsxfun( @eq , s , 1 ) ,2); end
  Tid = [ find( oT(:) ) ; Tid ];
  T = [ M.tri( oT ,:) ; size(oP,1) + T ];
  P = [ oP ; P ];
  
  [P,~,b] = unique( P , 'rows' , 'stable' ); T = reshape( b( T ) ,size(T) );

  nid = unique( T(:) );
  map = zeros( numel(nid) ,1);
  map( nid ) = 1:numel(nid);
  P = P( nid ,:);
  T = reshape( map( T ) ,size(T) );

  %[P,~,b] = unique( P , 'rows' , 'stable' ); T = reshape( b( T ) ,size(T) );
  [Tid,ord] = sort( Tid ); T = T( ord ,:);
  %
  
  
  for f = fieldnames( M ).', f = f{1};
    if strncmp( f , 'xyz' , 3 )
      C.(f) = bsxfun( @times , 1-P(:,3) , double( M.(f)( P(:,1) ,:,:,:,:,:) ) ) + bsxfun( @times , P(:,3) , double( M.(f)( P(:,2) ,:,:,:,:,:) ) );
    elseif strcmp( f , 'tri' )
      C.(f) = T;
    elseif strncmp( f , 'tri' , 3 )
      C.(f) = M.(f)( Tid ,:,:,:,:,:,:);
    end
  end
  if KP
    fn = fieldnames(C); fn = sort( fn( strncmp( fn , 'xyzParentEdge' ,13) ) );
    for f = fn(end:-1:1).', C = renameStructField( C , f{1} , [ f{1} , '_' ] ); end
    P( P(:,3) == 0 ,2) = 0;
    C.xyzParentEdge = P;
  end
  C.celltype = M.celltype;
  C = mergestruct( C , M ,'<' );
  
  
  function e = E(a,b), e = M.tri( t , [a,b] ); end
  function addCELLS( t , n )
    nt = sum( t ); if ~nt, return; end
    T = [ T ;  size( P ,1) + reshape( 1:(nt*n) ,nt,n) ];
    Tid = [ Tid ; find( t(:) ) ];
  end  
  function addWEDGES( t )
    nt = sum( t ); if ~nt, return; end
    W = [ W ;  size( P ,1) + reshape( 1:(nt*6) ,nt,6) ];
    Wid = [ Wid ; find( t(:) ) ];
  end  
  
end


function W = MeshClip_wedge2tets( W , XYZ )
% Split each triangular prism (wedge, celltype 13) into 3 tetrahedra.
%
% The diagonal of every quadrilateral face is chosen through that face's
% MINIMUM-INDEX vertex.  Because that choice depends only on the (global) node
% indices -- and two wedges that share a quad face share those very indices --
% both neighbours pick the same diagonal, so the tetrahedral mesh stays
% conforming ( no cracks / T-junctions ), and -- like a Delaunay tie-break by
% point id -- it needs no neighbour information, only a unique global node
% numbering.  This is the same principle the previous Python/VTK helper relied
% on ( vtkOrderedTriangulator ); here we use the fixed prism table instead of a
% general triangulation.
%
% References:
%   [1] W. J. Schroeder, B. Geveci, M. Malaterre, "Compatible Triangulations of
%       Spatial Decompositions", Proc. IEEE Visualization 2004 (VIS'04),
%       pp. 211-218.  ( general principle; basis of VTK's ordered triangulation,
%       i.e. what the previous Python/VTK helper used ).
%   [2] J. Dompierre, P. Labbe, M.-G. Vallet, R. Camarero, "How to Subdivide
%       Pyramids, Prisms and Hexahedra into Tetrahedra", Proc. 8th International
%       Meshing Roundtable (IMR'99), Lake Tahoe CA, 1999, pp. 195-204.
%       ( the minimum-index prism subdivision table implemented below ).
%   [3] R. M. J. Kramer, "Cutting Tetrahedra by Node Identifiers", Sandia
%       National Laboratories, report SAND-2015-3830, May 2015.  ( the
%       tetrahedron-clipping specialisation matching this file's use case ).
%
% Wedge node convention ( as produced by addWEDGES ):
%   nodes [1 2 3] and [4 5 6] are the two triangular faces,
%   vertical edges are 1-4, 2-5, 3-6.

  w   = W.tri;                       % [nW x 6] wedge connectivity
  wid = W.triWid(:);                 % [nW x 1] parent-cell id of each wedge
  nW  = size( w ,1 );

  if nW == 0
    W.tri = zeros(0,4); W.triWid = zeros(0,1); W.celltype = 10;
    return;
  end

  % 1) canonicalise: relabel so the smallest-index vertex sits in column 1.
  %    One relabelling per possible position of the minimum; each is a symmetry
  %    of the prism ( it maps {1,2,3}/{4,5,6} triangles and the 1-4,2-5,3-6
  %    vertical edges onto themselves ), so the [a b c | d e f] structure holds.
  PERMS = [ 1 2 3 4 5 6 ; ...
            2 3 1 5 6 4 ; ...
            3 1 2 6 4 5 ; ...
            4 5 6 1 2 3 ; ...
            5 6 4 2 3 1 ; ...
            6 4 5 3 1 2 ];
  [~,p] = min( w ,[],2 );                                  % position of the min
  w = w( sub2ind( size(w) , repmat( (1:nW).' ,1,6 ) , PERMS(p,:) ) );

  a = w(:,1); b = w(:,2); c = w(:,3);   % a = min index  -> bottom triangle a,b,c
  d = w(:,4); e = w(:,5); f = w(:,6);   %                   top triangle    d,e,f

  % Faces (a,b,e,d) and (c,a,d,f) both contain the global minimum a, so their
  % diagonals are forced through a ( a-e and a-f ).  This also rules out the two
  % invalid "cyclic" diagonal configurations.  Only the far face (b,c,f,e) is
  % still free: split it through its own minimum-index vertex.
  [~,mq] = min( [ b c e f ] ,[],2 );    % 1->b 2->c 3->e 4->f
  bf = ( mq == 1 | mq == 4 );           % min in {b,f} -> diagonal b-f, else c-e

  % 2) emit the 3 tets of every wedge ( row i, nW+i and 2*nW+i share wedge i ).
  T = zeros( 3*nW , 4 );
  T( 1:nW ,: ) = [ a d e f ];                                        % common tet
  T( nW   + find(~bf) ,: ) = [ a(~bf) b(~bf) c(~bf) e(~bf) ];        % diagonal c-e
  T( 2*nW + find(~bf) ,: ) = [ a(~bf) c(~bf) f(~bf) e(~bf) ];
  T( nW   + find( bf) ,: ) = [ a( bf) b( bf) c( bf) f( bf) ];        % diagonal b-f
  T( 2*nW + find( bf) ,: ) = [ a( bf) b( bf) f( bf) e( bf) ];
  Tid = repmat( wid , 3 , 1 );

  % 3) orient every tet to positive volume ( swapping two nodes flips the sign
  %    but keeps the same faces, so conformity is preserved ).
  if nargin > 1 && ~isempty( XYZ )
    v1 = XYZ( T(:,1) ,:); v2 = XYZ( T(:,2) ,:);
    v3 = XYZ( T(:,3) ,:); v4 = XYZ( T(:,4) ,:);
    neg = dot( v2-v1 , cross( v3-v1 , v4-v1 , 2 ) , 2 ) < 0;
    T( neg ,[3 4] ) = T( neg ,[4 3] );
  end

  W.tri      = T;
  W.triWid   = Tid;
  W.celltype = 10;

end


function str = renameStructField( str , oldFieldName , newFieldName )
%RENAMESTRUCTFIELD  rename a struct field, preserving field order (MathWorks FEX).
%   Local copy so the 'kp'/'keepparent' option works without an external dependency.
  if ~strcmp( oldFieldName , newFieldName )
    allNames = fieldnames( str );
    isOverwriting = ~isempty( find( strcmp( allNames , newFieldName ) , 1 ) );
    matchingIndex = find( strcmp( allNames , oldFieldName ) );
    if ~isempty( matchingIndex )
      allNames{ matchingIndex(1) } = newFieldName;
      [ str.(newFieldName) ] = deal( str.(oldFieldName) );
      str = rmfield( str , oldFieldName );
      if ~isOverwriting
        str = orderfields( str , allNames );
      end
    end
  end
end
