function [C,d] = MeshQuery( A , B , varargin )
%MESHQUERY  Transfer the fields of mesh A onto query points / mesh B.
%
%   C = MeshQuery( A , B )
%   [C,d] = MeshQuery( A , B , mode , THRES , outW )
%
%   For every query location, MeshQuery finds the closest element of A and
%   barycentrically interpolates A's per-vertex fields (xyz*) onto it; A's
%   per-face fields (tri*) are copied from the closest element. The result C
%   keeps B's geometry (C.xyz / C.tri) and gains the transferred A.* fields.
%
%   INPUTS
%     A     mesh struct with .xyz and (optionally) .tri, plus any .xyz*/.tri*
%           fields to transfer. Supported A cell types:
%             point cloud / vertices  -> nearest point (knnsearch)
%             edges (celltype 3)      -> closest segment  (distanceFrom)
%             triangles (celltype 5)  -> closest triangle (distanceFrom)
%             tetrahedra (celltype 10)-> containing tetra (tsearchn)
%           (A point cloud with no .tri is treated as vertex cells.)
%     B     query: an N-by-3 numeric array of points, OR a mesh struct (both
%           its vertices AND its face centers are queried and transferred).
%     mode  'NaN' (default): query points farther than THRES from A get NaN in
%                 the interpolated NUMERIC fields (non-numeric take the nearest
%                 vertex value). 'closest': every query is snapped to A (for
%                 tetrahedra, outside points are projected onto A's boundary).
%     THRES distance tolerance for the 'NaN' mode (default 1e-10).
%     outW  logical; if true, C also carries .xyzQueryWeights (the raw
%           barycentric weights per query point).
%
%   OUTPUTS
%     C     B's geometry with A's fields interpolated onto it.
%     d     distance from each query point to A (NaN for tetrahedral A, whose
%           containment test does not produce a distance).
%
%   NOTE: for edge/triangle A the barycentric weights come from distanceFrom;
%   any inaccuracy there propagates to the interpolated fields.
%
%   See also meshMapPoints, distanceFrom, vtkClosestElement, tsearchn.
if 0
  
  A = tetgen( sphereMesh(2) ,'a',0);
  A = MeshReOrderNodes( A , randperm(size(A.xyz,1)) );
  A.xyzX = A.xyz;
  
  XYZ = ndmat( [0 1],[0 1],[0 1] );
  XYZ = [1 1 1];
  XYZ = rand( 100 , 3 )*10;
  
  C = MeshQuery( A , XYZ , 'closest' );
  
  CC = MeshQuery( A , XYZ , 'NaN' , true );
  w = all( isnan(CC.xyzQueryWeights) ,2);
  CP = XYZ;
  [~,CP(w,:)]=vtkClosestElement( MeshBoundary(A) , XYZ( w ,:) );
  maxnorm( C.xyzX - CP )
  
  %%
end


  Omode = 'NaN';
  outW  = false;
  THRES = 1e-10;
  while numel(varargin)
    v = varargin{1}; varargin(1) = [];
    if ischar( v ) && isrow( v )
      Omode = v;
    elseif islogical( v ) && isscalar( v )
      outW  = v;
    elseif isnumeric( v ) && isscalar( v )
      THRES = v;
    else
      error('invalid varargin');
    end
  end  
  

  if isnumeric( B )
    B(:,end+1:3) = 0;
    B = struct( 'xyz' , B , 'tri' , ( 1:size(B,1) ).' );
  end
  
  if ~isfield( A , 'tri' ) || isempty( A.tri )
    A.tri = ( 1:size( A.xyz ,1) ).';   %point cloud -> vertex cells: nearest-point transfer
  end
  celltype = meshCelltype( A );
  if numel( celltype ) ~= 1
    error('MeshQuery:homogeneous','MeshQuery only handles homogeneous meshes (a single celltype).');
  end
  
  C = Mesh(B,0);
  
  XYZ = double( C.xyz );
  FC  = meshFacesCenter( C );
  if isequal( XYZ , FC ), FC = []; end
  
  
  AA  = Mesh(A,0);
  Eout = false( size(XYZ,1) ,1);
  Fout = false( size(FC ,1) ,1);
  d    = NaN( size(XYZ,1) ,1);   %always defined (tsearchn/celltype 10 leaves it NaN)
  if false
  elseif celltype == 1
    %[ Eid , ~ , d ] = vtkClosestPoint( struct('xyz',double(resize( AA.xyz ,[],3)),'tri',zeros(0,3)) , double( resize( XYZ ,[],3) ) );
    [ Eid , d ] = knnsearch( resize( AA.xyz ,[],3) , resize( XYZ ,[],3) );
    W = ones( size( Eid ) );
    switch lower( Omode )
      case 'nan'
        Eout = d > THRES;
    end

    if ~isempty( FC )
    Fout = true( size( FC ,1) , 1 );
    end
    
  elseif celltype == 3
    
    %[ Eid , ~ , d , W ] = distancePoint2Segments( XYZ , { AA.xyz , AA.tri } ); W = [ 1 - W , W ];

    %[Eid , cp , d ] = ClosestSegment( double( AA.xyz ) , double( AA.tri ) , double( XYZ ) );

    [ d , cp , Eid , W ] = distanceFrom( XYZ , AA );

    switch lower( Omode )
      case 'nan'
        Eout = d > THRES;
    end
    
    if ~isempty( FC )
    [ Fid , ~ , dFC  ] = distancePoint2Segments( FC , { AA.xyz , AA.tri } );

    switch lower( Omode )
      case 'nan'
        Fout = dFC > THRES;
    end
    end

  
  
  elseif celltype == 5
    %[ Eid , ~ , d , W ] = vtkClosestElement( AA , XYZ );
    %AA.xyz( : ,end+1:3) = 1;
    [ d , ~ , Eid , W ] = distanceFrom( XYZ , AA , true );
    switch lower( Omode )
      case 'nan'
        Eout = abs( d ) > THRES;
    end
    
    if ~isempty( FC )
    [ Fid , ~ , dFC     ] = vtkClosestElement( AA , FC  );

    switch lower( Omode )
      case 'nan'
        Fout = dFC > THRES;
    end
    end
    
    
  elseif celltype == 10
    [Eid,W] = tsearchn( AA.xyz , AA.tri , XYZ  );
    switch lower( Omode )
      case 'nan'
        Eout = isnan( Eid );
      case 'closest'
        S = MeshBoundary( MeshGenerateIDs( AA ) );
        w = find( isnan( Eid ) );
        [Eid_,~,~,W_] = vtkClosestElement( S , XYZ(w,:) );
        Eid( w ,:) = S.triID( Eid_ );
        
        W( w ,:) = 0;
        for r = 1:numel(Eid_)
          E_ = Eid_(r);
          
          tri   = S.tri( E_ ,:);
          tetra = AA.tri( S.triID(E_) ,:);
          
          for c = 1:3
            W( w(r) , tetra == tri(c) ) = W_(r,c);
          end
        end        
    end
    
    if ~isempty( FC )
    [Fid,W] = tsearchn( AA.xyz , AA.tri , FC  );
    switch lower( Omode )
      case 'nan'
        Fout = isnan( Fid );
      case 'closest'
        S = MeshBoundary( MeshGenerateIDs( AA ) );
        w = isnan( Fid );
        Fid_ = vtkClosestElement( S , FC(w,:) );
        Fid( w ,:) = S.triID( Fid_ );
    end
    end
    
  else
    error('not implemented for this celltype' );
  end

  if outW, WW = W; end
  
  if isempty( FC )
    Fid  = Eid;
    Fout = Eout;
  end
  
  switch lower( Omode )
    case 'nan'
      Eid( Eout ,:) = 1;
      W( Eout ,:)   = 1;
      Fid( Fout ,:) = 1;
  end
  
  for f = fieldnames( A ).', f = f{1};
    if strcmp( f ,'xyz' ), continue; end
    if strcmp( f ,'tri' ), continue; end
    if strncmp( f , 'xyz' , 3)
      if ~isnumeric( A.(f) )
        %non-numeric field: cannot interpolate -> value at the highest-weight
        %(nearest) incident vertex
        [~,b] = max(W,[],2); Wb = sparse( 1:size(W,1) , b , true , size(W,1) , size(W,2) );
        C.(f) = A.(f)( sum( double( A.tri( Eid ,:) ) .* Wb ,2) );
        continue;
      end
      w = isfinite( Eid );
      C.(f) = zeros( numel(Eid) ,size(A.(f),2));
      for c = 1:size( A.tri ,2)
        C.(f)(w ,:) = C.(f)(w ,:) + bsxfun( @times , double( A.(f)( A.tri( Eid(w) ,c) ,:) ) , W(w,c) );
      end
      if any( Eout )
        if ~isfloat( C.(f) ), C.(f) = double( C.(f) ); end
        C.(f)( Eout ,:,:,:,:) = NaN;
      end
    end
    if strncmp( f , 'tri' , 3)
      C.(f) = A.(f)( Fid ,:,:,:,:,:);

      if any( Fout )
        if ~isfloat( C.(f) ), C.(f) = double( C.(f) ); end
        C.(f)( Fout ,:,:,:,:) = NaN;
      end
    end
  end

  if outW, C.xyzQueryWeights = WW; end
  
end

