function [ p , d_tp ] = transform( p , varargin )

  if isempty(p)
    if nargout>1, d_tp = []; end
    return
  end



  if iscell( p )
    for c = 1:numel(p)
      p{c} = transform( p{c} ,varargin{:} );
    end
    return;
  end

%   if isstruct( p ) && isfield( p ,'xyz' )
%     p.xyz = transform( p.xyz ,varargin{:} );
%     if isfield( p ,'lmk' )
%       p.lmk = transform( p.lmk ,varargin{:} );
%     end
%     if isfield( p ,'bspmCoords' )
%       p.bspmCoords = transform( p.bspmCoords ,varargin{:} );
%     end
%     return;
%   end

  if isstruct( p ) && isfield( p , 'xyz' )
                                        p.xyz           = transform( p.xyz           , varargin{:} );
    if isfield( p , 'lmk' ),            p.lmk           = transform( p.lmk           , varargin{:} );  end
    if isfield( p , 'bspmCoords' ),     p.bspmCoords    = transform( p.bspmCoords    , varargin{:} );  end
    if isfield( p , 'markersCoords' ),  p.markersCoords = transform( p.markersCoords , varargin{:} );  end
    if isfield( p , 'xyzNORMALS' )
      R = maketransform( varargin{:} );
      R = R(1:3,1:3);
      R = R / cbrt( det(R) );
      p.xyzNORMALS = p.xyzNORMALS * R.';
    end
    if isfield( p , 'triNORMALS' )
      R = maketransform( varargin{:} );
      R = R(1:3,1:3);
      R = R / cbrt( det(R) );
      p.triNORMALS = p.triNORMALS * R.';
    end
    return;
  end


  if numel( varargin ) == 1 && isa( varargin{1} ,'function_handle')
    if isnumeric( p )
      p = feval( varargin{1} , p );
    end
    if isstruct( p ) && isfield( p , 'xyz' )
      p.xyz = feval( varargin{1} , p.xyz );
      if isfield( p , 'lmk' )
        p.lmk = feval( varargin{1} , p.lmk );
      end
      if isfield( p , 'bspmCoords' )
        p.bspmCoords = feval( varargin{1} , p.bspmCoords );
      end
      if isfield( p , 'markersCoords' )
        p.markersCoords = feval( varargin{1} , p.markersCoords );
      end
    end
    return;
  end

  %"procrustes type" transformation for points
  if nargout < 2 && numel( varargin ) == 2 && ischar( varargin{1} ) && ...
     isnumeric( varargin{2} ) && size( varargin{2} ,2) == size( p ,2) && size( varargin{2} ,1) > 3 && size( varargin{2} ,1) == size( p ,1) 

    R = MatchPoints( varargin{2} , p , varargin{1} );
    p = transform( p , R );
    return;
  end




  %fast classical 3d points transform by homogeneous matrix
  if numel(varargin) == 1 && nargout < 2 && isnumeric( varargin{1} ) && ( isequal( size( varargin{1} ),[4 4] )  ||  isequal( size( varargin{1} ),[3 3] ) ) && ...
      all( varargin{1}( end , 1:end-1 ) == 0 )  &&  varargin{1}(end,end) == 1
    T = varargin{1};
    
    if isa(p,'struct') && isfield( p , 'xyz' )
%       p.xyz(:,(end+1):3) = 0;
      p.xyz = bsxfun( @plus , p.xyz * T(1:end-1,1:end-1).' , T(1:end-1,end).' );
      if isfield( p , 'lmk' )
%         p.lmk(:,(end+1):3) = 0;
        p.lmk = bsxfun( @plus , p.lmk * T(1:end-1,1:end-1).' , T(1:end-1,end).' );
      end
      if isfield( p , 'bspmCoords' )
        p.bspmCoords = bsxfun( @plus , p.bspmCoords * T(1:end-1,1:end-1).' , T(1:end-1,end).' );
      end
      if isfield( p , 'markersCoords' )
        p.markersCoords = bsxfun( @plus , p.markersCoords * T(1:end-1,1:end-1).' , T(1:end-1,end).' );
      end
      
% xyzSSM_MEAN, xyzSSM_MODES      
      if isfield( p , 'xyzSSM_MEAN' )
        p.xyzSSM_MEAN = bsxfun( @plus , p.xyzSSM_MEAN * T(1:end-1,1:end-1).' , T(1:end-1,end).' );
      end
      if isfield( p , 'xyzSSM_MODES' )
        for m = 1:size( p.xyzSSM_MODES ,3)
          p.xyzSSM_MODES(:,:,m) = p.xyzSSM_MODES(:,:,m) * T(1:end-1,1:end-1).';
        end
      end
      return
      
    elseif isa(p,'struct') && isfield( p , 'vertices' )
%       p.vertices(:,(end+1):3) = 0;
      p.vertices = bsxfun( @plus , p.vertices * T(1:end-1,1:end-1).' , T(1:end-1,end).' );
      return

    elseif isa( p ,'double') && isequal( size(p) , [4 4] ) && isequal( size(T) , [4 4] )
      p = T * p;
      return;
      
    elseif isa( p ,'double') && isequal( size(p) , [3 3] ) && isequal( size(T) , [3 3] )
      p = T * p;
      return;

    elseif isa(p,'numeric') && size(p,2) == 3  && isequal( size(T) , [4 4] ) %&& size(p,1) > 3
      p = bsxfun( @plus , p * T(1:3,1:3).' , T(1:3,4).' );
      return;
      
    elseif isa(p,'numeric') && size(p,2) == 2  && isequal( size(T) , [3 3] ) %&& size(p,1) > 3
      p = bsxfun( @plus , p * T(1:2,1:2).' , T(1:2,3).' );
      return;

    elseif isa(p,'cell')
      for pp = 1:numel(p)
        p{pp} = transform( p{pp} , T );
      end
      return;
      
    end
    
  end

  if isa( p ,'double') && isequal( size(p) , [4 4] ) && isequal( p(4,:) , [0 0 0 1] )
    T = maketransform( varargin{:} );
    p = T * p;
    return;
  end
  
  
  if isa(p,'CubicMeshClass')
    T = maketransform( varargin{:} );
    
    if maxnorm( T(1:3,1:3)*T(1:3,1:3).' - eye(3) ) > 1e-8  ||...
       maxnorm( T(1:3,1:3).'*T(1:3,1:3) - eye(3) ) > 1e-8  ||...
       abs( det( T(1:3,1:3) ) - 1 ) > 1e-8  ||...
       maxnorm( T(4,:) - [0 0 0 1] ) > 1e-8
     error('the transformation should be a rigid body transform in homogeneous coordinates.');
    end
    
    p = p.rotateMesh( T(1:3,1:3) , [0;0;0] );
    p = p.translateMesh( T(1:3,4) );
  
    return;
  end
  
  if isa(p,'struct') && isfield( p , 'xyz' )
    p.xyz(:,(end+1):3) = 0;
    if nargin == 2 && isnumeric( varargin{1} )
      T = varargin{1};
      p.xyz = bsxfun( @plus , p.xyz * T(1:3,1:3).' , T(1:3,4).' );
    else
      p.xyz = transform( p.xyz , varargin{:} );
    end
    return;
  end

  if isa(p,'struct') && isfield( p , 'SpatialTransform' )
    
    R = p.SpatialTransform;
    R = maketransform( varargin{:} ) * R;
    p.SpatialTransform = R;
    
    return;
  end
  

  
  if isa(p,'struct') && isfield( p , 'Format' ) && isequal( p.Format , 'DICOM' )
    if ~isfield( p , 'ImageOrientationPatient' ) || ~isfield( p , 'ImagePositionPatient' )
      error('not proper dicom.');
    end
    
    if isfield( p , 'xSpatialTransform' )
      R = p.xSpatialTransform;
    else
      R = reshape( p.ImageOrientationPatient , 3 , 2 );
      R(:,3)= cross( R(:,1), R(:,2) );
      for cc = 1:3, for it = 1:5, R(:,cc) = R(:,cc)/sqrt( R(:,cc).' * R(:,cc) ); end; end

      R = [ R , p.ImagePositionPatient(:) ; 0 0 0 1 ];
    end
    
    R = maketransform( varargin{:} ) * R;
    if isfield( p , 'xSpatialTransform' )
      p.xSpatialTransform = R;
    end
    
    p.ImagePositionPatient(:)    = R(1:3,4  );
    p.ImageOrientationPatient(:) = R(1:3,1:2);
    
    return;
  end
  
  if isa(p,'cell') && numel(varargin)==1 && iscell(varargin{1})
    Ms = varargin{1};

    szP = size(p );
    szM = size(Ms);
    szP((end+1):max(numel(szP),numel(szM))) = 1;
    szM((end+1):max(numel(szP),numel(szM))) = 1;
    r = szP./szM;
    if any( mod(r,1) ), error('Non-singleton dimensions of M must match each of P.'); end 
    if any( szM( r ~= 1 ) ~= 1 ), error('only singletons can be expanded.'); end
    
    Ms = repmat( Ms , [r 1] );
    
    for m = 1:numel(Ms)
      p{m} = transform( p{m} , Ms{m} );
    end
    
    return;
  end
  
  if isa(p,'cell')
    M = maketransform( varargin{:} );
    p = cellfun( @(c)transform( c , M ) , p , 'UniformOutput' , false );
    return;
  end
  
  if size( p ,1) > 20 && size( p ,2) > 20 && nargout < 2
    sz_in  = size(p);
    sz_out = sz_in;

    OPTS = {};
    if iscell( varargin{end} )
      OPTS = varargin{end}; varargin(end) = [];
    end
    
    to_remove = [];
    for o = 1:numel(OPTS)
      if isnumeric( OPTS{o} ) && numel( OPTS{o} ) == 2
        to_remove = [ to_remove , o ];
        sz_out = OPTS{o};
      end
    end
    OPTS( to_remove ) = [];

    
    M = maketransform( varargin{:} );
    M = minv( M );
    
    xy = ndmat( 0:sz_out(2)-1 , 0:sz_out(1)-1 , 0 , 1 ) * M.';
    xy = bsxfun( @rdivide , xy(:,1:3) , xy(:,4) );
    
    p = permute( p ,[2,1,4,3,5:10] );
    if islogical( p )
      p = ~~InterpPointsOn3DGrid( p , 0:sz_in(2)-1 , 0:sz_in(1)-1 , 0 , xy ,'*nearest','closest',OPTS{:});
    else
      p =   InterpPointsOn3DGrid( p , 0:sz_in(2)-1 , 0:sz_in(1)-1 , 0 , xy ,'*linear','closest',OPTS{:});
    end
    p = ipermute( reshape( p , [ sz_out([2,1]) , size(p,2) ] ) ,[2,1,3]);
    
    return;
  end
    
  if isa(p,'struct')
    for f = fieldnames(p).', f = f{1};
      try, p.(f) = transform( p.(f) , varargin{:} ); end
    end
    return;
  end
  

%   try
%     T = maketransform( varargin{:} );
%     p = transform( p , T );
%     return;
%   end
  

  [ varargin,i,formatpoints ]= parseargs( varargin,'Format','formatpoints','pointsformat','fp','pf','$DEFS$',false );
  
  if ~formatpoints, [ varargin,formatpoints ]= parseargs(varargin,'Rows'                        ,'$FORCE$','rows'       ); end
  if ~formatpoints, [ varargin,formatpoints ]= parseargs(varargin,'RowsH'                       ,'$FORCE$','rowsh'      ); end
  if ~formatpoints, [ varargin,formatpoints ]= parseargs(varargin,'r2d','rows2d'                ,'$FORCE$','rows2d'     ); end
  if ~formatpoints, [ varargin,formatpoints ]= parseargs(varargin,'r2dh','rows2dh'              ,'$FORCE$','rows2dh'    ); end
  if ~formatpoints, [ varargin,formatpoints ]= parseargs(varargin,'col','Columns'               ,'$FORCE$','columns'    ); end
  if ~formatpoints, [ varargin,formatpoints ]= parseargs(varargin,'colsh','ColumnsH'            ,'$FORCE$','columnsh'   ); end
  if ~formatpoints, [ varargin,formatpoints ]= parseargs(varargin,'c2d','col2d','columns2d'     ,'$FORCE$','columns2d'  ); end
  if ~formatpoints, [ varargin,formatpoints ]= parseargs(varargin,'c2dh','cols2dh','columns2dh' ,'$FORCE$','columns2dh' ); end

  if ~formatpoints, formatpoints= 'rows'; end
  

  if isequal( formatpoints  , 'rows2d' ) && numel( varargin ) == 1 && nargout < 2
    if size(p,2) ~= 2, error('an nx2 matrix was expected'); end
    H = varargin{1}; if ~isequal( size( H ) , [3,3] ), error('a 3x3 matrix was expected'); end
    p(:,3) = 1;
    p = p * H.';
    p = bsxfun( @rdivide , p(:,1:2) , p(:,3) );
    return;
  end
  
  
  
%   dH = zeros(16,0);
  dH = eye(16);
  
  if numel( varargin ) == 1 && isnumeric( varargin{1} )
    
    H = varargin{1};
    if iscell( H )
      dH = H{2};
      H  = H{1};
    end

  else

    if nargout < 2
       H     = maketransform( varargin{:} );
    else
      [H,dH] = maketransform( varargin{:} );
    end

  end

  if isequal( H , eye(4) )   &&  nargout < 2
    return;
  end

  
  if size(H,1) < 4 || size(H,2) < 4
    H(4,4)= 1;
  end

  
  if strcmpi( formatpoints , 'rows' )
    if isequal( H(4,:) , [0 0 0 1] )

      if nargout > 1
%         d_tp = kron( [ p , ones(size(p,1),1) ] , eye(3,3) )*dH([1 2 3 5 6 7 9 10 11 13 14 15],:);
%         d_tp = d_tp( [ 1:3:end 2:3:end 3:3:end ] , : );
        
        d_tp = [ bsxfun( @plus , p * dH([ 1  5  9],:) , dH(13,:) ) ; ...
                 bsxfun( @plus , p * dH([ 2  6 10],:) , dH(14,:) ) ; ...
                 bsxfun( @plus , p * dH([ 3  7 11],:) , dH(15,:) ) ];
        
      end
      
      if ~isequal( H(1:3,1:3) , eye(3) )
        p = p*H(1:3,1:3).';
      end

      if ~isequal( H(1:3,4) , [0;0;0] )
        p = bsxfun( @plus , p , H(1:3,4).' );
      end

      return;

    end
  end
  
  
  if nargout > 1
    error('falta por implementar');
  end
  
  switch lower(formatpoints)
    case {'r','rows'}
      p= p(:,[1 2 3]).';
    case {'rh','rowsh'};
      p= p(:,[1 2 3 4]).';
      p= bsxfun( @rdivide , p ,  p(4,:) );
    case {'r2d','rows2d'};
      p= p(:,[1 2]).';
    case {'r2dh','rows2dh'};
      p([1 2 4],:)= p(:,[1 2 3]).';
      p= tp./repmat( p(4,:),4,1 );

    case {'c','col','columns'}
      p= p([1 2 3],:);
    case {'ch','colsh','columnsh'};
      p= p([1 2 3 4],:);
      p= tp./repmat( p(4,:),4,1 );
    case {'c2d','col2d','columns2d'};
      p= p([1 2],:);
    case {'c2dh','cols2dh','columns2dh'};
      p([1 2 4],:)= p([1 2 3],:);
      p= tp./repmat( p(4,:),4,1 );
  end
  p(4,:)=1;

  p= H*p;
  p= bsxfun( @rdivide , p ,  p(4,:) );
  p(4,:)=1;


  switch lower(formatpoints)
    case {'r','rows'}
      p= p([1 2 3],:).';
    case {'rh','rowsh'};
      p= p.';
    case {'r2d','rows2d'};
      try
        if abs( p(3,:) ) < 1e-10
          p= p([1 2],:).';
        else
          p= p([1 2 3],:).';
        end
      catch
          p= p([1 2],:).';
      end
    case {'r2dh','rows2dh'};
      if abs( p(3,:) ) < 1e-10
        p= p([1 2 4],:).';
      else
        p= p.';
      end

    case {'c','col','columns'}
      p= p([1 2 3],:);
    case {'ch','colsh','columnsh'};

    case {'c2d','col2d','columns2d'};
      if abs( p(3,:) ) < 1e-10
        p= p([1 2],:);
      else
        p= p([1 2 3],:);
      end
    case {'c2dh','cols2dh','columns2dh'};
      if abs( p(3,:) ) < 1e-10
        p= p([1 2 4],:);
      else
        
      end
  end

end
