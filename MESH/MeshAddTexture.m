function S = MeshAddTexture( M , varargin )

if 0
  M = read_OBJ( 'M.obj' );
  
  tic; Tn = MeshAddTexture( M , M.texture ,'n'); toc
  figure; plotMESH( Tn ,'EdgeColor','none','rgb','RGB','ac');
  
  tic; Ti = MeshAddTexture( M , M.texture ,'i'); toc
  figure; plotMESH( Ti ,'EdgeColor','none','rgb','RGB','ac');
  
  tic; Td = MeshAddTexture( M , M.texture ,'d'); toc
  figure; plotMESH( Td ,'EdgeColor','none','rgb','RGB','ac');
  
  tic; Ts = MeshAddTexture( M , M.texture ,'s'); toc
  figure; plotMESH( Ts ,'EdgeColor','none','rgb','RGB','ac');
end

  UV = [];
  T  = [];
  PLOT = false;
  mode = 'subdivide';
  FACTOR = 1;
  minEL = 0;


  try, UV = M.xyzUV; end
  try, T  = M.texture; end
  for v = varargin(:).', v = v{1};
    if ischar( v )
      switch lower(v)
        case {'n','nearest'},       mode = 'nearest';
        case {'i','interpolated'},  mode = 'interpolated';
        case {'d','delaunay'},      mode = 'delaunay';
        case {'s','subdivide'},     mode = 'subdivide';
        case {'plot'},              PLOT = true;
        otherwise, error('Unknown texturization mode.');
      end
      continue;
    end
    if isnumeric( v ) && ismatrix( v ) && size( v ,2) == 2 && size( v ,1) == size( M.xyz ,1)
      UV = v;
      continue;
    end
    if isnumeric( v ) && numel( v ) > 1
      T = v;
      continue;
    end
    if numel( v ) == 1 && ( islogical( v ) || v == 0 || v == 1 )
      PLOT = ~~v;
      continue;
    end
    if numel( v ) == 1 && v > 0 && v < 1
      FACTOR = v;
      continue;
    end
    if numel( v ) == 1 && v < 0
      minEL = -v;
      continue;
    end
    error('invalid argument.');
  end
  if isempty( UV ), error('UV coordinates must be provided.'); end
  if isempty( T  ), error('An image for texturizing must be provided.'); end
  if FACTOR ~= 1
    T = imresize( T , FACTOR );
  end

  U = bsxfun( @times , [ UV(:,1) , 1-UV(:,2) ] , [size(T,2),size(T,1)] ) + 0.5;
  
  switch mode
    case 'subdivide'

      S = Mesh( M );
      S.xyzXYZ = S.xyz;
      S.xyz = U;

      d = max( minEL , 1.5 );
      for it = 1:6
        if size( S.tri ,1) > 1e6, break; end
        w = meshQuality( S , 'maxl' ) > d; if ~any(w), break; end
        S = MeshSubdivide( S , w );
      end
    
      UV = S.xyz; S.xyz = S.xyzXYZ; S = rmfield( S , 'xyzXYZ');
%       if isfield( S ,'xyzUV')
%         S.xyzUV = UV;
%         S.xyzUV = S.xyzUV - 0.5;
%         S.xyzUV = bsxfun( @rdivide , S.xyzUV , [size(T,2),size(T,1)] );
%         S.xyzUV(:,2) = 1 - S.xyzUV(:,2);
%       end
    
      I = griddedInterpolant( { 1:size(T,2) , 1:size(T,1) } , double( permute(T(:,:,:),[2,1,3]) ) ,'linear','nearest');
      S.xyzRGB = I( UV );
      S.xyzRGB = reshape( S.xyzRGB , size( S.xyzRGB ,1) , [] );
      switch class( T )
        case 'uint8', S.xyzRGB = S.xyzRGB / 255;
      end

    case 'nearest'

      U = struct('xyz',U,'tri',double(M.tri));
      S = flip( meshSeparate( U ,'order',@meshSurface) ,1);

      state = warning( 'off' , 'MATLAB:delaunayTriangulation:ConsConsSplitWarnId' );
      onCLEAN = onCleanup( @()warning(state) );
      for s = 1:numel(S), %disp(s);
        V = S{s};
        bb = meshBB( V );
        x = ( floor( bb(1,1) ) - 0.5 ):( ceil( bb(2,1) ) + 0.5 );
        y = ( floor( bb(1,2) ) - 0.5 ):( ceil( bb(2,2) ) + 0.5 );
        
        Q = fastAppend( { Mesh( ndmat( x([1,end]) , y          ) , reshape( 1:2*numel(y) ,2,[]).' ) ,...
                          Mesh( ndmat( x          , y([1,end]) ) , reshape( 1:2*numel(x) ,[],2)   ) } );
        Q = fastAppend( { Q , MeshWireframe( V ) } );
        Q = MeshTidy( Q ,0,true);

      
        Q = delaunayTriangulation( Q.xyz(:,1) , Q.xyz(:,2) , Q.tri );
      
        Q = Mesh( Q );
        try
        CH = Mesh( MeshTidy( struct( 'xyz',V.xyz,'tri', convhulln( V.xyz ) ) ) ,'delaunay');
        Q = MeshRemoveFaces( Q , isnan( tsearchn( CH.xyz , CH.tri , meshFacesCenter( Q ) ) ) );
        end
        Q = MeshRemoveFaces( Q , isnan( tsearchn(  V.xyz ,  V.tri , meshFacesCenter( Q ) ) ) );
        S{s} = Q;
      end
      clearvars( 'onCLEAN' );
      
      S = fastAppend( S );
      I = griddedInterpolant( { 1:size(T,2) , 1:size(T,1) } , double( permute(T(:,:,:),[2,1,3]) ) ,'nearest','nearest');
      S.triRGB = I( meshFacesCenter( S ) );
      S.triRGB = reshape( S.triRGB , size( S.triRGB ,1) , [] );
      S.triRGB = cast( S.triRGB , class( T ) );
      S.xyz = meshMapPoints( S.xyz , U , M );

    case 'interpolated'

      U = struct('xyz',U,'tri',double(M.tri));
      S = flip( meshSeparate( U ,'order',@meshSurface) ,1);
      state = warning( 'off' , 'MATLAB:delaunayTriangulation:ConsConsSplitWarnId' );
      onCLEAN = onCleanup( @()warning(state) );
      for s = 1:numel(S), %disp(s);
        V = S{s};
        bb = meshBB( V );
        x = floor( bb(1,1) ):ceil( bb(2,1) );
        y = floor( bb(1,2) ):ceil( bb(2,2) );
        
        Q = fastAppend( { Mesh( ndmat( x([1,end]) , y          ) , reshape( 1:2*numel(y) ,2,[]).' ) ,...
                          Mesh( ndmat( x          , y([1,end]) ) , reshape( 1:2*numel(x) ,[],2)   ) } );
        Q = fastAppend( { Q , MeshWireframe( V ) } );
        Q = MeshTidy( Q ,0,true);
      
        Q = delaunayTriangulation( Q.xyz(:,1) , Q.xyz(:,2) , Q.tri );
      
        Q = Mesh( Q );
        CH = Mesh( MeshTidy( struct( 'xyz',V.xyz,'tri', convhulln( V.xyz ) ) ) ,'delaunay');
        Q = MeshRemoveFaces( Q , isnan( tsearchn( CH.xyz , CH.tri , meshFacesCenter( Q ) ) ) );
        Q = MeshRemoveFaces( Q , isnan( tsearchn(  V.xyz ,  V.tri , meshFacesCenter( Q ) ) ) );
        S{s} = Q;
      end
      clearvars( 'onCLEAN' );

      S = fastAppend( S );
      I = griddedInterpolant( { 1:size(T,2) , 1:size(T,1) } , double( permute(T(:,:,:),[2,1,3]) ) ,'linear','nearest');
      S.xyzRGB = I( S.xyz );
      S.xyzRGB = reshape( S.xyzRGB , size( S.xyzRGB ,1) , [] );
      switch class( T )
        case 'uint8', S.xyzRGB = S.xyzRGB / 255;
      end
      S.xyz = meshMapPoints( S.xyz , U , M );

    case 'delaunay'

      U = struct('xyz',U,'tri',double(M.tri));
      S = flip( meshSeparate( U ,'order',@meshSurface) ,1);
      for s = 1:numel(S), %disp(s);
        V = S{s};
        bb = meshBB( V );
        x = floor( bb(1,1) ):ceil( bb(2,1) );
        y = floor( bb(1,2) ):ceil( bb(2,2) );

        xy = ndmat( x , y );
        xy( isnan( tsearchn( V.xyz , V.tri , xy ) ) ,:) = [];

        Q = MeshWireframe( V );
        while 1
          w = meshQuality( Q ,'length') > 1.5;
          if ~any(w), break; end
          Q = MeshSubdivide( Q ,w);
        end
        xy( ismember( xy , Q.xyz ,'rows' ) ,:) = [];
        Q.xyz = [ Q.xyz ; xy ];
        
        state = warning( 'off' , 'MATLAB:delaunayTriangulation:ConsConsSplitWarnId' );
        onCLEAN = onCleanup( @()warning(state) );
        Q = delaunayTriangulation( Q.xyz(:,1) , Q.xyz(:,2) , Q.tri );
        clearvars( 'onCLEAN' );
      
        Q = Mesh( Q );
        Q = MeshRemoveFaces( Q , isnan( tsearchn( V.xyz , V.tri , meshFacesCenter( Q ) ) ) );
        S{s} = Q;
      end
      S = fastAppend( S );
      I = griddedInterpolant( { 1:size(T,2) , 1:size(T,1) } , double( permute(T(:,:,:),[2,1,3]) ) ,'linear','nearest');
      S.xyzRGB = I( S.xyz );
      S.xyzRGB = reshape( S.xyzRGB , size( S.xyzRGB ,1) , [] );
      switch class( T )
        case 'uint8', S.xyzRGB = S.xyzRGB / 255;
      end
      S.xyz = meshMapPoints( S.xyz , U , M );

    otherwise, error('Invalid mode.');
  end


  if PLOT
    S = plotMESH( S ,'FaceColor','flat','CData',reshape( S.triRGB ,[],1,3),'ne'); %headlight
  end

end


function S = fastAppend( S )

  offsets = cumsum( [ 0 ; cellfun( @(M)size(M.xyz,1) ,S(:) ) ] );
  
  xyz = cell2mat( cellfun(@(M)M.xyz , S(:),'un',0) );
  tri = cell2mat( arrayfun( @(s)S{s}.tri+offsets(s) ,1:numel(S) ,'un',0).' );

  S = struct( 'xyz' , xyz , 'tri' , tri );

end