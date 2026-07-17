function M = xyplaneMesh( Xs , Ys , Zs , mode )

  if nargin < 4
    mode = '/';
  end

  if nargin < 3
    Zs = 0;
  elseif ischar( Zs )
    mode = Zs;
    Zs = 0;
  end
  if numel( Zs ) ~= 1, error('Zs should be a single number. This function only return a plane.'); end

  if iscell( Xs ), Xs = linspace( Xs{:} ); end
  if iscell( Ys ), Ys = linspace( Ys{:} ); end
  
  M.xyz = ndmat( Xs , Ys );
  
  switch lower( mode )
    case {'del','d','delaunay'}
      M.tri = delaunayn( M.xyz );

    case {'r','random'}
      r = 2*( rand( size( M.xyz ) ) - 0.5 )/1000;
      r = reshape( r , [ numel(Xs) , numel(Ys) , 2 ] );
      r([1,end],:,:) = 0;
      r(:,[1,end],:) = 0;
      r = reshape( r ,size(M.xyz));
      M.xyz = M.xyz + r;
      
      M.tri = delaunayn( M.xyz + r );
    
    case {'/'}
      NX = numel( Xs );
      NY = numel( Ys );
      
      M.tri = [ (1:NX-1).' , (2:NX).'      , (NX+2:2*NX).' ;
                (1:NX-1).' , (NX+2:2*NX).' , (NX+1:2*NX-1).' ];
      M.tri = bsxfun( @plus , M.tri , vec( 0:NX:NX*(NY-2) ,3) );
      M.tri = reshape( permute( M.tri ,[2 3 1]) ,3,[]).';
      M.tri = sortrows( M.tri );
      
    case {'\'}
      NX = numel( Xs );
      NY = numel( Ys );
      
      M.tri = [ (1:NX-1).' , (2:NX).'      , (NX+1:2*NX-1).' ;
                (2:NX).'   , (NX+2:2*NX).' , (NX+1:2*NX-1).' ];
      M.tri = bsxfun( @plus , M.tri , vec( 0:NX:NX*(NY-2) ,3) );
      M.tri = reshape( permute( M.tri ,[2 3 1]) ,3,[]).';
      M.tri = sortrows( M.tri );
    
    case {'x'}
      error('not implemented yet.');
    
    otherwise
      error('invalid mode');
  end
  
  if isa( Zs , 'function_handle' )
    try
      M.xyz(:,3) = Zs( M.xyz(:,1) , M.xyz(:,2) );
    catch
      for z = 1:size( M.xyz ,1)
        M.xyz(z,3) = Zs( M.xyz(z,1) , M.xyz(z,2) );
      end
    end
  elseif iscell( Zs ) && isequal( size( Zs{1} ) , [4 4] )
    M.xyz(:,3) = 0;
    M.xyz = transform( M.xyz , Zs{1} );
  else
    M.xyz(:,end+1) = Zs;
  end

end
