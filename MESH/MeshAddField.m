function M = MeshAddField( M , varargin )

  nF = 0; try, nF = size( M.tri ,1); end
  nV = 0; try, nV = size( M.xyz ,1); end


  while ~isempty( varargin )
    Fname = varargin{1}; varargin(1) = [];
    val   = varargin{1}; varargin(1) = [];

    if Fname(1) == '-'
      Fname(1) = [];
      if isfield( M , Fname ), continue; end
    end

    if isa( val , 'function_handle' )
      val = val( M );
    end

    if ischar( val )
      if size( val ,1) > 1, error('a string was expected'); end
      val = { val };
    end
  
    
    %if isvector( val ), val = val(:); end
    n  = size( val , 1 );
    
    parent = '';
    if false
    elseif strncmp( Fname , 'xyz' , 3 )  &&  iscell( val ) && numel( val ) == 1 && isnumeric( val{1} ) && isvector( val{1} ) && numel( val{1} ) > 1 && all( val{1} > 0 ) && all( val{1} <= nV )
      v = zeros( nF ,1); v( val{1} ) = 1; val = v; n = nF;
      parent = 'tri'; Fname = Fname(4:end);
    elseif strncmp( Fname , 'tri' , 3 )  &&  iscell( val ) && numel( val ) == 1 && isnumeric( val{1} ) && isvector( val{1} ) && numel( val{1} ) > 1 && all( val{1} > 0 ) && all( val{1} <= nF )
      v = zeros( nF ,1); v( val{1} ) = 1; val = v; n = nF;
      parent = 'tri'; Fname = Fname(4:end);
    elseif strncmp( Fname , 'xyz' , 3 )  &&  ( nV == n || n == 1 )
      parent = 'xyz'; Fname = Fname(4:end);
    elseif strncmp( Fname , 'tri' , 3 )  &&  ( nF == n || n == 1 )
      parent = 'tri'; Fname = Fname(4:end);
    elseif n == nF  &&  n ~= nV
      parent = 'tri';
    elseif n ~= nF  &&  n == nV
      parent = 'xyz';
    elseif n == nF  &&  n == nV
      error('cannot determine if it is a field for VERTICES or for FACES.');
    elseif n ~= nF  &&  n ~= nV
      error('It cannot be a field for VERTICES or for FACES.');
    end
    
    if strcmp( parent , 'xyz' )
      if n == 1
        val = repmat( val , [ nV , 1 ] );
      end
      if size( val ,1) ~= nV
        error('It cannot be a field for VERTICES.');
      end
    end
    if strcmp( parent , 'tri' )
      if n == 1
        val = repmat( val , [ nF , 1 ] );
      end
      if size( val ,1) ~= nF
        error('It cannot be a field for FACES.');
      end
    end
      
    M.([ parent , Fname ]) = val;

  end

end
