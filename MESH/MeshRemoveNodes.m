function M = MeshRemoveNodes( M , w )

  if isa( w , 'function_handle' )

    try
      w = feval( w , M );
    catch
      try
        w = feval( w , M.xyz );
      catch
        error('invalid function to be evaluated on Mesh or on nodes coordinates.');
      end
    end
    
  elseif ischar( w )

    NOT = false;
    w = strrep(w,' ','');
    if w(1) == '~', NOT = true; w(1) = []; end

    if ~isfield( M , ['xyz',w] )
      error('There is no field %s for nodes.',w);
    end
    
    w = M.(['xyz',w]);
    if NOT, w = ~w; end

  elseif isnumeric( w ) && isfloat( w ) && size( w ,2) == size( M.xyz ,2) && size( w ,1) > 1
    
    w = ismember( M.xyz , w , 'rows' );
    
  elseif iscell( w ) && numel(w) == 1 && isnumeric( w{1} ) && isfloat( w{1} ) && size( w{1} ,2) == size( M.xyz ,2) && size( w{1} ,1) > 1

    w = w{1};
    w = ismember( M.xyz , w , 'rows' );
    w = {w};
    
  elseif isfield( M , 'xyzID' ) && isnumeric( w ) && numel(w) == size(w,1) && all( w < 0 )
    
    w = ismember( M.xyzID , -w );
  
  end


  NN = size( M.xyz , 1 );
  
  if islogical( w )
    
    if ~isvector( w ) || numel( w ) ~= NN
      error('Incorrect logical indexing');
    end
    if ~any( w ), return; end
    w = find(w);

  elseif iscell( w ) && isnumeric( w{1} )
    
    w = setdiff( 1:NN , w{1} );
    
  elseif iscell( w ) && islogical( w{1} )
    
    w = setdiff( 1:NN , find( w{1} ) );
    
  elseif isnumeric( w )
  
    if isempty( w ), return; end
    if any( w < 0 ) || any( mod( w , 1 ) )
      error('Indices must either be real positive integers.');
    end

    if max( w ) > NN
      error('Indices must be smaller than the number of nodes.');
    end
    
  else
    
    error('incorrect argument');
    
  end
  
  w = setdiff( 1:NN , w );
    
  Fs = fieldnames( M );
  
  for f = 1:numel( Fs )
    if ~strncmp( Fs{f} , 'xyz' , 3 ), continue; end
    M.(Fs{f}) = M.(Fs{f})( w ,:,:,:,:,: );
  end
  
  map = zeros( NN , 1 );
  map( w ) = 1:numel(w);
  
  if isfield( M , 'tri' )
    w = all( ismember( M.tri , w ) , 2 );
    for f = 1:numel( Fs )
      if ~strncmp( Fs{f} , 'tri' , 3 ), continue; end
      M.(Fs{f}) = M.(Fs{f})( w ,:,:,:,:,: );
    end
      
    try
      M.tri = feval( class( M.tri ) , reshape( map( M.tri ) ,size(M.tri) ) );
    catch
      M.tri = reshape( map( M.tri ) ,size(M.tri) );
    end  
  end
  
end



