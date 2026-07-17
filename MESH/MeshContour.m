function C = MeshContour( M , V , levels )

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
  V = double(V);
  levels = double( levels );

  mm = min( V(isfinite(V(:))) );
  MM = max( V(isfinite(V(:))) );
  
  C = MeshZeroContour( M , V*0+1 );
  for l = sort( unique( levels(:) ).' )
    if l < mm, continue; end
    if l > MM, continue; end
    CC = MeshZeroContour( M , V - l );
    if isempty( CC.tri ), continue; end
    C = MeshAppend( C , CC );
  end

end
