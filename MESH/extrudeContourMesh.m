function M = extrudeContourMesh( C , h , topo )

  if nargin < 2 || isempty( h ), h = 1; end
  if nargin < 3 || isempty( topo ), topo = 'square'; end

  switch lower(topo)
    case {'s','4','square'}, topo = 'square';
    otherwise, error('not implemented yet');
  end

  C0 = C;
  if isstruct( C )
    C = meshSeparate( C );
    C = cellfun( @(c)mesh2contours(c) , C ,'un',0);
  end

  if ~iscell( C ), C = {C}; end
  for c = 1:numel(C)
    [Z,iZ] = getPlane( C{c} ,'+z');
    C{c} = transform( C{c} , iZ );

    isclosed = isequal( C{c}(1,:) , C{c}(end,:) );
    n = size( C{c} ,1);

    if 0
    elseif isclosed && strcmp( topo , 'square' )
      n = n-1;
      C{c}(end,:) = [];
      C{c} = struct( 'xyz' , [ C{c} ; bsxfun( @plus , C{c} , [0,0,h] ) ] );
      C{c}.tri = [ bsxfun( @plus , [ 1 , 2 , n+2 , n+1 ] , (0:n-2).' ) ; n , 1 , n+1 , 2*n ];
      C{c}.celltype = 9;
    else
      error('not implemented yet');
    end

    C{c} = transform( C{c} , Z );
  end
  M = MeshAppend( C );

end
