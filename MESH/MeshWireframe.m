function M = MeshWireframe( M , mode )

  if nargin < 2, mode = 's'; end
  if isnumeric( mode ) && isscalar( mode )
    switch mode
      case 3, mode = 't';
      case 2, mode = 's';
      otherwise, error('Invalid mode. Only 2 or 3 are allowed.');
    end        
  end
  if ~ischar( mode )
    error('Invalid mode.');
  end
  switch lower(mode)
    case {'s','seg','segment','segments'}
    case {'t','tri','triangle','triangles'}
    otherwise, error('Invalid mode. Only ''Segment'' or ''Triangle'' are allowed.');
  end

  M.tri = meshEdges( M.tri );

  for f = fieldnames( M ).'
    if strcmp( f{1} , 'tri' ), continue; end
    if ~strncmp( f{1} , 'tri' , 3 ), continue; end
    M = rmfield( M , f{1} );
  end  
  
  
  switch lower(mode)
    case {'s','seg','segment','segments'}
      M.celltype = 3;
    case {'t','tri','triangle','triangles'}
      M.celltype = 5;
      nV = size( M.xyz ,1);

      M.xyz = [ M.xyz ; meshFacesCenter( M ) ];
      M.tri(:,3) = nV + ( 1:size(M.tri,1) );
      for f = fieldnames( M ).'
        if strncmp( f{1} , 'xyz' , 3 )
          M.(f{1})( (end+1):size( M.xyz ,1) ,:,:,:,:,:,:) = NaN;
        end
      end  
    otherwise, error('Invalid mode. Only ''Segment'' or ''Triangle'' are allowed.');
  end


end
