function [M,C] = MeshRemoveFaces( M , w , varargin )
% - if specified, tidy up the mesh

  if nargout > 1
    M_ = MeshRemoveFaces( M , w , varargin{:} );
    if iscell( w ), w = w{1};
    else,           w = {w};
    end
    C = MeshRemoveFaces( M , w , varargin{:} );
    M = M_;
    return;
  end

  tidyAfter = false;  %remove
  if numel( varargin ) && isscalar( varargin{1} ) && ( isnumeric( varargin{1} ) || islogical( varargin{1} ) )
    tidyAfter = ~~varargin{1};
    varargin(1) = [];
  end
  if numel( varargin ), error('Extra arguments remain in varargin.'); end


  if 0
  elseif isstring( w )
    w = ~cellfun( 'isempty' , regexp( M.triPIECE , char( w ) ,'once' ) );
  elseif iscell(w) && isstring( w{1} )
    w =  cellfun( 'isempty' , regexp( M.triPIECE , char( w{1} ) ,'once' ) );

    %w = ~strcmpi( M.triPIECE , char( w{1} ) );

  elseif isa( w , 'function_handle' )
    
    w = feval( w , M );
  

  elseif isstruct( w ) && isfield( w ,'xyz') && isfield( w ,'tri') && size( M.tri ,2) == size( w.tri ,2)

    %[nid,~,d] = vtkClosestPoint( struct('xyz',double(M.xyz)) , double( w.xyz ) );
    [nid,d] = knnsearch( double(M.xyz) , double( w.xyz ) );
    nid( d > 1e-5 ) = 0;
    nid( end+1:max(w.tri(:)) ) = 0;
    w.tri = nid( w.tri );
    w = ismember( sort( M.tri ,2) , sort( w.tri ,2) ,'rows');

  elseif iscell( w ) && numel(w) == 1 && isstruct( w{1} ) && isfield( w{1} ,'xyz') && isfield( w{1} ,'tri') && size( M.tri ,2) == size( w{1}.tri ,2)

    w = w{1};
    [nid,~,d] = vtkClosestPoint( struct('xyz',double(M.xyz)) , double( w.xyz ) );
    nid( d > 1e-5 ) = 0;
    nid( end+1:max(w.tri(:)) ) = 0;
    w.tri = nid( w.tri );
    w = ismember( sort( M.tri ,2) , sort( w.tri ,2) ,'rows');
    w = {w};

  elseif ischar( w ) && strncmp( w , 'PIECE==' ,7)
    w = w(8:end);
    w = ~builtin( 'cellfun' , 'isempty' , regexpi( M.triPIECE , w ,'once') );

  elseif ischar( w ) && strncmp( w , 'PIECE~=' ,7)
    w = w(8:end);
    w =  builtin( 'cellfun' , 'isempty' , regexpi( M.triPIECE , w ,'once') );

  elseif iscell( w ) && numel(w) == 1 && ischar( w{1} ) && strncmp( w{1} , 'PIECE==' ,7)
    w = w{1};
    w = w(8:end);
    w = ~builtin( 'cellfun' , 'isempty' , regexpi( M.triPIECE , w ,'once') );
    w = {w};

  elseif iscell( w ) && numel(w) == 1 && ischar( w{1} ) && strncmp( w{1} , 'PIECE~=' ,7)
    w = w{1};
    w = w(8:end);
    w =  builtin( 'cellfun' , 'isempty' , regexpi( M.triPIECE , w ,'once') );
    w = {w};

  elseif ischar( w )

    NOT = false;
    w = strrep(w,' ','');
    if w(1) == '~', NOT = true; w(1) = []; end
    
    if ~isfield( M , ['tri',w] )
      error('There is no field %s for faces.',w);
    end
    
    w = M.(['tri',w]);
    if isnumeric(w) && all( ismember([0,1],unique(w)) )
      w = ~~w;
    end
    if NOT, w = ~w; end
    
  elseif isnumeric( w ) && size( w ,2) == size( M.tri ,2) && size( w ,1) > 1
    
    w = ismember( sort( M.tri ,2) , sort( w ,2) ,'rows' );
    
  elseif iscell( w ) && numel(w) == 1 && isnumeric( w{1} ) && size( w{1} ,2) == size( M.tri ,2) && size( w{1} ,1) > 1

    w = w{1};
    w = ismember( sort( M.tri ,2) , sort( w ,2) ,'rows' );
    w = {w};
    
  elseif isfield( M , 'triID' ) && isnumeric( w ) && numel(w) == size(w,1) && all( w < 0 )
    
    w = ismember( M.triID , -w );
    
  elseif iscell( w ) && numel(w) == 1 && isfield( M , 'triID' ) && isnumeric( w{1} ) && numel(w{1}) == size(w{1},1) && all( w{1} < 0 )
    
    w = w{1};
    w = ismember( M.triID , -w );
    w = {w};
    
  end
  
  NN = size( M.tri , 1 );
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
      error('Indices must be smaller than the number of faces.');
    end
    
  else
    
    error('incorrect argument');
    
  end
  
  w = setdiff( 1:NN , w );
  
  Fs = fieldnames( M );
  
  for f = fieldnames( M ).', f=f{1};
    if ~strncmp( f , 'tri' , 3 ), continue; end
    sz = size( M.(f) ); sz(1) = numel(w);
    M.(f) = reshape( M.(f)( w ,:,:,:,:,:,:) , sz );
  end
  if isfield( M , 'celltype' ) && ~isscalar( M.celltype )
    M.celltype = M.celltype( w ,:);
  end
  
  if tidyAfter
    M = MeshTidy( M ,NaN,false);
  end
  
end
