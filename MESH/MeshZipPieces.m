function M = MeshZipPieces( MA , MB , varargin )

  if nargin < 2, MB = []; end
  if ischar( MB )
    if strcmpi( MB , 'relax' )
      
      isB = meshBoundaryNodes( MA );
      M = MeshZipPieces( MA );
      isB( meshBoundaryNodes( M ) ) = false;
      M.xyz( isB ,:) = NaN;
      M = MeshRelax( M );
      return;
    
    elseif strcmpi( MB , 'split' )

      M = MeshZipPieces( MA , [] , 'split' );
      return;
    
    elseif strcmpi( MB , 'fair' )

      M = MeshZipPieces( MA , [] , 'fair' );
      return;
    
    end
  end

  SPLIT = false;
  FAIR  = false;
  while numel( varargin )
    key = varargin{1}; varargin(1) = [];
    if strcmpi( key , 'split' ), SPLIT = true; end
    if strcmpi( key , 'fair' ),  FAIR  = true; end
  end

  PIECES = 0;
  try, PIECES = [ PIECES ; unique( MA.triPIECE(:) ) ]; end
  try, PIECES = [ PIECES ; unique( MB.triPIECE(:) ) ]; end

  A = MA;
  if isa( A , 'struct' )
    ct = meshCelltype( A );
    if ~isscalar( ct ), error('Mesh M should be homogeneous'); end
    if ct == 5
      A = MeshBoundary( A );
      ct = meshCelltype( A );
    end
    if ct ~= 3, error('Not a set of lines'); end
    if ~isfield( A , 'tri' ) || isempty( A.tri )
      error('Not a mesh');
    end
    A = mesh2contours( A );
    if ~isfield( MA , 'triPIECE' ), MA = MeshAddField( MA ,'triPIECE', max(PIECES)+1 ); PIECES = [ PIECES ; max(PIECES)+1 ]; end
  end

  B = MB;
  if ~isempty( MB ) && isa( B , 'struct' )
    ct = meshCelltype( B );
    if ~isscalar( ct ), error('Mesh B should be homogeneous'); end
    if ct == 5
      B = MeshBoundary( B );
      ct = meshCelltype( B );
    end
    if ct ~= 3, error('Not a set of lines'); end
    if ~isfield( B , 'tri' ) || isempty( B.tri )
      error('Not a mesh');
    end
    B = mesh2contours( B );
    if ~isfield( MB , 'triPIECE' ), MB = MeshAddField( MB ,'triPIECE', max(PIECES)+1 ); PIECES = [ PIECES ; max(PIECES)+1 ]; end
  end

  if ~iscell( A ) && any( isnan( A(:) ) )
    A = nans2split( A );
  else
    A = { A };
  end
  if ~isempty( MB ) && ~iscell( B ) && any( isnan( B(:) ) )
    B = nans2split( B );
  else
    B = { B };
  end

  if ~isempty( MB )

    J = {};
    for a = 1:numel(A)
      for b = 1:numel(B)
        J{end+1} = zipMesh( A{a} , B{b} , true );
      end
    end

  else

    J = {};
    for a = 1:numel(A)
      for b = a+1:numel(A)
        J{end+1} = zipMesh( A{a} , A{b} , true );
      end
    end

  end

  J = J{ argmin( cellfun( @(m)meshSurface(m) , J ) ) };

  if SPLIT

    J.xyzID = ( 1:size(J.xyz,1) ).';
    B = meshSeparate( MeshBoundary( J ) );
    J.xyzS( B{1}.xyzID ,1) =  1;
    J.xyzS( B{2}.xyzID ,1) = -1;

    J = MeshClip( J , J.xyzS , 'both' );
    J.xyz( abs( J.xyzS ) < 1 ,:) = NaN;
    J = MeshRelax( J );

    J = struct( 'xyz', J.xyz , 'tri' ,J.tri );
  
  end
  if FAIR
    B = MeshBoundary( J );
    J = g4Remesh( J , median( meshEdges(B,B) )/2 , 100 , 'HardBoundary',J,'quiet','Speed',0.5);
    %J = stickToMesh( J , MeshRelax(J) );

    M = MeshWeld( MA , J ); f = MA.xyz;
    if ~isempty( MB )
      M = MeshWeld( M , MB ); f = [ f ; MB.xyz ];
    end
    f = knnsearch( M.xyz , f );

    M = MeshFairing( M , f );

    J = MeshRemoveFaces( M , MA ,true);
    if ~isempty( MB )
      J = MeshRemoveFaces( J , MB ,true);
    end

    J = struct( 'xyz', J.xyz , 'tri' ,J.tri );
  end



  try, J = MeshAddField( J  ,'triPIECE', max(PIECES)+1 ); end

  M = MeshWeld( MA , J );
  if ~isempty( MB ) && isa( MB , 'struct' )
    M = MeshWeld( M , MB );
  end

end
