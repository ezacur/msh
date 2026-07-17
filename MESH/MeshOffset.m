function MO = MeshOffset( M , d )

  M = Mesh( M );
  M = MeshFixCellOrientation( M );

  bb = meshBB( M );
  bb(1,:) = bb(1,:) - 2*d;
  bb(2,:) = bb(2,:) + 2*d;

  dd = d/10;
  DM = I3D( [] , 'X' , bb(1,1):dd:bb(2,1) ,...
                 'Y' , bb(1,2):dd:bb(2,2) ,...
                 'Z' , bb(1,3):dd:bb(2,3) );

	dm = zeros( size( DM ,1:3) ); dm = dm(:);
% 	[~,~,dm] = vtkClosestPoint( M , DM.XYZ );
%   dm( dm > d*4 ) = NaN;
  w = find( isfinite( dm ) );
  
  vtkClosestElement( M ); CLEANOUT = onCleanup( @()vtkClosestElement( [] , [] ) );
  CHUNKSIZE = 100000;
  for i = 1:CHUNKSIZE:numel(w)
    ww = w( i:min( i+CHUNKSIZE-1 , numel(w)) );
    fprintf( '%3.1f %%\n' , ww(end) / numel(w) * 100 );
    [~,~,dm(ww)] = vtkClosestElement( DM.XYZ(ww,:) );
  end

  
  
  dm( ~isfinite( dm ) ) = d*10;
  DM.data(:) = dm;
  
  MO = Mesh( isosurface( DM , d ) );
  MO = meshSeparate( MO , 'sort' , @(m)min(m.xyz(:,3)) ,'select',1,'combine');

  
end
