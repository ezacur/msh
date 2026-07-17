function Z = MeshZipSmooth( A , B , voxelSize )

  A = MeshFixCellOrientation( A ); a = MeshBoundary( A ); if numel( meshSeparate(a) ) > 1, error('Only Single boundary meshes are allowed in A.'); end; a = mesh2contours( a );
  B = MeshFixCellOrientation( B ); b = MeshBoundary( B ); if numel( meshSeparate(b) ) > 1, error('Only Single boundary meshes are allowed in B.'); end; b = mesh2contours( b );
  
  Z = zipMesh( a , b ,true);
  
  [R,iR] = getPlane( Z );
  Ar = transform( A ,iR);
  Br = transform( B ,iR);
  Z  = transform( Z ,iR);
  
  DM = bsxfun( @plus , meshBB(Z,1.1) , [-1;1]*4*voxelSize );
  DM = I3D('X',DM(1,1):voxelSize:DM(2,1),'Y',DM(1,2):voxelSize:DM(2,2),'Z',DM(1,3):voxelSize:DM(2,3));
  
  M = MeshFixCellOrientation( MeshTidy( MeshAppend( Ar , Br ) ,0,true) );
  Z = MeshFixCellOrientation( MeshTidy( MeshAppend( M  , Z  ) ,0,true) );
  
  DM = MeshToSignedDistance( M , DM );
  DZ = MeshToSignedDistance( Z , DM );
  
  if isempty( which( 'approximatedHS' ) )
    addpath( 'C:\Dropbox\mTools_at_Apr2022\Tools\' );
  end
  to_11 = @(x)approximatedHS( x ,'cos',{0,0.25,-1,2});
  w = DM.data ~= DZ.data;
  try,    w = dilate3d( w );
  catch,  w = imdilate( w , strel( 'cube',3) );
  end
  w = imfill( w ,'holes' );
  DM.data(w) = NaN;
  
  DM.data = nonans( to_11( DM.data ) , 'mingradient' ,'replicate','initial', to_11( DZ.data) );
  
  U = Mesh( isosurface( DM , 0 ) );
  
  U  = transform( U  ,R);
  
  % plotMESH( { A1 , A2 , U } );
  U = MeshRemoveNodes( U , distanceFrom( U.xyz , A ) < voxelSize*4 );
  U = MeshRemoveNodes( U , distanceFrom( U.xyz , B ) < voxelSize*4 );



  % plotMESH( { A1 , A2 , U } )
  
  u = meshSeparate( MeshBoundary( U ) );
  if numel( u ) ~= 2, warning('2 boundaries are expected in U.'); end

  A = fun( @(v)zipMesh( mesh2contours( v ) , a ,true) , u ,'un',0);
  id = argmin( fun( @meshSurface , A  ) );
  u(id) = []; A = A{id};

  B = fun( @(v)zipMesh( mesh2contours( v ) , b ,true) , u ,'un',0);
  id = argmin( fun( @meshSurface , B  ) );
  u(id) = []; B = B{id};
  
  Z = MeshAppend( A , B , U );
  Z = MeshTidy( Z ,0,true);

end

function D = MeshToSignedDistance( M , res )

  if isempty( which( 'meshIsInterior' ) )
    addpath( 'C:\Dropbox\mTools_at_Apr2022\MESH\' );
  end

  if nargin < 2, res = 50; end
  
  if isnumeric( res )
    if numel(res) < 3 , res(end+1:3) = res(end); end

    Cs = {[],[],[]};
    BB = meshBB( M );
    for d = 1:3
      BB(:,d) = ( BB(:,d) - mean( BB(:,d) ) ) * 1.2 + mean( BB(:,d) );
      if res(d) > 0
        Cs{d} = linspace( BB(1,d) , BB(2,d) , res(d) );
      elseif res(d) < 0
        n = floor( - diff( BB(:,d) )/2/res(1) )+1;
        Cs{d} = ( -n:n ) * -res(d) + mean( BB(:,d) );
      end
    end
    
    D = I3D( NaN( [ numel(Cs{1}) , numel(Cs{2}) , numel(Cs{3}) ] ) ,'X',Cs{1},'Y',Cs{2},'Z',Cs{3} );
    
  elseif isa( res , 'I3D' )
    D = res;
    D.data(:) = NaN;
  end
  
  D.data(:) = distanceFrom( D.XYZ , M );
  w = meshIsInterior( M , D.XYZ );
  D.data(w) = -D.data(w);

end


