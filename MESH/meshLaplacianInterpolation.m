function interpolator = meshLaplacianInterpolation( M , C )

  if isstruct( M ) && isfield( M ,'xyz' ) && isfield( M ,'tri' )
    [~,~, M ] = tuftedLaplacianFree( M.xyz , M.tri );
  end


  if islogical( C ), C = find( C ); end
  k = numel( C );
  n = size( M ,1);
  
  U = 1:n;
  U( C ) = [];
  ids = [ C , U ];
  [ ~ , o ] = sort( ids );

  % reshuffle rows & columns of lap matrix
  M = M( ids , ids );
  B = M( : , 1:k );
  A = M( : , (k+1):end );
  
  % Convert to sparse for quicker computation
  interpolator = [ eye(k,k) ; - A \ B ];
  interpolator = interpolator( o ,:);

end

function [ A , B , laplacianFullMatrix ] = tuftedLaplacianFree( V , F )

  cmd0 = fileparts( mfilename('fullpath') );

  % get the path
  cmd0 = fullfile( cmd0 , 'meshLaplacian.exe' );

  [D,CLEAN] = sanddir( '#_****' );
  
  write_OBJ( Mesh(V,F) , fullfile( D , 'obj.obj') );
  
  cmd = sprintf( 'cd /d "%s" & "%s"  obj.obj  --writeLaplacian  --writeMass' , D , cmd0 );
  [~,~] = system(cmd);

  nV = size(V,1);

  A = load( fullfile( D , 'tufted_laplacian.spmat' ) );
  A = sparse( A(:,1) , A(:,2) , A(:,3) , nV , nV );
  
  B = load( fullfile( D , 'tufted_lumped_mass.spmat' ) );
  B = sparse( B(:,1) , B(:,2) , B(:,3) , nV , nV );
  
  if nargout > 2, laplacianFullMatrix = B\A; end
end
