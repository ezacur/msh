function M = MeshSubdivide_reduced( M , W )
  nT = size( M.tri , 1);   %number of faces
  nP = size( M.xyz , 1);   %number of points
  
  T  = M.tri;
  F  = ( 1:nT ).';         %face indexes
  
  allE = sort( [ T(:,[1 2]) ; T(:,[2 3]) ; T(:,[1 3]) ] ,2); %no remove the repeated
  E = unique( allE( [ W ; W + nT ; W + 2*nT ] , : ) ,'rows');
  [~,ET] = ismember( allE , E(:,1:2) , 'rows' );
  
  ET = reshape( ET , [nT,3] );
  LET = ~~ET;  %logical ET
  W = find( any( LET ,2) );
  P = ( (nP+1):(nP+size(E,1)) ).';
  for f = fieldnames( M ).'
    if ~strncmp( f{1} , 'xyz' , 3 ), continue; end
    M.(f{1}) = [ M.(f{1}) ; ( M.(f{1})( E(:,1) ,:,:,:,:,:) + M.(f{1})( E(:,2) ,:,:,:,:,:) )/2 ];
  end

  w = find( LET(:,1) & LET(:,2) & LET(:,3) );
  T = [ T ; [     T(w, 1 )   , P( ET(w, 1 ) ) , P( ET(w, 3 ) ) ] ;...
            [ P( ET(w, 1 ) ) ,     T(w, 2 )   , P( ET(w, 2 ) ) ] ;...
            [ P( ET(w, 3 ) ) , P( ET(w, 2 ) ) ,     T(w, 3 )   ] ;...
            [ P( ET(w, 1 ) ) , P( ET(w, 2 ) ) , P( ET(w, 3 ) ) ] ];
  F = [ F ; w ; w ; w ; w ];

  w = find( LET(:,1) & ~LET(:,2) & ~LET(:,3) );
  T = [ T ; [     T(w, 1 )   , P( ET(w, 1 ) ) ,     T(w, 3 )   ] ;...
            [ P( ET(w, 1 ) ) ,     T(w, 2 )   ,     T(w, 3 )   ] ];
  F = [ F ; w ; w ];

  w = find( ~LET(:,1) & LET(:,2) & ~LET(:,3) );
  T = [ T ; [     T(w, 1 )   ,     T(w, 2 )   , P( ET(w, 2 ) ) ] ;...
            [     T(w, 1 )   , P( ET(w, 2 ) ) ,     T(w, 3 )   ] ];
  F = [ F ; w ; w ];
  
  w = find( ~LET(:,1) & ~LET(:,2) & LET(:,3) );
  T = [ T ; [     T(w, 1 )   ,     T(w, 2 )   , P( ET(w, 3 ) ) ] ;...
            [     T(w, 2 )   ,     T(w, 3 )   , P( ET(w, 3 ) ) ] ];
  F = [ F ; w ; w ];
  
  w = find( LET(:,1) & LET(:,2) & ~LET(:,3) );
  T = [ T ; [ P( ET(w, 1 ) ) ,     T(w, 2 )   , P( ET(w, 2 ) ) ] ;...
            [     T(w, 1 )   , P( ET(w, 1 ) ) ,     T(w, 3 )   ] ;...
            [ P( ET(w, 1 ) ) , P( ET(w, 2 ) ) ,     T(w, 3 )   ] ];
  F = [ F ; w ; w ; w ];
  
  w = find( LET(:,1) & ~LET(:,2) & LET(:,3) );
  T = [ T ; [     T(w, 1 )   , P( ET(w, 1 ) ) , P( ET(w, 3 ) ) ] ;...
            [ P( ET(w, 1 ) ) ,     T(w, 3 )   , P( ET(w, 3 ) ) ] ;...
            [ P( ET(w, 1 ) ) ,     T(w, 2 )   ,     T(w, 3 )   ] ];
  F = [ F ; w ; w ; w ];

  w = find( ~LET(:,1) & LET(:,2) & LET(:,3) );
  T = [ T ; [ P( ET(w, 3 ) ) , P( ET(w, 2 ) ) ,     T(w, 3 )   ] ;...
            [     T(w, 1 )   , P( ET(w, 2 ) ) , P( ET(w, 3 ) ) ] ;...
            [     T(w, 1 )   ,     T(w, 2 )   , P( ET(w, 2 ) ) ] ];
  F = [ F ; w ; w ; w ];
  
  M.tri      = T;
  F(W)       = [];
  M.tri(W,:) = [];            %remove the original faces
  [F,ord] = sort( F );      %reorder the new faces in their "original position"
  M.tri = M.tri(ord,:);
  for f = fieldnames( M ).'
    if strcmp( f{1} , 'tri' ), continue; end
    if ~strncmp( f{1} , 'tri' , 3 ), continue; end
    M.(f{1}) = M.(f{1})(F,:,:,:,:,:,:);
  end

end
