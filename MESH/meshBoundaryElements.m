function b = meshBoundaryElements( M , its )

  if nargin < 2, its = 1; end

  M = struct( 'tri' , M.tri );
%   M = struct( 'xyz' , M.xyz , 'tri' , M.tri );
%   M.xyzID = ( 1:size( M.xyz ,1) ).';
  M.triID = ( 1:size( M.tri ,1) ).';
  
  b = false( size(M.tri,1) , 1);
  for i = 1:its
    B = MeshBoundary( M );
    b_ = any( ismember( M.tri , B.tri ) , 2);

    b( M.triID( b_ ) ) = true;
    M.tri( b_ ,:) = [];
    M.triID( b_ ,:) = [];
  end

end
