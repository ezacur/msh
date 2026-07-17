function bb = meshBB( M , sc )

  bb = [ min(M.xyz,[],1) ; max(M.xyz,[],1) ];
  
  if nargin > 1

    cs = @(x,s) ( x - mean(x) ) * s + mean( x );
    for c = 1:size( bb ,2)
      bb(:,c) = cs( bb(:,c) , sc );
    end
  
  end


end