function [Surface,Centroid] = meshSurface( M , centroid )
% [S,C] = meshSurface( M )
% 
% S is the total surface area of the triangles
% C is the centroid of the shell (not the contained volume)
% 


  if nargin < 2, centroid = false; end


  celltype = meshCelltype( M );
  if celltype ~= 5, error('only triangle meshes are allowed up to now.'); end
  
  Pa = M.xyz( M.tri(:,1) ,:);
  Pb = M.xyz( M.tri(:,2) ,:);
  Pc = M.xyz( M.tri(:,3) ,:);

  L1 = Pb - Pa;   L1(:,end+1:3) = 0;
  L2 = Pa - Pc;   L2(:,end+1:3) = 0;

  L2xL1 = cross( L2 , L1 ,2);
  
  doubleArea = sqrt( sum( L2xL1.^2 ,2) );
  Surface = sum( doubleArea )/2;
  
  if nargout > 1 || centroid
    Centroid = sum( bsxfun( @times , Pa + Pb + Pc , doubleArea ) ,1)/( 6 * Surface );
  end
  if centroid
    [ Centroid , Surface ] = deal( Surface , Centroid );
  end

end
