function [X,onB] = stickToMesh( M , X0 , BOUNDARYtreatment )

  if nargin < 3, BOUNDARYtreatment = 'none'; end

  if isstruct( X0 ) && isfield( X0 , 'xyz' )
    X = X0;
    X.xyz = stickToMesh( M , X0.xyz , BOUNDARYtreatment );
    return;
  end



  if ischar( BOUNDARYtreatment )
  elseif islogical( BOUNDARYtreatment ) && isscalar( BOUNDARYtreatment ) && ~BOUNDARYtreatment
    BOUNDARYtreatment = 'none';
  elseif islogical( BOUNDARYtreatment ) && isscalar( BOUNDARYtreatment ) && ~~BOUNDARYtreatment
    BOUNDARYtreatment = 'original';
  elseif isnumeric( BOUNDARYtreatment ) && isscalar( BOUNDARYtreatment ) && isnan( BOUNDARYtreatment )
    BOUNDARYtreatment = 'NaN';
  elseif isnumeric( BOUNDARYtreatment ) && isscalar( BOUNDARYtreatment ) && ~~BOUNDARYtreatment
    BOUNDARYtreatment = 'original';
  else
    error('Invalid BOUNDARYtreatment.');
  end

  switch lower( BOUNDARYtreatment )
    case {'original','o'},  BOUNDARYcheck = true;
    case {'nan','n'},       BOUNDARYcheck = true;
    case {'none'},          BOUNDARYcheck = false;
    otherwise, error('Invalid BOUNDARYtreatment.');
  end

  M = struct( 'xyz' , double(M.xyz) , 'tri' , double(M.tri) );
  M = MeshTidy( M ,0,true );

  [d,X] = distanceFrom( double(X0) , M , BOUNDARYcheck );

  switch lower( BOUNDARYtreatment )
    case {'original','o'},  onB = d < 0; X( onB ,:) = X0( onB ,:);
    case {'nan','n'},       onB = d < 0; X( onB ,:) = NaN;
    case {'none'},
  end
  
end


% function [X,onB] = stickToMesh( M , X , Nits , Prct , avoidBOUNDARIES )
% 
%   if nargin < 5, avoidBOUNDARIES = false; end
% 
%   M = struct('xyz',double(M.xyz),'tri',double(M.tri));
%   vtkClosestElement([],[]);
%   vtkClosestElement( M );
%   CLEANUP = onCleanup( @()vtkClosestElement([],[]) );
%   
%   X = double( X );
%   onB = false( size(X,1) , 1);
%   
%   if avoidBOUNDARIES
%     Bedges = MeshBoundary( M.tri );
%     Bedges = sort( Bedges ,2);
%     u = unique( Bedges );
%     Bedges( end + (1:numel(u)) , 2 ) = u;
%     Bedges = sortrows( Bedges );
%   
%     isBtri = any( ismember(  M.tri , u ) ,2);
%   end
%   
%   for it = 1:Nits+1
%     if ~avoidBOUNDARIES
%       [~,F] = vtkClosestElement( X );
%     else
%       w = ~onB;
%       [ eid(w) , F(w,:) , ~ , bc(w,:) ] = vtkClosestElement( X(w,:) );
%       e = isBtri( eid ); e( onB ) = false;
%       if any(e)
%         b = bc > 1e-5;
%         b( ~e ,:) = false;
%         b( all(b,2) ,:) = false;
%         b = M.tri(eid,:) .* b;
%         b = sort( b ,2);
%         b = b(:,2:3);
%   
%         e = find( ~~b(:,2) );
%         e = e( ismember( b(e,:) , Bedges ,'rows' ) );
%         onB( e ) = true;
%       end
%       F( onB ,:) = X( onB ,:);
%     end
%     if it > Nits,     X = F;
%     else,             X = X + ( F - X ) * Prct;
%     end
%   end
% 
% end
