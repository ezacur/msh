function [hF_,hB] = backFaceCulling( hF , varargin )
%BACKFACECULLING  Hide the back-facing faces of a patch (view-dependent culling).
%
%   [hF,hB] = backFaceCulling( hF )                back faces -> INVISIBLE patch
%   [hF,hB] = backFaceCulling( hF , 'Prop',v,... ) back faces styled instead
%
%   Sugar for  backFaceCullingSplit( hF , 'Visible','off' ) : hF keeps only its
%   camera-facing faces, tracking the camera; the invisible back patch hB costs
%   nothing per camera event (its Faces refresh only if you make it visible).
%   Deleting hB undoes the culling.
%
%   See also backFaceCullingSplit, silhouette.

  if nargin == 1
    varargin = {'Visible','off'};
  end

  hB = backFaceCullingSplit( hF , varargin{:} );

  if nargout, hF_ = hF; end

end
