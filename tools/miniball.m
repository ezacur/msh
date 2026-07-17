function [c,r] = miniball( x , onlyR )
%MINIBALL  Smallest enclosing ball (or axis-aligned ellipsoid) of a point set.
%
%   [C,R] = miniball( X )                  minimum enclosing BALL
%   [R,C] = miniball( X , true )           radius FIRST (see ONLYR below)
%   [C,R] = miniball( X , 'ellipse' )      minimum-volume AXIS-ALIGNED ellipsoid
%
%   X is an N-by-ND array of N points (one per row) in ND dimensions. The ball
%   mode works in ANY dimension; it returns the center C (1-by-ND) and the
%   radius R (scalar) of the unique smallest ball that contains every point.
%   The heavy lifting is done by the compiled MEX miniball_mx, a wrapper around
%   the Seb library (B. Gaertner's move-to-front "smallest enclosing ball"):
%   robust and exact to machine precision, with ND+1 points on the boundary.
%
%   ONLYR (2nd arg) selects the mode and is STRICTLY validated -- it must be
%   one of true / false / 0 / 1 / 'ellipse' / 'e'; anything else errors.
%     true or 1  : SWAP the outputs so the RADIUS comes out first, i.e.
%                  R = miniball(X,true) hands you the radius in the first output.
%     false or 0 : the default plain ball.
%     'ellipse'/'e' : the ellipsoid mode below.
%
%   'ELLIPSE' / 'E' mode (3-D ONLY): returns the minimum-VOLUME enclosing
%   ellipsoid restricted to be AXIS-ALIGNED. R is then a 1-by-3 vector of
%   semi-axes and C the center. It is found by scaling each axis by SC, taking
%   the enclosing ball of the scaled cloud, and minimizing the ellipsoid volume
%   R^3*prod(SC) over SC with FMINSEARCH (optimized in LOG(SC) so the scales
%   stay positive and the search cannot diverge). NOTE this is NOT the general
%   (arbitrarily oriented) Loewner-John MVEE; FMINSEARCH is a LOCAL optimizer
%   seeded at SC=1. It is a fast approximation that gives a reasonable volume,
%   not the true minimum-volume ellipsoid.
%
%   Non-finite rows of X are dropped before the computation. If no finite point
%   remains, C = NaN(1,ND) and R = NaN.
%
%   Example:
%     X = randn(500,3)*[3 0 0;0 1 0;0 0 1];   % an elongated cloud
%     [C,R]  = miniball( X );                  % tightest ball
%     [Ce,Re]= miniball( X , 'ellipse' );      % tightest axis-aligned ellipsoid
%
% See also miniball_mx, bwlargest, sphereMesh, transform.
if 0

  [C,R] = miniball( X , 'ellipse' );

  plot3d( X ,'okr','eq');
  hplotMESH( transform( sphereMesh(4) ,'s',R,'t',C) ,'[0.5]','nice')


  %%
end


  if nargin < 2, onlyR = false; end

  %-- validate ONLYR: true/false/0/1 (radius as 1st output) or 'ellipse'/'e' ----
  ELLIPSE = ( ischar(onlyR) || isstring(onlyR) ) && ( strcmpi(onlyR,'ellipse') || strcmpi(onlyR,'e') );
  if ~ELLIPSE
    if ~( ( islogical(onlyR) || isnumeric(onlyR) ) && isscalar(onlyR) && ( isequal(onlyR,0) || isequal(onlyR,1) ) )
      error( 'miniball:onlyR' , 'ONLYR must be true, false, 0 or 1, or the string ''ellipse'' / ''e''.' );
    end
    onlyR = logical( onlyR );
  end

  nd = size( x ,2);

  w = all( isfinite(x) ,2);
  x = x(w,:);

  if isempty( x )
    c = NaN(1,nd);
    r = NaN;
    return;
  end
  x = double(x);
  
  if ELLIPSE
    if nd ~= 3, error('miniball:ellipse','''ellipse'' mode requires 3-D points.'); end
    %optimize the per-axis scales in LOG space so they stay strictly POSITIVE.
    %The raw objective r^3*prod(sc) is UNBOUNDED BELOW for negative sc (odd #of
    %sign flips -> prod(sc)<0 -> fminsearch runs away to sc->-Inf). With
    %sc = exp(p) the volume is always > 0 and has a genuine interior minimum
    %(sc->0 blows r up, sc->Inf blows the volume up), so the search can't leak.
    p  = fminsearch( @(p)miniball_V( x , exp(p) ) , zeros(1,nd-1) );
    sc = exp( p );
    [c,r] = miniball_mx( bsxfun( @rdivide , x , [sc,1] ) );
    sc(nd) = 1;
    r = r  * sc;
    c = c .* sc;
    return;
  end



  [c,r] = miniball_mx( x );
  if onlyR
    [r,c] = deal(c,r);
  end

end
function V = miniball_V( x , sc )
  [~,r] = miniball_mx( bsxfun( @rdivide , x , [sc,1] ) );
  V = r^3 * prod( sc );
end

