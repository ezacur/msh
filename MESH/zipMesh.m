function M = zipMesh( A , B , surfaceTest )
% ZIPMESH   Triangulate ("zip") the band between two contours.
%
%   M = ZIPMESH( A , B )
%   M = ZIPMESH( A , B , surfaceTest )
%
%   Builds a triangle strip that stitches contour A to contour B, i.e. the
%   lateral surface that joins both polylines.  It is the core primitive used
%   by MeshZipPieces to bridge the boundaries of two mesh pieces.
%
% INPUTS
%   A , B        Contour vertices, one per row ( [nA x D] and [nB x D], D = 2
%                or 3 ), given in order along the contour.  Must be finite
%                ( no NaN ).
%
%                A contour is treated as CLOSED when its first and last rows
%                are identical ( the repeated closing vertex is dropped
%                internally ); otherwise it is treated as OPEN.  Any of the
%                four closed/open combinations is accepted.
%
%   surfaceTest  Logical, default FALSE.  When TRUE the strip is also built
%                with B reversed ( flip(B,1) ) and the version with the smaller
%                total edge length is returned, which resolves the relative
%                orientation of the two contours automatically.
%
% OUTPUT
%   M            Triangle mesh ( celltype 5 ) struct with fields:
%                  M.xyz = [ A ; B ]  ( contour vertices, closing dup. removed )
%                  M.tri              ( [nT x 3] triangles, indices into M.xyz )
%                Triangle count: nA+nB ( both closed ), nA+nB-2 ( both open ).
%
% ALGORITHM
%   Minimum-cost monotone path through the pairwise distance matrix ipd(A,B)
%   ( cost = sum of squared "rung" lengths ) -- a minimal-area-ish
%   triangulation of the band.  For closed contours the seam is anchored at
%   the closest pair of vertices.  Output triangles are consistently oriented.
%
% NOTES
%   * Heuristic: not guaranteed to be globally minimal-area nor free of
%     self-intersections for strongly non-convex / non-planar contours.
%   * Open contours are matched endpoint-to-endpoint ( use surfaceTest, or
%     pre-flip B, if the natural correspondence is reversed ).
%
% CLOSED-CLOSED SEAM  ( note to self, if a tube ever comes out twisted )
%   For two closed contours the strip is a tube whose seam is anchored at the
%   single closest pair of vertices ( the min(D(:)) + circshift step in the
%   Aclose && Bclose branch ), then unrolled with one monotone DP.  That is the
%   optimum only CONDITIONED on that seam, not the global optimum of the cyclic
%   problem.  In practice it is already optimal: an exhaustive seam sweep gave
%   0% area gain over many contour pairs -- ellipses rotated / offset, clustered
%   vertices, offset bumps ( verified 2026-07 ).  So normally: leave it alone.
%
%   If one day a closed-closed zip does come out twisted / higher-area, the
%   exact fix is to SEARCH the seam: fix A's cut at vertex 1 and try every
%   rotation of B, keeping the one with the smallest meshSurface:
%
%       Sopt = Inf; Bopt = [];
%       for b = 1:size(B,1)
%         Bb = B( [b:end,1:b] ,:);
%         Z  = zipMesh_fixedSeam( A , Bb , true );   % variant, see IMPORTANT
%         s  = meshSurface( Z );
%         if s < Sopt, Sopt = s; Bopt = Bb; end
%       end
%       Z = zipMesh_fixedSeam( A , Bopt , true );
%
%   IMPORTANT: that loop needs a variant of THIS function with the internal
%   re-anchor removed ( delete the min(D(:)) + circshift lines in the
%   Aclose && Bclose branch; keep the padding/wrap so the tube still closes ).
%   Calling the stock zipMesh in the loop is a NO-OP -- it re-anchors on its own,
%   so every b returns the identical mesh.  And do NOT try to dodge the closed
%   detection with  A(end,:) = A(end,:) + eps(...) : it works, but the repeated
%   seam vertex is then NOT dropped, leaving two seam nodes 1 ULP apart and
%   unwelded -> a hairline crack ( boundary becomes 1 loop instead of 2 ).
%
% See also MeshZipPieces, mesh2contours, ipd, MeshWeld.

  if nargin < 3, surfaceTest = false; end
  if ~isscalar( surfaceTest ), error('surfaceTest should be true/false.'); end

  if any( ~isfinite( A(:) ) ), error('Only simple contours are allowed in A.'); end
  if any( ~isfinite( B(:) ) ), error('Only simple contours are allowed in B.'); end

  if surfaceTest
    M  = zipMesh( A ,       B     , false ); s  = sum( meshQuality( MeshWireframe( M  ) , 'size') );
    Mf = zipMesh( A , flip( B ,1) , false ); sf = sum( meshQuality( MeshWireframe( Mf ) , 'size') );
    if s > sf
      M = Mf;
    end

    return;
  end


  Aclose = false;
  if isequal( A(1,:) , A(end,:) ), A(end,:) = []; Aclose = true; end

  Bclose = false;
  if isequal( B(1,:) , B(end,:) ), B(end,:) = []; Bclose = true; end


  D = ipd( A , B );
  if 0
  elseif  Aclose &&  Bclose
    m = min( D(:) ); [a,b] = find( D == m ,1);
    A = circshift( A , -a+1 ,1); D = circshift( D , -a+1 ,1);
    B = circshift( B , -b+1 ,1); D = circshift( D , -b+1 ,2);
  elseif ~Aclose &&  Bclose
    b = argmin( D(1,:) );
    B = circshift( B , -b+1 ,1); D = circshift( D , -b+1 ,2);
  elseif  Aclose && ~Bclose
    a = argmin( D(:,1) );
    A = circshift( A , -a+1 ,1); D = circshift( D , -a+1 ,1);
  end
  if 0
     plot3d( A ,'r'); text( A(:,1) , A(:,2) , A(:,3) , fun(@num2str,1:size(A,1)),'Color','r');
    hplot3d( B ,'b'); text( B(:,1) , B(:,2) , B(:,3) , fun(@num2str,1:size(B,1)),'Color','b'); ze
    maxnorm( D , ipd(A,B) )
    %%
  end
  if 0
  elseif  Aclose &&  Bclose
    D = [ D ; D(1,:) ];
    D = [ D , D(:,1) ];
  elseif ~Aclose &&  Bclose
    D = [ D , D(:,1) ];
  elseif  Aclose && ~Bclose
    D = [ D ; D(1,:) ];
  end
    
  D = D.^2;
  for i = 1:size(D,1)
    for j = 1:size(D,2)
      if i == 1 && j == 1, continue; end
      p1 = Inf; try, p1 = D(i-1,j); end
      p2 = Inf; try, p2 = D(i,j-1); end
      try, D(i,j) = min( p1 , p2 ) + D(i,j); end
    end
  end
  
  E = [ size(D,1) , size(D,2) ];
  while ~isequal( E(end,:) , [1,1] ) 
    p1 = E(end,:) - [1,0]; v1 = Inf; try, v1 = D(p1(1),p1(2)); end
    p2 = E(end,:) - [0,1]; v2 = Inf; try, v2 = D(p2(1),p2(2)); end
    if v1 < v2, E = [ E ; p1 ];
    else      , E = [ E ; p2 ];
    end
  end
  E( E(:,1) > size(A,1) ,1) = 1;
  E( E(:,2) > size(B,1) ,2) = 1;
   
  %%
  E(:,2) = E(:,2) + size(A,1);
  
  M = struct('xyz',[ A ; B ] ,'tri',[] );
  for e = 2:size(E,1)
    M.tri = [ M.tri ; unique( [ E(e-1,1) , E(e,1) , E(e,2) , E(e-1,2) ] ,'stable') ];
  end

end
