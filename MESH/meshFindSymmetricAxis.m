function R = meshFindSymmetricAxis( M , varargin )
%MESHFINDSYMMETRICAXIS  Find the best mirror-symmetry plane of a surface mesh.
%
%   R = meshFindSymmetricAxis( M )
%   R = meshFindSymmetricAxis( M , OptimizeOptions... )
%
%   Locates the plane about which the mesh M is most nearly mirror-symmetric.
%   The plane is found by minimizing, over rigid poses, the sum of squared
%   distances between every REFLECTED vertex and the mesh surface:
%
%       E(p) = sum_i  dist( reflect(v_i) , surface(M) )^2
%
%   where reflect() mirrors each vertex across the candidate plane. E is 0 for a
%   perfectly symmetric mesh and grows with asymmetry, so its minimizer is the
%   best-fitting symmetry plane.
%
%   INPUT
%     M    surface-mesh struct with fields .xyz (#V-by-3 vertex coordinates)
%          and .tri (#F-by-3 triangle indices). Both are cast to double.
%     ...  extra name/value options, forwarded verbatim to Optimize (to tweak
%          tolerances, line search or verbosity of the local refinement).
%
%   OUTPUT
%     R    4-by-4 homogeneous transform giving the POSE of the symmetry plane:
%          it maps a reference frame whose XY-plane is the mirror into the
%          coordinates of M. Hence, expressed in the frame of M:
%             plane normal         = R(1:3,3)
%             a point on the plane = R(1:3,4)
%             the mirror operator  = R * diag([1 1 -1 1]) * minv(R)
%
%   METHOD
%     A plane has 3 degrees of freedom, so the pose is parametrized minimally by
%     p = [ tz ; rx ; ry ]: a translation along z and rotations (in DEGREES)
%     about x then y (see the local subfunction RxRyTz; the rotation block is
%     Ry*Rx). The rigid DOF that do NOT move the z=0 plane -- roll about z and
%     x/y translation -- are deliberately omitted. The minimum is sought in two
%     stages (both derivative-free to the caller; a numerical Jacobian is used
%     internally):
%       1) ExhaustiveSearch : coarse adaptive-grid sweep that brackets the basin
%          (seeded at p = 0, initial radius 10, grown as needed);
%       2) Optimize         : conjugate + coordinate-descent local refinement.
%
%   NOTES / CAVEATS
%     * The search is LOCAL and greedy, seeded at the XY-plane through the origin
%       (normal +z, zero offset). The grid radius grows, but for a plane far from
%       that seed pre-orient M (e.g. by PCA) or pass a better start via Optimize.
%     * The distance metric is one-directional (reflected VERTICES -> surface)
%       and vertex-sampled; adequate for near-symmetric meshes.
%     * Every energy evaluation re-poses the whole mesh and rebuilds the
%       closest-point locator (vtkClosestElement); this is the cost bottleneck
%       on large meshes.
%
%   See also vtkClosestElement, ExhaustiveSearch, Optimize, minv, maketransform.

  M = struct( 'tri' , double( M.tri ) , 'xyz' , double( M.xyz ) );
  MR = M;
  

  p = [0;0;0];
  p = ExhaustiveSearch( @(z)ENER(z) , p , 10 , 3 , 'maxIterations', 50 ,'verbose');
  p = Optimize( @(z)ENER(z) , p , 'methods',{'conjugate','coordinate',5},'ls',{'quadratic','golden','quadratic'},struct('COMPUTE_NUMERICAL_JACOBIAN',{{'f'}}),varargin{:});
  
  R = minv( RxRyTz( p ) );
  
  function E = ENER( p )
    R = RxRyTz( p );
    
    MR.xyz = bsxfun( @plus , M.xyz * R(1:3,1:3).' , R(1:3,4).' );
    XYZ = MR.xyz;
    XYZ(:,3) = -XYZ(:,3);
    
    [~,~,d] = vtkClosestElement( MR , XYZ );
    E = sum( d.^2 );
  end
end

function R = RxRyTz( zxy )
% 4x4 homogeneous transform from the 3-DOF plane pose. NOTE the argument order:
%   zxy = [ tz ; rx ; ry ]  ->  z-translation, then rotations (deg) about x, y.
% The rotation block equals Ry*Rx. Equivalent to:
  %R = maketransform( 'rx' , p(2) , 'ry' , p(3) , 't' , [0,0,p(1)] );
  z = zxy(1); if abs(z) < eps(1), z = 0; end
  x = zxy(2); if abs(x) < eps(1), x = 0; end
  y = zxy(3); if abs(y) < eps(1), y = 0; end


  cx = cosd( x );
  sx = sind( x );
  cy = cosd( y );
  sy = sind( y );
  R = [ cy , sx*sy , cx*sy , 0 ;...
         0 ,    cx ,   -sx , 0 ;...
       -sy , cy*sx , cx*cy , z ;...
         0 ,     0 ,     0 , 1 ];

end
