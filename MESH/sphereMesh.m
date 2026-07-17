function S = sphereMesh( its , basePolyhedron )
%SPHEREMESH  Triangulated unit sphere by recursive subdivision of a polyhedron.
%
%   S = sphereMesh( its , basePolyhedron )
%       Returns a mesh struct S with fields .xyz (N x 3 vertices) and .tri
%       (M x 3 triangles) approximating the UNIT sphere. Built by subdividing
%       a base polyhedron 'its' times (1-to-4 split) and projecting the
%       vertices onto the sphere.
%
%   S = sphereMesh( M )    (M a struct) does NOT build anything: it just
%       PROJECTS the vertices of the existing mesh M onto the unit sphere
%       (M.xyz normalized to norm 1; M.tri untouched). Used internally.
%
%   its  - number of subdivision steps. Default 5. It may be:
%       * 0        -> the raw base polyhedron, NOT projected (so 'ico2' and
%                     'tetra', whose base coords are not unit-norm, come back
%                     un-normalized; 'ico'/'octa' are already on the sphere).
%       * N > 0    -> subdivide N times, then project ONCE at the end.
%       * N < 0    -> subdivide |N| times, projecting after EVERY step
%                     (rounder intermediate meshes; slightly different vertex
%                      positions than +N — same connectivity).
%       * a VECTOR -> its phases are applied in sequence (subdivisions
%                     accumulate); e.g. [2 -1] = 3 total subdivisions, the
%                     last one projected per-step.
%       For the icosahedron: V = 10*4^its + 2, F = 20*4^its
%       (its 0..5 -> 12/42/162/642/2562/10242 vertices).
%
%   basePolyhedron - starting solid (default 'icosahedron'). Case-insensitive,
%       with aliases:
%       'icosahedron' | 'ico' | 'i'   12 vertices, 20 faces (default; unit)
%       'icosahedron2'| 'ico2'| 'i2'  golden-ratio icosahedron (NOT unit at its=0)
%       'octahedron'  | 'octa'| 'o'   6 vertices, 8 faces (unit)
%       'tetrahedron' | 'tetra'|'t'   4 vertices, 4 faces (NOT unit at its=0)
%       'rand50' | 'rand500'          convex hull of 50/500 random sphere points
%
%   EXAMPLES:
%       S = sphereMesh();               % default: its=5 icosphere (10242 verts)
%       S = sphereMesh( 3 , 'octa' );   % octahedron subdivided 3 times
%       S = sphereMesh( M );            % snap an existing mesh onto the sphere
%
% See also MeshSubdivide, convhull.

  if nargin == 1 && isstruct( its )
    S = its;
    S.xyz = bsxfun( @rdivide , S.xyz , sqrt( sum( S.xyz.^2 ,2)  ) );
    return;
  end


  if nargin < 2, basePolyhedron = 'icosahedron'; end


  if nargin < 1, its = []; end
  if isempty( its ), its = 5; end

  switch lower( basePolyhedron )
    case {'rand50'}
      S.xyz = randn( 50 ,3);
      S.xyz = bsxfun( @rdivide , S.xyz , sqrt( sum( S.xyz.^2 ,2)  ) );
      S.tri = convhull( S.xyz );

    case {'rand500'}
      S.xyz = randn( 500 ,3);
      S.xyz = bsxfun( @rdivide , S.xyz , sqrt( sum( S.xyz.^2 ,2)  ) );
      S.tri = convhull( S.xyz );

    case {'icosahedron','icosahedrom','ico','i'}
      S = struct( 'xyz' , [0,0,1;0.894427,0,0.4472135;0.276393,0.850651,0.4472135;-0.723607,0.525731,0.4472135;...
                           -0.723607,-0.525731,0.4472135;0.276393,-0.850651,0.4472135;0.723607,0.525731,-0.4472135;...
                           -0.276393,0.850651,-0.4472135;-0.894427,0,-0.4472135;-0.276393,-0.850651,-0.4472135;...
                           0.723607,-0.525731,-0.4472135;0,0,-1] ,...
                  'tri' , [3,1,2;4,1,3;5,1,4;6,1,5;2,1,6;3,2,7;8,3,7;4,3,8;9,4,8;5,4,9;10,5,9;6,5,10;11,6,10;7,2,11;2,6,11;7,12,8;8,12,9;9,12,10;10,12,11;11,12,7]);

    case {'icosahedron2','icosahedrom2','ico2','i2'}
      p = (1+sqrt(5))/2;
      S = struct( 'xyz' , [  1 ,  p ,  0 ; -1 ,  p ,  0 ;  1 , -p ,  0 ; -1 , -p ,  0 ;
                             p ,  0 ,  1 ; -p ,  0 ,  1 ;  p ,  0 , -1 ; -p ,  0 , -1 ;
                             0 ,  1 ,  p ;  0 , -1 ,  p ;  0 ,  1 , -p ;  0 , -1 , -p ] ,...
                  'tri' , [1,2,9;1,5,7;1,7,11;1,9,5;1,11,2;2,6,9;2,8,6;2,11,8;3,4,12;3,5,10;3,7,5;3,10,4;3,12,7;4,6,8;4,8,12;4,10,6;5,9,10;6,10,9;7,12,11;8,11,12] );

    case {'octahedron','octahedrom','octa','o'}
      S = struct( 'xyz' , [-1,0,0;0,-1,0;0,0,-1;0,0,1;0,1,0;1,0,0] , ...
                  'tri' , [1,2,4;1,3,2;1,4,5;1,5,3;2,3,6;2,6,4;3,5,6;4,6,5] );

    case {'tetrahedron','tetrahedrom','tetra','t'}
      S = struct( 'xyz' , [  1 ,  1 ,  1 ; -1 , -1 ,  1 ;  
                            -1 ,  1 , -1 ;  1 , -1 , -1 ] ,...
                  'tri' , [1,2,4;1,3,2;1,4,3;2:4] );
      
    otherwise
      error('Unknown basePolyhedron.');
  end
 %if ~isfield( S ,'tri' ), S.tri = convhull( S.xyz ); end

  for i = 1:numel(its)
    for it = 1:abs( its(i) )
      S = MeshSubdivide( S ); S = struct('xyz',S.xyz,'tri',S.tri);
      if its(i) < 0, S = sphereMesh( S ); end
    end
    if its(i) > 0, S = sphereMesh( S ); end
  end

end
