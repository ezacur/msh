function M = fillContoursMesh( C , delta )
%FILLCONTOURSMESH  Triangle mesh filling the region bounded by closed contours.
%
%   M = fillContoursMesh( C )            fill the contour(s) C (default delta=Inf)
%   M = fillContoursMesh( C , delta )    with the interior density set by delta
%
%   C describes one or more (nearly) planar boundary curves, as either:
%     * a cell  {V1,V2,...}   each Vk an Np-by-2 or Np-by-3 list of the ordered
%                             loop vertices (all curves same dimension);
%     * a numeric  V          one curve, or several separated by NaN rows;
%     * a segment mesh struct (celltype 3)  -> its polylines (mesh2contours).
%   With SEVERAL curves the nested/inner loops become HOLES (the constrained
%   Delaunay keeps only what isInterior marks inside the loop set), so e.g. a
%   disc with two inner circles comes out as a disc with two holes.
%
%   DELTA sets how the interior is filled:
%     Inf   (DEFAULT)   NO interior nodes: just triangulate the polygon from its
%                       own boundary vertices (a constrained Delaunay fill).
%     > 0               add a near-equilateral triangular grid of spacing delta
%                       inside, for a roughly UNIFORM interior mesh.
%     []   (empty)      auto: delta = median boundary edge length.
%     < 0               delta = max(|delta|, median edge length) -- a FLOOR on
%                       triangle size, so the interior is never refined finer
%                       than the contour's own sampling.
%     NaN               FAN / centered triangulation: add the centroid and join
%                       every boundary segment to it. SINGLE, EXPLICITLY CLOSED
%                       curve only (first vertex repeated as the last); errors
%                       otherwise. Cheapest, best for star-shaped contours.
%   For the finite / Inf / empty / negative modes the interior is smoothed with
%   MeshRelax; for NaN it is returned as-is.
%
%   NON-PLANAR curves are handled: the filling is done in 2-D and mapped back to
%   3-D. Several curves (or a single OPEN curve) are projected onto their
%   best-fit plane (getPlane) and any off-plane deviation is restored by a
%   thin-plate-spline lift; a single CLOSED curve is developed into the plane
%   (MeshFlatten) and warped back through an interpolating spline. Either way
%   the boundary vertices are snapped back to their EXACT input coordinates, so
%   M reproduces C on its border and only the new interior nodes are inferred.
%
%   M is a mesh struct (M.xyz nodes, M.tri triangles = celltype 5). The first
%   size(boundary) rows of M.xyz are the contour vertices, in input order.
%
%   Run  edit fillContoursMesh  and evaluate the leading `if 0` block for a
%   worked 3-curve example.
%
% See also mesh2contours, MeshFillSmallestHole, MeshRelax, MeshFlatten,
%          getPlane, InterpolatingSplines, delaunayTriangulation.

if 0

%Boundary 1
ns=150;
t=linspace(0,2*pi,ns);
t=t(1:end-1);
r=6+2.*sin(5*t);
[x,y] = pol2cart(t,r);
z=1/10*x.^2;
V1=[x(:) y(:) z(:)];

%Boundary 2
ns=100;
t=linspace(0,2*pi,ns);
t=t(1:end-1);
[x,y] = pol2cart(t,ones(size(t)));
z=zeros(size(x));
V2=[x(:) y(:)+4 z(:)];

%Boundary 3
ns=75;
t=linspace(0,2*pi,ns);
t=t(1:end-1);
[x,y] = pol2cart(t,2*ones(size(t)));
z=zeros(size(x));
V3=[x(:) y(:)-0.5 z(:)*0];

%Create Euler angles to set directions
% E=[0.25*pi -0.25*pi 0];
% [R,~]=euler2DCM(E); %The true directions for X, Y and Z axis
% 
% V1=(R*V1')'; %Rotate polygon
% V2=(R*V2')'; %Rotate polygon
% V3=(R*V3')'; %Rotate polygon

M = fillContoursMesh( {V1,V2,V3} , 0.2 );
M = fillContoursMesh( V1 , 2 );



plotMESH( M );
hplot3d( {V1,V2,V3} , '.-r' );
  
%%  
end
  if nargin < 2
    delta = Inf;
  end

  if isstruct( C )
    C.celltype = meshCelltype( C );
    if C.celltype ~= 3, error('only LINE meshes are allowed'); end

    C = mesh2contours( C );
    C = nans2split( C );
  elseif isnumeric( C )
    C = nans2split( C );
  end
  
  if numel( C ) == 1
    if ~isequal( C{1}(1,:) , C{1}(end,:) )
      C{2,1} = [];
    end
  end
  
  
  if numel( C ) > 1

    C( cellfun('isempty',C) ) = [];
    
    nsd = cellfun( @(c)size(c,2) , C );
    if ~all( nsd == nsd(1) ), error('mixed Number of Spatial Dimensions'); end
    nsd = nsd(1);

    if isnan( delta )
      error('NaN delta corresponds to a centered triangulation which is only valid for single curves closed.');
    end
    
    if nsd == 3
      [Z,iZ] = getPlane( cell2mat( C(:) ) );
    elseif nsd == 2
      Z = eye(4); iZ = eye(4);
    else
      error('Number of Spatial Dimensions must be 2 or 3.');
    end

    %clean and open all curves and collect coordinates
    X = [];
    F = [];
    for c = 1:numel(C)
      C{c}( all( ~diff( C{c} , 1 , 1 ) ,2) ,:) = [];
      if isequal( C{c}(1,:) , C{c}(end,:) )
        C{c}(end,:) = [];
      end
      F = [ F ; size(X,1) + [ ( 1:size(C{c},1) ).' , [ 2:size(C{c},1) , 1 ].' ] ];
      X = [ X ; C{c} ];
    end
    X(:,end+1:3) = 0;
    nX = size( X ,1);

    X0 = X;
    X = transform( X , iZ );
  
  elseif ~isnan( delta )
    
    X0 = C{1};
    X  = MeshFlatten( X0 );
    X(end,:) = [];
    X(:,end+1:3) = 0;
    
    nX = size( X ,1);
    F = [ 1:nX ; 2:nX , 1 ].';
    Z = [];
    
  end
  
  if isnan( delta )
    M.xyz = C{1};
    
    M.tri = [ 1:size( M.xyz ,1)-1 ; 2:size( M.xyz ,1) ].';
    
    M.xyz( end+1 ,:) = mean( M.xyz ,1);
    M.tri(:,3) = size( M.xyz ,1);
    
    return;
  end
  
  
  if isempty( delta )
    delta = median( sqrt( sum( diff(X,1,1).^2 ,2) ) );
  end
  if delta < 0
    delta = max( -delta , median( sqrt( sum( diff(X,1,1).^2 ,2) ) ) );
  end
  if isfinite( delta )
    xs = [ min( X(:,1) ) , max( X(:,1) ) ];
    xs = mean(xs) + ( -ceil( diff(xs)/2/delta + 1 ):ceil( diff(xs)/2/delta + 1 ) )*delta;
    
    delta = delta * sin(2*pi/3);
    ys = [ min( X(:,2) ) , max( X(:,2) ) ];
    ys = mean(ys) + ( -ceil( diff(ys)/2/delta + 2 ):ceil( diff(ys)/2/delta + 2 ) )*delta;
    
    Y = ndmat( ys(1:end-1) , xs ); Y = Y(:,[2,1]);
    Y(1:2:end,1) = Y(1:2:end,1) + mean(diff(xs))/2;
    Y(:,end+1:3) = 0;
  else
    Y = [];
  end
  Y = [ X ; Y ];
%   Y = unique( Y , 'rows','stable' );
  
  
  DT = delaunayTriangulation( double( Y(:,1) ) , double( Y(:,2) ) , double( F ) );
  M.xyz = DT.Points; M.xyz(:,end+1:3) = 0;
  M.tri = DT.ConnectivityList;
  M.tri = M.tri( isInterior(DT) ,:);
  
%   M.xyz( nX+1:end ,:) = NaN;
  M = MeshRelax( M );
  
  if size( M.xyz ,1) > nX  &&  any( abs( X(:,3) ) > 1e-10 )
    M.xyz( nX+1:end ,:) = InterpolatingSplines( M.xyz( 1:nX ,:) , X , M.xyz( nX+1:end ,:) , 'rlogr' );
  end
%   M.xyz( 1:nX ,:) = Y( 1:nX ,:);
  
  if ~isempty( Z )
    M = transform( M , Z );
  else
    M.xyz = InterpolatingSplines( M.xyz( 1:nX ,:) , X0( 1:nX ,:) , M.xyz , 'r' );
  end
  M.xyz( 1:nX ,:) = X0( 1:nX ,:);

end

