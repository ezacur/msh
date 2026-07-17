function [eID, cp, d, bc, F] = bvhClosestElement( M , P , B , Dmax )
%bvhClosestElement  Closest mesh element (and point on it) to query points.
%
%   [e,cp,d,bc,F] = bvhClosestElement( M , P )          builds the blob, queries
%   [e,cp,d,bc,F] = bvhClosestElement( M , P , B )      reuses a prebuilt BVH
%   [e,cp,d,bc,F] = bvhClosestElement( {M,B} , P )      same, mesh+blob bundled:
%                                    "on the FIRST argument, find the closest
%                                    element to the SECOND one"
%   [ ... ]       = bvhClosestElement( ... , Dmax )     search radius
%
%   OUTPUTS (nP rows each):
%     e    row of M.tri of the closest element
%     cp   closest point ON that element (nP x 3)
%     d    Euclidean distance |p - cp|
%     bc   barycentric coordinates of cp within its element, ROBUST (cross-
%          product form, exact even for slivers of aspect 1e10; clamped to
%          [0,1] and renormalized) -- NaN-padded to size(M.tri,2) columns
%     F    classification of WHERE cp lies (struct):
%            .type        0 = no answer
%                         1 = on a VERTEX of the element
%                         2 = on an EDGE (interior of it)
%                         3 = on a FACE interior (triangle / tet face)
%                         4 = INSIDE the element (tet interior)
%            .onBoundary  true if cp lies ON the OPEN boundary of the mesh
%                         (triangle surfaces: a MeshBoundary edge or one of its
%                         vertices; wireframes: a free end). Always false for
%                         closed meshes, tets, point clouds and mixed meshes.
%          The active bc pattern (bc > tol) tells WHICH vertex/edge it is.
%
%   Dmax (scalar or nP-VECTOR, default Inf): SEARCH RADIUS. The best-so-far
%   bound is seeded with Dmax, so everything farther prunes from the very root
%   -- points beyond Dmax cost ~one node visit. A vector seeds a PER-POINT
%   upper bound (e.g. the nearest-vertex heuristic, see bench_vertexSeed);
%   inflate a tight bound by (1+1e-9) if the element AT the bound must still
%   be found. NO ANSWER convention (beyond Dmax):
%     e = 0 , d = Inf , cp = NaN , bc = NaN , F.type = 0
%   (non-finite query points give the same with d = NaN).
%
%   THREADS: the MEX parallelizes over the points following maxNumCompThreads
%   -- '-singleCompThread' sessions (or maxNumCompThreads(1)) stay serial, and
%   maxNumCompThreads(k) caps the MEX at k threads.
%
%   BLOB + FRAME: B is SELF-CONTAINED (owns .X/.Tri in build-frame space and a
%   global similarity .frame). Queries run in build space -- P and Dmax are
%   transformed in, cp and d are transformed out; e/bc/F are invariant. M is
%   used for a cheap staleness spot-check (mismatch -> warning + rebuild).
%   2-D meshes (2-column vertices) are supported end to end; dimensionality
%   always comes from the column count, never from the values.
%
%   Requires the compiled MEXes (BVH_mx, bvhClosestElement_mx).
%
% See also BVH, bvhIntersectRay, plotBVH, tsearchn, MeshBoundary.

  %bundled form: bvhClosestElement( {M,B} , P , Dmax )
  if iscell( M )
    if nargin > 2, Dmax = B; end          %shift: 3rd arg is Dmax in this form
    B = M{2};  M = M{1};
  elseif nargin < 3
    B = [];
  end
  if ~exist('Dmax','var') || isempty( Dmax ), Dmax = Inf; end
  if ~( isnumeric(Dmax) && isreal(Dmax) && isvector(Dmax) && all( Dmax(:) >= 0 ) )
    error('bvhClosestElement:Dmax', ...
          'Dmax must be a nonnegative scalar or an nP-vector (per-point bound seed).');
  end
  Dmax = double( Dmax(:) );
  if exist( 'bvhClosestElement_mx' ,'file') ~= 3
    error('bvhClosestElement:mex','bvhClosestElement_mx is not compiled (mex COMPFLAGS="$COMPFLAGS /openmp" bvhClosestElement_mx.cpp).');
  end

  P = double( P );  P(:,end+1:3) = 0;
  nP = size( P ,1);
  if ~isscalar( Dmax ) && numel( Dmax ) ~= nP
    error('bvhClosestElement:Dmax','per-point Dmax must have nP elements (%d vs %d).', numel(Dmax) , nP );
  end

  if isempty( B ), B = BVH( M ); end

  %staleness spot-check: if frame(B.X) no longer matches M.xyz (mesh edited
  %without updating B), warn and rebuild. Cost: 4 vertices, ~microseconds.
  Xw = double( M.xyz ); Xw(:,end+1:3) = 0;
  ok = size( Xw ,1) == size( B.X ,1) && isequal( size( double(M.tri) ) , size( B.Tri ) );
  if ok && size( Xw ,1) > 0
    ii = unique( round( linspace( 1 , size(Xw,1) , 4 ) ) );
    Yw = B.X(ii,:) * B.frame(1:3,1:3).' + B.frame(1:3,4).';
    ok = max(max( abs( Yw - Xw(ii,:) ) )) <= 1e-6 * max( 1 , max(max( abs( Xw(ii,:) ) )) );
  end
  if ~ok
    %ERROR on purpose (not a silent rebuild: B is a value, so an in-call
    %rebuild would be discarded and repeated on EVERY subsequent call).
    %Recover explicitly:  B = BVH(B, M)  (refit, same connectivity)
    %or  B = BVH(M)  (rebuild).
    error('bvhClosestElement:staleBVH', ...
          'B does not match M (stale or foreign blob). Refit it (B = BVH(B,M)) or rebuild it (B = BVH(M)).');
  end

  %global frame: query in BUILD space, un-transform the outputs at the end.
  %For a similarity frame  world = A*xf + t  with  A.'*A = s^2*I :
  %  pf = (p - t)*inv(A).' ,  inv(A) = A.'/s^2 ;  distances scale by s.
  Fr   = B.frame;
  hasF = ~isequal( Fr , eye(4) );
  if hasF
    Af = Fr(1:3,1:3);  tf = Fr(1:3,4).';
    s2 = trace( Af.'*Af )/3;  fscale = sqrt( s2 );
    P  = ( P - tf ) * ( Af / s2 );
    DmaxEff = Dmax / fscale;
  else
    DmaxEff = Dmax;
  end

  [ eID , cp , d ] = bvhClosestElement_mx( P , B , maxNumCompThreads , DmaxEff );

  X   = B.X;
  Tri = B.Tri;

  %barycentric coordinates of cp within its element -- ROBUST forms (see
  %vtkClosestElement/calcular_barycentric): triangle weights as ratios of
  %CROSS-product areas (the Gram form d00*d11-d01^2 cancels catastrophically
  %on slivers), then clamped to >= 0 and renormalized to sum 1.
  if nargout > 3
    bc = NaN( nP , size( Tri ,2) );
    k  = sum( Tri > 0 ,2);
    ke = zeros( nP ,1);
    w  = eID > 0;  ke(w) = k( eID(w) );

    w = find( ke == 1 );
    if ~isempty( w ), bc(w,1) = 1; end

    w = find( ke == 2 );
    if ~isempty( w )
      A  = X( Tri(eID(w),1) ,:);  Bv = X( Tri(eID(w),2) ,:);
      ab = Bv - A;
      t  = sum( ( cp(w,:) - A ).*ab ,2) ./ max( sum( ab.^2 ,2) , realmin );
      t  = min( max( t ,0) ,1);
      bc(w,1:2) = [ 1-t , t ];
    end

    w = find( ke == 3 );
    if ~isempty( w )
      A  = X( Tri(eID(w),1) ,:);
      v0 = X( Tri(eID(w),2) ,:) - A;
      v1 = X( Tri(eID(w),3) ,:) - A;
      v2 = cp(w,:) - A;
      n  = cross( v0 , v1 , 2 );
      nn = max( sum( n.*n ,2) , realmin );
      g  = sum( n .* cross( v0 , v2 , 2 ) ,2) ./ nn;   %weight of vertex 3
      b  = sum( n .* cross( v2 , v1 , 2 ) ,2) ./ nn;   %weight of vertex 2
      W  = [ 1-b-g , b , g ];
      W  = max( W , 0 );  W = W ./ sum( W , 2 );
      bc(w,1:3) = W;
    end

    w = find( ke == 4 );
    if ~isempty( w )
      A  = X( Tri(eID(w),1) ,:);  Bv = X( Tri(eID(w),2) ,:);
      C  = X( Tri(eID(w),3) ,:);  D  = X( Tri(eID(w),4) ,:);
      d0 = detr( Bv-A , C-A , D-A );
      d0( d0 == 0 ) = realmin;
      W  = [ detr( Bv-cp(w,:) , C-cp(w,:) , D-cp(w,:) ) ./ d0 , ...
             detr( cp(w,:)-A  , C-A       , D-A       ) ./ d0 , ...
             detr( Bv-A       , cp(w,:)-A , D-A       ) ./ d0 ];
      W  = [ W , 1 - sum( W ,2) ];
      W  = max( W , 0 );  W = W ./ sum( W , 2 );
      bc(w,1:4) = W;
    end
  end

  %classification of the closest-point FEATURE (+ open-boundary flag)
  if nargout > 4
    tolF = 1e-9;
    nz   = sum( bc > tolF , 2 );                 %active barycentric weights
    F = struct();
    F.type = zeros( nP ,1);
    w = eID > 0;
    F.type(w) = min( nz(w) , 4 );                %1 vtx, 2 edge, 3 face, 4 inside
    F.onBoundary = false( nP ,1);

    kk = size( Tri ,2);
    if kk == 3 && all( k == 3 )                  %pure triangle surface
      Bed = MeshBoundary( Tri );
      if ~isempty( Bed )
        Bed  = sort( Bed ,2);
        Bvx  = unique( Bed );
        %on a vertex: the single active node must be a boundary vertex
        wv = find( F.type == 1 );
        if ~isempty( wv )
          [~,imax] = max( bc(wv,1:3) ,[],2);
          nd = Tri( sub2ind( size(Tri) , eID(wv) , imax ) );
          F.onBoundary(wv) = ismember( nd , Bvx );
        end
        %on an edge: the two active nodes must form a boundary edge
        we = find( F.type == 2 );
        if ~isempty( we )
          E2 = zeros( numel(we) ,2);
          for q = 1:numel( we )
            act = find( bc(we(q),1:3) > tolF );
            E2(q,:) = Tri( eID(we(q)) , act );
          end
          F.onBoundary(we) = ismember( sort(E2,2) , Bed , 'rows' );
        end
      end
    elseif kk == 2                               %wireframe: free ends
      fe = MeshBoundary( Tri );                  %degree-1 node ids
      if ~isempty( fe )
        wv = find( F.type == 1 );
        if ~isempty( wv )
          [~,imax] = max( bc(wv,1:2) ,[],2);
          nd = Tri( sub2ind( size(Tri) , eID(wv) , imax ) );
          F.onBoundary(wv) = ismember( nd , fe(:) );
        end
      end
    end
    %tets, point clouds and mixed meshes: onBoundary stays false (documented)
  end

  %un-transform: from build-frame space back to world (e/bc/F are invariant;
  %misses propagate: NaN cp stays NaN, Inf d stays Inf)
  if hasF
    cp = cp * Af.' + tf;
    d  = d  * fscale;
  end

end

function v = detr( a , b , c )
%rowwise determinant det([a;b;c]) = a . (b x c)
  v = sum( a .* cross( b , c ,2) ,2);
end
