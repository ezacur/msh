function M = MeshSubdivide( M , varargin )
%MESHSUBDIVIDE  Selective, conforming subdivision of segment/triangle/tetra meshes.
%
%   MM = MeshSubdivide( M )              subdivide EVERY face (default scheme)
%   MM = MeshSubdivide( M , W )          subdivide only faces W (indices or a logical
%                                        mask); neighbours are split the MINIMUM
%                                        needed to keep the mesh CONFORMING (no
%                                        T-junctions)
%   MM = MeshSubdivide( M , 'scheme' )   pick the subdivision scheme (see SCHEMES)
%   MM = MeshSubdivide( M , W , 'scheme' )
%   MM = MeshSubdivide( M , ... , 'kp' ) keep each new node's parent edge in
%                                        M.xyzParentEdge (alias 'keepparent')
%   MM = MeshSubdivide( M , -EL )        refine adaptively until every edge <= EL
%   MM = MeshSubdivide( M , -EL , 'scheme' )  adaptive refinement WITH a scheme
%                                        (e.g. -EL,'pn' refines onto the curved
%                                        PN geometry until edges <= EL)
%   MM = MeshSubdivide( {M,Nits} , ... ) apply the subdivision Nits times
%
% M is a mesh struct: .xyz (N-by-3 nodes) and .tri (faces). The arity of .tri
% selects the family via meshCelltype: segments (3), triangles (5), tetrahedra
% (10). Per-node fields (prefix 'xyz*') are interpolated at the new nodes;
% per-cell fields ('tri*') are inherited by the child cells. Where the new nodes
% land is what the scheme decides.
%
% SCHEMES
%   segments (celltype 3)
%     'default'      split each segment at its midpoint (the trivial linear
%                    stationary scheme).                                      [8]
%     'pn'           curved Point-Normal SEGMENTS: same connectivity as 'default'
%                    but the new node is lifted onto the cubic PN curve built
%                    from the endpoint positions + vertex normals (M.xyzNORMALS,
%                    else meshNormals(M,'quadratic')) -- the boundary cubic of
%                    the PN triangles [4]: 4th-order accurate on smooth curves
%                    (radial error 3*th^4/8 on a circle of half-angle-per-segment
%                    th, vs th^2/2 of the linear midpoint), invariant to the
%                    normals' sign, and straight (z=NaN-marked) pieces stay
%                    straight. If a normal field is carried, the new nodes get
%                    the PN quadratic normals (marks stay unit-xy + NaN).     [4]
%     'pn<k>'        e.g. 'pn3': 1 -> k split placing the k-1 new nodes on the
%                    SAME first-generation PN cubic at t = j/k. Unlike applying
%                    'pn' repeatedly, no interpolation error accumulates (every
%                    node lies on the generation-1 curve) and it is cheaper.
%                    'pn2' == 'pn'.                                           [4]
%     'cornercutting'  Chaikin corner cutting: every segment keeps its middle
%                    half (cuts at 1/4 and 3/4) and the interior corners are
%                    bridged. APPROXIMATING, C1: the limit is the quadratic
%                    B-spline of the original nodes (closed curves SHRINK,
%                    ~(3/8)*th^2 per pass). Interior (valence-2) nodes are
%                    REMOVED; open ends / branchings are kept and attached. [6,9]
%     '4points'      Dyn-Levin-Gregory INTERPOLATING 4-point: midpoint nodes
%                    from the (-1 9 9 -1)/16 stencil over the 4 chain nodes
%                    around each segment. 4th order on smooth UNIFORMLY sampled
%                    curves -- the same error constant 3*th^4/8 as 'pn', with
%                    NO normals needed -- but it degrades on non-uniform
%                    sampling (measured 47x worse than 'pn'); one-sided rules
%                    at open ends / branchings.                              [7]
%   Picking a segment scheme: interpolating + smooth-accurate -> 'pn'/'pn<k>'
%   (needs normals; robust to non-uniform sampling) or '4points' (no normals;
%   wants uniform sampling); corner rounding / smoothing -> 'cornercutting';
%   plain refinement -> 'default'. Every segment scheme accepts a subset W
%   except 'cornercutting' (it moves/removes ORIGINAL nodes, so a partial
%   application would crack the selection boundary; the interpolating ones
%   never move original nodes and split cleanly anywhere).
%   triangles (celltype 5)
%     'default'/'linear4'  conforming RED refinement: marked triangles -> 4,
%                    neighbours -> 2, at edge midpoints (piecewise linear).     [1]
%     'triangular'   as 'default' but forces every touched triangle fully to 4. [1]
%     'linear3'      1 -> 3 split at the face centroid (barycentric insertion). [8]
%     'linear9'      1 -> 9 TRIADIC split: every edge at its thirds plus the
%                    face centroid (piecewise linear). SELECTIVE: neighbours of
%                    the marked faces are fan-triangulated (3 or 5 children),
%                    conforming in a single pass (no closure loop needed).    [8]
%     'pn9'          curved 1 -> 9: the 'linear9' connectivity with every new
%                    node on the GENERATION-1 PN patch: edge nodes on the
%                    boundary cubics at t = 1/3 and 2/3, the centroid at the
%                    patch barycenter (== the mean of the 6 interior edge
%                    control points; b111 cancels). ONE pass matches {M,2} of
%                    'pn' but with no accumulation of interpolation error.
%                    Selective, like 'pn'. ('pn9' on a SEGMENT mesh is the
%                    1 -> 9 'pn<k>' split instead.)                           [4]
%     'loop'         Loop's approximating C2 subdivision (all faces).           [2]
%     'loop_matrix'  Loop, also returning the sparse operator in M.LoopMatrix.  [2]
%     'butterfly'    modified-Butterfly INTERPOLATING subdivision (all faces).  [3]
%     'sqrt3'        Kobbelt sqrt(3)-subdivision: centroid insertion + flip of
%                    every original interior edge + valence-mask relaxation of
%                    the original vertices. APPROXIMATING, C2 away from
%                    extraordinary vertices, and the slowest-growing scheme:
%                    1 -> 3 per pass (two passes = a realigned 3-adic 1 -> 9).
%                    Original edges vanish; centroids get valence 6; original
%                    vertices keep theirs. Boundary edges are not flipped and
%                    boundary vertices are kept (paper's exact boundary rule
%                    not implemented). Needs an edge-manifold consistently
%                    wound mesh (checked). Whole-mesh only.                 [11]
%     'pn'           curved Point-Normal (PN) triangles: linear connectivity, but
%                    new edge nodes are lifted onto the cubic PN Bezier from the
%                    vertex normals (M.xyzNORMALS, else meshNormals).           [4]
%     'fakebutterfly'  Loop subdivision followed by a thin-plate-spline
%                    reprojection through a FarthestPointSampling subset of the
%                    ORIGINAL nodes (so the result interpolates them).     [2,10]
%     'safebutterfly'  'butterfly', falling back to 'fakebutterfly' on error. [3]
%   tetrahedra (celltype 10)
%     'default'      conforming RED refinement 1 -> 8, with consistent diagonals
%                    (tie-broken by vertex index) across all partial-edge cases. [5]
%     'linear4'      1 -> 4 split at the cell centroid (barycentric insertion). [8]
%
%   NON-MANIFOLD triangle meshes (an edge shared by 3+ faces): the linear / PN
%   family ('default', 'triangular', 'linear3', 'linear4', 'linear9', 'pn',
%   'pn9') accepts them BY CONSTRUCTION -- the shared edge gets ONE set of
%   split nodes used by every incident face, so seams and multi-region
%   boundaries stay conforming. The smoothing / stencil schemes need a real
%   manifold and ERROR out: 'loop', 'loop_matrix' and 'butterfly' with
%   MeshSubdivide:nonManifold (before, 'loop' silently used only 2 of the 3+
%   wings -- order-dependent -- and 'butterfly' crashed inside accumarray),
%   'sqrt3' with MeshSubdivide:sqrt3NonManifold.
%
% OUTPUT: the subdivided mesh. Extra fields may appear: M.xyzParentEdge ('kp'),
% M.LoopMatrix ('loop_matrix').
%
% REFERENCES
%   [1] Bank, Sherman, Weiser, "Refinement algorithms and data structures for
%       regular local mesh refinement", 1983 (red/green refinement).
%   [2] C. Loop, "Smooth Subdivision Surfaces Based on Triangles", MSc thesis,
%       University of Utah, 1987.
%   [3] Zorin, Schroder, Sweldens, "Interpolating Subdivision for Meshes with
%       Arbitrary Topology", SIGGRAPH 1996 (modified Butterfly). The stencil
%       cases mirror VTK's vtkButterflySubdivisionFilter.
%   [4] Vlachos, Peters, Boyd, Mitchell, "Curved PN Triangles", ACM I3D 2001.
%   [5] J. Bey, "Tetrahedral grid refinement", Computing 55(4):355-378, 1995.
%   [6] G. Chaikin, "An algorithm for high-speed curve generation", Computer
%       Graphics and Image Processing 3(4):346-349, 1974.  (corner cutting)
%   [7] N. Dyn, D. Levin, J. A. Gregory, "A 4-point interpolatory subdivision
%       scheme for curve design", Computer Aided Geometric Design 4(4):257-268,
%       1987.  (the interpolating 4-point scheme)
%   [8] N. Dyn, D. Levin, "Subdivision schemes in geometric modelling", Acta
%       Numerica 11:73-144, 2002.  (survey; covers the plain midpoint /
%       barycentric linear refinements as the trivial stationary schemes)
%   [9] R. F. Riesenfeld, "On Chaikin's algorithm", Computer Graphics and Image
%       Processing 4(3):304-310, 1975.  (proves Chaikin's limit curve is the
%       uniform quadratic B-spline -- the basis of the shrink/limit claims)
%   [10] F. L. Bookstein, "Principal warps: thin-plate splines and the
%       decomposition of deformations", IEEE PAMI 11(6):567-585, 1989.  (the
%       thin-plate-spline interpolation behind the 'fakebutterfly' reprojection)
%   [11] L. Kobbelt, "sqrt(3)-Subdivision", SIGGRAPH 2000, pp. 103-112.  (the
%       centroid-insertion + edge-flip + valence-relaxation scheme)
%
% See also Mesh, meshCelltype, meshNormals, MeshQuality, MeshSubdividePN.

if 0

  M = Mesh(1:3,1:3);
  MeshSubdivide( M );
  MeshSubdivide( M , 1:10 );
  MeshSubdivide( M , 'linear' );
  MeshSubdivide( M , 'butterfly' );
  MeshSubdivide( M , 1:10 , 'butterfly' );
  MeshSubdivide( M , 'linear' , 1:10 ,'kp');
  
end


if 0
  M.xyz = randn(2000,3); M.tri = delaunayn( M.xyz ); M.celltype = 10; M.triL = (1:size(M.tri,1)).';
  MM = MeshSubdivide( M , 1:11:size(M.tri,1) );
%   plMESH( MM )
  w = MeshQuality( MM , 'volume' ) < 0; unique( MM.triCASE(w) )
  %%
end
if 0
  p1=1;p2=2;p3=3;p4=4;p12=5;p13=6;p14=7;p23=8;p24=9;p34=10;
  T=[];F=[];w=[];

  V = struct('tri',T,'xyz',[0,0,0;2,0,0;1,1.6,0;1,0.6,1.4;1,0,0;0.5,0.8,0;0.5,0.3,0.7;1.5,0.8,0;1.5,0.3,0.7;1,1.1,0.7],'celltype',10);
  MeshQuality( V , 'volume' )
  plMESH(V);
  %%  
end
if 0
%   M.xyz = rand(20,2); M.tri = delaunayn( M.xyz ); M.triL = (1:size(M.tri,1)).';
  W = [ 1 8 10 20 25 ]; W( W>size(M.tri,1) ) = []; cm = rand(size(M.tri,1),3)/2+0.5;
  subplot(121); plotMESH(  M , 'td','L' ); colormap(cm); colorbar
  MM = MeshSubdivide( M , W );
  subplot(122); plotMESH( MM , 'td','L' ); colormap(cm); colorbar
  hplot3d( getv( FacesCenter( M ) , W , ':' ) ,'*r' )
  hplotMESH( MeshBoundary(MM) , '-','EdgeColor','r','LineWidth',3)
  %%
end

  %adaptive mode: a NEGATIVE scalar among the arguments is the edge-length
  %target EL; every OTHER argument (a scheme like 'pn', 'kp') is forwarded to
  %each pass, so  MeshSubdivide( M , -EL , 'pn' )  refines CURVED until all
  %edges are <= EL.
  wneg = find( cellfun( @(a) isnumeric(a) && isscalar(a) && a < 0 , varargin ) , 1 );
  if ~isempty( wneg )
    EL = -varargin{wneg};
    varargin( wneg ) = [];
    for it = 1:100
      w = any( meshQuality( M , 'lengths' ) > EL ,2);
      if ~any( w ), break; end
      M = MeshSubdivide( M , w , varargin{:} );
    end
    return;
  end

  if isstruct( M ) && isempty( M.tri )
    return;
  end

  W       = Inf;
  SubType = '';
  KP      = false;

  for v = 1:numel(varargin)
    if ischar( varargin{v} ) && ...
       ( strcmpi( varargin{v} ,'kp' ) || strcmpi( varargin{v} ,'keepparent' ) || strcmpi( varargin{v} ,'keepparentedge' ) )
      KP = true; continue;
    end
    if ischar( varargin{v} ),                                              SubType = varargin{v}; continue; end
    if isnumeric( varargin{v} ),                                           W = varargin{v}; continue; end
    if islogical( varargin{v} ),                                           W = varargin{v}; continue; end
  end

  if strcmpi( SubType , 'fakebutterfly' )
    if isempty( W ), W = 1000; end

    T = M.xyz;
    [~,ids] = FarthestPointSampling( T , 1 , 0 , W );
    M = MeshSubdivide( M , 'loop' );
    M.xyz = InterpolatingSplines( T(ids,:) , M.xyz(ids,:) , M.xyz , 'r' );
    
    return;
  end
  if strcmpi( SubType , 'safebutterfly' )
    try
      M = MeshSubdivide( M , 'butterfly' );
    catch
      M = MeshSubdivide( M , 'fakebutterfly' , W );
    end
    return;
  end

  
  if iscell( M )
    Nits = M{2};
    M    = M{1};
    for it = 1:Nits
      M = MeshSubdivide( M , W , SubType );
    end
    return;
  end
  
  if ischar( W ) && isempty( W ), W = Inf; end
  nT  = size( M.tri , 1);   %number of faces
  if isinf( W )
    W = 1:nT;
  elseif islogical( W )
    if numel(W) ~= nT, error('a number of triangles logical were expected'); end
    W = find(W);
    if any( W > nT ),   error('Index exceeds number of faces.');        end
    if isempty( W ), return; end
  else
    W = unique( W(:) ,'sorted');
    if any( mod(W,1) ), error('indexes must be integers or logicals.'); end
    if any( W < 1 ),    error('indexes must be positive integers.');    end
    if any( W > nT ),   error('Index exceeds number of faces.');        end
    if isempty( W ), return; end
  end
  W = W(:);

  M.celltype = meshCelltype( M );
  PNK = regexp( lower( SubType ) , '^pn(\d+)$' , 'tokens' , 'once' );   %'pn3','pn4',... = 1->k PN split
  if ~isempty( PNK )
    PNK = str2double( PNK{1} );
    if PNK < 2
      error( 'MeshSubdivide:pnK' , '''pn%d'': the split factor must be at least 2.' , PNK );
    elseif PNK == 2
      SubType = 'pn';  PNK = [];        %'pn2' IS 'pn'
    elseif M.celltype == 5 && PNK == 9
      SubType = 'pn9';  PNK = [];       %triangles: the 1 -> 9 triadic PN split
    elseif M.celltype ~= 3
      error( 'MeshSubdivide:pnK' , ...
        '''pn%d'': ''pn<k>'' is for segment meshes (for triangles only ''pn9'', the 1 -> 9 triadic PN split, exists).' , PNK );
    else
      SubType = 'pnk';
    end
  else
    PNK = [];
  end
  if strcmpi( SubType , 'pn' ) && M.celltype ~= 5 && M.celltype ~= 3
    error( 'MeshSubdivide:pnTriOnly' , 'PN subdivision is only valid for segment or triangular meshes.' );
  end
  T  = M.tri;
  F  = ( 1:nT ).';         %face indexes
  nP = size( M.xyz , 1);   %number of points
  
  if isempty( SubType ), SubType = 'default'; end

  if 0
  elseif M.celltype == 3   %segments case
    switch lower(SubType)
      case {'default','pn','pnk'},
        % segment split. 'default': each selected segment -> two halves at the
        % chord midpoint [linear]. 'pn': same connectivity, but the new node is
        % lifted onto the cubic PN curve built from the endpoint positions +
        % vertex normals (M.xyzNORMALS, else meshNormals(M,'quadratic')) -- the
        % boundary cubic of the PN triangles [4] at t=1/2, 4th-order accurate
        % on smooth curves. 'pn<k>' ('pn3','pn4',...): 1 -> k split placing the
        % k-1 new nodes on that SAME cubic at t = j/k; every new node lies on
        % the FIRST-generation curve, so unlike iterating 'pn' k times there is
        % no accumulation of interpolation error (and it is cheaper).
        if isempty( W ), return; end

        isPN = strcmpi( SubType , 'pn' ) || ~isempty( PNK );

        %PN: capture the ORIGINAL vertex normals before new nodes are appended
        if isPN
          if isfield( M , 'xyzNORMALS' ) && ~isempty( M.xyzNORMALS )
            NRM = M.xyzNORMALS;
            if size( NRM ,1) ~= nP
              error( 'MeshSubdivide:xyzNORMALS' , ...
                'M.xyzNORMALS has %d rows but the mesh has %d nodes: stale field from another mesh?' , size(NRM,1) , nP );
            end
          else
            NRM = meshNormals( M , 'quadratic' );    %the curve-accurate vertex normals
          end
        end

        %new points on the edges (at the midpoint, or at every t = j/k for
        %'pn<k>'), every xyz* field interpolated linearly at its t
        E   = T( W ,:);
        m   = size( E ,1);
        if isempty( PNK ), TS = 0.5; else, TS = ( 1:PNK-1 ) / PNK; end
        if isempty( PNK )
          for f = fieldnames( M ).', if ~strncmp( f{1} , 'xyz' , 3 ), continue; end
            M.(f{1}) = [ M.(f{1}) ; ( M.(f{1})( E(:,1) ,:,:,:,:) + M.(f{1})( E(:,2) ,:,:,:,:) )/2 ];
          end
        else
          for f = fieldnames( M ).', if ~strncmp( f{1} , 'xyz' , 3 ), continue; end
            F1 = M.(f{1})( E(:,1) ,:,:,:,:);
            F2 = M.(f{1})( E(:,2) ,:,:,:,:);
            for t = TS
              M.(f{1}) = [ M.(f{1}) ; (1-t)*F1 + t*F2 ];
            end
          end
        end
        P = nP + ( 1:numel(TS)*m ).';          %new point ids: one m-block per t
        if KP
          fn = fieldnames(M); fn = sort( fn( strncmp( fn , 'xyzParentEdge' ,13) ) );
          for f = fn(end:-1:1).', M = renameStructField( M , f{1} , [ f{1} , '_' ] ); end
          PE = [ (1:nP).' , zeros(nP,2) ];
          for t = TS
            PE = [ PE ; [ double(E) , t*ones(m,1) ] ];
          end
          M.xyzParentEdge = PE;
        end

        %PN: lift the new nodes onto the cubic PN curve. A term with a
        %NON-FINITE normal (the z=NaN marks of straight pieces, NaN normals off
        %degenerate cells) is DROPPED, not propagated: the true normal of a
        %straight piece is perpendicular to its chord, so its lift is exactly
        %zero anyway (and the Bezier of the resulting uniform collinear control
        %polygon is EXACTLY the linear t = j/k point).
        if isPN
          P1 = M.xyz( E(:,1) ,:);   P2 = M.xyz( E(:,2) ,:);
          N1 = NRM( E(:,1) ,:);     N2 = NRM( E(:,2) ,:);
          dd = P2 - P1;
          if isempty( PNK )
            %closed form at t = 1/2 (== the Bezier below evaluated there)
            t1 = sum(  dd .* N1 ,2) .* N1;   t1( ~isfinite( t1 ) ) = 0;
            t2 = sum( -dd .* N2 ,2) .* N2;   t2( ~isfinite( t2 ) ) = 0;
            M.xyz( P ,:) = ( P1 + P2 )/2 - ( t1 + t2 )/8;
          else
            %Vlachos edge control points b1 = (2P1+P2)/3 - ((P2-P1).N1)N1/3 and
            %b2 = (P1+2P2)/3 - ((P1-P2).N2)N2/3, then plain cubic Bezier at t
            c1 = sum(  dd .* N1 ,2) .* N1 / 3;   c1( ~isfinite( c1 ) ) = 0;
            c2 = sum( -dd .* N2 ,2) .* N2 / 3;   c2( ~isfinite( c2 ) ) = 0;
            b1 = ( 2*P1 + P2 )/3 - c1;
            b2 = ( P1 + 2*P2 )/3 - c2;
            for j = 1:numel( TS )
              t = TS(j);
              M.xyz( P( (j-1)*m + (1:m) ) ,:) = (1-t)^3*P1 + 3*(1-t)^2*t*b1 + 3*(1-t)*t^2*b2 + t^3*P2;
            end
          end

          %PN quadratic normal field at each t (only when a normal field is
          %carried). The final normalization counts NaN as 0, so the z=NaN
          %marks of straight pieces come out with a UNIT finite part, like
          %meshNormals returns them.
          if isfield( M , 'xyzNORMALS' ) && ~isempty( M.xyzNORMALS )
            NS   = N1 + N2;
            v    = 2 * sum( dd .* NS ,2 ) ./ sum( dd .* dd ,2 );
            n110 = NS - v .* dd;
            w = ~( sum( n110.^2 ,2 ) > eps );   n110(w,:) = NS(w,:);   %degenerate/marked -> linear
            n110 = rowUnit( n110 );
            for j = 1:numel( TS )
              t  = TS(j);
              Nm = (1-t)^2*N1 + 2*(1-t)*t*n110 + t^2*N2;
              w = ~( sum( Nm.^2 ,2 ) > eps );   Nm(w,:) = NS(w,:);     %antipodal -> linear
              Z  = Nm;  Z( isnan( Z ) ) = 0;
              nn = sqrt( sum( Z.^2 ,2) );  k = nn > 0;
              Nm(k,:) = Nm(k,:) ./ nn(k);
              M.xyzNORMALS( P( (j-1)*m + (1:m) ) ,:) = Nm;
            end
          end
        end

        %old and new faces: the chain  a - q1 - ... - q_{k-1} - b  per segment
        if isempty( PNK )
          T = [   T ; ...
                [ T( W ,1) , P        ] ;...
                [ P        , T( W ,2) ] ;...
              ];
          F = [ F ; W ; W ];      %original ids of the new faces
        else
          Q = reshape( P , m , [] );              %Q(i,j): node of segment i at t=j/k
          T = [ T ; [ T( W ,1) , Q(:,1) ] ];      F = [ F ; W ];
          for j = 1:PNK-2
            T = [ T ; [ Q(:,j) , Q(:,j+1) ] ];    F = [ F ; W ];
          end
          T = [ T ; [ Q(:,end) , T( W ,2) ] ];    F = [ F ; W ];
        end

      case {'cornercutting','chaikin','cc'}
        % Chaikin corner cutting [6]: every segment [a,b] is replaced by its
        % middle half [qa,qb], qa = (3a+b)/4, qb = (a+3b)/4, and every interior
        % (valence-2) node is CUT AWAY, its corner bridged by a new segment
        % joining the two cut points around it. APPROXIMATING, C1: iterating
        % converges to the uniform quadratic B-spline of the original nodes
        % (closed curves SHRINK: ~(3/8)*th^2 radial deficit per pass on a
        % circle of half-angle-per-segment th). Nodes of valence ~= 2 (open
        % ends, branchings, isolated points) are KEPT and joined to their
        % adjacent cut points: open curves keep their endpoints, branched
        % networks stay attached. This is the only scheme that REMOVES nodes
        % from the output. Every xyz* field gets the same 3/4-1/4 weights; the
        % corner bridges inherit the tri* fields of the segment ARRIVING at
        % the corner (or the lower-index one if the two orientations do not
        % chain). Whole-mesh only.
        if ~isequal( W(:).' , 1:nT ), error('cornercutting subdivision is only valid when all faces are considered.'); end

        val  = accumarray( double( T(:) ) , 1 , [nP,1] );
        keep = val ~= 2;                             %kept nodes: ends/branchings/isolated

        %the two cut points of every segment (all xyz* fields, 3/4-1/4)
        for f = fieldnames( M ).', if ~strncmp( f{1} , 'xyz' , 3 ), continue; end
          A_ = M.(f{1})( T(:,1) ,:,:,:,:);
          B_ = M.(f{1})( T(:,2) ,:,:,:,:);
          M.(f{1}) = [ M.(f{1}) ; 0.75*A_ + 0.25*B_ ; 0.25*A_ + 0.75*B_ ];
        end
        QA = nP +      ( 1:nT ).';                   %cut point near T(:,1)
        QB = nP + nT + ( 1:nT ).';                   %cut point near T(:,2)

        %children: the middle half of every segment...
        Tn = [ QA , QB ];                            Fn = ( 1:nT ).';
        %...the stubs that keep the valence~=2 nodes attached...
        w  = keep( T(:,1) );
        Tn = [ Tn ; [ double(T(w,1)) , QA(w) ] ];    Fn = [ Fn ; find(w) ];
        w  = keep( T(:,2) );
        Tn = [ Tn ; [ QB(w) , double(T(w,2)) ] ];    Fn = [ Fn ; find(w) ];
        %...and the corner bridges at every cut node: from the cut point of the
        %segment ARRIVING at it to that of the segment LEAVING it (flow order;
        %if the two orientations do not chain, lower segment index first)
        inc = [ [ double(T(:,1)) , (1:nT).' , QA , 2*ones(nT,1) ] ;...  %the segment LEAVES this node
                [ double(T(:,2)) , (1:nT).' , QB ,   ones(nT,1) ] ];    %the segment ARRIVES here (sorts first)
        inc = inc( ~keep( inc(:,1) ) ,:);
        inc = sortrows( inc , [1 4 2] );             %node | arriving-first | segment
        r1  = inc( 1:2:end ,:);   r2 = inc( 2:2:end ,:);
        Tn  = [ Tn ; [ r1(:,3) , r2(:,3) ] ];        Fn = [ Fn ; r1(:,2) ];

        %drop the cut original nodes and renumber
        ids = [ find( keep ) ; QA ; QB ];
        map = zeros( nP + 2*nT , 1 );   map( ids ) = 1:numel( ids );
        for f = fieldnames( M ).', if ~strncmp( f{1} , 'xyz' , 3 ), continue; end
          M.(f{1}) = M.(f{1})( ids ,:,:,:,:);
        end
        if KP
          fn = fieldnames(M); fn = sort( fn( strncmp( fn , 'xyzParentEdge' ,13) ) );
          for f = fn(end:-1:1).', M = renameStructField( M , f{1} , [ f{1} , '_' ] ); end
          M.xyzParentEdge = [ [ find(keep) , zeros(nnz(keep),2) ] ;...
                              [ double(T) , 0.25*ones(nT,1) ]     ;...
                              [ double(T) , 0.75*ones(nT,1) ]     ];
        end
        T = [ T*0 ; map( Tn ) ];      %the original faces are dummies, removed below
        F = [ F ; Fn ];

      case {'4points','4point','dlg'}
        % Dyn-Levin-Gregory interpolating 4-point scheme [7]: same connectivity
        % as 'default' (original nodes KEPT, midpoint nodes inserted) but each
        % new node uses the 4-point stencil  q = ( 9(Pa+Pb) - (Pl+Pr) )/16 ,
        % Pl / Pr being the chain neighbours BEYOND a and beyond b. C1 limit,
        % reproduces cubics: 4th order on smooth UNIFORMLY-sampled curves --
        % radial error 3*th^4/8 on a circle, the SAME constant as 'pn', but
        % needing NO normals. On non-uniform sampling it degrades (measured
        % 47x worse than 'pn' on a random-step circle: the stencil assumes a
        % uniform parametrization) -- prefer 'pn' there. Where the neighbour
        % beyond one side does not exist or is ambiguous (open END / BRANCHING,
        % valence ~= 2) the one-sided quadratic rule (3Pe + 6Po - Pn)/8 is
        % used; with no neighbour on either side, the plain midpoint. Straight
        % polylines reproduce exactly. xyz* fields other than xyz are
        % interpolated LINEARLY at the midpoint (data should not overshoot);
        % xyzNORMALS, if present, is linear too. ACCEPTS a subset W: the scheme
        % is interpolating (original nodes never move), so splitting only some
        % segments stays consistent -- the stencil reads its neighbours from
        % the whole mesh whether they are selected or not.
        if isempty( W ), return; end

        %the chain neighbour beyond each endpoint: the OTHER segment at a
        %valence-2 node, and its far node (0 = no unique neighbour)
        val = accumarray( double( T(:) ) , 1 , [nP,1] );
        inc = [ [ double(T(:,1)) , (1:nT).' ] ; [ double(T(:,2)) , (1:nT).' ] ];
        inc = inc( val( inc(:,1) ) == 2 ,:);
        inc = sortrows( inc , [1 2] );
        OS  = sparse( [ inc(1:2:end,1) ; inc(2:2:end,1) ] , [ inc(1:2:end,2) ; inc(2:2:end,2) ] , ...
                      [ inc(2:2:end,2) ; inc(1:2:end,2) ] , nP , nT );  %(node,segment) -> the other segment there
        sa = full( OS( double(T(:,1)) + nP*( (1:nT).' - 1 ) ) );
        sb = full( OS( double(T(:,2)) + nP*( (1:nT).' - 1 ) ) );
        wl = sa > 0;   L = zeros( nT ,1);
        L(wl) = double( T(sa(wl),1) ) + double( T(sa(wl),2) ) - double( T(wl,1) );
        wr = sb > 0;   R = zeros( nT ,1);
        R(wr) = double( T(sb(wr),1) ) + double( T(sb(wr),2) ) - double( T(wr,2) );
        wl = wl( W );   L = L( W );          %restrict to the SELECTED segments
        wr = wr( W );   R = R( W );          %(the lookup itself used the whole mesh)

        %interpolate the xyz* fields at the midpoints (linear), as in 'default'
        E = T( W ,:);
        for f = fieldnames( M ).', if ~strncmp( f{1} , 'xyz' , 3 ), continue; end
          M.(f{1}) = [ M.(f{1}) ; ( M.(f{1})( E(:,1) ,:,:,:,:) + M.(f{1})( E(:,2) ,:,:,:,:) )/2 ];
        end
        P = nP + ( 1:numel(W) ).';
        if KP
          fn = fieldnames(M); fn = sort( fn( strncmp( fn , 'xyzParentEdge' ,13) ) );
          for f = fn(end:-1:1).', M = renameStructField( M , f{1} , [ f{1} , '_' ] ); end
          PE = (1:nP).'; PE(:,2) = 0;
          PE = [ PE ; double(E) ];
          PE( 1:nP ,3) = 0;
          PE( nP+1:end ,3) = 0.5;
          M.xyzParentEdge = PE;
        end

        %...but the GEOMETRY gets the 4-point stencil
        Pa = M.xyz( E(:,1) ,:);   Pb = M.xyz( E(:,2) ,:);
        Q  = ( Pa + Pb )/2;                          %no neighbour on either side
        w  = wl & wr;
        Q(w,:) = ( 9*( Pa(w,:) + Pb(w,:) ) - ( M.xyz(L(w),:) + M.xyz(R(w),:) ) )/16;
        w  = ~wl & wr;                               %open end / branching at a
        Q(w,:) = ( 3*Pa(w,:) + 6*Pb(w,:) - M.xyz(R(w),:) )/8;
        w  = wl & ~wr;                               %open end / branching at b
        Q(w,:) = ( 3*Pb(w,:) + 6*Pa(w,:) - M.xyz(L(w),:) )/8;
        M.xyz( P ,:) = Q;

        T = [   T ; ...
              [ T( W ,1) , P        ] ;...
              [ P        , T( W ,2) ] ;...
            ];
        F = [ F ; W ; W ];

      otherwise
        error( 'MeshSubdivide:segScheme' , ...
               'unknown SubType ''%s'' for segments (use ''default'', ''pn'', ''pn<k>'', ''cornercutting'' or ''4points'').' , SubType );

    end
    
  elseif M.celltype == 5          %% triangular mesh
    switch lower(SubType)
      case {'default','triangular','linear4','pn'}
        % RED (1->4) refinement at edge midpoints. 'default'/'linear4' keep the
        % split minimal & conforming (neighbours -> 2); 'triangular' forces every
        % touched triangle to 4; 'pn' additionally lifts the new nodes onto the
        % curved PN surface.  Refs [1] (red refinement), [4] (Vlachos PN).
        if isempty( W ), return; end

        allE = sort( [ T(:,[1 2]) ; T(:,[2 3]) ; T(:,[1 3]) ] ,2); %no remove the repeated

        if strcmpi( SubType , 'triangular' )   %% case 1, in case of cuadrilateral faces, split the whole triengle into 4

          while 1
            E = allE( [ W ; W + nT ; W + 2*nT ] , : );
            E = unique( E , 'rows' );
            TtD = reshape( ismember( allE , E , 'rows' ) , nT , 3 );
            Wp = W;
            W = find( sum( TtD ,2) > 1 );
            if isequal( Wp , W ), break; end
          end

        elseif strcmpi( SubType , 'default' ) ||...
               strcmpi( SubType , 'linear4' ) ||...
               strcmpi( SubType , 'pn' )   %% case 2, keep to minimum the number of faces

          E = allE( [ W ; W + nT ; W + 2*nT ] , : );
          E = unique( E , 'rows' );

        end

        
        %triangles containing edges
        [~,ET] = ismember2ROWS( allE , E(:,1:2) );
        ET = reshape( ET , [nT,3] );
        LET = ~~ET;  %logical ET

        %original faces, to be removed at the end
        W = find( any( LET ,2) );

        %indexes of the new points
        P = ( (nP+1):(nP+size(E,1)) ).';

        %PN: capture the ORIGINAL vertex normals before new nodes are appended
        if strcmpi( SubType , 'pn' )
          if isfield( M , 'xyzNORMALS' ) && ~isempty( M.xyzNORMALS )
            NRM = M.xyzNORMALS;
            if size( NRM ,1) ~= nP
              error( 'MeshSubdivide:xyzNORMALS' , ...
                'M.xyzNORMALS has %d rows but the mesh has %d nodes: stale field from another mesh?' , size(NRM,1) , nP );
            end
          else
            NRM = meshNormals( M , 'angle' );
          end
        end

        %middle points on edges
        for f = fieldnames( M ).'
          if ~strncmp( f{1} , 'xyz' , 3 ), continue; end
          M.(f{1}) = [ M.(f{1}) ; ( M.(f{1})( E(:,1) ,:,:,:,:,:) + M.(f{1})( E(:,2) ,:,:,:,:,:) )/2 ];
        end
        if KP
          fn = fieldnames(M); fn = sort( fn( strncmp( fn , 'xyzParentEdge' ,13) ) );
          for f = fn(end:-1:1).', M = renameStructField( M , f{1} , [ f{1} , '_' ] ); end
          PE = (1:nP).'; PE(:,2) = 0;
          PE = [ PE ; double(E) ];
          PE( 1:nP ,3) = 0;
          PE( nP+1:end ,3) = 0.5;
          M.xyzParentEdge = PE;
        end

        %PN: lift each new edge node onto the cubic PN Bezier boundary at t=1/2:
        %  Pmid = (P1+P2)/2 - ( ((P2-P1).N1)N1 + ((P1-P2).N2)N2 )/8   (Vlachos et al. [4])
        %  (this is the exact t=1/2 point of the PN boundary curve b300..b030).
        if strcmpi( SubType , 'pn' )
          P1 = M.xyz( E(:,1) ,:);   P2 = M.xyz( E(:,2) ,:);
          N1 = NRM( E(:,1) ,:);     N2 = NRM( E(:,2) ,:);
          dd  = P2 - P1;
          w12 = sum(  dd .* N1 , 2 );        %(P2-P1).N1
          w21 = sum( -dd .* N2 , 2 );        %(P1-P2).N2
          M.xyz( P ,:) = ( P1 + P2 )/2 - ( w12 .* N1 + w21 .* N2 )/8;

          %PN quadratic normal field: N(1/2)=unit(N1/4 + n110/2 + N2/4), with the
          %geometry-coupled control normal n110=unit(N1+N2 - v*(P2-P1)),
          %v=2<P2-P1,N1+N2>/<P2-P1,P2-P1>. Only when a normal field is carried.
          if isfield( M , 'xyzNORMALS' ) && ~isempty( M.xyzNORMALS )
            NS   = N1 + N2;
            v    = 2 * sum( dd .* NS ,2 ) ./ sum( dd .* dd ,2 );
            n110 = NS - v .* dd;
            w = sum( n110.^2 ,2 ) <= eps;   n110(w,:) = NS(w,:);   %degenerate -> linear
            n110 = rowUnit( n110 );
            Nm   = 0.25*N1 + 0.5*n110 + 0.25*N2;
            w = sum( Nm.^2 ,2 ) <= eps;     Nm(w,:)   = NS(w,:);   %antipodal -> linear
            M.xyzNORMALS( P ,:) = rowUnit( Nm );
          end
        end

        %first, triangles to be divided into 4.
        w = find( LET(:,1) & LET(:,2) & LET(:,3) );
        T = [ T ; ...
              [     T(w, 1 )   , P( ET(w, 1 ) ) , P( ET(w, 3 ) ) ] ;...
              [ P( ET(w, 1 ) ) ,     T(w, 2 )   , P( ET(w, 2 ) ) ] ;...
              [ P( ET(w, 3 ) ) , P( ET(w, 2 ) ) ,     T(w, 3 )   ] ;...
              [ P( ET(w, 1 ) ) , P( ET(w, 2 ) ) , P( ET(w, 3 ) ) ] ;...
            ];
        F = [ F ; w ; w ; w ; w ];

        %triangles to be divided at edge 1-2
        w = find( LET(:,1) & ~LET(:,2) & ~LET(:,3) );
        T = [ T ; ...
              [     T(w, 1 )   , P( ET(w, 1 ) ) ,     T(w, 3 )   ] ;...
              [ P( ET(w, 1 ) ) ,     T(w, 2 )   ,     T(w, 3 )   ] ;...
            ];
        F = [ F ; w ; w ];

        %triangles to be divided at edge 2-3
        w = find( ~LET(:,1) & LET(:,2) & ~LET(:,3) );
        T = [ T ; ...
              [     T(w, 1 )   ,     T(w, 2 )   , P( ET(w, 2 ) ) ] ;...
              [     T(w, 1 )   , P( ET(w, 2 ) ) ,     T(w, 3 )   ] ;...
            ];
        F = [ F ; w ; w ];

        %triangles to be divided at edge 1-3
        w = find( ~LET(:,1) & ~LET(:,2) & LET(:,3) );
        T = [ T ; ...
              [     T(w, 1 )   ,     T(w, 2 )   , P( ET(w, 3 ) ) ] ;...
              [     T(w, 2 )   ,     T(w, 3 )   , P( ET(w, 3 ) ) ] ;...
            ];
        F = [ F ; w ; w ];

        %triangles to be divided at edge 1-2 && 2-3
        w = find( LET(:,1) & LET(:,2) & ~LET(:,3) );
        T = [ T ; ...
              [ P( ET(w, 1 ) ) ,     T(w, 2 )   , P( ET(w, 2 ) ) ] ;...
              [     T(w, 1 )   , P( ET(w, 1 ) ) ,     T(w, 3 )   ] ;...
              [ P( ET(w, 1 ) ) , P( ET(w, 2 ) ) ,     T(w, 3 )   ] ;...
            ];
        F = [ F ; w ; w ; w ];

        %triangles to be divided at edge 1-2 && 1-3
        w = find( LET(:,1) & ~LET(:,2) & LET(:,3) );
        T = [ T ; ...
              [     T(w, 1 )   , P( ET(w, 1 ) ) , P( ET(w, 3 ) ) ] ;...
              [ P( ET(w, 1 ) ) ,     T(w, 3 )   , P( ET(w, 3 ) ) ] ;...
              [ P( ET(w, 1 ) ) ,     T(w, 2 )   ,     T(w, 3 )   ] ;...
            ];
        F = [ F ; w ; w ; w ];

        %triangles to be divided at edge 2-3 && 1-3
        w = find( ~LET(:,1) & LET(:,2) & LET(:,3) );
        T = [ T ; ...
              [ P( ET(w, 3 ) ) , P( ET(w, 2 ) ) ,     T(w, 3 )   ] ;...
              [     T(w, 1 )   , P( ET(w, 2 ) ) , P( ET(w, 3 ) ) ] ;...
              [     T(w, 1 )   ,     T(w, 2 )   , P( ET(w, 2 ) ) ] ;...
            ];
        F = [ F ; w ; w ; w ];
    
      case {'linear3'}
        % 1 -> 3 split: insert the face centroid and fan to the 3 vertices. [linear]
        if isempty( W ), return; end

        %middle points on faces
        for f = fieldnames( M ).', if ~strncmp( f{1} , 'xyz' , 3 ), continue; end
          M.(f{1}) = [ M.(f{1}) ; ( M.(f{1})( T(W,1) ,:) + M.(f{1})( T(W,2) ,:) + M.(f{1})( T(W,3) ,:) )/3 ];
        end
        P = nP + ( 1:numel(W) ); P = P(:);    %indexes of the new points

        T = [ T ; ...
              [  P         ,  T(W, 2 )  ,  T(W, 3 )  ] ;...
              [  T(W, 1 )  ,  P         ,  T(W, 3 )  ] ;...
              [  T(W, 1 )  ,  T(W, 2 )  ,  P         ] ;...
            ];
        F = [ F ; W ; W ; W ];
        
      case {'linear9','pn9'}
        % TRIADIC 1 -> 9 refinement: every edge of the marked faces is split at
        % its THIRDS (two nodes per edge) and each fully-split face gains its
        % CENTROID, giving 9 children. SELECTIVE and conforming in a single
        % pass: a neighbour sharing 1 (or 2) split edges is fan-triangulated
        % into 3 (or 5) children, which induces NO new edge splits -- so there
        % is no closure loop -- and the two nodes of a shared edge are built
        % from the edge in CANONICAL (sorted) orientation, so both sides agree.
        % 'linear9' places the new nodes linearly [8]. 'pn9' lifts them onto
        % the GENERATION-1 PN patch [4]: edge nodes on the boundary cubics at
        % t = 1/3 and 2/3, and the centroid at the patch barycenter -- which
        % collapses to the MEAN of the 6 interior edge control points, since in
        % patch(1/3,1/3,1/3) = V/9 + 6*E/9 + 2*b111/9, with b111 = E + (E-V)/2,
        % the corner average V cancels. Hence ONE 'pn9' pass replaces {M,2} of
        % 'pn' with NO accumulation of interpolation error: every node lies on
        % the original patch. Non-finite normal TERMS are dropped, as in 'pn';
        % a carried M.xyzNORMALS gets the quadratic normal patch at the same
        % parameters. 'kp' is NOT supported here (a centroid has no parent
        % edge). ('pn9' on a SEGMENT mesh is the 1 -> 9 'pn<k>' split instead.)
        if isempty( W ), return; end

        isPN9 = strcmpi( SubType , 'pn9' );
        if isPN9      %capture the ORIGINAL vertex normals before any append
          if isfield( M , 'xyzNORMALS' ) && ~isempty( M.xyzNORMALS )
            NRM = M.xyzNORMALS;
            if size( NRM ,1) ~= nP
              error( 'MeshSubdivide:xyzNORMALS' , ...
                'M.xyzNORMALS has %d rows but the mesh has %d nodes: stale field from another mesh?' , size(NRM,1) , nP );
            end
          else
            NRM = meshNormals( M , 'angle' );
          end
        end

        allE = sort( [ T(:,[1 2]) ; T(:,[2 3]) ; T(:,[1 3]) ] ,2);
        E  = allE( [ W ; W + nT ; W + 2*nT ] , : );
        E  = unique( E , 'rows' );
        nE = size( E ,1);

        [~,ET] = ismember2ROWS( allE , E(:,1:2) );
        ET  = reshape( ET , [nT,3] );
        LET = ~~ET;
        W   = find( any( LET ,2) );          %every touched face gets replaced

        %two nodes per split edge, at 1/3 and 2/3 FROM THE LOWER vertex (so
        %both adjacent faces see the same pair); every xyz* field interpolated
        PA = ( (nP+1):(nP+nE) ).';
        PB = PA + nE;
        for f = fieldnames( M ).', if ~strncmp( f{1} , 'xyz' , 3 ), continue; end
          F1 = M.(f{1})( E(:,1) ,:,:,:,:,:);
          F2 = M.(f{1})( E(:,2) ,:,:,:,:,:);
          M.(f{1}) = [ M.(f{1}) ; (2/3)*F1 + (1/3)*F2 ; (1/3)*F1 + (2/3)*F2 ];
        end

        %the faces with their 3 edges split get 9 children and a CENTROID node
        w9 = find( all( LET ,2) );
        PC = ( nP + 2*nE ) + ( 1:numel(w9) ).';
        for f = fieldnames( M ).', if ~strncmp( f{1} , 'xyz' , 3 ), continue; end
          M.(f{1}) = [ M.(f{1}) ; ( M.(f{1})( T(w9,1) ,:,:,:,:,:) + M.(f{1})( T(w9,2) ,:,:,:,:,:) + M.(f{1})( T(w9,3) ,:,:,:,:,:) )/3 ];
        end

        %'pn9': lift the new nodes onto the generation-1 PN patch
        if isPN9
          Pu = M.xyz( E(:,1) ,:);   Pv = M.xyz( E(:,2) ,:);
          Nu = NRM( E(:,1) ,:);     Nv = NRM( E(:,2) ,:);
          dd = Pv - Pu;
          cu = sum(  dd .* Nu ,2) .* Nu / 3;   cu( ~isfinite( cu ) ) = 0;
          cv = sum( -dd .* Nv ,2) .* Nv / 3;   cv( ~isfinite( cv ) ) = 0;
          g1 = ( 2*Pu + Pv )/3 - cu;           %the two interior Bezier control
          g2 = ( Pu + 2*Pv )/3 - cv;           %points of each edge cubic
          M.xyz( PA ,:) = ( 8*Pu + 12*g1 +  6*g2 +   Pv )/27;    %cubic at t=1/3
          M.xyz( PB ,:) = (   Pu +  6*g1 + 12*g2 + 8*Pv )/27;    %cubic at t=2/3
          if ~isempty( w9 )
            i1 = ET(w9,1);   i2 = ET(w9,2);   i3 = ET(w9,3);
            M.xyz( PC ,:) = ( g1(i1,:) + g2(i1,:) + g1(i2,:) + g2(i2,:) + g1(i3,:) + g2(i3,:) )/6;
          end

          %quadratic normal patch for a carried normal field (NaN counts as 0
          %in the final normalization, as everywhere else)
          if isfield( M , 'xyzNORMALS' ) && ~isempty( M.xyzNORMALS )
            NS   = Nu + Nv;
            v    = 2 * sum( dd .* NS ,2 ) ./ sum( dd .* dd ,2 );
            n110 = NS - v .* dd;
            w = ~( sum( n110.^2 ,2 ) > eps );   n110(w,:) = NS(w,:);
            n110 = rowUnit( n110 );
            NA = ( 4*Nu + 4*n110 +   Nv )/9;   %the quadratic field at t=1/3
            NB = (   Nu + 4*n110 + 4*Nv )/9;   %and at t=2/3
            if isempty( w9 ),  NC = zeros( 0 , size(NRM,2) );
            else               %the quadratic patch at the barycenter
              NC = ( NRM( T(w9,1) ,:) + NRM( T(w9,2) ,:) + NRM( T(w9,3) ,:) ...
                   + 2*n110(i1,:) + 2*n110(i2,:) + 2*n110(i3,:) )/9;
            end
            NN = [ NA ; NB ; NC ];
            Z  = NN;  Z( isnan( Z ) ) = 0;
            nn = sqrt( sum( Z.^2 ,2) );  k = nn > 0;
            NN(k,:) = NN(k,:) ./ nn(k);
            M.xyzNORMALS( [ PA ; PB ; PC ] ,:) = NN;
          end
        end

        %per-face, per-slot pair of edge nodes ORDERED along the face
        %traversal (slots run p1->p2, p2->p3, p1->p3; the stored pair runs
        %from the LOWER vertex id, so swap where the slot runs against it)
        XY = [ 1 2 ; 2 3 ; 1 3 ];
        A = zeros( nT , 3 );   B = zeros( nT , 3 );
        for k = 1:3
          wk = LET(:,k);
          ak = zeros( nT ,1);   bk = zeros( nT ,1);
          ak(wk) = PA( ET(wk,k) );
          bk(wk) = PB( ET(wk,k) );
          sw = wk & ( double( T(:,XY(k,1)) ) > double( T(:,XY(k,2)) ) );
          tmp = ak(sw);  ak(sw) = bk(sw);  bk(sw) = tmp;
          A(:,k) = ak;   B(:,k) = bk;
        end

        %fully split -> the 9 children of the triadic lattice (all the child
        %orientations below were checked to preserve the parent's)
        w = w9;
        T = [ T ; ...
              [ double(T(w,1)) , A(w,1) , A(w,3) ] ;...
              [ A(w,1) , B(w,1) , PC ] ;...
              [ A(w,1) , PC     , A(w,3) ] ;...
              [ B(w,1) , double(T(w,2)) , A(w,2) ] ;...
              [ B(w,1) , A(w,2) , PC ] ;...
              [ PC     , A(w,2) , B(w,2) ] ;...
              [ PC     , B(w,2) , B(w,3) ] ;...
              [ A(w,3) , PC     , B(w,3) ] ;...
              [ B(w,3) , B(w,2) , double(T(w,3)) ] ];
        F = [ F ; repmat( w , 9 , 1 ) ];

        %ONE split edge -> fan of 3 from the opposite vertex
        w = find( LET(:,1) & ~LET(:,2) & ~LET(:,3) );
        T = [ T ; [T(w,3),T(w,1),A(w,1)] ; [T(w,3),A(w,1),B(w,1)] ; [T(w,3),B(w,1),T(w,2)] ];
        F = [ F ; w ; w ; w ];
        w = find( ~LET(:,1) & LET(:,2) & ~LET(:,3) );
        T = [ T ; [T(w,1),T(w,2),A(w,2)] ; [T(w,1),A(w,2),B(w,2)] ; [T(w,1),B(w,2),T(w,3)] ];
        F = [ F ; w ; w ; w ];
        w = find( ~LET(:,1) & ~LET(:,2) & LET(:,3) );
        T = [ T ; [T(w,2),T(w,3),B(w,3)] ; [T(w,2),B(w,3),A(w,3)] ; [T(w,2),A(w,3),T(w,1)] ];
        F = [ F ; w ; w ; w ];

        %TWO split edges -> 5 children: the corner triangle at the shared
        %vertex, a strip of two, and the far quad split towards the far vertex
        w = find( LET(:,1) & LET(:,2) & ~LET(:,3) );        %shared vertex p2
        T = [ T ; [B(w,1),double(T(w,2)),A(w,2)] ; [A(w,1),B(w,1),A(w,2)] ; [A(w,1),A(w,2),B(w,2)] ;...
                  [A(w,1),B(w,2),double(T(w,3))] ; [double(T(w,1)),A(w,1),double(T(w,3))] ];
        F = [ F ; w ; w ; w ; w ; w ];
        w = find( ~LET(:,1) & LET(:,2) & LET(:,3) );        %shared vertex p3
        T = [ T ; [B(w,2),double(T(w,3)),B(w,3)] ; [A(w,2),B(w,2),B(w,3)] ; [A(w,2),B(w,3),A(w,3)] ;...
                  [A(w,2),A(w,3),double(T(w,1))] ; [double(T(w,2)),A(w,2),double(T(w,1))] ];
        F = [ F ; w ; w ; w ; w ; w ];
        w = find( LET(:,1) & ~LET(:,2) & LET(:,3) );        %shared vertex p1
        T = [ T ; [A(w,3),double(T(w,1)),A(w,1)] ; [B(w,3),A(w,3),A(w,1)] ; [B(w,3),A(w,1),B(w,1)] ;...
                  [B(w,3),B(w,1),double(T(w,2))] ; [double(T(w,3)),B(w,3),double(T(w,2))] ];
        F = [ F ; w ; w ; w ; w ; w ];

      case {'loop_matrix'}
        % Loop approximating subdivision, assembling the sparse operator G with
        % Xnew = G*X and returning it in M.LoopMatrix.  Ref [2].
        if ~isequal( W(:).' , 1:nT ), error('LOOP subdivision is only valid when all faces are considered.'); end
        
        E = [ T(:,[1 2 3]) ; T(:,[2 3 1]) ; T(:,[1 3 2]) ];
        E(:,5) = repmat( ( 1:nT ).' , [3,1] );
        E(:,1:2) = sort( E(:,1:2) ,2);      allE = E(:,1:2);
        E = sortROWS( E , [1 2] );
        w = find( all( ~diff( E(:,1:2) , 1 , 1 ) ,2) );
        if any( diff( w ) == 1 )   %3+ equal keys in a row = a non-manifold edge
          error( 'MeshSubdivide:nonManifold' , ...
                 'this scheme''s stencils need an edge-manifold mesh: some edge is shared by 3 or more triangles (the linear / PN schemes do accept that).' );
        end
        E( w  ,4) = E( w+1 ,3);
        E( w+1,:) = [];
        E = sortROWS( E , 5 );
        nE = size( E ,1);
        
        B = E(:,4) == 0;
        B = uniqueROWS( [ E( B  ,1) ; E( B ,2) ] , 1 );
        
        %triangles containing edges
        [~,ET] = ismember2ROWS( allE , E(:,1:2) );
        ET = reshape( ET , [nT,3] );
        allE = [ allE ; allE(:,[2 1]) ];
        allE = uniqueROWS( allE , [1 2] );
        aEB = ismembc( allE , B );

        %indexes of the new points
        P = nP + ( 1:nE ); P = P(:);
        T = [ T ; ...
              [     T(:, 1 )   , P( ET(:, 1 ) ) , P( ET(:, 3 ) ) ] ;...
              [ P( ET(:, 1 ) ) ,     T(:, 2 )   , P( ET(:, 2 ) ) ] ;...
              [ P( ET(:, 3 ) ) , P( ET(:, 2 ) ) ,     T(:, 3 )   ] ;...
              [ P( ET(:, 1 ) ) , P( ET(:, 2 ) ) , P( ET(:, 3 ) ) ] ;...
            ];
        F = [ F ; W ; W ; W ; W ];
        
        %interpolate fields at middle points on edges
        for f = fieldnames( M ).'
          if ~strncmp( f{1} , 'xyz' , 3 ), continue; end
          if strcmp( f{1} , 'xyz' ), continue; end
          M.(f{1}) = [ M.(f{1}) ; ( M.(f{1})( E(:,1) ,:) + M.(f{1})( E(:,2) ,:) )/2 ];
        end

        G = NaN( 100 , 3 ); GN = 0;
        
        %for boundaries: middle points on edges
        w  = all( ismembc( E(:,1:2) , B ) ,2);
        ww = find( w );
        add2G( nP+ww , E(w,1:2) , 1/2 );

        %for internals: rule 3/8 , 1/8
        w  = ~w;
        ww = find( w );
        c1 = 3/8; c2 = 1/8;
        add2G( nP+ww , E(w,1:2) , c1 );
        add2G( nP+ww , E(w,3:4) , c2 );


        %correction of original ("even") on the boundary
        w = all( aEB ,2);
        add2G( allE(w,1) , allE(w,2) , 1/8 );
        
        %correction of original ("even") internal nodes 
        w = ~aEB(:,1);
        K = accumarray( allE(:,1) , 1 );
        beta = ( 5/8 - ( 3/8 + cos( 2*pi ./ K )/4 ).^2 ) ./ K;
        add2G( allE(w,1) , allE(w,2) , beta( allE(w,1) ) );

        
        %perform the interpolation
        ww = ( 1:nP ).';
        d = accumarray( G(1:GN,1) , G(1:GN,3) , [ nP + nE , 1 ] ); d = 1-d( ww );
        add2G( ww , ww , d );

        G = G( 1:GN ,:);
        
        G = sparse( G(:,1) , G(:,2) , G(:,3) , nP + nE , nP );
        M.xyz = G * double( M.xyz );
        M.LoopMatrix = G;
        
      case {'loop'}
        % Loop approximating C2 subdivision (in-place): edge nodes by the 3/8-1/8
        % rule, even nodes relaxed with the beta(K) mask, boundary rules applied. Ref [2].
        if ~isequal( W(:).' , 1:nT ),
%           M1 = MeshSubdivide( M , 'linear4' , W , 'kp' );
%           M2 = MeshSubdivide( M , 'loop_matrix' );
          error('LOOP subdivision is only valid when all faces are considered.'); end
        
        E = [ T(:,[1 2 3]) ; T(:,[2 3 1]) ; T(:,[1 3 2]) ];
        E(:,5) = repmat( ( 1:nT ).' , [3,1] );
        E(:,1:2) = sort( E(:,1:2) ,2);      allE = E(:,1:2);
        E = sortROWS( E , [1 2] );
        w = find( all( ~diff( E(:,1:2) , 1 , 1 ) ,2) );
        if any( diff( w ) == 1 )   %3+ equal keys in a row = a non-manifold edge
          error( 'MeshSubdivide:nonManifold' , ...
                 'this scheme''s stencils need an edge-manifold mesh: some edge is shared by 3 or more triangles (the linear / PN schemes do accept that).' );
        end
        E( w  ,4) = E( w+1 ,3);
        E( w+1,:) = [];
        E = sortROWS( E , 5 );
        nE = size( E ,1);
        
        B = E(:,4) == 0;
        B = uniqueROWS( [ E( B  ,1) ; E( B ,2) ] , 1 );
        
        %triangles containing edges
        [~,ET] = ismember2ROWS( allE , E(:,1:2) );
        ET = reshape( ET , [nT,3] );
        allE = [ allE ; allE(:,[2 1]) ];
        allE = uniqueROWS( allE , [1 2] );
        aEB = ismembc( allE , B );

        %indexes of the new points
        P = nP + ( 1:nE ); P = P(:);
        T = [ T ; ...
              [     T(:, 1 )   , P( ET(:, 1 ) ) , P( ET(:, 3 ) ) ] ;...
              [ P( ET(:, 1 ) ) ,     T(:, 2 )   , P( ET(:, 2 ) ) ] ;...
              [ P( ET(:, 3 ) ) , P( ET(:, 2 ) ) ,     T(:, 3 )   ] ;...
              [ P( ET(:, 1 ) ) , P( ET(:, 2 ) ) , P( ET(:, 3 ) ) ] ;...
            ];
        F = [ F ; W ; W ; W ; W ];
        
        %interpolate fields at middle points on edges
        for f = fieldnames( M ).'
          if ~strncmp( f{1} , 'xyz' , 3 ), continue; end
          if strcmp( f{1} , 'xyz' ), continue; end
          M.(f{1}) = [ M.(f{1}) ; ( M.(f{1})( E(:,1) ,:,:,:,:,:,:,:) + M.(f{1})( E(:,2) ,:,:,:,:,:,:,:) )/2 ];
        end

        M.xyz( nP + nE ,1) = 0;
        
        %for boundaries: middle points on edges
        w  = all( ismembc( E(:,1:2) , B ) ,2);
        ww = find( w );
        M.xyz( nP+ww ,: ) = ( M.xyz( E(w,1) ,:,:,:,:,:,:,:) + M.xyz( E(w,2) ,:,:,:,:,:,:,:) )/2;

        %for internals: rule 3/8 , 1/8
        w  = ~w; Ew = E(w,:);
        ww = find( w );
        c1 = 3/8; c2 = 1/8;
        M.xyz( nP+ww ,: ) = c1 * ( M.xyz( Ew(:,1) ,:,:,:,:,:,:,:) + M.xyz( Ew(:,2) ,:,:,:,:,:,:,:) ) +...
                            c2 * ( M.xyz( Ew(:,3) ,:,:,:,:,:,:,:) + M.xyz( Ew(:,4) ,:,:,:,:,:,:,:) );


        G = NaN( 10 , 3 ); GN = 0;

        %correction of original ("even") on the boundary
        w = all( aEB ,2); %aEw = aE(w,:);
        add2G( allE(w,1) , allE(w,2) , 1/8 );
        
        %correction of original ("even") internal nodes 
        w = ~aEB(:,1);
        K = accumarray( allE(:,1) , 1 );
        beta = ( 5/8 - ( 3/8 + cos( 2*pi ./ K )/4 ).^2 ) ./ K;
        add2G( allE(w,1) , allE(w,2) , beta( allE(w,1) ) );

        
        %perform the interpolation
        ww = ( 1:nP ).';
        d = accumarray( G( 1:GN ,1) , G( 1:GN ,3) , [ nP , 1 ] ); d = 1-d;
        add2G( ww , ww , d );
        
        G = G( 1:GN ,:);
        
        G = sparse( G(:,1) , G(:,2) , G(:,3) , nP , nP );
        M.xyz(ww,:) = G * double( M.xyz(ww,:) );

      case {'butterfly'}
        % Modified-Butterfly INTERPOLATING subdivision: original nodes are kept,
        % new edge nodes use the 8-point / extraordinary stencils. Cases mirror
        % VTK's vtkButterflySubdivisionFilter.  Ref [3].
        if ~isequal( W(:).' , 1:nT ), error('BUTTERFLY subdivision is only valid when all faces are considered.'); end
        
        E = [ T(:,[1 2 3]) ; T(:,[2 3 1]) ; T(:,[1 3 2]) ];
        E(:,5) = repmat( ( 1:nT ).' , [3,1] );
        E(:,1:2) = sort( E(:,1:2) ,2);      allE = E(:,1:2);
        E = sortROWS( E , [1 2] );
        w = find( all( ~diff( E(:,1:2) , 1 , 1 ) ,2) );
        if any( diff( w ) == 1 )   %3+ equal keys in a row = a non-manifold edge
          error( 'MeshSubdivide:nonManifold' , ...
                 'this scheme''s stencils need an edge-manifold mesh: some edge is shared by 3 or more triangles (the linear / PN schemes do accept that).' );
        end
        E( w  ,4) = E( w+1 ,3);
        E( w+1,:) = [];
        E = sortROWS( E , 5 );
        nE = size( E ,1);
        E(:,5) = ( 1:nE ).';
        
        K = accumarray( T(:) , 1 , [nP,1] );
        
        %triangles containing edges
        [~,ET] = ismember2ROWS( allE , E(:,1:2) );
        ET = reshape( ET , [nT,3] );
        aEB = E( ~E(:,4) ,1:2);

        %indexes of the new points
        P = nP + ( 1:nE ); P = P(:);
        T = [ T ; ...
              [     T(:, 1 )   , P( ET(:, 1 ) ) , P( ET(:, 3 ) ) ] ;...
              [ P( ET(:, 1 ) ) ,     T(:, 2 )   , P( ET(:, 2 ) ) ] ;...
              [ P( ET(:, 3 ) ) , P( ET(:, 2 ) ) ,     T(:, 3 )   ] ;...
              [ P( ET(:, 1 ) ) , P( ET(:, 2 ) ) , P( ET(:, 3 ) ) ] ;...
            ];
        F = [ F ; W ; W ; W ; W ];
        
        %interpolate fields at middle points on edges
        for f = fieldnames( M ).'
          if ~strncmp( f{1} , 'xyz' , 3 ), continue; end
          if strcmp( f{1} , 'xyz' ), continue; end
          M.(f{1}) = [ M.(f{1}) ; ( M.(f{1})( E(:,1) ,:) + M.(f{1})( E(:,2) ,:) )/2 ];
        end

        %butterfly stencil
        Ec = E(:,1) + 1i*E(:,2);
        [~,c] = ismember2ROWS( sort( E(:,[1,3]) , 2 ) , Ec ); cc = [ 6 7 ]; w=~~c; E(w,cc) = E(c(w),3:4); w = E(:,cc(1)) == E(:,1) | E(:,cc(1)) == E(:,2); E(w,cc(1)) = 0;  w = E(:,cc(2)) == E(:,1) | E(:,cc(2)) == E(:,2); E(w,cc(2)) = 0; w=all(~E(:,cc),2); E(w,cc(1))=-E(w,4);
        [~,c] = ismember2ROWS( sort( E(:,[1,4]) , 2 ) , Ec ); cc = [ 8 9 ]; w=~~c; E(w,cc) = E(c(w),3:4); w = E(:,cc(1)) == E(:,1) | E(:,cc(1)) == E(:,2); E(w,cc(1)) = 0;  w = E(:,cc(2)) == E(:,1) | E(:,cc(2)) == E(:,2); E(w,cc(2)) = 0; w=all(~E(:,cc),2); E(w,cc(1))=-E(w,3);
        [~,c] = ismember2ROWS( sort( E(:,[2,3]) , 2 ) , Ec ); cc = [10 11]; w=~~c; E(w,cc) = E(c(w),3:4); w = E(:,cc(1)) == E(:,1) | E(:,cc(1)) == E(:,2); E(w,cc(1)) = 0;  w = E(:,cc(2)) == E(:,1) | E(:,cc(2)) == E(:,2); E(w,cc(2)) = 0; w=all(~E(:,cc),2); E(w,cc(1))=-E(w,4);
        [~,c] = ismember2ROWS( sort( E(:,[2,4]) , 2 ) , Ec ); cc = [12 13]; w=~~c; E(w,cc) = E(c(w),3:4); w = E(:,cc(1)) == E(:,1) | E(:,cc(1)) == E(:,2); E(w,cc(1)) = 0;  w = E(:,cc(2)) == E(:,1) | E(:,cc(2)) == E(:,2); E(w,cc(2)) = 0; w=all(~E(:,cc),2); E(w,cc(1))=-E(w,3);
        
        E( ~E(:) ) = -Inf;
        E(:, 6:13) = sort( E(:,6:13) , 2 ,'descend');
        E( ~isfinite( E(:) ) ) = 0;
        E(:,10:13) = [];


        B = unique( E( ~E(:,4),1:2) );
        
        G = NaN( 100 ,3); GN = 0;
        %TYPES = zeros( nE , 1 );
        
        %cases from:
        %vtkButterflySubdivisionFilter.cxx

        
        %case boundary, edges belonging to only one triangle
        w = ~E(:,4);
        Ew = E(w,:); eid = Ew(:,5);
        for n = 1:numel(eid)
          p1 = Ew(n,1); p2 = Ew(n,2);
          
          R1 = aEB( any( aEB == p1 ,2) ,:);
          R1( any( R1 == p2 ,2) ,:) = [];
          R1 = R1(1,:); R1( R1 == p1 ) = [];
          
          R2 = aEB( any( aEB == p2 ,2) ,:);
          R2( any( R2 == p1 ,2) ,:) = [];
          R2 = R2(1,:); R2( R2 == p2 ) = [];
          
          add2G( eid(n) , [ R1 ; R2 ] , -1/16 );
          %add2G( eid(n) , R2 , -1/16 );
        end
        add2G( eid , Ew(:,1:2) , 9/16 );
        E(w,:) = []; %TYPES(eid) = 1;
        
        %boundary-boundary case (ears case)
        w = all( ismember( E(:,1:2) , B ) ,2);
        Ew = E(w,:); eid = Ew(:,5);
        add2G( eid , Ew(:,1:2) ,  1/2  );
        E(w,:) = []; %TYPES(eid) = 2;
        
        %regular-regular interior
        w = K(E(:,1)) == 6 & K(E(:,2)) == 6;
        Ew = E(w,:); eid = Ew(:,5);
        add2G( eid , Ew(:,1:2) ,  1/2  );
        add2G( eid , Ew(:,3:4) ,  1/8  );
        add2G( eid , Ew(:,6:9) , -1/16 );
        E(w,:) = []; %TYPES(eid) = 2;


        %extraordinary-regular
        w = K(E(:,1)) ~= 6 & K(E(:,2)) == 6;
        E = [ E(w,:) ; E(~w,:) ];
        w = K(E(:,1)) == 6 & K(E(:,2)) ~= 6;
        E = [ E(w,[2 1 3:end]) ; E(~w,:) ];
        
        
        w = K(E(:,1)) ~= 6 & K(E(:,2)) == 6;
        Ew = E(w,:); eid = Ew(:,5);
        
        TT = M.tri; TT = TT( any( ismember( TT , E(:,1:2) ) ,2) ,:);
        RS = meshRings( TT );
        
        for n = 1:numel(eid)
          R = CIRCULARshift( RS{ Ew(n,1) }.' , Ew(n,2) );
          k = numel(R);
          switch k
            case 0,     s = [1/2,1/2,1/8,1/8,-1/16,-1/16,-1/16,-1/16]; R =  Ew(n,[1:4 6:9]);
            case 3,     s = [ 5/12 , -1/12 , -1/12 ];
                        s = [ s , 1-sum(s) ]; R = [ R , Ew(n,1) ];
            case 4,     s = [ 3/8 , 0 , -1/8 , 0 ];
                        s = [ s , 1-sum(s) ]; R = [ R , Ew(n,1) ];
            otherwise,  s = ( 1/4 + cos( 2 * pi * (0:k-1) / k ) + 1/2 * cos( 4 * pi * (0:k-1) / k ) ) / k;
                        s = [ s , 1-sum(s) ]; R = [ R , Ew(n,1) ];
          end
          add2G( eid(n) , R , s );
        end
        E(w,:) = []; %TYPES(eid) = 3;
        
        
        %extraordinary-extraordinary
        w = K(E(:,1)) ~= 6  &  K(E(:,2)) ~= 6;
        Ew = E(w,:); eid = Ew(:,5);
        for n = 1:numel(eid)
          R = CIRCULARshift( RS{ Ew(n,1) }.' , Ew(n,2) );
          k = numel(R);
          switch k
            case 0,     s = [1/2,1/2,1/8,1/8,-1/16,-1/16,-1/16,-1/16]; R = Ew(n,[1:4 6:9]);
            case 3,     s = [ 5/12 , -1/12 , -1/12 ];
                        s = [ s , 1-sum(s) ]; R = [ R , Ew(n,1) ];
            case 4,     s = [ 3/8 , 0 , -1/8 , 0 ];
                        s = [ s , 1-sum(s) ]; R = [ R , Ew(n,1) ];
            otherwise,  s = ( 1/4 + cos( 2 * pi * (0:k-1) / k ) + 1/2 * cos( 4 * pi * (0:k-1) / k ) ) / k;
                        s = [ s , 1-sum(s) ]; R = [ R , Ew(n,1) ];
          end
          add2G( eid(n) , R , s/2 );
%           if ~isempty(find(isnan(G(1:GN,2))))
%             1;
%           end
          
          R = CIRCULARshift( RS{ Ew(n,2) }.' , Ew(n,1) );
          k = numel(R);
          switch k
            case 0,     s = [1/2,1/2,1/8,1/8,-1/16,-1/16,-1/16,-1/16]; R = Ew(n,[1:4 6:9]);
            case 3,     s = [ 5/12 , -1/12 , -1/12 ];
                        s = [ s , 1-sum(s) ]; R = [ R , Ew(n,2) ];
            case 4,     s = [ 3/8 , 0 , -1/8 , 0 ];
                        s = [ s , 1-sum(s) ]; R = [ R , Ew(n,2) ];
            otherwise,  s = ( 1/4 + cos( 2 * pi * (0:k-1) / k ) + 1/2 * cos( 4 * pi * (0:k-1) / k ) ) / k;
                        s = [ s , 1-sum(s) ]; R = [ R , Ew(n,2) ];
          end
          add2G( eid(n) , R , s/2 );
%           if ~isempty(find(isnan(G(1:GN,2))))
%             1;
%           end
%           GN
        end
        E(w,:) = []; %TYPES(eid) = 4;

        if size( E ,1), warning('there are still Edges not processed'); end
        
        G = G( 1:GN ,:);
        
        G = sparse( G(:,1) , abs(G(:,2)) , G(:,3) , nE , nP );
        M.xyz = [ M.xyz ; G * double( M.xyz ) ];

        %TYPES = [ NaN( nP , 1 ) ; TYPES ];
        %setappdata(0,'butterfly',TYPES);
        
      case {'sqrt3'}
        % Kobbelt sqrt(3)-subdivision [11]: insert the CENTROID of every face,
        % split 1 -> 3, then FLIP every original interior edge -- the pair
        % {(p,q,ma),(q,p,mb)} becomes {(p,mb,ma),(q,ma,mb)}, orientation
        % preserving -- so the original edges VANISH and the new edges join
        % adjacent centroids; finally RELAX the original vertices with the
        % valence mask  p <- (1-a(n))*p + a(n)*mean(old neighbours),
        % a(n) = (4-2cos(2pi/n))/9. APPROXIMATING, C2 away from extraordinary
        % vertices, and the slowest-growing refinement: 1 -> 3 per pass (two
        % passes give a REALIGNED 3-adic 1 -> 9); centroids come out with
        % valence 6, original vertices KEEP their valence. Boundaries: boundary
        % edges are not flipped and boundary vertices are NOT moved (the exact
        % boundary rule of the paper needs an every-second-pass treatment, not
        % implemented); a single triangle thus reduces exactly to 'linear3'.
        % Needs an edge-manifold, consistently wound mesh (checked, unlike
        % 'loop'/'butterfly' which silently assume it). Whole-mesh only.
        if ~isequal( W(:).' , 1:nT ), error('SQRT3 subdivision is only valid when all faces are considered.'); end

        %centroids of every face (all xyz* fields), from the ORIGINAL positions
        for f = fieldnames( M ).', if ~strncmp( f{1} , 'xyz' , 3 ), continue; end
          M.(f{1}) = [ M.(f{1}) ; ( M.(f{1})( T(:,1) ,:,:,:,:) + M.(f{1})( T(:,2) ,:,:,:,:) + M.(f{1})( T(:,3) ,:,:,:,:) )/3 ];
        end
        P = nP + ( 1:nT ).';                   %centroid node of each face

        %directed edges p->q with their face; sorting by the undirected key
        %puts the two copies of each interior edge together
        DE = double( [ T(:,[1 2]) ; T(:,[2 3]) ; T(:,[3 1]) ] );
        FA = repmat( ( 1:nT ).' , 3 , 1 );
        K  = sort( DE , 2 );
        [ K , ord ] = sortrows( K );
        DE = DE( ord ,:);   FA = FA( ord );
        w  = [ ~any( diff( K ,1,1) ,2) ; false ];        %row i pairs with row i+1
        pr = [ false ; w(1:end-1) ];
        if any( w & pr )
          error( 'MeshSubdivide:sqrt3NonManifold' , ...
                 'sqrt3 needs an edge-manifold mesh: some edge is shared by 3 or more triangles.' );
        end
        i1 = find( w );          i2 = i1 + 1;            %the interior pairs
        ib = find( ~w & ~pr );                           %the boundary edges

        p = DE(i1,1);   q = DE(i1,2);
        a = FA(i1);     b = FA(i2);
        if any( DE(i2,1) ~= q | DE(i2,2) ~= p )
          error( 'MeshSubdivide:sqrt3Winding' , ...
                 'sqrt3 needs consistent winding: some interior edge runs in the SAME direction in its two triangles.' );
        end

        %flip every interior edge (orientation preserving), keep the fan
        %triangle of every boundary edge
        Tn = [ [ p , P(b) , P(a) ] ;...
               [ q , P(a) , P(b) ] ;...
               [ DE(ib,:) , P( FA(ib) ) ] ];
        Fn = [ a ; b ; FA(ib) ];

        %relax the interior ORIGINAL vertices with the valence mask, from the
        %ORIGINAL neighbour positions (boundary vertices stay put)
        UE = K( [ i1 ; ib ] ,:);                         %each undirected edge once
        A2 = sparse( [ UE(:,1) ; UE(:,2) ] , [ UE(:,2) ; UE(:,1) ] , 1 , nP , nP );
        n  = full( sum( A2 ,2) );
        isB = false( nP ,1);   isB( K(ib,:) ) = true;
        wv = find( ~isB & n > 0 );
        al = ( 4 - 2*cos( 2*pi ./ n(wv) ) ) / 9;
        XYZ0 = double( M.xyz( 1:nP ,:) );
        NBavg = bsxfun( @rdivide , A2( wv ,:) * XYZ0 , n( wv ) );
        M.xyz( wv ,:) = bsxfun( @times , 1 - al , XYZ0( wv ,:) ) + bsxfun( @times , al , NBavg );

        T = [ T ; Tn ];
        F = [ F ; Fn ];

      case {'cubic hermite'}
        if ~isequal( W(:).' , 1:nT ), error('CUBIC HERMITE subdivision is only valid when all faces are considered.'); end
        
        error('not implemented yet');
        
      otherwise, error('unknown SubType for triangles.');
    end
    
  elseif M.celltype == 10          %% tetrahedral mesh
    switch lower(SubType)
      case {'default'},
        % conforming RED refinement of tetrahedra: 1 -> 8 when all 6 edges split,
        % with every partial-edge pattern handled and diagonals chosen by vertex
        % index so shared faces of adjacent tets agree (no T-junctions).  Ref [5].
        if isempty( W ), return; end

        allE = sort( [ T(:,[1 2]) ; T(:,[1 3]) ; T(:,[1 4]) ; T(:,[2 3]) ; T(:,[2 4]) ; T(:,[3 4]) ] ,2);

        while 1
          E = allE( [ W ; W + nT ; W + 2*nT ; W + 3*nT ; W + 4*nT ; W + 5*nT ] , : );
          E = unique( E , 'rows' );
          ET = reshape( ismember( allE , E , 'rows' ) , nT , 6 );
          Wp = W;

          W = false;
          W = W | sum( ET ,2) > 3;
          W = W | all( bsxfun(@eq, ET , [0 0 1 0 1 1] ) ,2);
          W = W | all( bsxfun(@eq, ET , [0 0 1 1 0 1] ) ,2);
          W = W | all( bsxfun(@eq, ET , [0 0 1 1 1 0] ) ,2);
          W = W | all( bsxfun(@eq, ET , [0 1 0 0 1 1] ) ,2);
          W = W | all( bsxfun(@eq, ET , [0 1 0 1 0 1] ) ,2);
          W = W | all( bsxfun(@eq, ET , [0 1 0 1 1 0] ) ,2);
          W = W | all( bsxfun(@eq, ET , [0 1 1 0 1 0] ) ,2);
          W = W | all( bsxfun(@eq, ET , [0 1 1 1 0 0] ) ,2);
          W = W | all( bsxfun(@eq, ET , [1 0 0 0 1 1] ) ,2);
          W = W | all( bsxfun(@eq, ET , [1 0 0 1 0 1] ) ,2);
          W = W | all( bsxfun(@eq, ET , [1 0 0 1 1 0] ) ,2);
          W = W | all( bsxfun(@eq, ET , [1 0 1 0 0 1] ) ,2);
          W = W | all( bsxfun(@eq, ET , [1 0 1 1 0 0] ) ,2);
          W = W | all( bsxfun(@eq, ET , [1 1 0 0 0 1] ) ,2);
          W = W | all( bsxfun(@eq, ET , [1 1 0 0 1 0] ) ,2);
          W = W | all( bsxfun(@eq, ET , [1 1 1 0 0 0] ) ,2);

          W = find( W );
          if isequal( Wp , W ), break; end
        end

        %tetras containing edges
%         ET = zeros( nT , 6 );
%         [ ~ , ET(:,1) ] = ismember( sort( T(:,[1 2]) ,2) , E , 'rows' );
%         [ ~ , ET(:,2) ] = ismember( sort( T(:,[1 3]) ,2) , E , 'rows' );
%         [ ~ , ET(:,3) ] = ismember( sort( T(:,[1 4]) ,2) , E , 'rows' );
%         [ ~ , ET(:,4) ] = ismember( sort( T(:,[2 3]) ,2) , E , 'rows' );
%         [ ~ , ET(:,5) ] = ismember( sort( T(:,[2 4]) ,2) , E , 'rows' );
%         [ ~ , ET(:,6) ] = ismember( sort( T(:,[3 4]) ,2) , E , 'rows' );
        
        [~,ET] = ismember2ROWS( allE , E(:,1:2) );
        ET = reshape( ET , [nT,6] );

        LET = ~~ET;  %logical ET
        %original faces to be removed
        W = find( any( LET ,2) ); %size(W)


        %indexes of the new points
        P = ( (nP+1):(nP+size(E,1)) ).';

        %middle points on edges
        for f = fieldnames( M ).'
          if ~strncmp( f{1} , 'xyz' , 3 ), continue; end
          M.(f{1}) = [ M.(f{1}) ; ( M.(f{1})( E(:,1) ,:) + M.(f{1})( E(:,2) ,:) )/2 ];
        end
        if KP
          fn = fieldnames(M); fn = sort( fn( strncmp( fn , 'xyzParentEdge' ,13) ) );
          for f = fn(end:-1:1).', M = renameStructField( M , f{1} , [ f{1} , '_' ] ); end
          PE = (1:nP).'; PE(:,2) = 0;
          PE = [ PE ; double(E) ];
          PE( 1:nP ,3) = 0;
          PE( nP+1:end ,3) = 0.5;
          M.xyzParentEdge = PE;
        end
        
        %C = zeros( nT , 1);

        %first, tetras to be divided into 8.
        w = find( all(  bsxfun(@eq, LET , ~~[1 1 1 1 1 1] ) ,2) );
        if ~isempty(w)
        p1  = T(w,1); p2  = T(w,2); p3  = T(w,3); p4  = T(w,4); try, p12 = P( ET(w, 1 ) ); end; try, p13 = P( ET(w, 2 ) ); end; try, p14 = P( ET(w, 3 ) ); end; try, p23 = P( ET(w, 4 ) ); end; try, p24 = P( ET(w, 5 ) ); end; try, p34 = P( ET(w, 6 ) ); end; ET(w,:) = 0;
        T = [ T ;    p1   ,  p12  ,  p13  ,  p14  ]; F = [ F ; w ]; %C = [ C ; w*0 + 8.1 ];
        T = [ T ;    p12  ,  p2   ,  p23  ,  p24  ]; F = [ F ; w ]; %C = [ C ; w*0 + 8.2 ];
        T = [ T ;    p13  ,  p23  ,  p3   ,  p34  ]; F = [ F ; w ]; %C = [ C ; w*0 + 8.3 ];
        T = [ T ;    p14  ,  p24  ,  p34  ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 8.4 ];
        T = [ T ;    p12  ,  p13  ,  p14  ,  p24  ]; F = [ F ; w ]; %C = [ C ; w*0 + 8.5 ];
        T = [ T ;    p12  ,  p13  ,  p24  ,  p23  ]; F = [ F ; w ]; %C = [ C ; w*0 + 8.6 ];
        T = [ T ;    p13  ,  p14  ,  p24  ,  p34  ]; F = [ F ; w ]; %C = [ C ; w*0 + 8.7 ];
        T = [ T ;    p13  ,  p23  ,  p34  ,  p24  ]; F = [ F ; w ]; %C = [ C ; w*0 + 8.8 ];
        end

        %tetras to be divided at edge 1-2
        w = find( all(  bsxfun(@eq, LET , ~~[1 0 0 0 0 0] ) ,2) );
        if ~isempty(w)
        p1  = T(w,1); p2  = T(w,2); p3  = T(w,3); p4  = T(w,4); try, p12 = P( ET(w, 1 ) ); end; try, p13 = P( ET(w, 2 ) ); end; try, p14 = P( ET(w, 3 ) ); end; try, p23 = P( ET(w, 4 ) ); end; try, p24 = P( ET(w, 5 ) ); end; try, p34 = P( ET(w, 6 ) ); end; ET(w,:) = 0;
        T = [ T ;    p1   ,  p12  ,  p3   ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 12.1 ];
        T = [ T ;    p12  ,  p2   ,  p3   ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 12.2 ];
        end

        %tetras to be divided at edge 1-3
        w = find( all(  bsxfun(@eq, LET , ~~[0 1 0 0 0 0] ) ,2) );
        if ~isempty(w)
        p1  = T(w,1); p2  = T(w,2); p3  = T(w,3); p4  = T(w,4); try, p12 = P( ET(w, 1 ) ); end; try, p13 = P( ET(w, 2 ) ); end; try, p14 = P( ET(w, 3 ) ); end; try, p23 = P( ET(w, 4 ) ); end; try, p24 = P( ET(w, 5 ) ); end; try, p34 = P( ET(w, 6 ) ); end; ET(w,:) = 0;
        T = [ T ;    p1   ,  p2   ,  p13  ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 13.1 ];
        T = [ T ;    p13  ,  p2   ,  p3   ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 13.2 ];
        end

        %tetras to be divided at edge 1-4
        w = find( all(  bsxfun(@eq, LET , ~~[0 0 1 0 0 0] ) ,2) );
        if ~isempty(w)
        p1  = T(w,1); p2  = T(w,2); p3  = T(w,3); p4  = T(w,4); try, p12 = P( ET(w, 1 ) ); end; try, p13 = P( ET(w, 2 ) ); end; try, p14 = P( ET(w, 3 ) ); end; try, p23 = P( ET(w, 4 ) ); end; try, p24 = P( ET(w, 5 ) ); end; try, p34 = P( ET(w, 6 ) ); end; ET(w,:) = 0;
        T = [ T ;    p1   ,  p2   ,  p3   ,  p14  ]; F = [ F ; w ]; %C = [ C ; w*0 + 14.1 ];
        T = [ T ;    p14  ,  p2   ,  p3   ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 14.2 ];
        end

        %tetras to be divided at edge 2-3
        w = find( all(  bsxfun(@eq, LET , ~~[0 0 0 1 0 0] ) ,2) );
        if ~isempty(w)
        p1  = T(w,1); p2  = T(w,2); p3  = T(w,3); p4  = T(w,4); try, p12 = P( ET(w, 1 ) ); end; try, p13 = P( ET(w, 2 ) ); end; try, p14 = P( ET(w, 3 ) ); end; try, p23 = P( ET(w, 4 ) ); end; try, p24 = P( ET(w, 5 ) ); end; try, p34 = P( ET(w, 6 ) ); end; ET(w,:) = 0;
        T = [ T ;    p1   ,  p2   ,  p23  ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 23.1 ];
        T = [ T ;    p1   ,  p23  ,  p3   ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 23.2 ];
        end

        %tetras to be divided at edge 2-4
        w = find( all(  bsxfun(@eq, LET , ~~[0 0 0 0 1 0] ) ,2) );
        if ~isempty(w)
        p1  = T(w,1); p2  = T(w,2); p3  = T(w,3); p4  = T(w,4); try, p12 = P( ET(w, 1 ) ); end; try, p13 = P( ET(w, 2 ) ); end; try, p14 = P( ET(w, 3 ) ); end; try, p23 = P( ET(w, 4 ) ); end; try, p24 = P( ET(w, 5 ) ); end; try, p34 = P( ET(w, 6 ) ); end; ET(w,:) = 0;
        T = [ T ;    p1   ,  p2   ,  p3   ,  p24  ]; F = [ F ; w ]; %C = [ C ; w*0 + 24.1 ];
        T = [ T ;    p1   ,  p24  ,  p3   ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 24.2 ];
        end

        %tetras to be divided at edge 3-4
        w = find( all(  bsxfun(@eq, LET , ~~[0 0 0 0 0 1] ) ,2) );
        if ~isempty(w)
        p1  = T(w,1); p2  = T(w,2); p3  = T(w,3); p4  = T(w,4); try, p12 = P( ET(w, 1 ) ); end; try, p13 = P( ET(w, 2 ) ); end; try, p14 = P( ET(w, 3 ) ); end; try, p23 = P( ET(w, 4 ) ); end; try, p24 = P( ET(w, 5 ) ); end; try, p34 = P( ET(w, 6 ) ); end; ET(w,:) = 0;
        T = [ T ;    p1   ,  p2   ,  p3   ,  p34  ]; F = [ F ; w ]; %C = [ C ; w*0 + 34.1 ];
        T = [ T ;    p1   ,  p2   ,  p34  ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 34.2 ];
        end

        %tetras to be divided at edge 1-3 & 2-4
        w = find( all(  bsxfun(@eq, LET , ~~[0 1 0 0 1 0] ) ,2) );
        if ~isempty(w)
        p1  = T(w,1); p2  = T(w,2); p3  = T(w,3); p4  = T(w,4); try, p12 = P( ET(w, 1 ) ); end; try, p13 = P( ET(w, 2 ) ); end; try, p14 = P( ET(w, 3 ) ); end; try, p23 = P( ET(w, 4 ) ); end; try, p24 = P( ET(w, 5 ) ); end; try, p34 = P( ET(w, 6 ) ); end; ET(w,:) = 0;
        T = [ T ;    p1   ,  p13  ,  p4   ,  p24  ]; F = [ F ; w ]; %C = [ C ; w*0 + 1324.1 ];
        T = [ T ;    p13  ,  p3   ,  p4   ,  p24  ]; F = [ F ; w ]; %C = [ C ; w*0 + 1324.2 ];
        T = [ T ;    p1   ,  p13  ,  p24  ,  p2   ]; F = [ F ; w ]; %C = [ C ; w*0 + 1324.3 ];
        T = [ T ;    p13  ,  p3   ,  p24  ,  p2   ]; F = [ F ; w ]; %C = [ C ; w*0 + 1324.4 ];
        end

        %tetras to be divided at edge 1-2 & 3-4
        w = find( all(  bsxfun(@eq, LET , ~~[1 0 0 0 0 1] ) ,2) );
        if ~isempty(w)
        p1  = T(w,1); p2  = T(w,2); p3  = T(w,3); p4  = T(w,4); try, p12 = P( ET(w, 1 ) ); end; try, p13 = P( ET(w, 2 ) ); end; try, p14 = P( ET(w, 3 ) ); end; try, p23 = P( ET(w, 4 ) ); end; try, p24 = P( ET(w, 5 ) ); end; try, p34 = P( ET(w, 6 ) ); end; ET(w,:) = 0;
        T = [ T ;    p1   ,  p12  ,  p3   ,  p34  ]; F = [ F ; w ]; %C = [ C ; w*0 + 1234.1 ];
        T = [ T ;    p1   ,  p12  ,  p34  ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 1234.2 ];
        T = [ T ;    p12  ,  p2   ,  p3   ,  p34  ]; F = [ F ; w ]; %C = [ C ; w*0 + 1234.3 ];
        T = [ T ;    p12  ,  p2   ,  p34  ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 1234.4 ];
        end

        %tetras to be divided at edge 1-4 & 2-3
        w = find( all(  bsxfun(@eq, LET , ~~[0 0 1 1 0 0] ) ,2) );
        if ~isempty(w)
        p1  = T(w,1); p2  = T(w,2); p3  = T(w,3); p4  = T(w,4); try, p12 = P( ET(w, 1 ) ); end; try, p13 = P( ET(w, 2 ) ); end; try, p14 = P( ET(w, 3 ) ); end; try, p23 = P( ET(w, 4 ) ); end; try, p24 = P( ET(w, 5 ) ); end; try, p34 = P( ET(w, 6 ) ); end; ET(w,:) = 0;
        T = [ T ;    p1   ,  p2   ,  p23  ,  p14  ]; F = [ F ; w ]; %C = [ C ; w*0 + 1423.1 ];
        T = [ T ;    p1   ,  p23  ,  p3   ,  p14  ]; F = [ F ; w ]; %C = [ C ; w*0 + 1423.2 ];
        T = [ T ;    p14  ,  p2   ,  p23  ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 1423.3 ];
        T = [ T ;    p14  ,  p23  ,  p3   ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 1423.4 ];
        end


        %tetras to be divided at edge 1-2 & 1-3 & 2-3
        w = find( all(  bsxfun(@eq, LET , ~~[1 1 0 1 0 0] ) ,2) );
        if ~isempty(w)
        p1  = T(w,1); p2  = T(w,2); p3  = T(w,3); p4  = T(w,4); try, p12 = P( ET(w, 1 ) ); end; try, p13 = P( ET(w, 2 ) ); end; try, p14 = P( ET(w, 3 ) ); end; try, p23 = P( ET(w, 4 ) ); end; try, p24 = P( ET(w, 5 ) ); end; try, p34 = P( ET(w, 6 ) ); end; ET(w,:) = 0;
        T = [ T ;    p1   ,  p12  ,  p13  ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 121323.1 ];
        T = [ T ;    p12  ,  p2   ,  p23  ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 121323.2 ];
        T = [ T ;    p13  ,  p23  ,  p3   ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 121323.3 ];
        T = [ T ;    p12  ,  p23  ,  p13  ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 121323.4 ];
        end

        %tetras to be divided at edge 2-3 & 2-4 & 3-4
        w = find( all(  bsxfun(@eq, LET , ~~[0 0 0 1 1 1] ) ,2) );
        if ~isempty(w)
        p1  = T(w,1); p2  = T(w,2); p3  = T(w,3); p4  = T(w,4); try, p12 = P( ET(w, 1 ) ); end; try, p13 = P( ET(w, 2 ) ); end; try, p14 = P( ET(w, 3 ) ); end; try, p23 = P( ET(w, 4 ) ); end; try, p24 = P( ET(w, 5 ) ); end; try, p34 = P( ET(w, 6 ) ); end; ET(w,:) = 0;
        T = [ T ;    p1   ,  p2   ,  p23  ,  p24  ]; F = [ F ; w ]; %C = [ C ; w*0 + 232434.1 ];
        T = [ T ;    p1   ,  p24  ,  p34  ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 232434.2 ];
        T = [ T ;    p1   ,  p23  ,  p3   ,  p34  ]; F = [ F ; w ]; %C = [ C ; w*0 + 232434.3 ];
        T = [ T ;    p1   ,  p24  ,  p23  ,  p34  ]; F = [ F ; w ]; %C = [ C ; w*0 + 232434.4 ];
        end

        %tetras to be divided at edge 1-3 & 1-4 & 3-4
        w = find( all(  bsxfun(@eq, LET , ~~[0 1 1 0 0 1] ) ,2) );
        if ~isempty(w)
        p1  = T(w,1); p2  = T(w,2); p3  = T(w,3); p4  = T(w,4); try, p12 = P( ET(w, 1 ) ); end; try, p13 = P( ET(w, 2 ) ); end; try, p14 = P( ET(w, 3 ) ); end; try, p23 = P( ET(w, 4 ) ); end; try, p24 = P( ET(w, 5 ) ); end; try, p34 = P( ET(w, 6 ) ); end; ET(w,:) = 0;
        T = [ T ;    p1   ,  p2   ,  p13  ,  p14  ]; F = [ F ; w ]; %C = [ C ; w*0 + 131434.1 ];
        T = [ T ;    p13  ,  p2   ,  p3   ,  p34  ]; F = [ F ; w ]; %C = [ C ; w*0 + 131434.2 ];
        T = [ T ;    p14  ,  p2   ,  p34  ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 131434.3 ];
        T = [ T ;    p14  ,  p2   ,  p13  ,  p34  ]; F = [ F ; w ]; %C = [ C ; w*0 + 131434.4 ];
        end

        %tetras to be divided at edge 1-2 & 1-4 & 2-4
        w = find( all(  bsxfun(@eq, LET , ~~[1 0 1 0 1 0] ) ,2) );
        if ~isempty(w)
        p1  = T(w,1); p2  = T(w,2); p3  = T(w,3); p4  = T(w,4); try, p12 = P( ET(w, 1 ) ); end; try, p13 = P( ET(w, 2 ) ); end; try, p14 = P( ET(w, 3 ) ); end; try, p23 = P( ET(w, 4 ) ); end; try, p24 = P( ET(w, 5 ) ); end; try, p34 = P( ET(w, 6 ) ); end; ET(w,:) = 0;
        T = [ T ;    p1   ,  p12  ,  p3   ,  p14  ]; F = [ F ; w ]; %C = [ C ; w*0 + 121424.1 ];
        T = [ T ;    p12  ,  p2   ,  p3   ,  p24  ]; F = [ F ; w ]; %C = [ C ; w*0 + 121424.2 ];
        T = [ T ;    p14  ,  p24  ,  p3   ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 121424.3 ];
        T = [ T ;    p12  ,  p24  ,  p3   ,  p14  ]; F = [ F ; w ]; %C = [ C ; w*0 + 121424.4 ];
        end



        %tetras to be divided at edge 2-4 & 3-4
        w = find( all(  bsxfun(@eq, LET , ~~[0 0 0 0 1 1] ) ,2) ); w( T(w,2) > T(w,3) ) = [];
        if ~isempty(w)
        p1  = T(w,1); p2  = T(w,2); p3  = T(w,3); p4  = T(w,4); try, p12 = P( ET(w, 1 ) ); end; try, p13 = P( ET(w, 2 ) ); end; try, p14 = P( ET(w, 3 ) ); end; try, p23 = P( ET(w, 4 ) ); end; try, p24 = P( ET(w, 5 ) ); end; try, p34 = P( ET(w, 6 ) ); end; ET(w,:) = 0;
        T = [ T ;    p1   ,  p24  ,  p34  ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 2434.1 ];
        T = [ T ;    p1   ,  p2   ,  p3   ,  p34  ]; F = [ F ; w ]; %C = [ C ; w*0 + 2434.2 ];
        T = [ T ;    p1   ,  p34  ,  p24  ,  p2   ]; F = [ F ; w ]; %C = [ C ; w*0 + 2434.3 ];
        end
        w = find( all(  bsxfun(@eq, LET , ~~[0 0 0 0 1 1] ) ,2) ); w( T(w,2) < T(w,3) ) = [];
        if ~isempty(w)
        p1  = T(w,1); p2  = T(w,2); p3  = T(w,3); p4  = T(w,4); try, p12 = P( ET(w, 1 ) ); end; try, p13 = P( ET(w, 2 ) ); end; try, p14 = P( ET(w, 3 ) ); end; try, p23 = P( ET(w, 4 ) ); end; try, p24 = P( ET(w, 5 ) ); end; try, p34 = P( ET(w, 6 ) ); end; ET(w,:) = 0;
        T = [ T ;    p1   ,  p24  ,  p34  ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 2434.4 ];
        T = [ T ;    p1   ,  p3   ,  p24  ,  p2   ]; F = [ F ; w ]; %C = [ C ; w*0 + 2434.5 ];
        T = [ T ;    p1   ,  p24  ,  p3   ,  p34  ]; F = [ F ; w ]; %C = [ C ; w*0 + 2434.6 ];
        end

        %tetras to be divided at edge 2-3 & 2-4
        w = find( all(  bsxfun(@eq, LET , ~~[0 0 0 1 1 0] ) ,2) ); w( T(w,3) > T(w,4) ) = [];
        if ~isempty(w)
        p1  = T(w,1); p2  = T(w,2); p3  = T(w,3); p4  = T(w,4); try, p12 = P( ET(w, 1 ) ); end; try, p13 = P( ET(w, 2 ) ); end; try, p14 = P( ET(w, 3 ) ); end; try, p23 = P( ET(w, 4 ) ); end; try, p24 = P( ET(w, 5 ) ); end; try, p34 = P( ET(w, 6 ) ); end; ET(w,:) = 0;
        T = [ T ;    p1   ,  p2   ,  p23  ,  p24  ]; F = [ F ; w ]; %C = [ C ; w*0 + 2324.1 ];
        T = [ T ;    p1   ,  p3   ,  p4   ,  p24  ]; F = [ F ; w ]; %C = [ C ; w*0 + 2324.2 ];
        T = [ T ;    p1   ,  p24  ,  p23  ,  p3   ]; F = [ F ; w ]; %C = [ C ; w*0 + 2324.3 ];
        end
        w = find( all(  bsxfun(@eq, LET , ~~[0 0 0 1 1 0] ) ,2) ); w( T(w,3) < T(w,4) ) = [];
        if ~isempty(w)
        p1  = T(w,1); p2  = T(w,2); p3  = T(w,3); p4  = T(w,4); try, p12 = P( ET(w, 1 ) ); end; try, p13 = P( ET(w, 2 ) ); end; try, p14 = P( ET(w, 3 ) ); end; try, p23 = P( ET(w, 4 ) ); end; try, p24 = P( ET(w, 5 ) ); end; try, p34 = P( ET(w, 6 ) ); end; ET(w,:) = 0;
        T = [ T ;    p1   ,  p2   ,  p23  ,  p24  ]; F = [ F ; w ]; %C = [ C ; w*0 + 2324.4 ];
        T = [ T ;    p1   ,  p4   ,  p24  ,  p23  ]; F = [ F ; w ]; %C = [ C ; w*0 + 2324.5 ];
        T = [ T ;    p1   ,  p23  ,  p3   ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 2324.6 ];
        end

        %tetras to be divided at edge 2-3 & 3-4
        w = find( all(  bsxfun(@eq, LET , ~~[0 0 0 1 0 1] ) ,2) ); w( T(w,2) > T(w,4) ) = [];
        if ~isempty(w)
        p1  = T(w,1); p2  = T(w,2); p3  = T(w,3); p4  = T(w,4); try, p12 = P( ET(w, 1 ) ); end; try, p13 = P( ET(w, 2 ) ); end; try, p14 = P( ET(w, 3 ) ); end; try, p23 = P( ET(w, 4 ) ); end; try, p24 = P( ET(w, 5 ) ); end; try, p34 = P( ET(w, 6 ) ); end; ET(w,:) = 0;
        T = [ T ;    p1   ,  p23  ,  p3   ,  p34  ]; F = [ F ; w ]; %C = [ C ; w*0 + 2334.1 ];
        T = [ T ;    p1   ,  p2   ,  p23  ,  p34  ]; F = [ F ; w ]; %C = [ C ; w*0 + 2334.2 ];
        T = [ T ;    p1   ,  p34  ,  p4   ,  p2   ]; F = [ F ; w ]; %C = [ C ; w*0 + 2334.3 ];
        end
        w = find( all(  bsxfun(@eq, LET , ~~[0 0 0 1 0 1] ) ,2) ); w( T(w,2) < T(w,4) ) = [];
        if ~isempty(w)
        p1  = T(w,1); p2  = T(w,2); p3  = T(w,3); p4  = T(w,4); try, p12 = P( ET(w, 1 ) ); end; try, p13 = P( ET(w, 2 ) ); end; try, p14 = P( ET(w, 3 ) ); end; try, p23 = P( ET(w, 4 ) ); end; try, p24 = P( ET(w, 5 ) ); end; try, p34 = P( ET(w, 6 ) ); end; ET(w,:) = 0;
        T = [ T ;    p1   ,  p23  ,  p3   ,  p34  ]; F = [ F ; w ]; %C = [ C ; w*0 + 2334.4 ];
        T = [ T ;    p1   ,  p4   ,  p23  ,  p34  ]; F = [ F ; w ]; %C = [ C ; w*0 + 2334.5 ];
        T = [ T ;    p1   ,  p23  ,  p4   ,  p2   ]; F = [ F ; w ]; %C = [ C ; w*0 + 2334.6 ];
        end

        %tetras to be divided at edge 1-2 & 1-3
        w = find( all(  bsxfun(@eq, LET , ~~[1 1 0 0 0 0] ) ,2) ); w( T(w,2) > T(w,3) ) = [];
        if ~isempty(w)
        p1  = T(w,1); p2  = T(w,2); p3  = T(w,3); p4  = T(w,4); try, p12 = P( ET(w, 1 ) ); end; try, p13 = P( ET(w, 2 ) ); end; try, p14 = P( ET(w, 3 ) ); end; try, p23 = P( ET(w, 4 ) ); end; try, p24 = P( ET(w, 5 ) ); end; try, p34 = P( ET(w, 6 ) ); end; ET(w,:) = 0;
        T = [ T ;    p1   ,  p12  ,  p13  ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 1213.1 ];
        T = [ T ;    p2   ,  p3   ,  p13  ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 1213.2 ];
        T = [ T ;    p13  ,  p12  ,  p2   ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 1213.3 ];
        end
        w = find( all(  bsxfun(@eq, LET , ~~[1 1 0 0 0 0] ) ,2) ); w( T(w,2) < T(w,3) ) = [];
        if ~isempty(w)
        p1  = T(w,1); p2  = T(w,2); p3  = T(w,3); p4  = T(w,4); try, p12 = P( ET(w, 1 ) ); end; try, p13 = P( ET(w, 2 ) ); end; try, p14 = P( ET(w, 3 ) ); end; try, p23 = P( ET(w, 4 ) ); end; try, p24 = P( ET(w, 5 ) ); end; try, p34 = P( ET(w, 6 ) ); end; ET(w,:) = 0;
        T = [ T ;    p1   ,  p12  ,  p13  ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 1213.4 ];
        T = [ T ;    p3   ,  p13  ,  p12  ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 1213.5 ];
        T = [ T ;    p12  ,  p2   ,  p3   ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 1213.6 ];
        end

        %tetras to be divided at edge 1-2 & 1-4
        w = find( all(  bsxfun(@eq, LET , ~~[1 0 1 0 0 0] ) ,2) ); w( T(w,2) > T(w,4) ) = [];
        if ~isempty(w)
        p1  = T(w,1); p2  = T(w,2); p3  = T(w,3); p4  = T(w,4); try, p12 = P( ET(w, 1 ) ); end; try, p13 = P( ET(w, 2 ) ); end; try, p14 = P( ET(w, 3 ) ); end; try, p23 = P( ET(w, 4 ) ); end; try, p24 = P( ET(w, 5 ) ); end; try, p34 = P( ET(w, 6 ) ); end; ET(w,:) = 0;
        T = [ T ;    p1   ,  p12  ,  p3   ,  p14  ]; F = [ F ; w ]; %C = [ C ; w*0 + 1214.1 ];
        T = [ T ;    p2   ,  p4   ,  p3   ,  p14  ]; F = [ F ; w ]; %C = [ C ; w*0 + 1214.2 ];
        T = [ T ;    p14  ,  p12  ,  p3   ,  p2   ]; F = [ F ; w ]; %C = [ C ; w*0 + 1214.3 ];
        end
        w = find( all(  bsxfun(@eq, LET , ~~[1 0 1 0 0 0] ) ,2) ); w( T(w,2) < T(w,4) ) = [];
        if ~isempty(w)
        p1  = T(w,1); p2  = T(w,2); p3  = T(w,3); p4  = T(w,4); try, p12 = P( ET(w, 1 ) ); end; try, p13 = P( ET(w, 2 ) ); end; try, p14 = P( ET(w, 3 ) ); end; try, p23 = P( ET(w, 4 ) ); end; try, p24 = P( ET(w, 5 ) ); end; try, p34 = P( ET(w, 6 ) ); end; ET(w,:) = 0;
        T = [ T ;    p1   ,  p12  ,  p3   ,  p14  ]; F = [ F ; w ]; %C = [ C ; w*0 + 1214.4 ];
        T = [ T ;    p4   ,  p14  ,  p3   ,  p12  ]; F = [ F ; w ]; %C = [ C ; w*0 + 1214.5 ];
        T = [ T ;    p12  ,  p2   ,  p3   ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 1214.6 ];
        end

        %tetras to be divided at edge 1-2 & 2-4
        w = find( all(  bsxfun(@eq, LET , ~~[1 0 0 0 1 0] ) ,2) ); w( T(w,1) > T(w,4) ) = [];
        if ~isempty(w)
        p1  = T(w,1); p2  = T(w,2); p3  = T(w,3); p4  = T(w,4); try, p12 = P( ET(w, 1 ) ); end; try, p13 = P( ET(w, 2 ) ); end; try, p14 = P( ET(w, 3 ) ); end; try, p23 = P( ET(w, 4 ) ); end; try, p24 = P( ET(w, 5 ) ); end; try, p34 = P( ET(w, 6 ) ); end; ET(w,:) = 0;
        T = [ T ;    p12  ,  p2   ,  p3   ,  p24  ]; F = [ F ; w ]; %C = [ C ; w*0 + 1224.1 ];
        T = [ T ;    p1   ,  p12  ,  p3   ,  p24  ]; F = [ F ; w ]; %C = [ C ; w*0 + 1224.2 ];
        T = [ T ;    p24  ,  p4   ,  p3   ,  p1   ]; F = [ F ; w ]; %C = [ C ; w*0 + 1224.3 ];
        end
        w = find( all(  bsxfun(@eq, LET , ~~[1 0 0 0 1 0] ) ,2) ); w( T(w,1) < T(w,4) ) = [];
        if ~isempty(w)
        p1  = T(w,1); p2  = T(w,2); p3  = T(w,3); p4  = T(w,4); try, p12 = P( ET(w, 1 ) ); end; try, p13 = P( ET(w, 2 ) ); end; try, p14 = P( ET(w, 3 ) ); end; try, p23 = P( ET(w, 4 ) ); end; try, p24 = P( ET(w, 5 ) ); end; try, p34 = P( ET(w, 6 ) ); end; ET(w,:) = 0;
        T = [ T ;    p12  ,  p2   ,  p3   ,  p24  ]; F = [ F ; w ]; %C = [ C ; w*0 + 1224.4 ];
        T = [ T ;    p4   ,  p12  ,  p3   ,  p24  ]; F = [ F ; w ]; %C = [ C ; w*0 + 1224.5 ];
        T = [ T ;    p12  ,  p4   ,  p3   ,  p1   ]; F = [ F ; w ]; %C = [ C ; w*0 + 1224.6 ];
        end

        %tetras to be divided at edge 1-2 & 2-3
        w = find( all(  bsxfun(@eq, LET , ~~[1 0 0 1 0 0] ) ,2) ); w( T(w,1) > T(w,3) ) = [];
        if ~isempty(w)
        p1  = T(w,1); p2  = T(w,2); p3  = T(w,3); p4  = T(w,4); try, p12 = P( ET(w, 1 ) ); end; try, p13 = P( ET(w, 2 ) ); end; try, p14 = P( ET(w, 3 ) ); end; try, p23 = P( ET(w, 4 ) ); end; try, p24 = P( ET(w, 5 ) ); end; try, p34 = P( ET(w, 6 ) ); end; ET(w,:) = 0;
        T = [ T ;    p12  ,  p2   ,  p23  ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 1223.1 ];
        T = [ T ;    p1   ,  p12  ,  p23  ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 1223.2 ];
        T = [ T ;    p23  ,  p3   ,  p1   ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 1223.3 ];
        end
        w = find( all(  bsxfun(@eq, LET , ~~[1 0 0 1 0 0] ) ,2) ); w( T(w,1) < T(w,3) ) = [];
        if ~isempty(w)
        p1  = T(w,1); p2  = T(w,2); p3  = T(w,3); p4  = T(w,4); try, p12 = P( ET(w, 1 ) ); end; try, p13 = P( ET(w, 2 ) ); end; try, p14 = P( ET(w, 3 ) ); end; try, p23 = P( ET(w, 4 ) ); end; try, p24 = P( ET(w, 5 ) ); end; try, p34 = P( ET(w, 6 ) ); end; ET(w,:) = 0;
        T = [ T ;    p12  ,  p2   ,  p23  ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 1223.4 ];
        T = [ T ;    p3   ,  p12  ,  p23  ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 1223.5 ];
        T = [ T ;    p12  ,  p3   ,  p1   ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 1223.6 ];
        end

        %tetras to be divided at edge 1-3 & 1-4
        w = find( all(  bsxfun(@eq, LET , ~~[0 1 1 0 0 0] ) ,2) ); w( T(w,3) > T(w,4) ) = [];
        if ~isempty(w)
        p1  = T(w,1); p2  = T(w,2); p3  = T(w,3); p4  = T(w,4); try, p12 = P( ET(w, 1 ) ); end; try, p13 = P( ET(w, 2 ) ); end; try, p14 = P( ET(w, 3 ) ); end; try, p23 = P( ET(w, 4 ) ); end; try, p24 = P( ET(w, 5 ) ); end; try, p34 = P( ET(w, 6 ) ); end; ET(w,:) = 0;
        T = [ T ;    p1   ,  p2   ,  p13  ,  p14  ]; F = [ F ; w ]; %C = [ C ; w*0 + 1314.1 ];
        T = [ T ;    p3   ,  p2   ,  p14  ,  p13  ]; F = [ F ; w ]; %C = [ C ; w*0 + 1314.2 ];
        T = [ T ;    p14  ,  p2   ,  p3   ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 1314.3 ];
        end
        w = find( all(  bsxfun(@eq, LET , ~~[0 1 1 0 0 0] ) ,2) ); w( T(w,3) < T(w,4) ) = [];
        if ~isempty(w)
        p1  = T(w,1); p2  = T(w,2); p3  = T(w,3); p4  = T(w,4); try, p12 = P( ET(w, 1 ) ); end; try, p13 = P( ET(w, 2 ) ); end; try, p14 = P( ET(w, 3 ) ); end; try, p23 = P( ET(w, 4 ) ); end; try, p24 = P( ET(w, 5 ) ); end; try, p34 = P( ET(w, 6 ) ); end; ET(w,:) = 0;
        T = [ T ;    p1   ,  p2   ,  p13  ,  p14  ]; F = [ F ; w ]; %C = [ C ; w*0 + 1314.4 ];
        T = [ T ;    p4   ,  p2   ,  p14  ,  p13  ]; F = [ F ; w ]; %C = [ C ; w*0 + 1314.5 ];
        T = [ T ;    p13  ,  p2   ,  p3   ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 1314.6 ];
        end

        %tetras to be divided at edge 1-4 & 3-4
        w = find( all(  bsxfun(@eq, LET , ~~[0 0 1 0 0 1] ) ,2) ); w( T(w,1) > T(w,3) ) = [];
        if ~isempty(w)
        p1  = T(w,1); p2  = T(w,2); p3  = T(w,3); p4  = T(w,4); try, p12 = P( ET(w, 1 ) ); end; try, p13 = P( ET(w, 2 ) ); end; try, p14 = P( ET(w, 3 ) ); end; try, p23 = P( ET(w, 4 ) ); end; try, p24 = P( ET(w, 5 ) ); end; try, p34 = P( ET(w, 6 ) ); end; ET(w,:) = 0;
        T = [ T ;    p14  ,  p2   ,  p34  ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 1434.1 ];
        T = [ T ;    p1   ,  p2   ,  p34  ,  p14  ]; F = [ F ; w ]; %C = [ C ; w*0 + 1434.2 ];
        T = [ T ;    p34  ,  p2   ,  p1   ,  p3   ]; F = [ F ; w ]; %C = [ C ; w*0 + 1434.3 ];
        end
        w = find( all(  bsxfun(@eq, LET , ~~[0 0 1 0 0 1] ) ,2) ); w( T(w,1) < T(w,3) ) = [];
        if ~isempty(w)
        p1  = T(w,1); p2  = T(w,2); p3  = T(w,3); p4  = T(w,4); try, p12 = P( ET(w, 1 ) ); end; try, p13 = P( ET(w, 2 ) ); end; try, p14 = P( ET(w, 3 ) ); end; try, p23 = P( ET(w, 4 ) ); end; try, p24 = P( ET(w, 5 ) ); end; try, p34 = P( ET(w, 6 ) ); end; ET(w,:) = 0;
        T = [ T ;    p14  ,  p2   ,  p34  ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 1434.4 ];
        T = [ T ;    p3   ,  p2   ,  p34  ,  p14  ]; F = [ F ; w ]; %C = [ C ; w*0 + 1434.5 ];
        T = [ T ;    p14  ,  p2   ,  p1   ,  p3   ]; F = [ F ; w ]; %C = [ C ; w*0 + 1434.6 ];
        end

        %tetras to be divided at edge 1-4 & 2-4
        w = find( all(  bsxfun(@eq, LET , ~~[0 0 1 0 1 0] ) ,2) ); w( T(w,1) > T(w,2) ) = [];
        if ~isempty(w)
        p1  = T(w,1); p2  = T(w,2); p3  = T(w,3); p4  = T(w,4); try, p12 = P( ET(w, 1 ) ); end; try, p13 = P( ET(w, 2 ) ); end; try, p14 = P( ET(w, 3 ) ); end; try, p23 = P( ET(w, 4 ) ); end; try, p24 = P( ET(w, 5 ) ); end; try, p34 = P( ET(w, 6 ) ); end; ET(w,:) = 0;
        T = [ T ;    p14  ,  p24  ,  p3   ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 1424.1 ];
        T = [ T ;    p1   ,  p24  ,  p3   ,  p14  ]; F = [ F ; w ]; %C = [ C ; w*0 + 1424.2 ];
        T = [ T ;    p24  ,  p1   ,  p3   ,  p2   ]; F = [ F ; w ]; %C = [ C ; w*0 + 1424.3 ];
        end
        w = find( all(  bsxfun(@eq, LET , ~~[0 0 1 0 1 0] ) ,2) ); w( T(w,1) < T(w,2) ) = [];
        if ~isempty(w)
        p1  = T(w,1); p2  = T(w,2); p3  = T(w,3); p4  = T(w,4); try, p12 = P( ET(w, 1 ) ); end; try, p13 = P( ET(w, 2 ) ); end; try, p14 = P( ET(w, 3 ) ); end; try, p23 = P( ET(w, 4 ) ); end; try, p24 = P( ET(w, 5 ) ); end; try, p34 = P( ET(w, 6 ) ); end; ET(w,:) = 0;
        T = [ T ;    p14  ,  p24  ,  p3   ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 1424.4 ];
        T = [ T ;    p2   ,  p24  ,  p3   ,  p14  ]; F = [ F ; w ]; %C = [ C ; w*0 + 1424.5 ];
        T = [ T ;    p14  ,  p1   ,  p3   ,  p2   ]; F = [ F ; w ]; %C = [ C ; w*0 + 1424.6 ];
        end

        %tetras to be divided at edge 1-3 & 3-4
        w = find( all(  bsxfun(@eq, LET , ~~[0 1 0 0 0 1] ) ,2) ); w( T(w,1) > T(w,4) ) = [];
        if ~isempty(w)
        p1  = T(w,1); p2  = T(w,2); p3  = T(w,3); p4  = T(w,4); try, p12 = P( ET(w, 1 ) ); end; try, p13 = P( ET(w, 2 ) ); end; try, p14 = P( ET(w, 3 ) ); end; try, p23 = P( ET(w, 4 ) ); end; try, p24 = P( ET(w, 5 ) ); end; try, p34 = P( ET(w, 6 ) ); end; ET(w,:) = 0;
        T = [ T ;    p13  ,  p2   ,  p3   ,  p34  ]; F = [ F ; w ]; %C = [ C ; w*0 + 1334.1 ];
        T = [ T ;    p1   ,  p2   ,  p13  ,  p34  ]; F = [ F ; w ]; %C = [ C ; w*0 + 1334.2 ];
        T = [ T ;    p34  ,  p2   ,  p4   ,  p1   ]; F = [ F ; w ]; %C = [ C ; w*0 + 1334.3 ];
        end
        w = find( all(  bsxfun(@eq, LET , ~~[0 1 0 0 0 1] ) ,2) ); w( T(w,1) < T(w,4) ) = [];
        if ~isempty(w)
        p1  = T(w,1); p2  = T(w,2); p3  = T(w,3); p4  = T(w,4); try, p12 = P( ET(w, 1 ) ); end; try, p13 = P( ET(w, 2 ) ); end; try, p14 = P( ET(w, 3 ) ); end; try, p23 = P( ET(w, 4 ) ); end; try, p24 = P( ET(w, 5 ) ); end; try, p34 = P( ET(w, 6 ) ); end; ET(w,:) = 0;
        T = [ T ;    p13  ,  p2   ,  p3   ,  p34  ]; F = [ F ; w ]; %C = [ C ; w*0 + 1334.4 ];
        T = [ T ;    p4   ,  p2   ,  p13  ,  p34  ]; F = [ F ; w ]; %C = [ C ; w*0 + 1334.5 ];
        T = [ T ;    p13  ,  p2   ,  p4   ,  p1   ]; F = [ F ; w ]; %C = [ C ; w*0 + 1334.6 ];
        end

        %tetras to be divided at edge 1-3 & 2-3
        w = find( all(  bsxfun(@eq, LET , ~~[0 1 0 1 0 0] ) ,2) ); w( T(w,1) > T(w,2) ) = [];
        if ~isempty(w)
        p1  = T(w,1); p2  = T(w,2); p3  = T(w,3); p4  = T(w,4); try, p12 = P( ET(w, 1 ) ); end; try, p13 = P( ET(w, 2 ) ); end; try, p14 = P( ET(w, 3 ) ); end; try, p23 = P( ET(w, 4 ) ); end; try, p24 = P( ET(w, 5 ) ); end; try, p34 = P( ET(w, 6 ) ); end; ET(w,:) = 0;
        T = [ T ;    p13  ,  p23  ,  p3   ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 1323.1 ];
        T = [ T ;    p1   ,  p23  ,  p13  ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 1323.2 ];
        T = [ T ;    p23  ,  p1   ,  p2   ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 1323.3 ];
        end
        w = find( all(  bsxfun(@eq, LET , ~~[0 1 0 1 0 0] ) ,2) ); w( T(w,1) < T(w,2) ) = [];
        if ~isempty(w)
        p1  = T(w,1); p2  = T(w,2); p3  = T(w,3); p4  = T(w,4); try, p12 = P( ET(w, 1 ) ); end; try, p13 = P( ET(w, 2 ) ); end; try, p14 = P( ET(w, 3 ) ); end; try, p23 = P( ET(w, 4 ) ); end; try, p24 = P( ET(w, 5 ) ); end; try, p34 = P( ET(w, 6 ) ); end; ET(w,:) = 0;
        T = [ T ;    p13  ,  p23  ,  p3   ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 1323.4 ];
        T = [ T ;    p4   ,  p23  ,  p2   ,  p13  ]; F = [ F ; w ]; %C = [ C ; w*0 + 1323.5 ];
        T = [ T ;    p13  ,  p1   ,  p2   ,  p4   ]; F = [ F ; w ]; %C = [ C ; w*0 + 1323.6 ];
        end


        ET( ~any( ET ,2) ,:) = [];
        ET = unique( ~~ET , 'rows' );
        if ~isempty( ET )
          error( 'MeshSubdivide:incompleteTet' , ...
                 'tetra subdivision incomplete: %d edge-split pattern(s) not handled.' , size(ET,1) );
        end
        
      case {'linear4'}
        % 1 -> 4 split: insert the cell centroid and fan to the 4 faces. [linear]
        if isempty( W ), return; end

        %indexes of the new points
        P = ( ( nP + 1):( nP + numel(W) ) ).';

        %middle points on faces
        for f = fieldnames( M ).', if ~strncmp( f{1} , 'xyz' , 3 ), continue; end
          M.(f{1}) = [ M.(f{1}) ; ( M.(f{1})( T(W,1) ,:) + M.(f{1})( T(W,2) ,:) + M.(f{1})( T(W,3) ,:) + M.(f{1})( T(W,4) ,:) )/4 ];
        end

        T = [ T ; ...
              [  P        , T(W, 2 ) , T(W, 3 ) , T(W, 4 ) ] ;...
              [  T(W, 1 ) , P        , T(W, 3 ) , T(W, 4 ) ] ;...
              [  T(W, 1 ) , T(W, 2 ) , P        , T(W, 4 ) ] ;...
              [  T(W, 1 ) , T(W, 2 ) , T(W, 3 ) , P        ] ;...
            ];
        F = [ F ; W ; W ; W ; W ];
        
      otherwise, error('Only ''default'' SubType is valid for tetrahedra.');
    end

  else
    error( 'MeshSubdivide:celltype' , ...
           'unsupported celltype %g (only 3=segments, 5=triangles, 10=tetrahedra).' , M.celltype );
  end

  M.tri      = T;
  F(W)       = [];
  M.tri(W,:) = [];            %remove the original faces
  %try, C(W,:)= []; end
  
  
  [F,ord] = sort( F );      %reorder the new faces in their "original position"
  M.tri = M.tri(ord,:);
  %try, C     = C(ord); end
  
  for f = fieldnames( M ).'
    if strcmp( f{1} , 'tri' ), continue; end
    if ~strncmp( f{1} , 'tri' , 3 ), continue; end
    M.(f{1}) = M.(f{1})(F,:,:,:,:,:,:);
  end
  %try, M.triCASE = C; end


  
  function add2G( aaa , bbb , vvv )
    nnn = numel(bbb);
    if ~nnn, return; end
    
    bbb = bbb(:);
    if numel(aaa) > 1 && numel(aaa) ~= nnn
      aaa = repmat( aaa , ceil( size(bbb)./size(aaa) ) );
      aaa = aaa(:); aaa = aaa(1:nnn);
    end
    if numel(vvv) ~= nnn && numel(vvv) > 1
      vvv = repmat( vvv , ceil( size(bbb)./size(vvv) ) );
      vvv = vvv(:); vvv = vvv(1:nnn);
    end

%     switch nnn
%       case 1,     iii = GN + ( 1 );
%       case 2,     iii = GN + ( 1:2 );
%       case 3,     iii = GN + ( 1:3 );
%       case 4,     iii = GN + ( 1:4 );
%       case 5,     iii = GN + ( 1:5 );
%       case 6,     iii = GN + ( 1:6 );
%       case 7,     iii = GN + ( 1:7 );
%       case 8,     iii = GN + ( 1:8 );
%       case 9,     iii = GN + ( 1:9 );
%       case 10,    iii = GN + ( 1:10 );
%       case 11,    iii = GN + ( 1:11 );
%       case 12,    iii = GN + ( 1:12 );
%       case 13,    iii = GN + ( 1:13 );
%       case 14,    iii = GN + ( 1:14 );
%       case 15,    iii = GN + ( 1:15 );
%       case 16,    iii = GN + ( 1:16 );
%       case 17,    iii = GN + ( 1:17 );
%       case 18,    iii = GN + ( 1:18 );
%       case 19,    iii = GN + ( 1:19 );
%       case 20,    iii = GN + ( 1:20 );
%       otherwise,  
        iii = ( GN + 1 ):( GN + nnn );
%     end
    
    sz_G = size(G,1);
    if sz_G <= iii(end)
      %G( ( sz_G+1 ):( iii(end)*2 ) ,:) = NaN;
      G( end+1:(iii(end)*2) ,:) = NaN;
    end
    
    G( iii ,1) = aaa(:);
    G( iii ,2) = bbb(:);
    G( iii ,3) = vvv(:);
    
    GN = GN + nnn;
  end
  
end

function A = uniqueROWS( A , cols )
  if isempty( A ), return; end
  A = sortROWS( A , cols );
  w = [ true ; any( diff( A , 1 , 1 ) ,2) ];
  A = A(w,:);
end
function A = sortROWS( A , cols )
  for c = cols(end:-1:1)
    [~,ord] = sort( A(:,c) );
    A = A(ord,:);
  end
end
function varargout = ismember2ROWS( A , B )
  try
    if size( A ,2) == 2
      A = A(:,1) + 1i*A(:,2);
    elseif size(A,2) > 2
      error('A has more than 2 columns');
    end
    if size( B ,2) == 2
      B = B(:,1) + 1i*B(:,2);
    elseif size(B,2) > 2
      error('B has more than 2 columns');
    end

    [varargout{1:nargout}] = ismember( A , B );
  catch
    [varargout{1:nargout}] = ismember( A , B , 'rows' );
  end
end
function R = CIRCULARshift( R , b )
  if R(1) ~= R(end)
    R = [];
  else
    id = find( R == b );
    R = R( [ id:end-1 , 1:id-1 ] );
  end
end
function N = rowUnit( N )
%normalize each row to unit length; zero-norm rows are left untouched.
  nn = sqrt( sum( N.^2 , 2 ) );
  k  = nn > 0;
  N(k,:) = N(k,:) ./ nn(k);
end
function str = renameStructField( str , oldFieldName , newFieldName )
%RENAMESTRUCTFIELD  rename a struct field, preserving field order (MathWorks FEX).
  if ~strcmp( oldFieldName , newFieldName )
    allNames = fieldnames( str );
    isOverwriting = ~isempty( find( strcmp( allNames , newFieldName ) , 1 ) );
    matchingIndex = find( strcmp( allNames , oldFieldName ) );
    if ~isempty( matchingIndex )
      allNames{ matchingIndex(1) } = newFieldName;
      [ str.(newFieldName) ] = deal( str.(oldFieldName) );
      str = rmfield( str , oldFieldName );
      if ~isOverwriting
        str = orderfields( str , allNames );
      end
    end
  end
end