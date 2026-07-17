function N = meshNormals( M , mode )
%MESHNORMALS  Face or vertex normals of a triangle/segment mesh (unit length).
%
%   N = meshNormals( M )            per-FACE normals (one row per M.tri row)
%   N = meshNormals( M , false )    same (also numeric 0)
%   N = meshNormals( M , k )        per-FACE normals SMOOTHED k times over the
%                                   edge-adjacency graph (k a positive integer).
%                                   CONTINUES from M.triNORMALS when present --
%                                   k passes then j passes == k+j passes,
%                                   bit-exact; without the field it starts from
%                                   the raw face normals
%   N = meshNormals( M , 'uniform' )  per-VERTEX normals, equal weight per face
%   N = meshNormals( M , 'area'    )  per-VERTEX normals, area-weighted        [2]
%   N = meshNormals( M , 'angle'   )  per-VERTEX normals, incident-angle-weighted
%                                     -- the "pseudonormal": the CORRECT normal
%                                     for the inside/outside signed-distance sign
%                                     test.                                  [1,3]
%   N = meshNormals( M , 'best'    )  per-VERTEX, robust variant (see BEST below)
%   N = meshNormals( M , 'reciprocal' ) per-VERTEX, reciprocal-weighted: the
%                                     ACCURACY option. Segments: 1/length --
%                                     the smooth-curve normal, EXACT on circles
%                                     under any sampling. Triangles: Max's
%                                     corner weights sin(angle)/(|e1||e2|) --
%                                     EXACT on any mesh inscribed in a sphere
%                                     (probe: 1.5e-6 vs 6..13 deg of the other
%                                     modes on an irregular random sphere)  [2]
%   N = meshNormals( M , 'quadratic'  ) per-VERTEX, the normal of the chord-length
%                                     PARABOLA through the vertex and its chain
%                                     neighbours (segment meshes only)
%   N = meshNormals( M , true )       == 'uniform'
%
% M is a mesh struct (.xyz nodes, .tri cells). Vertex modes aggregate the face
% normals with meshF2V; if M.triNORMALS exists it is used as the face normals
% -- TRUSTED blindly: remove or rebuild the field after editing the mesh (a
% wrong ROW COUNT errors as a stale field, but a same-size field over moved
% nodes cannot be detected). Triangle meshes (celltype 5) and segment meshes /
% polylines (celltype 3) are supported. Every returned row has norm exactly 1,
% counting NaN marker entries as 0 (see normalizeRows).
%
% NOISY meshes -- the usual recipe: smooth the FACE normals first, then build
% the vertex normals FROM them. EVERY mode picks up M.triNORMALS when present;
% in particular the numeric mode CONTINUES smoothing it, so consecutive calls
% ACCUMULATE (5 then 10 == 15, bit-exact) and progressive smoothing needs no
% recomputation:
%     M.triNORMALS = meshNormals( M );       %seed = raw face normals (RESET)
%     M.triNORMALS = meshNormals( M , 5 );   %5 neighbour-averaging passes
%     M.triNORMALS = meshNormals( M , 10 );  %10 more: 15 accumulated
%     Nv = meshNormals( M , 'area' );        %vertex normals from the smoothed field
% Re-running the whole block is safe (the seed line restarts it); re-running
% the k-line alone keeps smoothing. Since k passes diffuse WHATEVER lives in
% M.triNORMALS, storing any per-face direction field there (fibers, gradients)
% and calling meshNormals(M,k) smooths THAT field over the mesh -- directions
% only: rows are renormalized to unit length every pass.
% Verified: noisy sphere (2% radial noise) 7.7 -> 5.0 deg mean vertex error;
% noisy 2-D circle (1%) 4.8 -> 1.5 deg. Segment meshes: same recipe, same
% mechanics ('quadratic' is the exception: it refits the geometry and only
% takes the sign / fallbacks from the field). CAVEAT for genuinely-3-D curves:
% it helps only while the noise stays below the per-segment curvature
% deflection, sigma <~ h^2*kappa/2 (there, 3-D circle: faces 6.6 -> vertices
% 0.6 deg). Above it the LOCAL osculating planes see noise instead of
% curvature and NO normal smoothing can recover the truth (sigma = h/6: faces
% 46 deg, still 42 after k=25) -- smooth the GEOMETRY instead. 2-D and
% coplanar pieces are immune: their plane is global, not local.
%
% SEGMENTS  (celltype 3, per-face normals)
%   A line segment in 3-D has no single normal but a whole CIRCLE of them (the
%   plane perpendicular to it), so a normal is well defined only RELATIVE to a
%   plane. meshNormals splits the mesh into CONNECTED COMPONENTS (meshSeparate,
%   joining segments that share a node) and decides PER COMPONENT, using a
%   node-to-line / node-to-plane distance threshold TH = 1e-8:
%     * 2-D input (M.xyz has 2 columns): every normal is the in-plane
%       perpendicular, returned 2-D (N has 2 columns), sign following p1->p2.
%       There is no ambiguity (the ambient plane is fixed) and nothing is marked.
%     * a COPLANAR component (nodes within TH of a single plane, not collinear):
%       the in-plane perpendicular w.r.t. that plane (a normal 3-D vector).
%     * a genuinely 3-D component: each segment's perpendicular within the LOCAL
%       osculating plane fitted to it and its node-sharing neighbours. The
%       ORIENTATION of those local planes -- which the fit alone leaves
%       arbitrary: the same curve in another pose used to flip sides along its
%       length -- is made consistent by PROPAGATING it over the component from
%       the most planar neighbourhood (itself oriented '+z', like the coplanar
%       case, so nearly-coplanar curves agree with the coplanar branch).
%       Straight sub-runs INSIDE such a component (neighbourhoods that fix no
%       plane) INHERIT the plane of their curved neighbours: continuous
%       normals, no marks, no warnings.
%     * a SINGLE segment, or a STRAIGHT / collinear polyline (nodes within TH of
%       a line): the normal is UNDEFINED -- any perpendicular is equally valid.
%       One is chosen (parallel to the XY plane, or to the XZ plane if the piece
%       is vertical) and MARKED by setting its z component to NaN, so callers can
%       tell it was an arbitrary pick. (Only 3-D input reaches this: for a 2-D
%       mesh even a lone segment has a well-defined in-plane normal.)
%   CHIRALITY: in every branch the normal follows its OWN segment, p1->p2
%   rotated -90 deg within the (oriented) plane -- the right-hand rule, exactly
%   as in 2-D. Reversing one segment's connectivity flips THAT normal; the
%   propagation above orients the planes, never the normals, so it preserves
%   these per-segment flips (they are information, not noise).
%
%   VERTEX modes on segment meshes (the analogs keep each mode's MEANING):
%     'angle' == 'uniform': at a polyline vertex the pseudonormal [1,3] reduces
%         to the plain BISECTOR of the two incident normals (the Gauss-map arc
%         integral = the max-margin direction: margin cos45 at a right corner,
%         vs ~0.05 for the length weightings), so 'angle' simply maps to it.
%     'area' -> LENGTH-weighted (the measure analog == the normal of the
%         central-difference tangent P2-P0). NOTE it is the WORST choice for
%         recovering a smooth curve's normal (2x the bisector's error);
%     'reciprocal' (1/length) is the accuracy one: EXACT radial normals on a
%         circle whatever the (non-uniform) sampling, first-order exact
%         elsewhere -- the discrete curvature-vector direction.
%     'quadratic': tangent of the chord-length-parametrized parabola through
%         the vertex and its chain neighbours (closed form: u0/L0 + u1/L1, so it
%         EQUALS 'reciprocal' on 2-D / coplanar pieces), normal taken in the
%         3-point plane -- differs from 'reciprocal' only under torsion. Open
%         ends use the ONE-SIDED parabola through the first 3 chain points;
%         branching vertices (valence > 2), collinear stencils and degenerate
%         cases fall back to the bisector; the sign follows the incident face
%         normals (chain-order invariant).
%     'best': same machinery as triangles, seeded with the bisector (at
%         valence-2 vertices the bisector already IS the max-margin direction;
%         miniball only matters at branchings).
%   Vertex normals over NaN-MARKED pieces stay marked: unit finite part, NaN z.
%
% BEST  ('best')
%   Uses the angle-weighted pseudonormal [1,3] -- arithmetically the correct
%   normal for the inside/outside sign test -- EXCEPT at vertices where it falls
%   OUTSIDE the cone of the incident face normals (skewed / spiky fans, where its
%   margin collapses and the sign test turns unreliable). There it is replaced by
%   the MAX-MARGIN direction of that cone,
%       u* = argmax_{|u|=1}  min_i ( n_i . u ),
%   i.e. the axis of the smallest spherical cap enclosing the incident unit
%   normals. That axis equals the NORMALIZED CENTER of the minimum enclosing BALL
%   of the normal tips: for the ball's support points |n_i-c| = r, so
%   n_i.c = (1+|c|^2-r^2)/2 is EQUAL for all of them -- hence c/|c| is the point
%   of equal, maximal minimum dot. So u* is obtained EXACTLY (no iteration) as the
%   normalized minimum-enclosing-ball center via miniball [4]. A vertex whose cone
%   is EMPTY (normals spanning more than a half-space => cap radius > 90 deg) has
%   NO valid decision direction; there the angle-weighted normal is kept and one
%   warning reports how many.
%   REQUIRES the compiled miniball MEX (miniball_mx) for that fallback; smooth
%   meshes never reach it ( best == angle ).
%
% REFERENCES
%   [1] G. Thurmer, C. A. Wuthrich, "Computing vertex normals from polygonal
%       facets", Journal of Graphics Tools 3(1):43-46, 1998.  (angle weighting)
%   [2] N. Max, "Weights for computing vertex normals from facet normals",
%       Journal of Graphics Tools 4(2):1-6, 1999.  (area & alternative weights;
%       his sin/(|e1||e2|) corner weighting = the triangle 'reciprocal' mode,
%       exact on sphere-inscribed meshes)
%   [3] J. A. Baerentzen, H. Aanaes, "Signed Distance Computation Using the Angle
%       Weighted Pseudonormal", IEEE TVCG 11(3):243-253, 2005.  (proves the
%       angle-weighted vertex/edge pseudonormal yields the CORRECT inside/outside
%       sign -- the reason 'angle'/'best' are the right default for that test).
%   [4] B. Gartner, "Fast and Robust Smallest Enclosing Balls", ESA 1999,
%       LNCS 1643:325-338.  (the exact minimum-enclosing-ball whose normalized
%       center is the max-margin cap axis used by 'best'; = miniball_mx here).
%
% See also meshF2V, meshEsuE, meshEsuP, miniball, normalizeRows.

  if nargin < 2 || isempty( mode ), mode = false; end

  if ~ischar(mode) && isscalar( mode ) && isnumeric( mode ) && mode == 0, mode = false; end   %numeric 0 = face normals (NaN falls to the validation below)
  if islogical( mode ) && mode
    mode = 'uniform';
  end

  if ischar( mode )
    if isfield( M , 'triNORMALS' )
      N = M.triNORMALS;
      if size( N ,1) ~= size( M.tri ,1)
        error( 'meshNormals:triNORMALS' , ...
          'M.triNORMALS has %d rows but the mesh has %d faces: stale field from another mesh? Remove it, or rebuild it with meshNormals(M,false).' , ...
          size( N ,1) , size( M.tri ,1) );
      end
    else
      N = meshNormals( M ,false );
    end

    ct = meshCelltype( M );

    switch lower(mode)
      case {'u','uniform'},   N = meshF2V( M , N ,'sum'    );
      case {'g','angle'  }
        if ct == 3   %at a polyline vertex the pseudonormal IS the bisector (see help)
          N = meshF2V( M , N ,'sum'    );
        else
          N = meshF2V( M , N ,'angles' );
        end
      case {'a','area'   }
        if ct == 3   %the measure analog for curves: length weighting
          N = meshF2V( M , N ,'length' );
        else
          N = meshF2V( M , N ,'area'   );
        end
      case {'r','reciprocal','reciproco'}
        switch ct
          case 3      %1/length: the smooth-curve normal, exact on circles
            N = meshF2V( M , bsxfun( @rdivide , N , meshQuality( M ,'length' ) ) ,'sum' );
          case 5      %Max's corner weights [2]: exact on sphere-inscribed meshes
            N = triReciprocalVertexNormals( M , N );
          otherwise
            error( 'meshNormals:mode' , '''reciprocal'' is only defined for segment (celltype 3) and triangle (celltype 5) meshes.' );
        end
      case {'q','quadratic'}
        if ct ~= 3
          error( 'meshNormals:mode' , '''quadratic'' is only defined for segment meshes (celltype 3).' );
        end
        N = segQuadraticVertexNormals( M , N );
      case {'best'}
        % angle-weighted pseudonormal (the exact-arithmetic correct one for the
        % inside/outside sign test), EXCEPT where it falls outside the cone of
        % its incident face normals (skewed sharp spikes: its margin collapses
        % and the sign test becomes numerically unreliable). There it is
        % replaced by the MAXIMUM-MARGIN direction of the cone,
        %     u* = argmax_u min_i ( n_i . u )   (the spherical 1-center),
        % computed EXACTLY as the normalized center of the minimum enclosing
        % ball of the unit normal tips -- no iterative optimization. A vertex
        % whose cone is EMPTY (normals spanning more than a halfspace) admits
        % NO single decision direction at all: the angle-weighted normal is
        % kept there and a single warning reports how many.
        if ct == 3   %valence-2 bisector == max margin already; miniball acts at branchings
          X = meshF2V( M , N , 'sum'    );
        else
          X = meshF2V( M , N , 'angles' );
        end

        % vectorized screening: vertices with some incident face normal at
        % negative dot with their angle-weighted normal
        FI  = repmat( ( 1:size(M.tri,1) ).' , size(M.tri,2) , 1 );
        bad = find( accumarray( double( M.tri(:) ) , ...
                double( sum( N(FI,:) .* X(M.tri(:),:) ,2) < 0 ) , [size(X,1),1] , @max ) );

        if ~isempty( bad )
          L = meshEsuP( M , 0 );
          nbad = 0;
          for x = bad(:).'
            NS = N( L(:,x) ,:);
            c  = miniball( NS );  c = c(:).';       %max-margin direction (cap center)
            if all( NS * c.' > 0 )
              X(x,:) = c;
            else
              nbad = nbad + 1;                      %empty cone: keep angle-weighted
            end
          end
          if nbad
            warning( 'meshNormals:emptyNormalCone' , ...
              [ '%d vertices have an EMPTY normal cone (no single direction can make the ' , ...
                'sign test valid there); the angle-weighted normal was kept.' ] , nbad );
          end
        end
        N = X;

      otherwise
        error( 'meshNormals:mode' , [ 'invalid vertex weighting ''%s'' (use ''uniform'', ''angle'', ' , ...
               '''area'', ''best'', or -- segments only -- ''reciprocal'', ''quadratic'').' ] , mode );
    end
    
    N = normalizeRows( N );
    return;
  end
  if isnumeric( mode )
    % numeric mode = number of smoothing iterations of the FACE normals: it must
    % be a clean count. Anything else (negative, Inf, NaN, fractional, or
    % non-scalar) used to fall through SILENTLY to the raw face normals -- or to
    % hang, for Inf -- so now it is an error.
    if ~isscalar( mode ) || ~isfinite( mode ) || mode <= 0 || mode ~= round( mode )
      error( 'meshNormals:mode' , ...
        'a numeric mode (number of smoothing iterations) must be a finite positive integer scalar.' );
    end
    if isfield( M , 'triNORMALS' )
      N = M.triNORMALS;         %CONTINUE from the stored field: the iteration is
                                %memoryless, so k passes then j passes == k+j
                                %passes, bit-exact. meshNormals(M,false) resets.
      if size( N ,1) ~= size( M.tri ,1)
        error( 'meshNormals:triNORMALS' , ...
          'M.triNORMALS has %d rows but the mesh has %d faces: stale field from another mesh? Remove it, or rebuild it with meshNormals(M,false).' , ...
          size( N ,1) , size( M.tri ,1) );
      end
    else
      N = meshNormals( M ,false );
    end
    if size( M.tri ,1) == 0, return; end   %empty mesh: nothing to smooth

    NEIGS = double( meshEsuE( M , false , 'n' ) );
    NEIGS = NEIGS + speye( size(NEIGS) );
    NEIGS = NEIGS ./ full( sum( NEIGS ,2) );
    for it = 1:mode
      N0 = N;
      N = NEIGS * N;
      N = normalizeRows( N );
      w = any( ~isfinite( N ) ,2);
      N(w,:) = N0(w,:);
    end
    %NO final normalizeRows here: every pass already normalizes, and the output
    %must be EXACTLY the loop state so that continuing from M.triNORMALS is
    %bit-identical to having run all the passes in one call (an extra polish
    %pass here re-touched the rare +-1 ulp rows and broke that, by 1 ulp).
    return;
  end


  M.celltype = meshCelltype( M );

  switch M.celltype
    case 3

      if size( M.xyz ,2) < 3
        %----- 2-D segment mesh: the normal is simply the in-plane perpendicular
        %      (dx,dy) -> (dy,-dx). Unambiguous, 2-D, nothing to separate or mark.
        N = ( M.xyz( M.tri(:,2) ,:) - M.xyz( M.tri(:,1) ,:) ) * [0 -1;1 0];

      elseif size( M.tri ,1) == 0
        %----- EMPTY segment mesh: nothing to compute (and meshSeparate -> graph
        %      chokes on 0 cells with a cryptic error).
        N = zeros( 0 , 3 );

      else
        %----- 3-D segment mesh: decide PER CONNECTED COMPONENT (meshSeparate).
        %      Each component is classified (collinear -> ambiguous+NaN mark,
        %      coplanar -> in-plane, else -> local osculating) in segNormals3D.
        TH = 1e-8;                                       %node-to-line / node-to-plane tolerance
        M.triID4NORMALS__ = ( 1:size( M.tri ,1) ).';     %label the segments to scatter back
        MS = meshSeparate( M );                          %connected pieces; each keeps its label
        N  = NaN( size( M.tri ,1) , 3 );
        for i = 1:numel( MS )
          N( MS{i}.triID4NORMALS__ ,:) = segNormals3D( MS{i} , TH );
        end
      end

    case 5
      M.xyz(:,end+1:3) = 0;
      N = cross( M.xyz( M.tri(:,2) ,:) - M.xyz( M.tri(:,1) ,:) , M.xyz( M.tri(:,3) ,:) - M.xyz( M.tri(:,1) ,:) , 2 );
      
    otherwise
      error( 'meshNormals:celltype' , ...
        'normals are implemented for segment (celltype 3) and triangle (celltype 5) meshes only, got celltype %d (for tetrahedra, extract the boundary surface first).' , M.celltype );
  end

  N = normalizeRows( N );
  
end
function N = normalizeRows( N )
% rows are RE-normalized (up to 5 passes) until their norm is EXACTLY 1, with NaN
% entries counting as 0 towards the norm and staying in place: the NaN-marked
% "ambiguous" segment normals thus come out with a UNIT finite part (also after
% vertex aggregation, where several marks may have been summed), never all-NaN.
% A single division does not always land on 1 in floating point, so passes after
% the first touch ONLY the rows that are not exact yet -- essentially free (a
% rare row still oscillating at +-1 ulp after the 5 passes is accepted). A
% zero row (degenerate face) turns all-NaN on the first division, as always.
  Z = N;  Z( isnan( Z ) ) = 0;
  nn = sqrt( sum( Z.^2 ,2) );
  w  = find( nn ~= 1 & isfinite( nn ) );
  if isempty( w ), return; end
  N(w,:) = N(w,:) ./ nn( w );
  w  = w( nn( w ) > 0 );          %zero rows just became (and stay) all-NaN
  for it = 1:4
    Z  = N(w,:);  Z( isnan( Z ) ) = 0;
    nn = sqrt( sum( Z.^2 ,2) );
    k  = nn ~= 1;
    if ~any( k ), break; end
    w  = w( k );
    N(w,:) = N(w,:) ./ nn( k );
  end
end


function N = segNormals3D( P , TH )
% Per-face normals of ONE connected 3-D segment component P (P.xyz N-by-3,
% P.tri M-by-2). Classified by the max node distance to a best-fit line / plane:
%   collinear (single segment or straight polyline)  -> AMBIGUOUS: one normal
%       parallel to the XY plane (or XZ if the piece is vertical), z set to NaN.
%   coplanar (nodes within TH of one plane, not a line) -> in-plane perpendicular.
%   otherwise (genuinely 3-D) -> per-segment perpendicular within the LOCAL
%       osculating plane (segment + node-sharing neighbours), the plane
%       ORIENTATIONS made consistent by propagation over the component and
%       straight sub-runs inheriting their neighbours' plane (see below).
  nseg = size( P.tri ,1);
  V    = P.xyz( unique( P.tri ) ,:);

  %----- collinear? all nodes within TH of their best-fit line -> AMBIGUOUS
  A = bsxfun( @minus , V , mean( V ,1) );
  [~,~,W] = svd( A , 0 );
  d = W(:,1).';                                    %best-fit line direction
  perp = A - (A * d.') * d;                        %node offsets perpendicular to that line
  if all( sqrt( sum( perp.^2 ,2) ) <= TH )
    N = repmat( markedNormal( d ) , nseg , 1 );    %same marked normal for the whole straight piece
    return;
  end

  %----- coplanar? all nodes within TH of their best-fit plane -> IN-PLANE
  [Pl,iPl] = getPlane( V ,'+z');
  if all( distance2Plane( V , Pl ) <= TH )
    XY = transform( P.xyz , iPl );
    XY = XY( P.tri(:,2) ,:) - XY( P.tri(:,1) ,:);
    XY = XY * [0 -1 0;1 0 0;0 0 0];
    N  = XY * Pl(1:3,1:3).';
    return;
  end

  %----- genuinely 3-D -> per-segment perpendicular within the LOCAL osculating
  %      plane (svd fit of the segment + its node-sharing neighbours). The fit
  %      leaves each plane's SIGN (binormal direction) arbitrary -- raw svd sign
  %      is pose-dependent and used to flip sides along one same curve -- so the
  %      signs are made consistent by BREADTH-FIRST PROPAGATION over the segment
  %      adjacency graph, seeded at the most planar neighbourhood (largest 2nd
  %      singular value), the seed itself oriented with the same '+z' convention
  %      as the coplanar branch (nearly-coplanar curves thus agree with it). A
  %      neighbourhood that is itself COLLINEAR (a straight run inside the
  %      component) fixes no plane at all: it INHERITS the propagated plane of
  %      its neighbour -- continuous normals instead of the former arbitrary
  %      (and getPlane-warning-flooding) fit of a degenerate neighbourhood.
  %      The normal is then chord x binormal == the chord projected onto its
  %      plane and rotated -90 deg IN it: it keeps following p1->p2, so
  %      reversing one segment flips its normal (2-D-like right-hand chirality)
  %      -- propagation orients the planes, never the normals.
  S = meshEsuE( P , true );
  B = zeros( nseg , 3 );                           %local plane normal (binormal); [0 0 0] = fixes no plane
  Q = zeros( nseg , 1 );                           %2nd singular value: planarity strength (seed pick)
  for k = 1:nseg
    Wk = P.xyz( unique( P.tri( [ S{k} ; k ] ,:) ) ,:);
    A  = bsxfun( @minus , Wk , mean( Wk ,1) );
    [~,D,W] = svd( A , 0 );
    perp = A - ( A * W(:,1) ) * W(:,1).';          %node offsets perpendicular to the local best-fit line
    if all( sqrt( sum( perp.^2 ,2) ) <= TH ), continue; end   %straight run: inherit below
    B(k,:) = W(:,3).';
    Q(k)   = D(2,2);
  end

  chord = P.xyz( P.tri(:,2) ,:) - P.xyz( P.tri(:,1) ,:);

  if ~any( Q )     %pathological: NO neighbourhood fixes a plane, yet the piece as a
    N = zeros( nseg , 3 );                         %whole is neither collinear nor coplanar
    for k = 1:nseg                                 %(needs degeneracies): mark each segment
      N(k,:) = markedNormal( chord(k,:) );         %on its own direction, as a lone straight piece
    end
    return;
  end

  [~,seed] = max( Q );
  if B(seed,3) < 0, B(seed,:) = -B(seed,:); end    %'+z' seed, as the coplanar branch
  seen = false( nseg ,1 );  seen(seed) = true;
  queue = zeros( nseg ,1 );  queue(1) = seed;  qh = 1;  qt = 1;
  while qh <= qt                                   %breadth-first over the adjacency
    i = queue(qh);  qh = qh + 1;
    for j = S{i}.'
      if seen(j), continue; end
      seen(j) = true;
      if ~any( B(j,:) )                            %straight run: inherit the neighbour's plane
        B(j,:) = B(i,:);
      elseif B(j,:) * B(i,:).' < 0                 %flip to the already-visited side
        B(j,:) = -B(j,:);
      end
      qt = qt + 1;  queue(qt) = j;
    end
  end

  N = cross( chord , B , 2 );                      %projected chord rotated -90 deg in its plane
                                                   %(made norm-1 by the caller's normalizeRows)
  for k = find( ~seen ).'                          %unreachable segments: cannot happen out of
    N(k,:) = markedNormal( chord(k,:) );           %meshSeparate (connected), kept as a guard
  end
end


function n = markedNormal( d )
% Ambiguous straight piece of direction d: no plane fixes its normal, so pick one
% perpendicular parallel to the XY plane (or to the XZ plane when d is vertical)
% and MARK it by setting the z component to NaN (the x,y part stays unit).
  d = d / norm( d );
  if hypot( d(1) , d(2) ) > 1e-8
    n = [ d(2) , -d(1) ];  n = n / norm( n );      %perpendicular within the XY plane
    n = [ n , NaN ];
  else                                             %vertical: perpendicular within the XZ plane (y = 0)
    n = [ 1 , 0 , NaN ];
  end
end


function V = triReciprocalVertexNormals( M , NF )
% 'reciprocal' vertex normals of a TRIANGLE mesh: N. Max's weighting [2]. Each
% incident face contributes its (unit) normal NF weighted by
%     w = sin(corner angle) / ( |e1| |e2| )  ==  |e1 x e2| / ( |e1|^2 |e2|^2 )
% with e1,e2 the two edges leaving the vertex in that face -- equivalently, the
% raw corner cross product over the squared edge lengths. This is EXACT (radial
% normals) for ANY mesh inscribed in a sphere, whatever the triangulation: the
% surface mirror of the segment 1/length weighting being exact on circles.
% The weight is per CORNER, not per face (each vertex of a face gets a
% different one), hence this dedicated aggregation instead of meshF2V. Taking
% NF as input (instead of recomputing the crosses) keeps M.triNORMALS honored,
% like every other vertex mode; when NF are the plain face normals both forms
% coincide. Degenerate corners (zero edge / zero area) yield w = NaN|0 and
% poison their vertex, as in the other modes.
  T   = double( M.tri );
  xyz = M.xyz;  xyz(:,end+1:3) = 0;
  nV  = size( xyz ,1);
  V   = zeros( nV , size( NF ,2) );
  for r = 1:3
    v  = T(:,r);   a = T(:, mod(r,3)+1 );   b = T(:, mod(r+1,3)+1 );
    e1 = xyz(a,:) - xyz(v,:);
    e2 = xyz(b,:) - xyz(v,:);
    w  = sqrt( sum( cross( e1 , e2 ,2).^2 ,2) ) ./ ( sum( e1.^2 ,2) .* sum( e2.^2 ,2) );
    for c = 1:size( NF ,2)
      V(:,c) = V(:,c) + accumarray( v , w .* NF(:,c) , [nV,1] );
    end
  end
end


function V = segQuadraticVertexNormals( M , NF )
% 'quadratic' vertex normals of a segment mesh: the tangent at each vertex is
% that of the chord-length-parametrized PARABOLA through the vertex and its two
% chain neighbours. Its closed form at the vertex is  T = u0/L0 + u1/L1  (the
% 1/length sum of the unit chords -- which is why it EQUALS 'reciprocal' on
% 2-D / coplanar pieces), and the normal is the perpendicular to T within the
% 3-point plane, SIGNED to agree with the bisector of the incident face normals
% (that keeps the p1->p2 chirality and makes the result invariant to the chain
% direction). Open ends use the ONE-SIDED parabola through the first 3 chain
% points, with the tangent evaluated AT the end. Fallback = the bisector itself:
% branching vertices (no single curve through), collinear/degenerate stencils,
% lone-segment ends, and ambiguous sign (anti-aligned incident normals).
  xyz  = M.xyz;
  is2D = size( xyz ,2) < 3;
  V = meshF2V( M , NF ,'sum' );                    %bisector: fallback + sign reference
  if size( M.tri ,1) == 0, return; end             %no segments: all-bisector (-> NaN);
                                                   %meshEsuP('cell') on 0 cells returns
                                                   %a DOUBLE (accumarray fill), not a cell
  E = meshEsuP( M ,'cell' );                       %incident segments per vertex

  for v = 1:numel( E )
    switch numel( E{v} )
      case 2                                       %chain interior
        nb = M.tri( E{v} ,:).';  nb = double( nb( nb ~= v ) );
        if numel( nb ) ~= 2, continue; end         %a degenerate [v v] cell -> bisector
        P1 = xyz( v ,:);
        d0 = P1 - xyz( nb(1) ,:);   L0 = norm( d0 );
        d1 = xyz( nb(2) ,:) - P1;   L1 = norm( d1 );
        if L0 == 0 || L1 == 0, continue; end
        u0 = d0 / L0;   u1 = d1 / L1;
        T  = u0/L0 + u1/L1;                        %parabola tangent at the vertex

      case 1                                       %open end: one-sided parabola
        nb = M.tri( E{v} ,:);  nb = double( nb( nb ~= v ) );
        if numel( nb ) ~= 1, continue; end
        en = E{ nb };  en = en( en ~= E{v} );      %the chain's next segment...
        if numel( en ) ~= 1, continue; end         %...unless lone segment / branching
        nx = M.tri( en ,:);  nx = double( nx( nx ~= nb ) );
        if numel( nx ) ~= 1, continue; end
        P1 = xyz( v ,:);  P2 = xyz( nb ,:);  P3 = xyz( nx ,:);
        L1 = norm( P2 - P1 );   L2 = norm( P3 - P2 );
        if L1 == 0 || L2 == 0, continue; end
        a  = L1;   b  = L1 + L2;
        T  = ( ( P2 - P1 )*b^2 - ( P3 - P1 )*a^2 ) / ( a*b*(b - a) );   %Q'(0)
        u0 = ( P2 - P1 ) / L1;   u1 = ( P3 - P2 ) / L2;                 %span the 3-point plane

      otherwise                                    %branching / isolated -> bisector
        continue;
    end

    if is2D
      Nq = [ T(2) , -T(1) ];                       %-90 deg, as everywhere else
    else
      W  = cross( u0 , u1 );                       %the 3-point plane normal
      if norm( W ) <= 1e-8, continue; end          %collinear stencil -> bisector
      Nq = cross( W , T );                         %in-plane perpendicular to the tangent
    end
    s = Nq .* V( v ,:);
    s = sum( s( ~isnan( s ) ) );                   %NaN counts as 0 (marked bisector)
    if     s < 0, V( v ,:) = -Nq;
    elseif s > 0, V( v ,:) =  Nq;
    end                                            %s == 0: ambiguous sign -> keep bisector
  end
end
