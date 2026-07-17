function [M,w] = MeshFixCellOrientation( M , f )
%MESHFIXCELLORIENTATION  Make the cell orientation of a mesh CONSISTENT.
%
%   [M,w] = MeshFixCellOrientation( M )
%   [M,w] = MeshFixCellOrientation( M , f )
%
%   TRIANGLE meshes (celltype 5): each connected component is re-oriented so
%   that every shared edge is traversed in OPPOSITE directions by its two
%   faces (consistent orientation) and, if needed, the whole component is
%   flipped so that its SIGNED VOLUME is positive (outward normals). Every
%   component is treated independently -- a hollow sphere gets BOTH shells
%   oriented outward (by design: the hollow only exists once tetrahedralized).
%   Open sheets have no inside: they come out CONSISTENT, the global side
%   given by the sign of their (open) flux -- arbitrary but deterministic.
%
%   f -- optional list of reference faces (GLOBAL indices into M.tri):
%        f(k) > 0 pins face  f(k) with its CURRENT orientation;
%        f(k) < 0 pins face -f(k) FLIPPED.
%        A component containing a reference follows it and gets NO volume
%        correction (the reference wins over outward); the first reference
%        inside a component is the binding one. Components without any
%        reference are oriented outward as usual.
%
%   TETRAHEDRAL meshes (celltype 10): tets with negative volume get their
%   last two vertices swapped (f is not supported there).
%
%   w -- logical mask of the cells that were FLIPPED w.r.t. the input.
%
%   NOTES
%   - Pure MATLAB, no VTK / convex hull / auxiliary queries: one pass over the
%     edges builds the face-adjacency PARITY, a single CONNCOMP on a doubled
%     graph splits every orientable component into its two orientation
%     classes, and the signed volume picks the outward one.
%   - Derived fields (triNORMALS, xyzNORMALS, ...) are NOT updated: by
%     convention, recomputing them after the flip is the CALLER's job.
%   - Non-orientable components (Moebius-like) cannot be made consistent: a
%     warning is issued and there the result is arbitrary.
%   - IDEMPOTENT: re-running on an already consistent + outward mesh flips
%     nothing. This is the pure-MATLAB orienter meshIsInterior relies on
%     (addNormalsToM_tri) to guarantee OUTWARD normals for its sign tests.
%   - PERFORMANCE: this CONNCOMP formulation is ~7x faster than the old
%     seed-and-grow front (which seeded from a convhulln / VTK and propagated
%     face by face). It scales ~linearly: ~0.03 s @ 60k faces, ~0.1 s @ 200k,
%     ~0.33 s @ 600k. About half of that (~0.18 s @ 600k) is the MATLAB graph()
%     object build + conncomp overhead. FOR GIGANTIC MESHES (millions of faces)
%     where that half matters, replace the graph()+conncomp block (the parity
%     doubled graph) with a VECTORIZED UNION-FIND on the same edge-parity pairs
%     (fa,fb,opp): union node a with b (or b+nF) per pair, then read the roots --
%     this removes the graph-object overhead and keeps the exact same result.
%     Below a few hundred-k faces it is NOT worth it (the current cost is tiny).
%
% See also meshNormals, meshIsInterior, meshCelltype, MeshTidy.

  if nargin < 2, f = []; end

  switch meshCelltype( M )
    case 3
      error('not implemented yet.');

    case 5
      w = fixTri( double(M.xyz) , double(M.tri) , double(f(:)) );
      M.tri( w ,[2,3]) = M.tri( w ,[3,2]);

    case 10
      if ~isempty( f ), error('not implemented for fixed faces for this celltype'); end

      P1 = M.xyz( M.tri(:,1) ,:);
      P2 = M.xyz( M.tri(:,2) ,:);
      P3 = M.xyz( M.tri(:,3) ,:);
      P4 = M.xyz( M.tri(:,4) ,:);

      L1 = P2 - P1;
      L3 = P1 - P3;
      L4 = P4 - P1;

      A1 = cross( L3 , L1 , 2 );
      vs = dot( A1 , L4 ,2);

      w = vs < 0;

      if any(w)
        M.tri(w,[3,4]) = M.tri(w,[4,3]);
      end

    otherwise
      error('not implemented for this celltype.');
  end

end


function w = fixTri( xyz , tri , f )
% one pass over the edges -> face-adjacency parity -> one CONNCOMP -> outward.
  nF = size( tri ,1);
  nV = size( xyz ,1);
  if nF == 0, w = false(0,1); return; end

  % the 3 directed edges of every face, keyed UNDIRECTED (exact in double:
  % key < (nV+1)^2 << 2^53), keeping the traversal direction as a flag.
  E   = [ tri(:,[1,2]) ; tri(:,[2,3]) ; tri(:,[3,1]) ];
  FI  = repmat( (1:nF).' , 3 ,1);
  lh  = E(:,1) < E(:,2);                       %directed edge runs low->high
  key = min(E,[],2) * (nV+1) + max(E,[],2);
  [ key , ord ] = sort( key );
  FI = FI(ord);  lh = lh(ord);

  pr  = find( key(1:end-1) == key(2:end) );    %consecutive rows = same edge ->
  fa  = FI(pr);  fb = FI(pr+1);                %adjacent faces (a non-manifold
                                               %edge with k faces gets CHAINED)
  opp = lh(pr) ~= lh(pr+1);                    %opposite traversal = CONSISTENT

  % DOUBLED graph: node i = face i as-is, node i+nF = face i flipped.
  % Consistent neighbours link as-is<->as-is (and flipped<->flipped);
  % inconsistent ones link as-is<->flipped. Each orientable component then
  % splits into exactly TWO classes = its two global orientations.
  s   = [ fa             ; fa + nF        ];
  t   = [ fb + nF*(~opp) ; fb + nF*opp    ];
  cid = conncomp( graph( s , t , [] , 2*nF ) );

  % face components: the (unordered) pair of classes identifies the component
  [ ~ , ~ , comp ] = unique( min( cid(1:nF) , cid(nF+1:2*nF) ).' );
  nC = max( comp );

  % seed per component: the FIRST reference face (its sign = required parity)
  % or, without references, the lowest face id (parity settled by volume below)
  seed = accumarray( comp , (1:nF).' , [nC,1] , @min );
  spar = zeros( nC ,1);
  hasf = false( nC ,1);
  for k = 1:numel(f)
    c = comp( abs( f(k) ) );
    if ~hasf(c), hasf(c) = true; seed(c) = abs( f(k) ); spar(c) = f(k) < 0; end
  end

  K  = cid( seed + nF*spar ).';                %chosen class per component
  w0 = cid( 1:nF      ).' == K( comp );        %face as-is   lands in the class
  w1 = cid( nF+1:2*nF ).' == K( comp );        %face flipped lands in the class
  if any( w0 & w1 )
    warning( 'MeshFixCellOrientation:nonOrientable' , ...
             'non-orientable component(s) (Moebius-like): the orientation there is arbitrary.' );
  end
  w = w1 & ~w0;

  % outward: flip every NON-referenced component with negative signed volume.
  % centered first -- the volume of a CLOSED surface is translation-invariant,
  % and centering avoids the catastrophic cancellation of summing huge
  % origin-based tet volumes when the mesh lives far from the origin.
  xyz = xyz - mean( xyz ,1);
  dv  = dot( xyz( tri(:,1) ,:) , cross( xyz( tri(:,2) ,:) , xyz( tri(:,3) ,:) , 2 ) , 2 ) / 6;
  dv( w ) = -dv( w );
  neg = ( accumarray( comp , dv , [nC,1] ) < 0 ) & ~hasf;
  w( neg( comp ) ) = ~w( neg( comp ) );
end
