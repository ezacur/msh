function [X,C] = mesh2contours( B )
%MESH2CONTOURS  Trace the connected polylines (contours) of an edge mesh.
%
%   [X,C] = mesh2contours( B )
%
%   Walks a polyline mesh (celltype 3 -- a list of line segments) and returns
%   its connected chains, following each one edge by edge. Both a full mesh
%   struct and a bare edge list are accepted.
%
%   INPUT
%     B   either a mesh struct with fields .tri (#E-by-2 edge list) and .xyz
%         (#V-by-D vertex coordinates), or just a numeric #E-by-2 edge list.
%
%   OUTPUT
%     C   cell column, one entry per connected contour, each a row vector of
%         vertex indices in traversal order, SORTED by length (longest first).
%         A CLOSED loop repeats its first vertex at the end ([1 2 3 1]); an OPEN
%         chain does not ([1 2 3]).
%     X   * if B is a struct: the contour vertices B.xyz(C{i},:) stacked into one
%           array, consecutive contours separated by a row of NaN (so plot3(X(:,1),
%           X(:,2),X(:,3)) draws every contour with gaps between them).
%         * if B is a bare edge list: X is just C (no coordinates to gather).
%
%   NOTES
%     * Cost is linear in the number of edges. Traversal is via a node->incident
%       edge index, and each edge is consumed exactly once.
%     * Intended for simple polylines. At a junction (a vertex of valence >= 3)
%       one branch is followed and the rest are emitted as separate chains, so
%       the decomposition of a branching graph is not unique.
%     * An empty edge list yields C = {} and an empty X.
%
%   See also meshCelltype, MeshBoundary, meshEdges.

  asMESH = true;
  if ~isstruct( B )
    asMESH = false;
    B = struct( 'tri' , B , 'xyz' , zeros( max(B(:)) , 3 , 0 ) );
  end

  B.celltype = meshCelltype( B );
  if B.celltype ~= 3, error('only valid for polyline mesh type (celltype = 3).'); end

  E  = double( B.tri );
  nE = size( E ,1 );

  C = cell(0,1);

  if nE > 0
    nN = max( E(:) );

    % node -> incident edges, ordered exactly as the old find(T==node,1) scanned
    % the 2-by-nE layout T = [a b].' : endpoint 1 of edge k sits at linear index
    % 2k-1, endpoint 2 at 2k. Sorting (node,linIdx) reproduces that first-
    % occurrence order, so the traced sequences match a plain linear scan while
    % costing O(#edges) instead of O(#edges^2).
    node   = [ E(:,1)         ; E(:,2)       ];
    other  = [ E(:,2)         ; E(:,1)       ];
    edge   = [ (1:nE).'       ; (1:nE).'     ];
    linIdx = [ (2*(1:nE).'-1) ; (2*(1:nE).') ];
    [~,o]  = sortrows( [ node , linIdx ] );
    node = node(o); other = other(o); edge = edge(o);

    cnt   = accumarray( node , 1 , [nN,1] );   % degree of every node id 1..nN
    start = cumsum( [1;cnt] );                 % node v owns start(v):start(v+1)-1
    pos   = start(1:nN);                       % moving read pointer per node

    consumed  = false( nE ,1 );
    remaining = nE;
    R = zeros( 1 , nE+1 );                      % chain buffer (max length nE+1)

    while remaining > 0
      %pick a start: an open end (valence-1 node) if any, else any remaining node
      MB = MeshBoundary( E(~consumed,:) );
      if isempty( MB )
        l = E( find( ~consumed , 1 ) , 1 );
      else
        l = MB(1);
      end

      N = 1; R(1) = l;
      while true
        p = pos(l); lim = start(l+1)-1;
        while p <= lim && consumed( edge(p) ), p = p+1; end   %skip used edges
        pos(l) = p;
        if p > lim, break; end                                %no edge left at l
        k  = edge(p);
        nl = other(p);
        consumed(k) = true;  remaining = remaining - 1;
        N = N+1; R(N) = nl;
        l = nl;
      end

      C{end+1,1} = R(1:N);
    end
  end

  n = cellfun( 'prodofsize' , C );
  [~,ord] = sort( n , 'descend' );
  C = C(ord);

  if asMESH
    ncol = size( B.xyz ,2 );
    if isempty( C )
      X = zeros( 0 , ncol );
    else
      total = sum(n) + numel(C) - 1;            %contour blocks + one NaN row between
      X = NaN( total , ncol );
      w = 0;
      for c = 1:numel(C)
        if c > 1, w = w + 1; end                %leave the preallocated NaN row
        idx = C{c};
        X( w+(1:numel(idx)) ,:) = B.xyz( idx ,:);
        w = w + numel(idx);
      end
    end
  else
    X = C;
  end

end
