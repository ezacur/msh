function M = MeshReOrderNodes( M , X )
%MESHREORDERNODES  Renumber the nodes of a mesh (geometry is unchanged).
%
%  M = MeshReOrderNodes( M )         FLIP the node order (node k -> Nn+1-k)
%  M = MeshReOrderNodes( M , perm )  apply the PERMUTATION: new node i is the
%                                    old node perm(i)  (M.xyz -> M.xyz(perm,:)).
%                                    perm MUST be a permutation of 1:Nn -- a
%                                    partial vector is an ERROR (use the cells).
%  M = MeshReOrderNodes( M , {A} )   nodes A go FIRST (in A's order); the rest
%                                    keep their relative order. So
%                                        MeshReOrderNodes( M , {1:10:100} )
%                                    replaces the old idiom
%                                        MeshReOrderNodes( M , unique([1:10:100,1:Nn],'stable') )
%  M = MeshReOrderNodes( M , {A,B} ) nodes A go FIRST, nodes B go LAST, the
%                                    rest stay in between in their order
%                                    ( {[],B} just pushes B to the end ).
%  M = MeshReOrderNodes( M , X )     put FIRST the nodes matching the ROWS of
%                                    the coordinate matrix X (exact match), in
%                                    X's order; the rest keep their order.
%                                    X may also be a mesh struct (uses X.xyz).
%
%  Every xyz* field is reordered along with the coordinates (shapes preserved),
%  and M.tri is remapped accordingly; zero entries in tri (padding of mixed-cell
%  meshes) are preserved. The connectivity/geometry is IDENTICAL before and
%  after -- only the numbering changes.
%

% 
% using isomorphism from graph toolbox
% perm = isomorphism( simplify( graph( table( meshEdges(L1) ,'VariableNames',{'EndNodes'} ) ) ) , simplify( graph( table( meshEdges(L0) ,'VariableNames',{'EndNodes'} ) ) ) )'
% 
%
% using igraph toolbox  (setup by:  enablePYTHON( 'igraph' ) )
% perm = double( py.mesh_iso.find_permutation( int32( vec( meshEdges( L0 ) ).' - 1 ) , int32( vec( meshEdges( L1 ) ).' - 1 ) , int32( max( max( L0.tri(:) ) , max( L1.tri(:) ) ) ) ) ) + 1
% 
% >>>>>>>>>> mesh_iso.py >>>>>>>>>>
% import igraph as ig
% def find_permutation(e0_1d, e1_1d, num_nodes):
%     N = int(num_nodes)
%     e0_list = list(e0_1d)
%     e1_list = list(e1_1d)
%     n0 = len(e0_list) // 2
%     edges0 = [(e0_list[i], e0_list[i + n0]) for i in range(n0)]
%     n1 = len(e1_list) // 2
%     edges1 = [(e1_list[i], e1_list[i + n1]) for i in range(n1)]
%     g0 = ig.Graph(int(num_nodes), edges=edges0, directed=False)
%     g0.simplify()
%     p0 = g0.canonical_permutation()
%     g1 = ig.Graph(int(num_nodes), edges=edges1, directed=False)
%     g1.simplify()
%     p1 = g1.canonical_permutation()
%     M_B = [0]*N
%     for c in range(N):
%         M_B[p1[c]] = p0[c]
%     return M_B
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% 


  Nn = size( M.xyz , 1);

  if nargin < 2

    for f = fieldnames( M ).', f = f{1};
      if ~strncmp( f , 'xyz',3), continue; end
      M.(f) = flip( M.(f) , 1 );
    end

    w = M.tri ~= 0;                              %0 = padding (mixed cells): keep it
    M.tri(w) = ( Nn + 1 ) - M.tri(w);

  elseif isstruct( X )
    
    M = MeshReOrderNodes( M , X.xyz );

  elseif iscell( X ) || ( isnumeric( X ) && ( isvector( X ) || isempty( X ) ) )

    order = completePermutation( X , Nn , 'node' );
    for f = fieldnames( M ).', f = f{1};
      if ~strncmp( f , 'xyz',3), continue; end
      M.(f) = M.(f)( order ,:,:,:,:,:,:);        %extra colons: keep >2-dim fields' shape
    end


    T = M.tri(:); w = ~~T;
    T(w) = iperm( order , T(w) );
    M.tri = reshape( T , size( M.tri ) );

  elseif ~isvector( X )


    [~,b] = ismember( M.xyz , X , 'rows' );
    b( ~b ) = Inf;
    [~,order] = sort( b );

    for f = fieldnames( M ).', f = f{1};
      if ~strncmp( f , 'xyz',3), continue; end
      M.(f) = M.(f)( order ,:,:,:,:,:,:);        %extra colons: keep >2-dim fields' shape
    end

    T = M.tri(:); w = ~~T;                       %same 0-padding-safe remap as above
    T(w) = iperm( order , T(w) );                %(also keeps a 1-cell tri as a ROW)
    M.tri = reshape( T , size( M.tri ) );

  end


end


function order = completePermutation( X , N , what )
%a plain numeric vector MUST be a full permutation of 1:N (else error). The
%partial forms are OPT-IN via a cell: {first} or {first,last} lists of distinct
%ids -> completed to a full permutation, the unlisted ids keeping their relative
%order in between. (same helper as in MeshReOrderFaces)
  if ~iscell( X )
    if ~isnumeric( X ) || ~( isvector( X ) || isempty( X ) )
      error('a permutation of 1:%d, or a {first} / {first,last} cell of %s ids, was expected.',N,what);
    end
    order = X(:).';
    if numel( order ) ~= N || ~isequal( sort( order ) , 1:N )
      error('a PERMUTATION of 1:%d was expected (use {ids} to put some %ss first).',N,what);
    end
    return;
  end

  if numel( X ) > 2, error('at most a {first,last} pair of %s lists is allowed.',what); end
  A = []; B = [];
  if numel( X ) >= 1, A = X{1}(:).'; end
  if numel( X ) == 2, B = X{2}(:).'; end
  AB = [ A , B ];
  if ~isnumeric( AB ) || ~isreal( AB ) || any( ~isfinite( AB ) ) || any( AB ~= round( AB ) ) || ...
     any( AB < 1 ) || any( AB > N ) || ...
     numel( unique( A ) ) ~= numel( A ) || numel( unique( B ) ) ~= numel( B )
    error('%s ids must be DISTINCT integers between 1 and %d.',what,N);
  end
  if any( ismember( A , B ) )
    error('the {first} and {last} lists OVERLAP: a %s cannot be both first and last.',what);
  end
  rest = 1:N;  rest( ismember( rest , AB ) ) = [];
  order = [ A , rest , B ];
end