function M = MeshReOrderFaces( M , X )
%MESHREORDERFACES  Renumber the FACES (cells) of a mesh (geometry is unchanged).
%
%  M = MeshReOrderFaces( M )          FLIP the face order (face k -> Nf+1-k)
%  M = MeshReOrderFaces( M , perm )   apply the PERMUTATION: new face i is the
%                                     old face perm(i). perm MUST be a
%                                     permutation of 1:Nf -- a partial vector
%                                     is an ERROR (use the cells).
%  M = MeshReOrderFaces( M , {A} )    faces A go FIRST (in A's order); the rest
%                                     keep their relative order.
%  M = MeshReOrderFaces( M , {A,B} )  faces A go FIRST, faces B go LAST, the
%                                     rest stay in between in their order
%                                     ( {[],B} just pushes B to the end ).
%
%  Every tri* field (and a per-face celltype vector) is reordered along; nodes
%  are untouched. The mesh is IDENTICAL before and after -- only the face
%  numbering changes. Sibling of MeshReOrderNodes (same calling idioms).
%
% See also MeshReOrderNodes, MeshRemoveFaces.

  Nf = size( M.tri , 1);

  if nargin < 2
    order = Nf:-1:1;
  else
    order = completePermutation( X , Nf , 'face' );
  end

  for f = fieldnames( M ).', f = f{1};
    if ~strncmp( f , 'tri',3), continue; end
    M.(f) = M.(f)( order ,:,:,:,:,:,:);          %extra colons: keep >2-dim fields' shape
  end
  if isfield( M , 'celltype' ) && ~isscalar( M.celltype )
    M.celltype = M.celltype( order ,:);
  end

end


function order = completePermutation( X , N , what )
%a plain numeric vector MUST be a full permutation of 1:N (else error). The
%partial forms are OPT-IN via a cell: {first} or {first,last} lists of distinct
%ids -> completed to a full permutation, the unlisted ids keeping their relative
%order in between. (same helper as in MeshReOrderNodes)
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
