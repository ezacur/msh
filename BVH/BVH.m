function B = BVH( M , arg2 , varargin )
%BVH  BVH blob over the mesh elements: build / transform / refit.
%
%   B = BVH( M )              builds the blob over M.tri: binned-SAH BVH4,
%                                AABB node volumes by default. Celltypes:
%                                vertices (1 node/row), segments (2),
%                                triangles (3), TETRAHEDRA (4); 0-padded MIXED
%                                rows (zeros trailing). 4 nonzero nodes ALWAYS
%                                mean a tetrahedron.
%   B = BVH( M , leaf )       SAH leaf bounds: a PAIR [minLeaf maxLeaf]
%                                (adaptive: n <= minLeaf is always a leaf, the
%                                SAH may keep up to maxLeaf together when
%                                splitting does not pay; default [8 32]) or a
%                                SCALAR s == [s s] (fixed leaf size; Inf = one
%                                single leaf = brute force, for references).
%   B = BVH( M , leaf , VOLUME , ... )
%                                node volume type, all under the SAME SAH
%                                partition:
%                                  'aabb'   (default) boxes in the build frame
%                                  'sphere' bounding spheres
%                                  'obb'    per-slot oriented boxes (each slot
%                                           PCA-aligned to ITS contained nodes)
%                                  'kdop'   14-DOPs: AABB + the 4 diagonal
%                                           direction pairs (1,±1,±1)
%                                  'rss'    rectangle-swept spheres (PQP-style:
%                                           per-slot PCA rectangle + radius)
%                                  'lss'    line-swept spheres / capsules
%                                           (alias 'capsule'; per-slot PCA
%                                           segment + radius -- tubular
%                                           geometry's natural volume)
%                                extra option 'noframe': skip the centroid/PCA
%                                build frame (world-aligned blob; comparison
%                                piece -- costs precision far from the origin
%                                and tightness on diagonal geometry).
%
%   B = BVH( B , T )          TRANSFORM: similarities (rotation+translation+
%                                uniform scale, reflections included) FOLD into
%                                the global frame in O(1) -- no node is touched;
%                                queries transform instead. A NON-similarity T
%                                ERRORS (BVH:notSimilarity): refit against
%                                the transformed mesh instead --
%                                    try,   Bt = BVH( B , T  );
%                                    catch, Bt = BVH( B , Mt );
%                                    end
%                                T: 4x4 / 3x4 / 3x3 (3-D linear); for 2-D blobs
%                                also 3x3 homogeneous 2-D / 2x3 / 2x2 (lifted
%                                with z-scale = the 2-D scale).
%   B = BVH( B , M2 )         REFIT to a DEFORMED mesh (same connectivity):
%                                keeps the hierarchy AND the frame, recomputes
%                                the stored geometry and every node volume from
%                                its own element range (exact; compiled). O(n).
%
%   THE GLOBAL FRAME. The blob stores its geometry in a BUILD frame:
%   world = R*build + mu. At build time the frame is set to the centroid
%   (float node bounds keep full RELATIVE precision even for meshes far from
%   the origin) and, when the vertex cloud is anisotropic (principal-axis
%   ratio >= 2), to the PCA axes -- elongated DIAGONAL geometry gets tight
%   axis-aligned boxes in its OWN frame. Every later similarity composes into
%   the frame in O(1). Use plotBVH(B,M) to SEE all of this.
%
%   B fields (blob version 3, SELF-CONTAINED -- a plain value: cacheable,
%   save/load and parfor safe):
%     .bounds4          per-slot node volumes, 4-wide (single, conservatively
%                       rounded OUTWARD: the float bounds always contain the
%                       double geometry)
%     .node4            pool de nodos FUNDIDO (uint8): por nodo, sus bounds4 Y
%                       sus child4 contiguos, con stride multiplo de 64. Es lo
%                       que recorren las queries: una visita = un registro, en
%                       vez de dos lecturas en arrays distintos. Se construye
%                       aqui (antes lo rehacia cada llamada al mex y lo tiraba).
%                       Derivado de bounds4/child4: si los tocas a mano, refunde.
%     .child4,.srange   4-ary tree: children (>0 node, -1 leaf, 0 empty) and
%                       per-slot element ranges into .perm (int32)
%     .perm,.pkV,.pkS,.pkT,.pkE  element permutation + packed leaf data
%                       (vertices, spheres, node counts, original ids)
%     .X,.Tri           the geometry itself, in BUILD-frame space
%     .frame            4x4 global similarity: world = frame(buildSpace)
%     .eCenter,.eRadius per-element bounding spheres (build space)
%     .nsd              2|3 from size(M.xyz,2) AT BUILD -- NEVER inferred from
%                       the values (an all-z==0 3-D mesh stays 3-D and may
%                       deform out of plane)
%     .vol,.volume      1|2 and 'sphere'|'aabb'
%
%   Requires the compiled BVH_mx (mex BVH_mx.cpp).
%
% See also bvhClosestElement, bvhIntersectRay, plotBVH.

  %update modes: BVH( B , T ) affine / BVH( B , M ) refit-to-deformed
  if isstruct( M ) && isfield( M , 'isBVH' )
    if isstruct( arg2 )
      B = bvhRefit( M , arg2 );
    else
      B = bvhTransform( M , arg2 );
    end
    return;
  end

  if exist( 'BVH_mx' ,'file') ~= 3
    error('BVH:mex','BVH_mx is not compiled (run: mex BVH_mx.cpp).');
  end

  %options: volume type + 'noframe'
  volume  = 'aabb';
  noframe = false;
  for a = 1:numel( varargin )
    o = varargin{a};  if isempty(o), continue; end
    if ~( ischar(o) || isstring(o) ), error('BVH:opt','options must be strings.'); end
    switch lower( char(o) )
      case {'sphere','aabb','obb','kdop','rss','lss'}, volume = lower( char(o) );
      case 'capsule',                                  volume = 'lss';
      case 'noframe',                                  noframe = true;
      otherwise, error('BVH:opt','unknown option ''%s'' (use aabb|sphere|obb|kdop|rss|lss|noframe).', char(o));
    end
  end

  %SAH leaf bounds: pair [minLeaf maxLeaf]; scalar s == [s s]; Inf = one leaf
  sahLeaf = [ 8 , 32 ];   %post-PreTri4 las hojas son baratas: [8 32] gana a
                          %[2 16] un 14-26% en closest-point (near/mid/far) y
                          %queda a ~3% del optimo de rayos [16 64]
  if nargin >= 2 && ~isempty( arg2 )
    if numel( arg2 ) == 2, sahLeaf = double( arg2(:) ).';
    else,                  sahLeaf = double( arg2 ) * [ 1 , 1 ];
    end
  end
  sahLeaf = max( 1 , min( sahLeaf , 1e9 ) );     %clamp (Inf -> huge -> one leaf)

  nsd = size( M.xyz ,2);          %2 or 3, from the COLUMN COUNT only (never
                                  %from the values: all-z==0 stays 3-D)
  X = double( M.xyz ); X(:,end+1:3) = 0;
  T = double( M.tri );
  nE = size( T ,1);

  %build frame: ALWAYS center at the centroid (float bounds keep relative
  %precision far from the origin); PCA-align when anisotropy >= 2 (tight
  %boxes for elongated diagonal geometry). world = R*build + mu.
  mu = mean( X ,1);
  if any( ~isfinite( mu ) ) || noframe, mu = [0 0 0]; end
  R  = eye(3);
  Xc = X - mu;
  if size( X ,1) > 2 && ~noframe
    if nsd == 2
      C = Xc(:,1:2).' * Xc(:,1:2);
      [ V2 , D2 ] = eig( ( C + C.' )/2 );
      [ ev , ord ] = sort( sqrt( max( diag(D2) ,0) ) , 'descend' );  V2 = V2(:,ord);
      if ev(1) > 0 && ev(1) >= 2*ev(2)
        if det( V2 ) < 0, V2(:,2) = -V2(:,2); end
        R = blkdiag( V2 , 1 );
      end
    else
      C = Xc.' * Xc;
      [ V3 , D3 ] = eig( ( C + C.' )/2 );
      [ ev , ord ] = sort( sqrt( max( diag(D3) ,0) ) , 'descend' );  V3 = V3(:,ord);
      if ev(1) > 0 && ev(1) >= 2*ev(3)
        if det( V3 ) < 0, V3(:,3) = -V3(:,3); end
        R = V3;
      end
    end
  end
  X = Xc * R;                                    %BUILD-frame coordinates

  [ eC , eR ]   = elementSpheres( X , T );
  [ eLo , eHi ] = elementBoxes( X , T );

  volFlag = find( strcmp( volume , { 'sphere' , 'aabb' , 'obb' , 'kdop' , 'rss' , 'lss' } ) );
  [ eVv , eNv ] = elementVerts( X , T );
  [ b4 , c4 , r4 , pm , pk4 , pk4id , s4 ] = ...
      BVH_mx( eC , eR , eLo , eHi , sahLeaf(1) , sahLeaf(2) , volFlag , eVv , eNv );
  pm = double( pm );
  [ pkV , pkS , pkT ] = packLeafData( X , T , eC , eR , pm );

  B = struct( 'isBVH'   , true         ,...
              'version' , 3            ,...
              'vol'     , volFlag      ,...
              'volume'  , volume       ,...
              'bounds4' , b4           ,...
              'child4'  , c4           ,...
              'node4'   , fuseNodes( b4 , c4 ) ,...
              'srange'  , r4           ,...
              'perm'    , pm           ,...
              'pkV'     , pkV          ,...
              'pkS'     , pkS          ,...
              'pkT'     , pkT          ,...
              'pkE'     , int32( pm )  ,...
              'pk4'     , pk4          ,...
              'pk4id'   , pk4id        ,...
              's4'      , s4           ,...
              'eCenter' , eC           ,...
              'eRadius' , eR           ,...
              'nE'      , nE           ,...
              'leafSize', sahLeaf      ,...
              'frame'   , [ R , mu.' ; 0 0 0 1 ] ,...
              'nsd'     , nsd          ,...
              'X'       , X            ,...
              'Tri'     , int32( T )   );   %conectividad SIEMPRE int32 en el blob
  checkBlobIndices( B );    %una vez aqui, para que las consultas no lo repitan
end

%% ------------------------------------------------------------------ helpers
function [eC,eR] = elementSpheres( X , T )
%per-element bounding spheres: center = mean of its (nonzero) nodes, radius =
%max distance from center to them. Bounds the element (its convex hull).
  nE  = size( T ,1);
  cnt = sum( T > 0 ,2);
  eC  = zeros( nE ,3);
  for c = 1:size( T ,2)
    w = T(:,c) > 0;
    eC(w,:) = eC(w,:) + X( T(w,c) ,:);
  end
  eC = eC ./ max( cnt ,1);
  eR = zeros( nE ,1);
  for c = 1:size( T ,2)
    w = T(:,c) > 0;
    eR(w) = max( eR(w) , sqrt( sum( ( X( T(w,c) ,:) - eC(w,:) ).^2 ,2) ) );
  end
end

function [eLo,eHi] = elementBoxes( X , T )
%per-element AABBs (0-padded rows: zeros trailing; empty rows get a point box)
  nE  = size( T ,1);
  eLo = inf( nE ,3);  eHi = -inf( nE ,3);
  for c = 1:size( T ,2)
    w = T(:,c) > 0;
    V = X( T(w,c) ,:);
    eLo(w,:) = min( eLo(w,:) , V );
    eHi(w,:) = max( eHi(w,:) , V );
  end
  w = ~( T(:,1) > 0 );                 %empty rows: harmless degenerate box
  eLo(w,:) = 0;  eHi(w,:) = 0;
end

function [eV,eN] = elementVerts( X , T )
%element vertices in ORIGINAL order (12 doubles/element, zero-padded) + counts
%-- the geometry the builder needs for the vertex-based volumes (obb, kdop)
  nE = size( T ,1);
  eN = sum( T > 0 ,2);
  eV = zeros( 12 , nE );
  for c = 1:size( T ,2)
    w = T(:,c) > 0;
    eV( 3*c-2:3*c , w ) = X( T(w,c) ,:).';
  end
end

function [pkV,pkS,pkT] = packLeafData( X , T , eC , eR , pm )
%leaf element data packed at BUILD time, in perm order: raw vertices (12
%doubles/element, zero-padded), element spheres (4 doubles), node counts.
  nE  = size( T ,1);
  Tp  = T( pm ,:);
  pkT = int32( sum( Tp > 0 ,2) );
  pkV = zeros( 12 , nE );
  for c = 1:size( Tp ,2)
    w = Tp(:,c) > 0;
    pkV( 3*c-2:3*c , w ) = X( Tp(w,c) ,:).';
  end
  pkS = [ eC( pm ,:) , eR( pm ) ].';
end

%% ---------------------------------------------------------------- transform
function B = bvhTransform( B , T )
%SIMILARITY -> fold into the global frame, O(1), nothing else moves.
%NON-similarity (anisotropic/shear) -> ERROR: the caller decides, typically
%refitting against the transformed mesh (which keeps the hierarchy):
%    Mt = transform( M , T );
%    try,   Bt = BVH( B , T  );      %similarity: O(1) fold
%    catch, Bt = BVH( B , Mt );      %general affine: O(n) refit
%    end
  [ A , t ] = parseTransform( T , B );

  F  = [ A , t.' ; 0 0 0 1 ] * B.frame;
  Af = F(1:3,1:3);
  s2 = trace( Af.' * Af ) / 3;
  if s2 > 0 && norm( Af.'*Af - s2*eye(3) ,'fro') <= 1e-9 * 3 * s2
    B.frame = F;                                   %similarity: fold, O(1)
  else
    error( 'BVH:notSimilarity' , ...
           [ 'T is not a similarity (A''*A ~= s^2*I): it cannot fold into the ' ...
             'frame for distance queries. Refit against the transformed mesh ' ...
             'instead:  B = BVH( B , transformedMesh ).' ] );
  end
end

function [A,t] = parseTransform( T , B )
%accepted forms: 4x4 / 3x4 / 3x3 (3-D linear). For 2-D blobs (B.nsd == 2) also
%3x3 homogeneous 2-D (last row [0 0 1]) / 2x3 / 2x2, LIFTED with z-scale = the
%2-D scale so the padded z stays consistent under scaling.
  nsd = 3;  if isstruct( B ) && isfield( B ,'nsd'), nsd = B.nsd; end
  if nsd == 2 && ( isequal( size(T) , [2,2] ) || isequal( size(T) , [2,3] ) || ...
                 ( isequal( size(T) , [3,3] ) && isequal( T(3,:) , [0 0 1] ) ) )
    A2 = T(1:2,1:2);
    if size(T,2) > 2, t2 = T(1:2,end).'; else, t2 = [0 0]; end
    s  = sqrt( max( trace( A2.'*A2 )/2 , 0 ) );
    A  = blkdiag( A2 , s );
    t  = [ t2 , 0 ];
  elseif isequal( size(T) , [4,4] ), A = T(1:3,1:3); t = T(1:3,4).';
  elseif isequal( size(T) , [3,4] ), A = T(: ,1:3);  t = T(: ,4).';
  elseif isequal( size(T) , [3,3] ), A = T;          t = [0 0 0];
  else,  error('BVH:transform','T must be 4x4, 3x4 or 3x3 (2-D blobs also accept 3x3 homogeneous 2-D, 2x3, 2x2).');
  end
end

%% -------------------------------------------------------------------- refit
function B = bvhRefit( B , M )
%persistent-hierarchy refit: keep perm/child4/srange AND the frame; recompute
%ONLY the geometry from the CURRENT coordinates -- every slot from its own
%element range (exact for both volume types, no child-merge compounding).
  X = double( M.xyz ); X(:,end+1:3) = 0;
  T = double( M.tri );
  if size( T ,1) ~= B.nE
    error('BVH:refit','refit requires the SAME connectivity (%d vs %d elements). Rebuild instead.', size(T,1) , B.nE );
  end
  %world -> current build frame (the frame is KEPT: it is material, not spatial)
  if ~isequal( B.frame , eye(4) )
    Af = B.frame(1:3,1:3);
    s2 = trace( Af.'*Af )/3;
    X  = ( X - B.frame(1:3,4).' ) * ( Af / s2 );   %(X-t)*inv(A).' , inv(A)=A.'/s2
  end
  B.X   = X;
  B.Tri = int32( T );

  [ eC , eR ] = elementSpheres( X , T );
  B.eCenter = eC;
  B.eRadius = eR;
  [ eLo , eHi ] = elementBoxes( X , T );
  [ eVv , eNv ] = elementVerts( X , T );
  [ B.bounds4 , B.pk4 ] = BVH_mx( eC , eR , eLo , eHi , B.srange , B.child4 , B.pkE , B.vol , eVv , eNv );
  B.node4 = fuseNodes( B.bounds4 , B.child4 );   %el refit cambia bounds: refundir
  [ B.pkV , B.pkS ] = packLeafData( X , T , eC , eR , double( B.perm ) );
  checkBlobIndices( B );
end

%% -------------------------------------------- validacion de indices del blob
function checkBlobIndices( B )
%Comprueba ENTEROS los indices del blob: hijos, rangos de slot y lane ids del
%pool PreTri4. Se hace UNA VEZ, aqui (al construir y al refitear), vectorizado.
%
%   Antes lo hacia cada mex de consulta en CADA llamada, recorriendo 4*nB + 4*nN
%   entradas: 277 000 comprobaciones de rango y ~78 us en una malla de 87k caras,
%   mas que la consulta misma en lotes pequenos -- y revalidando un blob que es
%   inmutable. Los mexes ahora solo muestrean (O(1)); la garantia fuerte esta
%   aqui. Si manipulas bounds4/child4/srange/s4/pk4id a mano, vuelve a llamarla.
  nN = size( B.child4 ,2);
  nE = size( B.pkS ,2);
  nB = size( B.pk4 ,2);
  c  = int32( B.child4(:) );
  if any( c < -1 | c > int32( nN ) )
    error('BVH:corruptBlob','indice de hijo fuera de [-1,%d].', nN );
  end
  if any( int32( B.pk4id(:) ) < 0 | int32( B.pk4id(:) ) > int32( nE ) )
    error('BVH:corruptBlob','lane id de PreTri4 fuera de [0,%d].', nE );
  end
  %rangos de slot: solo los slots OCUPADOS (child4 ~= 0) tienen que ser validos
  occ = reshape( B.child4 ,4,[]) ~= 0;                 %4 x nN
  r   = double( B.srange );  lo = r(1:2:8,:);  hi = r(2:2:8,:);
  if any( lo(occ) < 1 ) || any( hi(occ) < lo(occ) ) || any( hi(occ) > nE )
    error('BVH:corruptBlob','rango de slot fuera de [1,%d] o invertido.', nE );
  end
  s  = double( B.s4 );  s0 = s(1:2:8,:);  sn = s(2:2:8,:);
  if any( sn(occ) < 0 )
    error('BVH:corruptBlob','contador de bloques PreTri4 negativo.');
  end
  w = occ & sn > 0;
  if any( s0(w) < 1 ) || any( s0(w) + sn(w) - 1 > nB )
    error('BVH:corruptBlob','rango de bloques PreTri4 fuera de [1,%d].', nB );
  end
end

%% ------------------------------------------------- registro de nodo fundido
function nd = fuseNodes( b4 , c4 )
%NODE4: por cada nodo, sus 4 volumenes de slot Y sus 4 hijos CONTIGUOS.
%
%   Las queries necesitan las dos cosas en cada visita, pero viven en dos
%   mxArrays distintos (dos direcciones sin relacion = dos fallos de cache
%   independientes por nodo, y el prefetcher no puede seguir dos streams).
%   Fundirlos deja una visita = un registro contiguo.
%
%   Antes esto lo hacia CADA llamada al mex (un vector temporal + memcpy sobre
%   los nN nodos, ~10 ns/nodo: 110 us en 87k caras, 580 us en 200k) y se tiraba
%   al volver. Es dato derivado del blob con la MISMA vida que el blob, asi que
%   se construye aqui una vez -- el mismo argumento que justifica cachear el BVH.
%
%   MEDIDO Y DESCARTADO: rellenar el stride a multiplo de 64 (128 en vez de 112
%   para aabb) para que cada nodo ocupe 2 lineas de cache exactas en vez de
%   cruzarlas (~2.75 de media). Da entre -5.8% y +4.5% segun el regimen, o sea
%   RUIDO, y cuesta un 12% mas de memoria: se deja el stride justo.
%
%   srange/s4 NO se funden aqui a proposito: solo se leen en las HOJAS, y
%   meterlos haria que cada nodo INTERNO (que son la mayoria) arrastrase 176 B
%   en vez de 112 sin usarlos.
  S      = size( b4 ,1);                  %floats de bounds por nodo (por volumen)
  nN     = size( b4 ,2);
  stride = 4*S + 16;
  nd = zeros( stride , nN , 'uint8' );
  nd( 1:4*S        , : ) = reshape( typecast( single( b4(:) ) ,'uint8') , 4*S , nN );
  nd( 4*S+(1:16)   , : ) = reshape( typecast( int32(  c4(:) ) ,'uint8') , 16  , nN );
  nd = nd(:);
end
