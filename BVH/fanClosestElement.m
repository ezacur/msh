function [eID, cp, d, bc, F] = fanClosestElement( M , P , Dmax )
%fanClosestElement  Localizador APROXIMADO por abanico del nodo mas cercano.
%
%   [e,cp,d,bc,F] = fanClosestElement( M , P )             point-BVH sobre los
%                                     nodos de M construido al vuelo (etapa 1)
%   [ ... ]       = fanClosestElement( {M,node}  , P )     UN nodo semilla para
%                                     TODOS los puntos (sin etapa 1)
%   [ ... ]       = fanClosestElement( {M,nodes} , P )     un nodo semilla POR
%                                     PUNTO (vector nP, sin etapa 1)
%   [ ... ]       = fanClosestElement( {M,Bnodes} , P )    reutiliza un
%                                     point-BVH sobre los nodos de M:
%                Bnodes = BVH( struct('xyz',M.xyz,'tri',(1:size(M.xyz,1)).') )
%   [ ... ]       = fanClosestElement( ... , Dmax )        radio de busqueda
%                                     (escalar o nP-vector)
%
%   Dos etapas: (1) el NODO mas cercano a cada punto -- por point-BVH, o dado
%   por el llamador -- y (2) la distancia EXACTA al ABANICO de elementos
%   incidentes a ese nodo (todos los celltypes: puntos/segmentos/triangulos/
%   TETS/mixtas 0-padded). NUNCA toca el arbol de elementos, ni siembra nada:
%   sembrar el arbol exacto no gana nada (medido, ver bench_seedCeiling).
%
%   GARANTIAS Y LIMITES (aproximado a conciencia):
%     - d es SIEMPRE una cota superior de la distancia real (el abanico es un
%       subconjunto real de la malla): d_fan >= d_exacta.
%     - con nodos de etapa 1 el elemento ganador es el correcto tipicamente el
%       95-100% de las veces; los fallos se concentran cerca de la superficie
%       y su error esta acotado por la escala LOCAL de la malla.
%     - puntos INTERIORES de mallas de tets pueden no dar d = 0 si el tet
%       contenedor no toca al nodo mas cercano.
%     - un nodo semilla sin elementos incidentes (vertice suelto) da miss.
%   Dmax corta la etapa 1 por la distancia AL NODO (nodo mas alla de Dmax ->
%   miss) y filtra el resultado final (d >= Dmax -> miss). Miss -> e=0, d=Inf,
%   cp/bc=NaN; punto no finito -> d=NaN. bc REGION-EXACTAS del MEX (ceros
%   exactos en aristas/vertices, >=0, suman 1); F.type/F.onBoundary con 5
%   salidas, como bvhClosestElement.
%
%   El abanico (EsuP CSR) se recalcula de M en CADA llamada (O(nnz) sort): la
%   funcion es deliberadamente COMPONIBLE y sin blob propio -- la unica
%   estructura reutilizable es el point-BVH estandar Bnodes. Para lotes
%   grandes el recalculo es despreciable; para bucles de pocas queries,
%   pasa los nodos ya calculados y no pagas nada.
%
%   Bnodes debe estar construido sobre TODOS los nodos de M en su orden
%   (tri = (1:nV).'): sus ids de elemento son entonces ids de nodo de M.
%
%   Requiere el MEX:  mex COMPFLAGS="$COMPFLAGS /openmp" -lut fanClosestElement_mx.cpp
%
% See also bvhClosestElement, BVH, bvhIntersectRay, bench_fanClosestElement.

  %target: M solo, o {M,semilla} donde semilla = node | nodes | Bnodes
  seedN = [];  Bn = [];
  if iscell( M )
    S2 = M{2};  M = M{1};
    if isstruct( S2 )
      Bn = S2;
    elseif isnumeric( S2 ) && ~isempty( S2 )
      seedN = double( S2(:) );
    else
      error('fanClosestElement:seed', ...
            'the 2nd cell element must be a seed node (scalar), an nP-vector of nodes, or a point-BVH blob.');
    end
  end
  if nargin < 3 || isempty( Dmax ), Dmax = Inf; end
  if ~( isnumeric(Dmax) && isreal(Dmax) && isvector(Dmax) && all( Dmax(:) >= 0 ) )
    error('fanClosestElement:Dmax','Dmax must be a nonnegative scalar or an nP-vector.');
  end
  Dmax = double( Dmax(:) );
  if exist( 'fanClosestElement_mx' ,'file') ~= 3
    error('fanClosestElement:mex','fanClosestElement_mx is not compiled (mex COMPFLAGS="$COMPFLAGS /openmp" -lut fanClosestElement_mx.cpp).');
  end

  V = double( M.xyz );  V(:,end+1:3) = 0;  nV = size( V ,1);
  T = double( M.tri );  nEl = size( T ,1);  kkT = size( T ,2);
  P = double( P );  P(:,end+1:3) = 0;  nP = size( P ,1);
  if ~isscalar( Dmax ) && numel( Dmax ) ~= nP
    error('fanClosestElement:Dmax','per-point Dmax must have nP elements (%d vs %d).', numel(Dmax) , nP );
  end

  %---- etapa 1: el nodo semilla de cada punto -----------------------------
  if ~isempty( seedN )
    if isscalar( seedN ), seedN = repmat( seedN , nP , 1 ); end
    if numel( seedN ) ~= nP
      error('fanClosestElement:seed','nodes must be a scalar or an nP-vector (%d vs %d).', numel(seedN) , nP );
    end
    if any( seedN < 1 | seedN > nV | seedN ~= round( seedN ) )
      error('fanClosestElement:seed','seed nodes must be integers in [1,%d].', nV );
    end
    nodes = seedN;
  else
    Mn = struct( 'xyz' , V , 'tri' , (1:nV).' );
    if isempty( Bn )
      Bn = BVH( Mn );
    else
      %Bnodes debe ser un point-BVH sobre TODOS los nodos de M, en su orden
      %(la comprobacion de coordenadas la hace bvhClosestElement por spot-check)
      if ~isfield( Bn ,'Tri') || size( Bn.Tri ,2) ~= 1 || size( Bn.Tri ,1) ~= nV ...
                              || ~isequal( Bn.Tri(:) , (1:nV).' )
        error('fanClosestElement:B', ...
              'Bnodes is not a point-BVH over the nodes of M (build it: BVH( struct(''xyz'',M.xyz,''tri'',(1:size(M.xyz,1)).'') )).');
      end
    end
    nodes = bvhClosestElement( {Mn,Bn} , P , Dmax );   %corta por distancia al NODO
  end

  %---- abanicos (EsuP como CSR), recalculados de M en cada llamada --------
  w   = T > 0;
  nid = T(w);
  eid = repmat( (1:nEl).' , 1 , kkT );  eid = eid(w);
  [ nid , o ] = sort( nid );
  fanEl    = int32( eid(o) );
  fanStart = int32( cumsum( [ 1 ; accumarray( nid , 1 , [ nV , 1 ] ) ] ) );

  %---- etapa 2: barrido exacto del abanico en el MEX ----------------------
  if nargout > 3
    [ eID , cp , d , bc4 ] = fanClosestElement_mx( P , V , int32(T) , int32(nodes) , ...
                                                   fanStart , fanEl , maxNumCompThreads , Dmax );
    bc = bc4( : , 1:kkT );
    bc( eID == 0 , : ) = NaN;
  else
    [ eID , cp , d ] = fanClosestElement_mx( P , V , int32(T) , int32(nodes) , ...
                                             fanStart , fanEl , maxNumCompThreads , Dmax );
  end

  %---- clasificacion del feature del cp (+ borde abierto), como la familia --
  if nargout > 4
    tolF = 1e-9;
    nz   = sum( bc > tolF , 2 );
    F = struct();
    F.type = zeros( nP ,1);
    wq = eID > 0;
    F.type(wq) = min( nz(wq) , 4 );              %1 vtx, 2 edge, 3 face, 4 inside
    F.onBoundary = false( nP ,1);

    k = sum( T > 0 ,2);
    if kkT == 3 && all( k == 3 )                 %pure triangle surface
      Bed = MeshBoundary( T );
      if ~isempty( Bed )
        Bed  = sort( Bed ,2);
        Bvx  = unique( Bed );
        wv = find( F.type == 1 );
        if ~isempty( wv )
          [~,imax] = max( bc(wv,1:3) ,[],2);
          nd = T( sub2ind( size(T) , eID(wv) , imax ) );
          F.onBoundary(wv) = ismember( nd , Bvx );
        end
        we = find( F.type == 2 );
        if ~isempty( we )
          E2 = zeros( numel(we) ,2);
          for q = 1:numel( we )
            act = find( bc(we(q),1:3) > tolF );
            E2(q,:) = T( eID(we(q)) , act );
          end
          F.onBoundary(we) = ismember( sort(E2,2) , Bed , 'rows' );
        end
      end
    elseif kkT == 2                              %wireframe: free ends
      fe = MeshBoundary( T );
      if ~isempty( fe )
        wv = find( F.type == 1 );
        if ~isempty( wv )
          [~,imax] = max( bc(wv,1:2) ,[],2);
          nd = T( sub2ind( size(T) , eID(wv) , imax ) );
          F.onBoundary(wv) = ismember( nd , fe(:) );
        end
      end
    end
    %tets, nubes de puntos y mixtas: onBoundary queda false (documentado)
  end
end
