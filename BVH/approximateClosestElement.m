function [eID, cp, d, bc, F] = approximateClosestElement( M , P , B , Dmax )
%APPROXIMATECLOSESTELEMENT  Localizador APROXIMADO por 1-ring-of-nearest-vertex.
%
%   [e,cp,d,bc,F] = approximateClosestElement( M , P )         construye blob
%   [ ... ]       = approximateClosestElement( M , P , B )     reutiliza blob
%   [ ... ]       = approximateClosestElement( {M,B} , P [, Dmax] )
%   [ ... ]       = approximateClosestElement( ...  , Dmax )   radio (escalar o nP-vector)
%   Ba            = approximateClosestElement( M )             solo construir el blob
%
%   Heuristica en dos etapas fusionadas en el MEX: (1) el VERTICE mas cercano
%   al punto (BVH de puntos sobre los vertices usados de la malla) y (2) la
%   distancia EXACTA a su abanico de elementos incidentes (EsuP, empacado CSR
%   en el blob). Todos los celltypes (puntos/segmentos/triangulos/TETS/mixtas
%   0-padded).
%
%   GARANTIAS Y LIMITES (aproximado a conciencia):
%     - d es SIEMPRE una cota superior de la distancia real (el abanico es un
%       subconjunto real de la malla): d_apx >= d_exacta.
%     - el elemento ganador es el correcto tipicamente el 95-100%% de las veces
%       (ver bench_approximate); los fallos se concentran cerca de la
%       superficie y su error esta acotado por la escala local de la malla.
%     - puntos INTERIORES de mallas de tets pueden no detectarse como d = 0 si
%       el tet contenedor no toca al vertice mas cercano.
%   Mismas convenciones que bvhClosestElement: Dmax escalar o nP-vector (aqui
%   corta por la distancia AL VERTICE, etapa 1); miss -> e=0, d=Inf, cp/bc=NaN;
%   punto no finito -> d=NaN; bc robustas; F.type/F.onBoundary con 5 salidas.
%
%   El blob es PROPIO (BVH de vertices + abanicos + elementos empacados), NO
%   intercambiable con el de BVH(M); coste de construccion similar al blob
%   exacto. Blob rancio (malla editada) -> ERROR, reconstruye con
%   Ba = approximateClosestElement(M).
%
%   Requiere el MEX:  mex COMPFLAGS="$COMPFLAGS /openmp" -lut approximateClosestElement_mx.cpp
%
% See also bvhClosestElement, BVH, bench_approximate.

  %forma constructor: Ba = approximateClosestElement( M )
  if nargin == 1
    eID = buildApprox( M );
    return;
  end

  %bundled form: approximateClosestElement( {M,B} , P , Dmax )
  if iscell( M )
    if nargin > 2, Dmax = B; end          %shift: 3rd arg is Dmax in this form
    B = M{2};  M = M{1};
  elseif nargin < 3
    B = [];
  end
  if ~exist('Dmax','var') || isempty( Dmax ), Dmax = Inf; end
  if ~( isnumeric(Dmax) && isreal(Dmax) && isvector(Dmax) && all( Dmax(:) >= 0 ) )
    error('approximateClosestElement:Dmax', ...
          'Dmax must be a nonnegative scalar or an nP-vector.');
  end
  Dmax = double( Dmax(:) );
  if exist( 'approximateClosestElement_mx' ,'file') ~= 3
    error('approximateClosestElement:mex','approximateClosestElement_mx is not compiled (mex COMPFLAGS="$COMPFLAGS /openmp" -lut approximateClosestElement_mx.cpp).');
  end

  P = double( P );  P(:,end+1:3) = 0;
  nP = size( P ,1);
  if ~isscalar( Dmax ) && numel( Dmax ) ~= nP
    error('approximateClosestElement:Dmax','per-point Dmax must have nP elements (%d vs %d).', numel(Dmax) , nP );
  end

  if isempty( B ), B = buildApprox( M ); end
  if ~isfield( B , 'approx' )
    error('approximateClosestElement:B', ...
          'B is not an approximate blob (build it with Ba = approximateClosestElement(M); the BVH(M) blob is a different animal).');
  end

  %staleness spot-check (mismo patron que bvhClosestElement)
  Xw = double( M.xyz ); Xw(:,end+1:3) = 0;
  ok = size( Xw ,1) == size( B.X ,1) && isequal( size( double(M.tri) ) , size( B.elTri ) );
  if ok && size( Xw ,1) > 0
    ii = unique( round( linspace( 1 , size(Xw,1) , 4 ) ) );
    Yw = B.X(ii,:) * B.frame(1:3,1:3).' + B.frame(1:3,4).';
    ok = max(max( abs( Yw - Xw(ii,:) ) )) <= 1e-6 * max( 1 , max(max( abs( Xw(ii,:) ) )) );
  end
  if ~ok
    error('approximateClosestElement:staleBVH', ...
          'B does not match M (stale or foreign blob). Rebuild it: Ba = approximateClosestElement(M).');
  end

  %global frame: query en espacio BUILD, des-transformar al final
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

  Tri = double( B.elTri );

  %baricentricas REGION-EXACTAS del MEX (igual que bvhClosestElement): aristas
  %y vertices con ceros exactos, reconstruccion a precision de maquina en
  %slivers, suma 1 y en [0,1] forzados; devuelve 4 columnas -> se recorta al
  %ancho de faceta de la malla
  if nargout > 3
    [ eID , cp , d , bc4 ] = approximateClosestElement_mx( P , B , maxNumCompThreads , DmaxEff );
    bc = bc4( : , 1:size( Tri ,2) );
    bc( eID == 0 , : ) = NaN;
  else
    [ eID , cp , d ] = approximateClosestElement_mx( P , B , maxNumCompThreads , DmaxEff );
  end

  %clasificacion del FEATURE del punto mas cercano (+ flag de borde abierto)
  if nargout > 4
    tolF = 1e-9;
    nz   = sum( bc > tolF , 2 );
    F = struct();
    F.type = zeros( nP ,1);
    w = eID > 0;
    F.type(w) = min( nz(w) , 4 );
    F.onBoundary = false( nP ,1);

    k  = sum( Tri > 0 ,2);                       %nonzero nodes per face
    kk = size( Tri ,2);
    if kk == 3 && all( k == 3 )
      Bed = MeshBoundary( Tri );
      if ~isempty( Bed )
        Bed  = sort( Bed ,2);
        Bvx  = unique( Bed );
        wv = find( F.type == 1 );
        if ~isempty( wv )
          [~,imax] = max( bc(wv,1:3) ,[],2);
          nd = Tri( sub2ind( size(Tri) , eID(wv) , imax ) );
          F.onBoundary(wv) = ismember( nd , Bvx );
        end
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
    elseif kk == 2
      fe = MeshBoundary( Tri );
      if ~isempty( fe )
        wv = find( F.type == 1 );
        if ~isempty( wv )
          [~,imax] = max( bc(wv,1:2) ,[],2);
          nd = Tri( sub2ind( size(Tri) , eID(wv) , imax ) );
          F.onBoundary(wv) = ismember( nd , fe(:) );
        end
      end
    end
  end

  %des-transformar: de espacio build a mundo (e/bc/F invariantes; los misses
  %propagan: cp NaN sigue NaN, d Inf sigue Inf)
  if hasF
    cp = cp * Af.' + tf;
    d  = d  * fscale;
  end

end

%------------------------------------------------------------------ blob
function Ba = buildApprox( M )
%blob aproximado: BVH de PUNTOS sobre los vertices USADOS (hojas grandes:
%los puntos son baratos de testear) + abanicos CSR + elementos empacados +
%kernels 4-wide precomputados: pt4 (puntos en bloques SoA para AVX) y, en
%mallas de PUROS triangulos, fan4 (los abanicos como bloques PreTri4)
  T = double( M.tri );
  if isempty( T )
    error('approximateClosestElement:build','the mesh has no elements.');
  end
  w  = T > 0;
  [ re , ~ ] = find( w );                    %elemento de cada par (elem,vertice)
  ve = T( w );                               %vertice de cada par
  uv = unique( ve );                         %vertices USADOS (orden ascendente)
  nV = numel( uv );

  Ba = BVH( struct( 'xyz' , M.xyz , 'tri' , uv(:) ) , [ 32 , 128 ] );  %hojas grandes:
  %con el kernel Pt4 los tests de punto son casi gratis y la travesia manda
  %(barrido medido: [32 128] optimo en near/mid/far)

  [ ~ , loc ] = ismember( ve , uv );         %fila del blob de cada par
  [ sl , o ] = sort( loc );
  fe  = re( o );
  cnt = accumarray( sl , 1 , [ nV , 1 ] );
  fanStart    = [ 0 ; cumsum( cnt ) ];
  Ba.fanStart = int32( fanStart );
  Ba.fanEl    = int32( fe );

  Xb  = Ba.X;                                %vertices en espacio BUILD
  nEl = size( T ,1);
  elV = zeros( 12 , nEl );
  for kk = 1:size( T ,2)
    idx = T(:,kk);  mrow = idx > 0;
    elV( 3*kk-2:3*kk , mrow ) = Xb( idx(mrow) ,:).';
  end
  Ba.elV    = elV;
  Ba.elT    = int32( sum( T > 0 , 2 ) );
  Ba.elTri  = int32( T );

  %pt4: los vertices empacados (orden perm, el de las hojas) en bloques SoA
  %de 4 [x0..x3 y0..y3 z0..z3] -- hojas con cargas AVX alineadas, lanes
  %filtradas por indice en el kernel (el padding replica el ultimo punto)
  P3  = Ba.pkS( 1:3 , : );                   %centros = los vertices, packed
  nBk = ceil( nV / 4 );
  P3  = [ P3 , repmat( P3(:,end) , 1 , 4*nBk - nV ) ];
  Ba.pt4 = reshape( permute( reshape( P3 , 3 , 4 , nBk ) , [2 1 3] ) , 12 , nBk );

  %fan4: SOLO mallas de puros triangulos -- cada abanico como bloques
  %PreTri4 (A, AB, AC en SoA de 4; lanes de relleno con id 0, filtradas)
  if size( T ,2) == 3 && all( T(:) > 0 )
    nb4       = ceil( cnt / 4 );
    f4s       = [ 0 ; cumsum( nb4 ) ];
    NB        = f4s(end);
    pos       = ( 0:numel(sl)-1 ).' - fanStart( sl );      %posicion dentro del abanico
    blk       = f4s( sl ) + floor( pos/4 );                %bloque (0-based)
    lane      = mod( pos , 4 );
    A  = Xb( T(fe,1) ,:);
    AB = Xb( T(fe,2) ,:) - A;
    AC = Xb( T(fe,3) ,:) - A;
    cols = [ A , AB , AC ];                                %n x 9
    fan4 = zeros( 36*NB , 1 );
    for g = 1:9
      fan4( blk*36 + (g-1)*4 + lane + 1 ) = cols(:,g);
    end
    f4id = zeros( 4*NB , 1 , 'int32' );
    f4id( blk*4 + lane + 1 ) = int32( fe );
    Ba.fan4      = reshape( fan4 , 36 , NB );
    Ba.fan4id    = reshape( f4id , 4 , NB );
    Ba.fan4Start = int32( f4s );
  end

  Ba.approx = 1;
end
