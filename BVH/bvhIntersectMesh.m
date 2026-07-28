function [ C , K , PR ] = bvhIntersectMesh( MA , MB )
%bvhIntersectMesh  Curva de interseccion de dos superficies trianguladas.
%
%   [ C , K , PR ] = bvhIntersectMesh( MA , MB )
%
%   Devuelve la curva donde las superficies MA y MB se cortan, como una malla de
%   SEGMENTOS (celltype 3 = VTK_LINE, 2 nodos por celda):
%     C.xyz   nC x 3   nodos de la curva
%     C.tri   nS x 2   segmentos
%     K       nC x 8   la CLAVE CANONICA de cada nodo (su procedencia):
%                        [ kind onB ia ib ja jb t s ]
%                        kind 1 = vertice ia
%                        kind 2 = punto en la arista canonica (ia<ib), parametro t
%                        kind 3 = cruce de la arista (ia<ib) de A con la (ja<jb) de B
%     PR      nS x 2   el par de caras ( caraA , caraB ) que produjo cada segmento
%
%   POR QUE CLAVES Y NO COORDENADAS. Cada punto de la curva se identifica por su
%   procedencia COMBINATORIA (que arista lo genero y con que parametro), no por su
%   posicion. Asi dos triangulos vecinos que comparten una arista producen claves
%   BIT-IDENTICAS y soldar es un `unique` por filas, sin tolerancia geometrica
%   ninguna. Es la tecnica de MESH\MeshZeroContour.m llevada a la interseccion, y
%   descansa en que los predicados de orientacion son EXACTOS (bvhKernels.h):
%   "que aristas cruzan" es entonces una decision combinatoria que depende solo de
%   (arista, triangulo) y nunca de por que triangulo se llego, asi que los dos
%   triangulos que comparten una arista coinciden SIEMPRE.
%
%   LA REGLA DE PARIDAD (zonas coplanares). Dos triangulos coplanares que se
%   solapan lo hacen en un AREA, no en una curva, y el kernel devuelve el borde de
%   ese solape. Sumar los bordes par a par sobre una region coincidente daria toda
%   la red interna de aristas -- medido: 1003 segmentos de ruido contra 30 utiles,
%   y la curva llena de bifurcaciones. La correccion es que la contribucion
%   coplanar es el borde de la UNION de los solapes: una arista INTERIOR a la
%   union la recorren DOS poligonos de solape y una del BORDE solo UNO, asi que se
%   cancelan las de multiplicidad PAR. La cancelacion se aplica SOLO a la
%   contribucion coplanar: en el caso transversal una multiplicidad 2 significa el
%   mismo segmento hallado por los dos triangulos que comparten una arista, y
%   cancelarlo borraria curva de verdad.
%
%   BLOB REUTILIZABLE. Cualquiera de los dos argumentos admite el bundle {M,B}.
%   El que se aprovecha es el del SEGUNDO: el recorrido desciende el arbol de MB
%   con los triangulos de MA, asi que un blob en el primero no tiene donde
%   emplearse y se ignora. El de MB se usa si sirve (mismo tamano, mismo frame y
%   cage 'aabb'); si no, se reconstruye. Construirlo es lo mas caro de la llamada
%   -- 67 de 90 ms en una malla de 87k caras -- asi que pasarlo hecho compensa
%   cuando se cortan varias mallas contra la misma.
%
%   AUTOINTERSECCION. Si las dos mallas son BIT-IDENTICAS se cambia de pregunta:
%   se buscan pares de triangulos NO ADYACENTES que se corten. Un triangulo corta
%   trivialmente consigo mismo y con quien comparte arista o vertice, y eso no es
%   un autocorte -- sin excluirlo, una esfera cerrada devolvia su red de aristas
%   entera (medido: 134 de 162 nodos con grado >2 en vez de curva vacia). La
%   deteccion es por igualdad EXACTA, no por semejanza: si MB es MA deformada, los
%   ids de vertice viven en el mismo espacio pero las superficies son distintas y
%   excluir vecinos falsearia el resultado.
%
%   Requiere triTriPairs_mx compilado (mex triTriPairs_mx.cpp).
%
% See also MeshZeroContour, bvhClosestElement, bvhIntersectRay, BVH.

  [ MA , BA ] = asMesh( MA , 'MA' );                              %#ok<ASGLU>
  [ MB , BB ] = asMesh( MB , 'MB' );
  C  = struct( 'xyz' , zeros(0,3) , 'tri' , zeros(0,2) , 'celltype' , 3 );
  K  = zeros(0,8);  PR = zeros(0,2);
  if isempty( MA.tri ) || isempty( MB.tri ), return; end

  %AUTOINTERSECCION: si las dos mallas son BIT-IDENTICAS, la pregunta es otra --
  %un triangulo corta trivialmente consigo mismo y con sus vecinos por la arista o
  %el vertice compartidos, y eso no es un autocorte. Se detecta por igualdad exacta
  %(no por "parecidas"): si B es A deformada, los ids de vertice viven en el mismo
  %espacio pero las superficies son distintas, y excluir vecinos seria un error.
  selfMode = isequal( MA.xyz , MB.xyz ) && isequal( MA.tri , MB.tri );

  P = candidatePairs( MA , MB , BB );
  if selfMode && ~isempty( P )
    P = P( P(:,1) < P(:,2) ,:);                    %cada par NO ordenado, una vez
    Ta = MA.tri( P(:,1) ,:);  Tb = MA.tri( P(:,2) ,:);
    sh = any( Ta(:,1) == Tb ,2) | any( Ta(:,2) == Tb ,2) | any( Ta(:,3) == Tb ,2);
    P  = P( ~sh ,:);                               %fuera los que comparten vertice
  end
  if isempty( P ), return; end

  [ Kr , XYZ , S , PRr , CODE ] = triTriPairs_mx( MA.xyz , MA.tri , ...
                                                 MB.xyz , MB.tri , P );
  if isempty( Kr ), return; end

  %--- soldadura EXACTA: por clave, y ademas por igualdad exacta de COORDENADAS
  %    para resolver el hueco de identidad entre mallas (el mismo punto es
  %    vertice de A y de B con ids distintos, luego dos claves). Es licito porque
  %    la coordenada de una clave de vertice es una coordenada de ENTRADA.
  [ Ku , ~ , g ] = unique( Kr , 'rows' , 'stable' );
  first = zeros( size(Ku,1) ,1);
  for i = 1:size(Ku,1), first(i) = find( g == i ,1); end
  XW = XYZ( first ,:);
  [ ~ , ia , ic ] = unique( XW , 'rows' , 'stable' );
  K  = Ku( ia ,:);   XW = XW( ia ,:);   g = ic( g );

  Sw = g( S );                                  %segmentos en indices soldados
  keep = Sw(:,1) ~= Sw(:,2);                    %fuera los degenerados
  Sw = Sw( keep ,:);  PRs = PRr( keep ,:);
  if isempty( Sw ), C.xyz = XW; return; end
  Sw = sort( Sw ,2);

  %--- de que rama viene cada segmento
  isCop = CODE( pairIndex( P , PRs ) ) == 3;

  %--- paridad SOLO en la contribucion coplanar
  Sc = Sw( isCop ,:);  Pc = PRs( isCop ,:);
  if ~isempty( Sc )
    [ u , ~ , h ] = unique( Sc , 'rows' );
    cnt = accumarray( h ,1);
    odd = find( mod( cnt ,2) == 1 );            %borde de la UNION
    keepC = ismember( h , odd );
    %un representante por segmento superviviente
    [ ~ , rep ] = unique( h( keepC ) , 'stable' );
    idx = find( keepC );  idx = idx( rep );
    Sc = Sc( idx ,:);  Pc = Pc( idx ,:);
  end
  %--- la contribucion transversal se deduplica sin mas
  St = Sw( ~isCop ,:);  Pt = PRs( ~isCop ,:);
  if ~isempty( St )
    [ St , it ] = unique( St , 'rows' );  Pt = Pt( it ,:);
  end

  Sall = [ St ; Sc ];  PR = [ Pt ; Pc ];
  if isempty( Sall ), C.xyz = XW; return; end
  [ Sall , iu ] = unique( Sall , 'rows' );  PR = PR( iu ,:);

  %--- quitar nodos que no participan en ningun segmento
  used = false( size(XW,1) ,1);  used( Sall(:) ) = true;
  ren  = zeros( size(XW,1) ,1);  ren( used ) = 1:sum(used);
  C.xyz = XW( used ,:);
  C.tri = ren( Sall );
  K     = K( used ,:);
end

%% ---------------------------------------------------------------- helpers
function [ M , B ] = asMesh( M , who )
%admite M o el bundle {M,B}. El blob que SE USA es el del SEGUNDO argumento de
%bvhIntersectMesh: el recorrido desciende el arbol de B con los triangulos de A,
%asi que un blob en el primer argumento no tiene donde emplearse (se ignora, y
%esta documentado arriba en lugar de fingir que sirve).
  B = [];
  if iscell( M )
    if numel( M ) >= 2, B = M{2}; end
    M = M{1};
  end
  if ~isstruct( M ) || ~isfield( M ,'xyz') || ~isfield( M ,'tri')
    error('bvhIntersectMesh:input','%s must be a struct with .xyz and .tri.', who );
  end
  M.xyz = double( M.xyz );  M.xyz(:,end+1:3) = 0;
  M.tri = double( M.tri );
  if ~isempty( M.tri ) && size( M.tri ,2) ~= 3
    error('bvhIntersectMesh:tri','%s: only TRIANGLE surfaces (nF x 3).', who );
  end
end

function ok = blobFits( B , M )
%¿sirve este blob para esta malla? Chequeo baratito, igual que el de staleness de
%bvhClosestElement: tamanos, y 4 vertices llevados del build frame al mundo.
  ok = false;
  if ~isstruct( B ), return; end
  for f = { 'node4','child4','srange','pkE','vol','frame','X','Tri' }
    if ~isfield( B , f{1} ), return; end
  end
  if double( B.vol ) ~= 2, return; end             %el recorrido solo hace 'aabb'
  X = M.xyz;
  if size( B.X ,1) ~= size( X ,1) || ~isequal( size(double(B.Tri)) , size(M.tri) )
    return;
  end
  if size( X ,1) > 0
    ii = unique( round( linspace( 1 , size(X,1) , 4 ) ) );
    Y  = B.X(ii,:) * B.frame(1:3,1:3).' + B.frame(1:3,4).';
    if max(max(abs( Y - X(ii,:) ))) > 1e-6 * max( 1 , max(max(abs( X(ii,:) ))) )
      return;
    end
  end
  ok = true;
end

function k = pairIndex( P , PRs )
%fila de P a la que corresponde cada par ( caraA , caraB )
  [ ~ , k ] = ismember( PRs , P , 'rows' );
  if any( k == 0 )
    error('bvhIntersectMesh:prov','internal: segment provenance not in the pair list.');
  end
end

function P = candidatePairs( A , B , Bb )
%pares candidatos, recorriendo el BVH de B con cada triangulo de A. Conservador:
%nunca descarta una interseccion real, solo deja pasar candidatos de mas.
%
%   Reutiliza el blob que venga en el bundle {MB,BB} si sirve (mismo tamano, mismo
%   frame y cage 'aabb'); si no, lo construye. Construirlo es lo mas caro de la
%   llamada -- 67 de 90 ms en una malla de 87k caras -- asi que pasarlo hecho
%   cuando se van a cortar varias mallas contra la misma B compensa.
%
%   El camino por FUERZA BRUTA sigue aqui como reserva y como REFERENCIA de
%   validacion. OJO: los dos son conservadores pero ninguno contiene al otro (la
%   fuerza bruta compara AABBs en MUNDO y el recorrido en el BUILD FRAME de B);
%   lo que test_bvhIntersectMesh comprueba es que los pares que SE CORTAN estan en
%   los dos y que la curva sale identica.
  if exist( 'bvhTriPairs_mx' ,'file') == 3
    if ~blobFits( Bb , B )
      Bb = BVH( struct('xyz',B.xyz,'tri',B.tri) );    %cage 'aabb' por defecto
    end
    P = bvhTriPairs_mx( toBuildFrame( A.xyz , Bb ) , A.tri , Bb );
    if isempty( P ), P = zeros(0,2); end
    return;
  end
  P = bruteForcePairs( A , B );
end

function P = bruteForcePairs( A , B )
%todas las parejas cuyas AABB de cara se solapan. O(nFA*nFB): solo para mallas
%pequenas y medianas, y para validar el recorrido.
  a = faceBoxes( A );  b = faceBoxes( B );
  nA = size( a ,1);
  if nA * size(b,1) > 4e8
    error('bvhIntersectMesh:tooBig', ...
          '%d x %d pares: la fuerza bruta no escala; compila bvhTriPairs_mx.' , ...
          nA , size(b,1) );
  end
  I = cell( nA ,1);
  for i = 1:nA
    w = a(i,1) <= b(:,4) & a(i,4) >= b(:,1) & ...
        a(i,2) <= b(:,5) & a(i,5) >= b(:,2) & ...
        a(i,3) <= b(:,6) & a(i,6) >= b(:,3);
    if any( w ), I{i} = [ repmat(i,sum(w),1) , find(w) ]; end
  end
  P = vertcat( I{:} );
  if isempty( P ), P = zeros(0,2); end
end

function X = toBuildFrame( X , B )
%el blob guarda su geometria en el BUILD frame: hay que llevar A hasta alli.
%Para una semejanza  mundo = A*xf + t  con  A.'*A = s^2 I :  xf = (x-t)*inv(A).'
  X = double( X );  X(:,end+1:3) = 0;
  if isequal( B.frame , eye(4) ), return; end
  Af = B.frame(1:3,1:3);
  s2 = trace( Af.'*Af )/3;
  X  = ( X - B.frame(1:3,4).' ) * ( Af / s2 );
end

function Bx = faceBoxes( M )
  X = M.xyz;  T = M.tri;  Bx = zeros( size(T,1) ,6);
  for c = 1:3
    V = [ X(T(:,1),c) , X(T(:,2),c) , X(T(:,3),c) ];
    Bx(:,c)   = min( V ,[],2);
    Bx(:,c+3) = max( V ,[],2);
  end
end
