function test_cache( mode )
%TEST_CACHE  Bateria de verificacion del patron de cache (estado CONSOLIDADO).
%
%   test_cache()          corre todo; cualquier fallo es un assert con nombre.
%   test_cache('record')  regenera los golden de los demos. SOLO tras un cambio
%                         deliberado de comportamiento, nunca para "que pase".
%
%   Correr ANTES y DESPUES de tocar cacheHandle.m, cacheProxy.m o cachedOwner.m.
%   Los demos (cartDemo/routeDemo) se comparan BYTE A BYTE contra su golden:
%   asi se detecta cualquier cambio de comportamiento, no solo los que revientan.
%
%   Uso:  addpath('C:\repos\msh','C:\repos\msh\dev','C:\repos\msh\dev\tests'); test_cache
%
% See also cachedOwner, cartDemo, routeDemo.

  if nargin < 1, mode = 'check'; end
  here = fileparts( mfilename( 'fullpath' ) );

  %% -------- 1) los demos, byte a byte contra su golden
  %el sidecar apunta EN QUE RELEASE se grabaron: un golden que difiere en otra
  %version de MATLAB puede ser solo formato de print, y el mensaje lo dice
  relf = fullfile( here , 'goldens_release.txt' );
  rel  = '?';  if exist( relf , 'file' ), rel = strtrim( fileread( relf ) ); end
  for d = { 'cartDemo' , 'routeDemo' , 'serieDemo' , 'factoDemo' , 'ivpDemo' , ...
            'collatzDemo' , 'linsysDemo' }
    out = evalc( d{1} );
    gf  = fullfile( here , [ d{1} '_expected.txt' ] );
    if strcmp( mode , 'record' )
      fid = fopen( gf , 'w' );  fwrite( fid , out );  fclose( fid );
      fprintf( 'GRABADO golden de %s\n' , d{1} );
    else
      assert( isequal( out , fileread( gf ) ) , ...
              'test_cache:goldenDiff' , ...
              [ 'la salida de %s cambio respecto al golden. OJO: goldens grabados en R%s y ' ...
                'esto corre en R%s -- si las releases difieren puede ser solo formato de ' ...
                'print: compara ANTES de tocar el motor.' ] , d{1} , rel , version('-release') );
      fprintf( 'OK  %-9s identico al golden\n' , d{1} );
    end
  end
  if strcmp( mode , 'record' )
    fid = fopen( relf , 'w' );  fwrite( fid , version('-release') );  fclose( fid );
    fprintf('goldens regenerados (R%s); vuelve a correr test_cache()\n', version('-release') );
    return;
  end

  %% -------- 2) guardas de Define: nombres que sombrean y politicas invalidas
  C = cart();  C.DEBUG = false;  C = C.Demo();
  expectError( @() C.Define( 'PRICE'  , @(c)1 ) , 'cachedOwner:shadow'         );  % propiedad
  expectError( @() C.Define( 'Add'    , @(c)1 ) , 'cachedOwner:shadow'         );  % metodo
  expectError( @() C.Define( 'CP'     , @(c)1 ) , 'cachedOwner:shadow'         );  % la puerta
  expectError( @() C.Define( 'total_' , @(c)1 ) , 'cachedOwner:reservedSuffix' );  % sufijo _
  expectError( @() C.Define( '3mal'   , @(c)1 ) , 'cachedOwner:name'           );  % no es varname
  expectError( @() C.Define( 'x' , @(c)1 , 'noexiste' , [] ) , 'cachedOwner:unknownEvent' );
  expectError( @() pDouble()                    , 'cachedOwner:twoIncremental' );  % 2 incrementales
  %miembros HIDDEN y PRIVADOS: properties()/methods() no los ven y se colaban
  for h = { 'CP_EVENTS' , 'DEFS' , 'CACHE' , 'cachedAccess' , 'displayCachedView' , 'loadobj' }
    expectError( @() C.Define( h{1} , @(c)1 ) , 'cachedOwner:shadow' );
  end
  %firmas AMBIGUAS: nargin -1 y -2 caian del lado absoluto en silencio
  expectError( @() C.Define( 'amb1' , @(c)1 , 'additem' , @(varargin)   1 ) , ...
               'cachedOwner:ambiguousHandler' );
  expectError( @() C.Define( 'amb2' , @(c)1 , 'additem' , @(v,varargin) 1 ) , ...
               'cachedOwner:ambiguousHandler' );
  %y la forma EXPLICITA de "incremental con args variables" si pasa
  Cx = C.Define( 'amb3' , @(c) 0 , 'additem' , @(v,c,varargin) v + 1 );
  v0 = Cx.amb3;                        % poblarla: fire solo recorre lo CACHEADO %#ok<NASGU>
  Cx = Cx.Add( 'x' , 1 , 1 );
  a3 = Cx.amb3;
  assert( a3 == 1 , 'test_cache:varargIncremental' , ...
          '@(v,c,varargin) debe valer como incremental y dio %g.' , a3 );
  fprintf( 'OK  Define    rechaza sombras (tambien Hidden/privadas), erratas y doble incremental\n' );

  %% -------- 3) fire valida los nombres (la errata NO puede dar cache rancia)
  T = pTypo();  v = T.dbl;                                                  %#ok<NASGU>
  expectError( @() T.BumpTypo( 5 ) , 'cachedOwner:unknownEvent' );
  T2 = T.BumpOK( 5 );
  d2 = T2.dbl;                 %capturar: assert evalua los args del mensaje SIEMPRE
  assert( d2 == 10 , 'test_cache:fireOK' , ...
          'el camino correcto de fire deberia dar dbl=10 y dio %g.' , d2 );
  fprintf( 'OK  fire      valida eventos; el camino correcto funciona\n' );

  %% -------- 4) camino de extension: subclase con subsref propio -> obj.access
  P = pSub();
  pd = P.dbl;
  assert( pd == 14 , 'test_cache:subclassAccess' , ...
          'la subclase con subsref propio deberia dar 14 y dio %g.' , pd );
  fprintf( 'OK  extension subclase con subsref propio llama al motor de la base\n' );

  %% -------- 5) indexar dentro de un valor que es objeto con subsref propio
  I = pInner();
  it = I.inner.total;
  assert( it == 0 , 'test_cache:innerSubsref' , ...
          'indexar dentro de un valor-objeto deberia dar 0 y dio %g.' , it );
  fprintf( 'OK  anidado   .miCP.total respeta el subsref del valor\n' );

  %% -------- 6) techo del log: 200 edits incrementales sin leer -> colapsa
  R = route().Demo();  v = R.len;                                           %#ok<NASGU>
  for i = 1:200, R = R.AddPoint( i , 0 ); end
  t = evalc( 'disp( R.CP )' );
  assert( contains( t , 'pendiente: 1 edit' ) , 'test_cache:logCap' , ...
          'tras 200 edits el log deberia ser 1 marcador; la tabla dice: %s' , strtrim(t) );
  rl = R.len;  rf = R.len_;    %len_ FUERZA recalculo: una sola vez, no en el mensaje
  assert( rl == rf , 'test_cache:logCapValue' , ...
          'el valor tras colapsar el log no coincide con el recalculo: %g vs %g.' , rl , rf );
  fprintf( 'OK  MAXLOG    200 edits -> 1 marcador, y el valor sigue correcto\n' );

  %% -------- 7) COW + replay + aislamiento entre copias (el nucleo)
  A = cart();  A.DEBUG = false;  A = A.Demo();  A = A.Add( 'cafe' , 3 , 2 );
  v = A.total;  B = A;  A = A.Add( 'pan' , 1 , 1 );                         %#ok<NASGU>
  at = A.total;  bt = B.total;
  assert( at == 7 && bt == 6 , 'test_cache:cowReplay' , ...
          'COW roto: A=%g (esperado 7), hermana B=%g (esperado 6).' , at , bt );
  A = A.SetPrice( 'cafe' , 4 );
  at = A.total;
  assert( at == 9 , 'test_cache:syncAbs' , ...
          'el sync absoluto deberia dar 9 y dio %g.' , at );
  fprintf( 'OK  nucleo    incremental 6->7, hermana intacta en 6, sync 9\n' );

  %% -------- 8) guarda de ciclos
  Z = cart();  Z.DEBUG = false;  Z = Z.Demo();  Z = Z.Define( 'bucle' , @(c) c.bucle );
  expectError( @() Z.bucle , 'cacheHandle:cycle' );
  fprintf( 'OK  ciclos    una CP que se pide a si misma corta con error\n' );

  %% -------- 9) proxy huerfano
  Q = cacheProxy();
  expectError( @() Q.total , 'cacheProxy:orphan' );
  fprintf( 'OK  proxy     huerfano da error con nombre\n' );

  %% -------- 10) save/load: DEFS viaja, CACHE revive, el log pendiente no sobrevive
  fn = fullfile( tempdir , 'test_cache_sl.mat' );
  cleaner = onCleanup( @() delete( fn ) );                                  %#ok<NASGU>
  S = cart();  S.DEBUG = false;  S = S.Demo();  S = S.Add( 'cafe' , 3 , 2 );
  v = S.total;  S = S.Add( 'pan' , 1 , 1 );          % total queda PENDIENTE %#ok<NASGU>
  W = S;                                             % copia que comparte handle
  save( fn , 'S' , 'W' );  clear S W
  L = load( fn );  S = L.S;  W = L.W;  S.DEBUG = true;
  t1 = evalc( 'a = S.total;' );  t2 = evalc( 'b = S.total;' );
  assert( a == 7 && b == 7 , 'test_cache:slValue' , ...
          'tras load el total deberia ser 7 y dio %g y %g.' , a , b );
  assert( contains( t1 , 'MISS' ) && contains( t2 , 'HIT' ) , 'test_cache:slRevive' , ...
          'tras load se esperaba MISS y luego HIT; trazas: [%s] [%s]' , strtrim(t1) , strtrim(t2) );
  wt = W.total;
  assert( wt == 7 , 'test_cache:slSibling' , ...
          'la copia que compartia handle deberia dar 7 y dio %g.' , wt );
  %el atajo INCREMENTAL sigue vivo tras el viaje: ejercita DEFS.inc/.incH, que
  %ahora se precomputan en Define y viajan en el .mat
  t0 = evalc( 'S = S.Add( ''te'' , 2 , 1 );' );      %DEBUG esta on: no ensuciar %#ok<NASGU>
  t3 = evalc( 'c = S.total;' );
  assert( c == 9 , 'test_cache:slIncremental' , ...
          'tras load el incremental deberia dar 9 y dio %g.' , c );
  assert( contains( t3 , 'incremental' ) , 'test_cache:slIncrementalPath' , ...
          'tras load el replay deberia usar el atajo incremental; traza: %s' , strtrim(t3) );
  %EVOLUCION: un .mat guardado antes de RENOMBRAR un evento trae CPs que
  %declaran el nombre viejo -> parecerian insensibles y servirian dato rancio
  %en silencio. loadobj avisa; aqui se llama al detector con inputs fabricados
  %(no se puede tener dos versiones de una clase en una sesion).
  defs = struct( 'val' , struct( 'compute' , @(o)1 , 'events' , struct('viejo',[]) , ...
                                 'inc' , '' , 'incH' , [] ) );
  lastwarn( '' );
  %evalc y NO warning('off'): un warning desactivado ni salta ni toca lastwarn
  evalc( 'st = cachedOwner.warnStaleEvents( ''pFake'' , defs , { ''nuevo'' } );' );
  [ ~ , wid ] = lastwarn;
  assert( isscalar( st ) && strcmp( wid , 'cachedOwner:staleEvents' ) , ...
          'test_cache:staleEvents' , ...
          'el detector deberia avisar de 1 CP rancia y dio %d avisos (id "%s").' , numel(st) , wid );
  st = cachedOwner.warnStaleEvents( 'pFake' , defs , { 'viejo' } );
  assert( isempty( st ) , 'test_cache:staleEventsClean' , ...
          'sin renombrado no deberia avisar y devolvio %d.' , numel(st) );
  fprintf( 'OK  save/load definiciones viajan, cache revive (MISS->HIT), incremental sigue, rename avisa\n' );

  %% -------- 11) la cache no engorda el .mat
  G = cart();  G.DEBUG = false;  G = G.Demo();  G = G.Add( 'x' , 1 , 1 );
  G = G.CP.total.set( zeros( 400 , 400 ) );          % 1.28 MB sembrados
  save( fn , 'G' );  d = dir( fn );
  assert( d.bytes < 50e3 , 'test_cache:matSize' , '.mat de %.0f KB: la cache viaja?!' , d.bytes/1024 );
  fprintf( 'OK  peso .mat %.1f KB con 1.28 MB cacheados: la cache no viaja\n' , d.bytes/1024 );

  %% -------- 12) RemoveCP revive el handle muerto (si no, cache apagada para siempre)
  M = cart();  M.DEBUG = false;  M = M.Demo();
  ws = warning( 'off' , 'MATLAB:structOnObject' );
  sm = struct( M );  delete( sm.CACHE );             % matar el handle a mano
  warning( ws );
  M = M.RemoveCP( 'ticket' );  M.DEBUG = true;
  m1 = evalc( 'x1 = M.total;' );  m2 = evalc( 'x2 = M.total;' );
  assert( x1 == 0 && x2 == 0 , 'test_cache:reviveValue' , ...
          'tras revivir el handle las lecturas deberian dar 0 y dieron %g y %g.' , x1 , x2 );
  assert( contains( m1 , 'MISS' ) && contains( m2 , 'HIT' ) , 'test_cache:revive' , ...
          'se esperaba MISS y luego HIT; trazas: [%s] [%s]' , strtrim(m1) , strtrim(m2) );
  fprintf( 'OK  RemoveCP  con handle muerto revive la cache (MISS->HIT), no la apaga\n' );

  %% -------- 13) aviso al cachear un objeto HANDLE (el COW no aisla ahi)
  %clave UNICA: el aviso solo salta la primera vez por clase.clave (memo persistente)
  [ ~ , u ] = fileparts( tempname );
  nm = matlab.lang.makeValidName( [ 'hbox' u ] );
  Hc = cart();  Hc.DEBUG = false;
  Hc = Hc.Define( nm , @(c) containers.Map( 'n' , 1 ) );
  lastwarn( '' );
  t = evalc( 'hv = Hc.( nm );' );                                         %#ok<NASGU>
  [ ~ , wid ] = lastwarn;
  assert( strcmp( wid , 'cachedOwner:handleValue' ) , 'test_cache:handleWarn' , ...
          'cachear un handle deberia avisar y llego "%s".' , wid );
  t2 = evalc( 'hv2 = Hc.( nm );' );                                       %#ok<NASGU>
  fprintf( 'OK  handle    avisa al cachear un objeto handle (aislamiento no garantizado)\n' );

  %% -------- 14) una EDICION no puede escalar con el nº de CPs cacheadas
  %fire clona el store y solo toca las claves AFECTADAS, asi que una CP
  %insensible y cacheada debe costar ~0 por edicion. Es el unico test de
  %RENDIMIENTO de la bateria y existe porque los goldens no lo verian: el
  %escalado puede volver en cualquier refactor sin cambiar una sola salida.
  %Guarda GRUESA y por RATIO (nunca en us) para que valga en otra maquina.
  %El umbral se calibro CON LA REGRESION PUESTA, no a ojo: en R2022a da 2.6x con
  %el fire actual y 6.3x devolviendo fire al que reconstruia el store campo a
  %campo. 4x parte esa horquilla, con ~1.5x de margen a cada lado. Si algun dia
  %falla por ruido y no por un cambio real, sube el umbral -- pero comprueba
  %antes que sigue cazando la version vieja.
  %TRES RONDAS INTERCALADAS Y LA MEDIANA, y no es adorno: de una sola toma este
  %ratio da entre 1.0x y 4.2x SOBRE EL MISMO CODIGO (medido), porque las dos
  %medidas cargan ~700 us de despacho fijo y el denominador -- decimas de ms --
  %se lleva todo el ruido de la maquina. El test saltaba por eso, no por una
  %regresion. Intercalar A/B en la misma sesion y quedarse con la mediana es lo
  %que funciona aqui.
  %RECALIBRADO (16-ago-2026) reintroduciendo la regresion ENTERA -- reconstruir
  %el store campo a campo Y leer obj.DEFS por clave --: da 12.3 ms con 64 CPs
  %(el "12 ms" de la auditoria que la encontro) y mediana 12.4x, contra 2.4x del
  %codigo actual. El limite de 4x parte esa horquilla con holgura por los dos
  %lados. OJO A LA SENSIBILIDAD: solo la mitad de aquella regresion (el rebuild,
  %sin tocar el hoisting) mide 3.8x y se colaria por debajo del limite; para
  %cazar tambien esa mitad habria que medir `fire` amortizado (N disparos dentro
  %de un metodo, sin el peaje del subsref), que es mas sensible pero ya no mide
  %el camino real.
  N = cart();  N.DEBUG = false;  N = N.Demo();
  w = N.total;  w = N.ticket;  w = N.topPrice;  w = N.currency;           %#ok<NASGU>
  M = cart();  M.DEBUG = false;  M = M.Demo();        % objeto APARTE: si creciera
  w = M.total;  w = M.ticket;  w = M.topPrice;  w = M.currency;           %#ok<NASGU>
  for i = 1:60, M = M.Define( sprintf( 'k%02d' , i ) , @(o) i ); end      % el mismo, las
  for i = 1:60, M.( sprintf( 'k%02d' , i ) ); end     % 60 claves nuevas caerian en el
                                                      % store COMPARTIDO y ensuciarian
                                                      % tambien la medida de 4 CPs
  rr = zeros( 1 , 3 );
  for j = 1:3
    tFew  = timeit( @() N.Add( 'x' , 1 , 1 ) );       % 4 CPs cacheadas
    tMany = timeit( @() M.Add( 'x' , 1 , 1 ) );       % 64 CPs cacheadas
    rr(j) = tMany / tFew;
  end
  rat = median( rr );
  assert( rat < 4 , 'test_cache:fireScaling' , ...
          [ 'una edicion con 64 CPs cacheadas cuesta %.1fx la de 4 (mediana de 3 rondas, ' ...
            'limite 4x; rondas: %s): el bucle de fire ha vuelto a escalar con el numero ' ...
            'de claves.' ] , rat , mat2str( round( rr , 2 ) ) );
  fprintf( 'OK  fire      edicion con 64 CPs cacheadas = %.1fx la de 4 (mediana de 3, limite 4x)\n' , rat );

  %% -------- 15) la clase es ESCALAR y se obliga (no hay cuarta puerta)
  S1 = cart();  S1.DEBUG = false;  S1 = S1.Demo();
  expectError( @() [ S1 , S1 ]            , 'cachedOwner:notScalar' );  % horzcat
  expectError( @() [ S1 ; S1 ]            , 'cachedOwner:notScalar' );  % vertcat
  expectError( @() cat( 2 , S1 , S1 )     , 'cachedOwner:notScalar' );  % cat
  expectError( @() repmat( S1 , 1 , 3 )   , 'cachedOwner:notScalar' );  % via ()
  expectError( @() S1( ones(1,3) )        , 'cachedOwner:notScalar' );  % paren-ref
  expectError( @() growIt( S1 )           , 'cachedOwner:notScalar' );  % A(3) = obj
  expectError( @() S1( 2 )                , 'cachedOwner:notScalar' );
  %el no-op SI pasa, y arrastra el resto de la cadena
  s1t = S1(1).total;
  assert( s1t == 0 , 'test_cache:scalarNoop' , ...
          'o(1).total deberia seguir la cadena y dar 0; dio %g.' , s1t );
  assert( isa( S1(1) , 'cart' ) , 'test_cache:scalarNoop2' );
  %`empty` es estatico y sobrevive a proposito (MATLAB lo usa por dentro): lo que
  %se exige es que el 0x0 sea INERTE, no que no exista
  S0 = cart.empty;
  assert( isempty( S0 ) , 'test_cache:emptyExists' , ...
          'cart.empty deberia existir y ser vacia; size %s.' , mat2str(size(S0)) );
  expectError( @() S0.total , 'cachedOwner:notScalar' );
  expectError( @() S0(1)    , 'cachedOwner:notScalar' );
  %el no-op tambien tiene que tragar un STATEMENT del plano de control: la
  %version que recursaba exigia 1 salida y moria en unassignedOutputs
  v0 = S1.total;                                                          %#ok<NASGU>
  evalc( 'S1(1).CP.total.delete' );
  t = evalc( 'disp( S1.CP )' );
  assert( contains( t , 'sin calcular' ) , 'test_cache:noopStmt' , ...
          'o(1).CP.x.delete deberia borrar el valor; la tabla dice: %s' , strtrim(t) );
  %y el proxy es escalar tambien: un array de proxies daba maxrhs criptico
  Pp = S1.CP;  Qq = [ Pp , Pp ];
  expectError( @() Qq.total , 'cacheProxy:notScalar' );
  fprintf( 'OK  escalar   concat/repmat/A(3)=o/o(k) bloqueados; o(1) no-op (statements incluidos); proxy escalar\n' );

  %% -------- 16) escribir una CP da un error que explica que hacer
  expectError( @() setCP( S1 ) , 'cachedOwner:readOnlyCP' );
  S1.DEBUG = true;  S1.DEBUG = false;        % las propiedades DE VERDAD se escriben
  fprintf( 'OK  readOnly  o.total = x remite a CP.set/Define; las propiedades normales OK\n' );

  %% -------- 17) la TABLA marca los dos fallos que dan dato malo
  %un warning se pierde (-batch, script largo, warning('off') ajeno); la tabla es
  %el sitio persistente donde alguien mira el estado, asi que el aviso vive ahi
  %tambien. Requisito extra: una tabla LIMPIA no cambia ni un byte -- de eso ya
  %responden los goldens del grupo 1, que la imprimen cuatro veces.
  N1 = cart();  N1.DEBUG = false;  N1 = N1.Demo();
  tc = evalc( 'disp( N1.CP )' );
  assert( ~contains( tc , 'AVISO' ) && ~contains( tc , '! ' ) , 'test_cache:tableClean' , ...
          'la tabla sin problemas no debe llevar marcador; salio: %s' , tc );
  wsx = warning( 'off' , 'cachedOwner:handleValue' );
  N1 = N1.CP.total.set( containers.Map( 'n' , 1 ) );      % valor HANDLE
  warning( wsx );
  th = evalc( 'disp( N1.CP )' );
  assert( contains( th , '! total' ) && contains( th , 'CON AVISO' ) && ...
          contains( th , 'containers.Map' ) , 'test_cache:tableHandle' , ...
          'la tabla deberia marcar el valor handle; salio: %s' , th );
  %la rama del evento huerfano, con inputs fabricados (ver cpNote)
  s1 = cachedOwner.cpNote( struct('viejo',[]) , { 'nuevo' } , false , [] );
  s2 = cachedOwner.cpNote( struct('viejo',[]) , { 'viejo' } , false , [] );
  s3 = cachedOwner.cpNote( struct()           , { 'x' }     , true  , containers.Map );
  assert( contains( s1 , 'YA NO ESTAN' ) && isempty( s2 ) && contains( s3 , 'handle' ) , ...
          'test_cache:cpNote' , 'cpNote: huerfano="%s" limpio="%s" handle="%s"' , s1 , s2 , s3 );
  fprintf( 'OK  tabla     marca valor handle y evento huerfano; la tabla limpia no cambia\n' );

  %% -------- 18) un handler INCREMENTAL no puede leer el dueno
  %el log guarda EDICIONES, no estados: en el replay de N edits el handler
  %recibe el objeto de AHORA, asi que uno que lo lea da un valor MAL -- y en
  %silencio, porque con UNA sola edicion pendiente acierta (comprobado antes de
  %la guarda: media incremental de [2 4 6] = 3.78 en vez de 4). Tres capas:
  %clasificador, guarda de Define (anonimo) y red del replay (opaco / .mat viejo).
  %--- el clasificador, caso a caso
  assert( strcmp( cachedOwner.ownerUse( @(v,o,d) v + o.X ) , 'yes' ) , ...
          'test_cache:ownerUseYes' , 'leer o.X deberia dar "yes".' );
  assert( strcmp( cachedOwner.ownerUse( @(v,~,d) v + d ) , 'no' ) , ...
          'test_cache:ownerUseTilde' , 'el argumento ~ deberia dar "no".' );
  assert( strcmp( cachedOwner.ownerUse( @(v,o,d) v + d ) , 'no' ) , ...
          'test_cache:ownerUseUnused' , 'no nombrarlo deberia dar "no".' );
  assert( strcmp( cachedOwner.ownerUse( @(v,o,d) v + numel('a o b') ) , 'no' ) , ...
          'test_cache:ownerUseLiteral' , 'un "o" dentro de un texto no es un uso.' );
  assert( strcmp( cachedOwner.ownerUse( @(v,o,d) v + d.o + otro(d) ) , 'no' ) , ...
          'test_cache:ownerUseToken' , 'ni d.o (campo) ni otro (otra palabra) son usos.' );
  assert( strcmp( cachedOwner.ownerUse( @pOwnerInc ) , 'unknown' ) , ...
          'test_cache:ownerUseOpaque' , 'un handler no anonimo no se puede inspeccionar.' );
  %--- Define rechaza el anonimo que lo lee, y admite la forma correcta
  R1 = route().Demo();
  expectError( @() R1.Define( 'avgX' , @(r) mean( r.XY(:,1) ) , ...
                              'addpoint' , @(v,r,d) ( v*(size(r.XY,1)-1) + d ) / size(r.XY,1) ) , ...
               'cachedOwner:incrementalUsesOwner' );
  R1 = R1.Define( 'nx' , @(r) size( r.XY ,1) , 'addpoint' , @(v,~,d) v + 1 );
  R1 = R1.AddPoint(0,0);  q0 = R1.nx;                                       %#ok<NASGU>
  R1 = R1.AddPoint(3,0);  R1 = R1.AddPoint(7,0);
  q1 = R1.nx;  q2 = R1.nx_;
  assert( q1 == 3 && q1 == q2 , 'test_cache:tildeIncremental' , ...
          '@(v,~,d) debe seguir siendo incremental y valer 3: replay %g, recompute %g.' , q1 , q2 );
  %--- la RED: handler opaco. Clave unica -> el aviso (memo persistente) salta
  %aunque la suite se corra dos veces en la misma sesion
  [ ~ , u ] = fileparts( tempname );
  cp = matlab.lang.makeValidName( [ 'seg' u ] );
  cpu = [ cp '_' ];
  R2 = route().Demo();
  R2 = R2.Define( cp , @(r) mean( vecnorm( diff( r.XY ) ,2,2) ) , 'addpoint' , @pOwnerInc );
  R2 = R2.AddPoint(0,0);  R2 = R2.AddPoint(3,0);
  g0 = R2.( cp );                            % poblarla (un tramo de 3)         %#ok<NASGU>
  R2.DEBUG = true;                           % a partir de aqui, todo por evalc
  e1 = evalc( 'R2 = R2.AddPoint(7,0); g1 = R2.( cp );' );   % UNA pendiente: el dueno vale
  evalc( 'g1r = R2.( cpu );' );
  assert( g1 == g1r && contains( e1 , 'incremental' ) , 'test_cache:oneEditKeepsShortcut' , ...
          'con 1 edicion el atajo debe usarse y acertar: replay %g, recompute %g.' , g1 , g1r );
  evalc( 'R2 = R2.AddPoint(11,0); R2 = R2.AddPoint(16,0);' );   % VARIAS pendientes
  lastwarn( '' );
  e2 = evalc( 'g2 = R2.( cp );' );
  [ ~ , wid ] = lastwarn;
  evalc( 'g2r = R2.( cpu );' );
  assert( g2 == g2r , 'test_cache:manyEditsDropShortcut' , ...
          'con varias ediciones el atajo debe descartarse: replay %g, recompute %g.' , g2 , g2r );
  assert( strcmp( wid , 'cachedOwner:incrementalOpaque' ) , 'test_cache:opaqueWarn' , ...
          'se esperaba el aviso del atajo descartado y llego "%s".' , wid );
  assert( contains( e2 , 'descartado' ) , 'test_cache:opaqueTrace' , ...
          'la traza deberia decir que descarta el atajo; dijo: %s' , strtrim(e2) );
  %--- la otra rama del aviso (handler que YA venia leyendo el dueno en un .mat)
  lastwarn( '' );
  evalc( 'cachedOwner.warnIncOwner( ''pFake'' , cp , ''yes'' , 3 );' );
  [ ~ , wid2 ] = lastwarn;
  assert( strcmp( wid2 , 'cachedOwner:incrementalUsesOwner' ) , 'test_cache:matWarn' , ...
          'la rama de .mat viejo deberia avisar y llego "%s".' , wid2 );
  fprintf( 'OK  dueno     Define rechaza el incremental que lee el objeto; con handler opaco el replay descarta el atajo\n' );

  %% -------- 19) CPs COMPUESTAS: una CP que sale de otra
  %antes, la compuesta no declaraba ningun evento -> el motor la daba por
  %insensible y SOBREVIVIA sirviendo un valor sacado de datos que ya no valian
  %(comprobado: base pasaba de 3 a 13 y la que salia de ella seguia diciendo 6).
  %Ahora la lectura queda apuntada como dependencia y la compuesta cae con ella.
  %--- desde FUERA de la clase: el caso original
  K = pComp();
  b0 = K.base;  d0 = K.viaAccess;                                          %#ok<NASGU>
  K2 = K.Define( 'outside' , @(o) o.base * 100 );     % compute definido AQUI
  o0 = K2.outside;
  assert( o0 == 1000 , 'test_cache:compValue' , 'la compuesta deberia valer 1000 y dio %g.' , o0 );
  t = evalc( 'disp( K2.CP )' );
  assert( contains( t , 'sale de: base' ) , 'test_cache:compTable' , ...
          'la tabla deberia decir de que sale la CP compuesta; dijo: %s' , t );
  K2 = K2.Bump();                                     % 'chg' invalida base
  o1 = K2.outside;  o1r = K2.outside_;
  assert( o1 == 2000 && o1 == o1r , 'test_cache:compStale' , ...
          'tras el evento la compuesta deberia valer 2000: dio %g (recompute %g).' , o1 , o1r );
  %--- TRANSITIVIDAD con el intermedio ya cacheado (el caso que obliga a
  %arrastrar las dependencias de la dependencia: sin eso, A->B->base con B
  %cacheado dejaba A rancia)
  K3 = pComp();
  K3 = K3.Define( 'mid' , @(o) o.base + 5 );
  K3 = K3.Define( 'top' , @(o) o.mid  * 2 );
  m0 = K3.mid;                                        % se cachea ANTES que top
  p0 = K3.top;
  assert( p0 == 30 , 'test_cache:chainValue' , 'la cadena deberia dar 30 y dio %g.' , p0 );
  K3 = K3.Bump();
  p1 = K3.top;  p1r = K3.top_;
  assert( p1 == 50 && p1 == p1r , 'test_cache:chainStale' , ...
          'la cadena entera deberia caer: top dio %g (recompute %g).' , p1 , p1r );
  %--- una dependencia rota GANA a la declaracion del handler
  K4 = pComp();  s0 = K4.stubborn;
  assert( s0 == 13 , 'test_cache:stubbornValue' , 'stubborn deberia valer 13 y dio %g.' , s0 );
  K4 = K4.Bump();
  s1 = K4.stubborn;
  assert( s1 == 23 , 'test_cache:stubbornDep' , ...
          [ 'stubborn declara un handler que se queda con el valor viejo, pero SALE de base: ' ...
            'deberia caer y dar 23, y dio %g.' ] , s1 );
  %--- DENTRO de la clase, el punto no resuelve la CP: error traducido
  expectError( @() K.dotInside , 'cachedOwner:cpFromInside' );
  %--- los dos gestos que cambian un valor SIN evento tiran a los dependientes
  K5 = pComp();  v5 = K5.base;  a5 = K5.viaAccess;                         %#ok<NASGU>
  K5.DEBUG = true;
  r5 = evalc( 'z = K5.base_;' );                                           %#ok<NASGU>
  assert( contains( r5 , 'tiradas las CPs' ) , 'test_cache:recompDrops' , ...
          'el recalculo forzado deberia tirar a viaAccess; traza: %s' , strtrim(r5) );
  K5.DEBUG = false;
  K6 = pComp();  v6 = K6.base;  a6 = K6.viaAccess;                         %#ok<NASGU>
  K6 = K6.CP.base.set( 100 );                         % sembrado a mano
  a6b = K6.viaAccess;
  assert( a6b == 102 , 'test_cache:seedDrops' , ...
          'tras sembrar base=100 la compuesta deberia valer 102 y dio %g.' , a6b );
  %--- las deps SOBREVIVEN al replay: compuesta CON handler incremental. El
  %refresco guardaba takeDeps() a secas -- y un handler incremental no lee
  %nada --, asi que la dependencia de base se perdia en el PRIMER replay y a
  %la siguiente caida de base la compuesta sobrevivia sirviendo rancio en
  %silencio (comprobado: servia 305 con la verdad en 905).
  K7 = pIncDep();
  c0 = K7.comp;
  assert( c0 == 300 , 'test_cache:incDepSeed' , 'comp deberia nacer en 300 y dio %g.' , c0 );
  K7 = K7.Bump( 5 );                                  % solo toca B: base sobrevive
  c1 = K7.comp;                                       % replay INCREMENTAL
  assert( c1 == 305 , 'test_cache:incDepReplay' , 'el replay deberia dar 305 y dio %g.' , c1 );
  t19 = evalc( 'disp( K7.CP )' );
  assert( contains( t19 , 'sale de: base' ) , 'test_cache:incDepKept' , ...
          'tras el replay la tabla deberia seguir diciendo "sale de: base"; dijo: %s' , strtrim(t19) );
  K7 = K7.ChgA( 9 );                                  % cae base -> comp cae CON ella
  c2 = K7.comp;  c2r = K7.comp_;
  assert( c2 == 905 && c2 == c2r , 'test_cache:incDepDrop' , ...
          'tras caer base la compuesta deberia dar 905: replay %g, recompute %g.' , c2 , c2r );
  %--- la TRAZA no puede contradecirse. Una compuesta pasa por la 1a pasada de
  %fire como INSENSIBLE (no declara eventos) y se anota como superviviente; si
  %la 2a pasada se la lleva y nadie la quita de esa lista, la traza la saca en
  %las DOS a la vez ("caen {x(compuesta)} | sobreviven {x}"). Lo destapo el demo
  %de collatz, que es el primer ejemplo con CPs compuestas de verdad.
  K8 = pComp();  v8 = K8.base;  a8 = K8.viaAccess;                         %#ok<NASGU>
  K8.DEBUG = true;
  t8 = evalc( 'K8 = K8.Bump();' );
  K8.DEBUG = false;
  cae  = regexp( t8 , 'caen \{([^}]*)\}'       , 'tokens' , 'once' );
  sobr = regexp( t8 , 'sobreviven \{([^}]*)\}' , 'tokens' , 'once' );
  assert( contains( cae{1} , 'viaAccess' ) && ~contains( sobr{1} , 'viaAccess' ) , ...
          'test_cache:traceContradiction' , ...
          [ 'la traza se contradice: viaAccess cae por dependencia pero sigue en la lista de ' ...
            'supervivientes. caen={%s} sobreviven={%s}' ] , cae{1} , sobr{1} );
  %--- REDEFINIR A tambien arrastra a B. Los eventos, el recalculo forzado y el
  %sembrado ya lo hacian; Define y RemoveCP no, y eran la misma clase de
  %agujero: cambiar lo que A significa sin avisar a quien salio de A.
  K9 = pComp();  b9 = K9.base;  v9 = K9.viaAccess;                         %#ok<NASGU>
  assert( v9 == 12 , 'test_cache:redefSeed' , 'viaAccess deberia nacer en 12 y dio %g.' , v9 );
  K9 = K9.Define( 'base' , @(o) 999 , 'chg' , [] );
  r9 = K9.viaAccess;
  assert( r9 == 1001 , 'test_cache:redefDrops' , ...
          [ 'redefinir base deberia tirar a viaAccess (que sale de ella): esperado 1001, dio %g. ' ...
            'Si sale 12, el valor viejo ha sobrevivido a un compute que ya no existe.' ] , r9 );
  %--- RemoveCP(A) tambien: antes servia el valor viejo Y reventaba al forzarlo
  KA = pComp();  bA = KA.base;  vA = KA.viaAccess;                         %#ok<NASGU>
  KA = KA.RemoveCP( 'base' );
  %ahora CAE (antes servia el valor viejo) y al recalcularse dice con claridad
  %que le falta la CP de la que salia, en vez del "Unrecognized field name" crudo
  expectError( @() KA.viaAccess , 'cachedOwner:noCP' );
  fprintf( 'OK  compuesta CP que sale de otra: cae con ella (evento, recompute, set, Define y RemoveCP); la traza no se contradice\n' );

  %% -------- 20) el OLVIDO TOTAL de eventos se avisa (guarda barata)
  %una CP que lee datos del objeto y no declara ni un evento sobrevive a TODO y
  %sirve rancio en silencio. Es lo unico juzgable sin un mapa propiedad->evento,
  %asi que es un aviso (heuristico) y no un error. Claves UNICAS: el aviso es
  %una vez por clase.clave y la suite puede correrse dos veces en una sesion.
  [ ~ , u2 ] = fileparts( tempname );
  q = @(k) matlab.lang.makeValidName( sprintf( 'ne%s%d' , u2 , k ) );
  E = cart();  E.DEBUG = false;  E = E.Demo();
  lastwarn( '' );
  ev1 = evalc( 'E1 = E.Define( q(1) , @(c) mean( c.PRICE ) );' );          %#ok<NASGU>
  [ ~ , w1 ] = lastwarn;
  assert( strcmp( w1 , 'cachedOwner:noEvents' ) , 'test_cache:noEventsWarn' , ...
          'una CP que lee datos sin declarar eventos deberia avisar y llego "%s".' , w1 );
  te = evalc( 'disp( E1.CP )' );
  assert( contains( te , 'no declara ningun evento' ) && contains( te , 'CON AVISO' ) , ...
          'test_cache:noEventsTable' , 'la tabla deberia marcarlo; salio: %s' , te );
  %y los tres casos en los que NO debe saltar
  lastwarn( '' );
  evalc( 'E.Define( q(2) , @(c) 42 );' );                     % no lee nada
  evalc( 'E.Define( q(3) , @(c) c.total * 2 );' );            % lee otra CP, no datos
  evalc( 'E.Define( q(4) , @(c) mean( c.PRICE ) , ''changeprice'' , [] );' );  % declara
  [ ~ , w2 ] = lastwarn;
  assert( isempty( w2 ) , 'test_cache:noEventsFalse' , ...
          'no deberia avisar de una constante, de leer otra CP ni de declarar evento; llego "%s".' , w2 );
  %el analizador, caso a caso
  assert( isequal( cachedOwner.propsRead( @(c) mean(c.PRICE) , {'PRICE','QTY'} ) , {'PRICE'} ) , ...
          'test_cache:propsReadHit' , 'deberia ver c.PRICE.' );
  assert( isempty( cachedOwner.propsRead( @(c) numel('c.PRICE') , {'PRICE'} ) ) , ...
          'test_cache:propsReadLiteral' , 'un c.PRICE dentro de un texto no es una lectura.' );
  assert( isempty( cachedOwner.propsRead( @(c) s.PRICE , {'PRICE'} ) ) , ...
          'test_cache:propsReadOther' , 's.PRICE no es del objeto.' );
  assert( isempty( cachedOwner.propsRead( @pOwnerInc , {'PRICE'} ) ) , ...
          'test_cache:propsReadOpaque' , 'un handler no anonimo no se puede inspeccionar.' );
  fprintf( 'OK  noEvents  avisa (y marca en la tabla) la CP que lee datos y no declara ningun evento\n' );

  %% -------- 21) un ATAJO ROTO deja de ser invisible
  %la cascada se traga el error del handler A PROPOSITO (que falle no puede
  %romper una lectura), y ese era el precio: la clase seguia dando el valor
  %correcto, por el camino caro, PARA SIEMPRE, sin una sola senal. Ahora avisa
  %una vez, lleva la cuenta y lo marca en la tabla. Claves unicas: el memo de
  %avisos vive en la sesion.
  [ ~ , u3 ] = fileparts( tempname );
  bn  = matlab.lang.makeValidName( [ 'brk' u3 ] );
  B1  = cart();  B1.DEBUG = false;  B1 = B1.Demo();  B1 = B1.Add( 'a' , 1 , 1 );
  B1  = B1.Define( bn , @(c) 42 , 'additem' , @(v,~,d) error('boom:bad','handler roto') );
  x0  = B1.( bn );                                                        %#ok<NASGU>
  B1  = B1.Add( 'b' , 2 , 1 );
  lastwarn( '' );
  tb1 = evalc( 'x1 = B1.( bn );' );                                       %#ok<NASGU>
  [ ~ , wb1 ] = lastwarn;
  assert( x1 == 42 , 'test_cache:brokenValue' , ...
          'el dato debe seguir bien pese al handler roto: dio %g.' , x1 );
  assert( strcmp( wb1 , 'cachedOwner:handlerFailed' ) , 'test_cache:brokenWarn' , ...
          'un handler que revienta deberia avisar y llego "%s".' , wb1 );
  %la SEGUNDA vez no repite el aviso, pero la cuenta sube y la tabla lo dice
  B1 = B1.Add( 'c' , 3 , 1 );
  lastwarn( '' );
  evalc( 'x2 = B1.( bn );' );
  [ ~ , wb2 ] = lastwarn;
  assert( isempty( wb2 ) , 'test_cache:brokenOnce' , ...
          'el aviso es una vez por CP y evento; se repitio con "%s".' , wb2 );
  tt = evalc( 'disp( B1.CP )' );
  assert( contains( tt , 'ATAJO ROTO' ) && contains( tt , 'additem (x2)' ) , ...
          'test_cache:brokenTable' , 'la tabla deberia marcar el atajo roto y su cuenta; salio: %s' , tt );
  %redefinir la CP borra el historial (esa anotacion era de la definicion vieja)
  B1 = B1.Define( bn , @(c) 42 );
  tt2 = evalc( 'disp( B1.CP )' );
  assert( ~contains( tt2 , 'ATAJO ROTO' ) , 'test_cache:brokenForget' , ...
          'al redefinir deberia limpiarse la marca; salio: %s' , tt2 );
  %la rama del SYNC ABSOLUTO avisa igual
  bn2 = matlab.lang.makeValidName( [ 'brs' u3 ] );
  B2  = cart();  B2.DEBUG = false;  B2 = B2.Demo();  B2 = B2.Add( 'a' , 1 , 1 );
  B2  = B2.Define( bn2 , @(c) 7 , 'changeprice' , @(v,o) error('boom:abs','absoluto roto') );
  y0  = B2.( bn2 );                                                       %#ok<NASGU>
  B2  = B2.SetPrice( 'a' , 5 );
  lastwarn( '' );
  evalc( 'y1 = B2.( bn2 );' );
  [ ~ , wb3 ] = lastwarn;
  assert( y1 == 7 && strcmp( wb3 , 'cachedOwner:handlerFailed' ) , 'test_cache:brokenAbs' , ...
          'el sync absoluto roto deberia avisar: valor %g, aviso "%s".' , y1 , wb3 );
  fprintf( 'OK  atajoroto handler que revienta: avisa una vez, cuenta las veces y lo marca la tabla\n' );

  %% -------- 22) indexar dentro del valor NO trunca las listas separadas por comas
  %pedir siempre UNA salida se comia el resto en silencio: {L.items.n} devolvia
  %un elemento en vez de tres, y la MISMA expresion pasando el valor por una
  %variable devolvia tres. Ahora el numero lo decide el VALOR (se le pregunta),
  %tanto en la lectura directa como por el plano de control.
  L = pList();
  assert( numel( { L.items.n } ) == 3 , 'test_cache:csStruct' , ...
          '{L.items.n} deberia dar 3 y dio %d.' , numel( { L.items.n } ) );
  assert( numel( { L.cellv{:} } ) == 3 , 'test_cache:csCell' , ...
          '{L.cellv{:}} deberia dar 3 y dio %d.' , numel( { L.cellv{:} } ) );
  [ c1 , c2 , c3 ] = L.items.n;
  assert( isequal( [ c1 c2 c3 ] , [1 2 3] ) , 'test_cache:csMultiLHS' , ...
          '[a,b,c] = L.items.n deberia dar 1 2 3 y dio %g %g %g.' , c1 , c2 , c3 );
  assert( isequal( [ L.items.n ] , [1 2 3] ) , 'test_cache:csConcat' , ...
          '[L.items.n] deberia dar la lista entera.' );
  %el mismo resultado que pasando el valor por una variable (era la asimetria)
  iv = L.items;
  assert( isequal( { L.items.n } , { iv.n } ) , 'test_cache:csSame' , ...
          'con y sin variable intermedia deberia salir lo mismo.' );
  %una CP normal sigue dando UNA salida, y el statement sigue mostrando valor
  assert( isequal( L.plain , 43 ) , 'test_cache:csPlain' , 'una CP escalar deberia dar 43.' );
  tp = evalc( 'L.plain' );
  assert( contains( tp , '43' ) , 'test_cache:csStatement' , ...
          'el statement deberia mostrar el valor; salio: %s' , strtrim(tp) );
  %el plano de control, en sus dos formas
  assert( numel( { L.CP.items.n } ) == 3 , 'test_cache:csProxyChain' , ...
          '{L.CP.items.n} deberia dar 3 y dio %d.' , numel( { L.CP.items.n } ) );
  PL = L.CP;
  assert( numel( { PL.items.n } ) == 3 && numel( { PL.cellv{:} } ) == 3 , ...
          'test_cache:csProxyLoose' , 'con el proxy en una variable tambien deberian ser 3.' );
  %...sin romper las OPERACIONES del proxy, que no son indexacion
  v0 = L.plain;                                                           %#ok<NASGU>
  evalc( 'L.CP.plain.delete' );                       % statement: CERO salidas
  t22 = evalc( 'disp( L.CP )' );
  assert( contains( t22 , 'plain' ) && contains( t22 , 'sin calcular' ) , ...
          'test_cache:csProxyOps' , '.delete deberia seguir funcionando; tabla: %s' , t22 );
  L2 = L.CP.plain.removeCP;
  assert( ~contains( evalc( 'disp( L2.CP )' ) , 'plain' ) , 'test_cache:csProxyRemove' , ...
          '.removeCP deberia seguir devolviendo el dueno nuevo.' );
  %y el recalculo forzado tampoco trunca (cuenta con el valor que ya hay)
  assert( numel( { L.items_.n } ) == 3 , 'test_cache:csForced' , ...
          '{L.items_.n} deberia dar 3 y dio %d.' , numel( { L.items_.n } ) );
  %tras una edicion, la lista crece con el valor
  L3 = L.Grow();
  assert( numel( { L3.items.n } ) == 4 , 'test_cache:csAfterEdit' , ...
          'tras crecer deberia dar 4 y dio %d.' , numel( { L3.items.n } ) );
  fprintf( 'OK  cs-list   indexar dentro del valor devuelve la lista entera (directo, proxy y forzado)\n' );

  %% -------- 23) decide el evento MAS ESPECIFICO, tambien frente al incremental
  %fire ya daba la ultima palabra al mas especifico para el [] (un [] especifico
  %no lo pisa el handler de un general), pero el ATAJO incremental se colaba solo
  %con estar en el lote: un incremental declarado en el evento GENERAL adelantaba
  %al absoluto del ESPECIFICO, que no corria nunca (daba 105 en vez de 999).
  Sp = pSpec();
  k0 = Sp.k;  m0 = Sp.m;                              % 100 y 100
  assert( k0 == 100 && m0 == 100 , 'test_cache:specSeed' , ...
          'las dos CPs deberian nacer en 100 y dieron %g y %g.' , k0 , m0 );
  Sp = Sp.Add( 5 );                                   % dispara ESPECIFICO + general
  k1 = Sp.k;  m1 = Sp.m;
  assert( k1 == 999 , 'test_cache:specWins' , ...
          [ 'con el especifico absoluto y el general incremental debe ganar el ESPECIFICO ' ...
            '(999) y dio %g.' ] , k1 );
  assert( m1 == 105 , 'test_cache:specInc' , ...
          'cuando el incremental ES el especifico debe usarse: esperado 105, dio %g.' , m1 );
  %y con el general a solas, cada una usa lo que declaro para el general
  S2 = pSpec();  v2 = S2.k;  S2 = S2.Bump();                              %#ok<NASGU>
  assert( S2.k == 101 , 'test_cache:specGenInc' , ...
          'el incremental del general debe correr cuando es el que gana: esperado 101, dio %g.' , S2.k );
  S3 = pSpec();  v3 = S3.m;  S3 = S3.Bump();                              %#ok<NASGU>
  assert( S3.m == 777 , 'test_cache:specGenAbs' , ...
          'el absoluto del general deberia dar 777 y dio %g.' , S3.m );
  %el atajo perdido no se cuela por la puerta de atras: dos ediciones seguidas
  %con el especifico ganando siguen dando el absoluto, no una cadena incremental
  S4 = pSpec();  v4 = S4.k;  S4 = S4.Add(1);  S4 = S4.Add(2);             %#ok<NASGU>
  assert( S4.k == 999 , 'test_cache:specChain' , ...
          'dos ediciones con el especifico ganando deberian dar 999 y dieron %g.' , S4.k );
  fprintf( 'OK  especifico el evento mas especifico decide tambien frente al atajo incremental\n' );

  %% -------- 24) VerifyCPs: caza al handler que MIENTE y al evento no declarado
  %son los dos fallos que ninguna guarda puede ver -- no fallan, sirven otra
  %cosa -- y los dos se manifiestan igual: la cache dice A y el compute dice B.
  V1 = cart();  V1.DEBUG = false;  V1 = V1.Demo();  V1 = V1.Add( 'cafe' , 3 , 2 );
  v = V1.total;  v = V1.ticket;  v = V1.topPrice;  v = V1.currency;       %#ok<NASGU>
  V1 = V1.Add( 'pan' , 1 , 1 );
  assert( isempty( V1.VerifyCPs() ) , 'test_cache:verifyClean' , ...
          'un cart sano no deberia tener ninguna CP rancia y VerifyCPs marco: %s' , ...
          strjoin( V1.VerifyCPs() , ', ' ) );
  %(1) un handler que MIENTE: no falla, devuelve otra cosa
  V2 = cart();  V2.DEBUG = false;  V2 = V2.Demo();  V2 = V2.Add( 'a' , 1 , 1 );
  V2 = V2.Define( 'liar' , @(c) sum( c.PRICE ) , 'additem' , @(v,~,d) v + 999 );
  v = V2.liar;  V2 = V2.Add( 'b' , 2 , 1 );                               %#ok<NASGU>
  [ bad2 , R2 ] = V2.VerifyCPs();
  assert( isequal( bad2 , { 'liar' } ) , 'test_cache:verifyLiar' , ...
          'deberia marcar solo "liar" y marco: %s' , strjoin( bad2 , ', ' ) );
  kk = find( strcmp( { R2.name } , 'liar' ) );
  assert( R2(kk).served == 1000 && R2(kk).truth == 3 , 'test_cache:verifyLiarVals' , ...
          'el informe deberia decir sirve 1000 / verdad 3 y dijo %g / %g.' , R2(kk).served , R2(kk).truth );
  %(2) el olvido PARCIAL de eventos: declara additem pero NO changeprice
  V3 = cart();  V3.DEBUG = false;  V3 = V3.Demo();  V3 = V3.Add( 'a' , 2 , 1 );
  V3 = V3.Define( 'avgP' , @(c) mean( c.PRICE ) , 'additem' , [] );
  v = V3.avgP;  V3 = V3.SetPrice( 'a' , 9 );                              %#ok<NASGU>
  assert( isequal( V3.VerifyCPs() , { 'avgP' } ) , 'test_cache:verifyForgot' , ...
          'deberia marcar la CP que no declaro changeprice.' );
  %NO TOCA LA CACHE: el pendiente sigue pendiente y se puede repetir
  V4 = cart();  V4.DEBUG = false;  V4 = V4.Demo();  V4 = V4.Add( 'a' , 1 , 1 );
  v = V4.total;  V4 = V4.Add( 'b' , 2 , 1 );                              %#ok<NASGU>
  t24a = evalc( 'disp( V4.CP )' );
  b1 = V4.VerifyCPs();  b2 = V4.VerifyCPs();
  t24b = evalc( 'disp( V4.CP )' );
  assert( isempty( b1 ) && isempty( b2 ) && isequal( t24a , t24b ) , 'test_cache:verifyPure' , ...
          'VerifyCPs no debe cambiar el estado de la cache ni al repetirse.' );
  %con sintaxis de FUNCION (nargout=0) imprime el informe; con punto NO puede,
  %porque subsref siempre pide una salida -- la misma regla que hace que
  %o.total muestre ans --, y entonces devuelve la lista
  V5 = cart();  V5.DEBUG = false;  V5 = V5.Demo();
  t24c = evalc( 'VerifyCPs( V5 )' );
  assert( contains( t24c , 'nada que verificar' ) , 'test_cache:verifyEmpty' , ...
          'sin valores cacheados deberia decirlo; dijo: %s' , strtrim(t24c) );
  assert( iscell( V5.VerifyCPs ) , 'test_cache:verifyDot' , ...
          'con sintaxis de punto deberia devolver la lista, no imprimir.' );
  t24d = evalc( 'VerifyCPs( V2 )' );
  assert( contains( t24d , '! liar' ) && contains( t24d , 'MIENTE' ) , 'test_cache:verifyReport' , ...
          'el informe deberia marcar liar y explicar; salio: %s' , t24d );
  fprintf( 'OK  VerifyCPs caza al handler que miente y al evento no declarado, sin tocar la cache\n' );

  %% -------- 25) refreshDefs: al cargar, el CODIGO gana al .mat
  %sin el gancho, un objeto cargado usa los compute y handlers DEL FICHERO para
  %siempre (la logica guardada ganaba a la nueva, y la unica salida era que cada
  %subclase escribiera su propio loadobj). pRefresh re-registra 'val' en su
  %refreshDefs; aqui se guarda la definicion "vieja" y se comprueba que tras
  %load sirve la del codigo -- y que una CP NO re-registrada conserva la suya.
  F5 = pRefresh();
  F5 = F5.Define( 'val'   , @(~) 1 , 'chg' , [] );    % la definicion "vieja" (el .mat)
  F5 = F5.Define( 'extra' , @(~) 7 );                 % CP ad hoc: nadie la re-registra
  pre5 = F5.val;
  assert( pre5 == 1 , 'test_cache:refreshPre' , 'antes de guardar val deberia dar 1 y dio %g.' , pre5 );
  save( fn , 'F5' );  L5 = load( fn );  G5 = L5.F5;
  v5 = G5.val;  e5 = G5.extra;
  assert( v5 == 2 , 'test_cache:refreshWins' , ...
          'tras load la definicion del CODIGO deberia ganar (2) y val dio %g.' , v5 );
  assert( e5 == 7 , 'test_cache:refreshKeeps' , ...
          'una CP no re-registrada deberia conservar la del fichero (7) y dio %g.' , e5 );
  fprintf( 'OK  refresh   refreshDefs re-registra al cargar: el codigo gana, lo ad hoc se conserva\n' );

  %% -------- 26) LA PRECONDICION TACITA: el compute tiene que ser PURO
  %VerifyCPs calcula "la verdad" RECALCULANDO desde cero, asi que una CP cuyo
  %valor dependa de la HISTORIA (de en que ORDEN se pidieron las cosas) y no
  %solo de los datos no es verificable: la herramienta canta discrepancias que
  %no lo son. Le paso a ivp y no era hipotetico -- la version obvia partia cada
  %tramo en n pasos, asi que la rejilla dependia del ancla: con DT=0.1, pedir
  %y(1) y luego y(1/3) daba 2.3197758575243266 y al reves 2.3197762151071655.
  %El demo se libraba por los pelos (4,5,6 son multiplos exactos de DT). La cura
  %fue anclar la rejilla al origen y ajustar a ella los tiempos pedidos.
  %Este grupo lo guarda con los tiempos QUE LO ROMPIAN, no con los del demo.
  f13 = @(t,y) cos(t) .* y;
  Pa  = ivp( f13 , 1 , 0.1 );  Pa = Pa.At( 1 );      Pa = Pa.At( 1/3 );   % desordenado
  Pb  = ivp( f13 , 1 , 0.1 );  Pb = Pb.At( 1/3 );    Pb = Pb.At( 1 );     % ordenado
  assert( isequaln( Pa.sol , Pb.sol ) , 'test_cache:ivpOrder' , ...
          [ 'el orden de las peticiones cambia la solucion (%.17g vs %.17g): la rejilla ha ' ...
            'vuelto a depender del ancla y `sol` ya no es funcion pura de los datos.' ] , ...
          Pa.sol(end,2) , Pb.sol(end,2) );
  assert( isempty( Pa.VerifyCPs() ) , 'test_cache:ivpVerifyOrder' , ...
          'VerifyCPs marca `sol` como rancia solo por el orden de las peticiones.' );
  %y el ajuste a la rejilla es visible y documentado, no un redondeo escondido
  [ ~ , ~ , ts ] = Pa.At( 0.75 );
  assert( abs( ts - 0.8 ) < 1e-12 , 'test_cache:ivpSnap' , ...
          'At debe devolver el tiempo AJUSTADO a la rejilla (0.75 -> 0.8 con DT=0.1) y dio %g.' , ts );
  fprintf( 'OK  pureza    el compute de ivp no depende del ORDEN de las peticiones (si no, VerifyCPs mentiria)\n' );

  %% -------- 27) ProfileCPs: medir si cachear compensa
  %La herramienta separa DOS preguntas que es facil confundir: (1) ¿deberia ser
  %una CP? -> trabajo contra HIT; (2) ¿conviene guardarla? -> MISS menos HIT por
  %lectura extra. Sin la primera, una CP CONSTANTE como currency (@(c)'EUR')
  %salia como "PAGA 4x", porque el camino MISS carga un par de cientos de us de
  %maquinaria pase lo que pase.
  P1 = cart();  P1.DEBUG = false;  P1 = P1.Demo();  P1 = P1.Add( 'a' , 1 , 1 );
  PR = P1.ProfileCPs();
  assert( numel( PR ) == 4 && all( isfield( PR , ...
            { 'name','work','compute','hit','bytes','gain','verdict' } ) ) , ...
          'test_cache:profileShape' , 'ProfileCPs deberia dar 4 entradas con los 7 campos.' );
  %el MISS contiene al compute pelado, asi que deberia costar MAS -- pero se
  %comparan DOS MEDIDAS distintas y el ruido puede invertirlas cuando las dos
  %son de microsegundos (visto en linsys: work 4.53 ms, MISS 3.50 ms). La
  %guarda es GRUESA a proposito: caza una inversion gorda, no el ruido.
  assert( all( [ PR.compute ] >= 0.5 * [ PR.work ] ) , 'test_cache:profileWork' , ...
          'el camino MISS cuesta MUCHO menos que el compute pelado que contiene: algo se mide mal.' );
  kc = strcmp( { PR.name } , 'currency' );          % @(c) 'EUR': trabajo ~0
  assert( contains( PR(kc).verdict , 'NO DEBERIA SER CP' ) , 'test_cache:profileConst' , ...
          'una CP constante deberia salir como "no deberia ser CP"; dijo: %s' , PR(kc).verdict );
  %y la CP cara del otro extremo tiene que salir claramente a favor
  P2 = collatz();  P2 = P2.Reserve( 20000 );
  QR = P2.ProfileCPs();
  kl = strcmp( { QR.name } , 'lens' );
  assert( contains( QR(kl).verdict , 'PAGA' ) && QR(kl).gain > 0 , ...
          'test_cache:profileHeavy' , ...
          'una tabla de 20000 cadenas de Collatz deberia salir como que PAGA; dijo: %s' , QR(kl).verdict );
  %con sintaxis de FUNCION imprime (misma regla que VerifyCPs)
  tp = evalc( 'ProfileCPs( P1 )' );
  assert( contains( tp , 'DOS PREGUNTAS DISTINTAS' ) && contains( tp , 'trabajo' ) , ...
          'test_cache:profileReport' , 'el informe deberia explicar las dos preguntas; salio: %s' , tp );
  fprintf( 'OK  ProfileCPs mide trabajo/MISS/HIT y distingue "no deberia ser CP" de "no conviene guardarla"\n' );

  %% -------- 28) el CONSEJERO (o.ADVISE): avisa de la CP que no se sostiene
  %Cuenta misses (gratis) y, cuando una CP lleva ADVISE_N recalculos, mide EN
  %CONDICIONES y opina. La medida NO se puede tomar dentro del propio MISS: ahi
  %se esta dentro del enter/onCleanup de esa clave y las tomas salen
  %contaminadas (comprobado: se tragaba casos que ProfileCPs caza con margen
  %7x). Por eso se apunta dentro y se emite en la siguiente lectura, ya sin nada
  %en vuelo. Clave UNICA: el memo de avisos vive en la sesion.
  [ ~ , u4 ] = fileparts( tempname );
  an = matlab.lang.makeValidName( [ 'adv' u4 ] );
  lastwarn( '' );
  for i = 1:6                                  % objetos NUEVOS: cada uno un miss
    Av = cart();  Av.DEBUG = false;  Av.ADVISE = true;
    Av = Av.Define( an , @(c) 1 );             % CP ridiculamente barata
    evalc( 'x = Av.( an ); x = Av.( an );' );  % miss + hit (el hit emite)
  end
  [ ~ , wa ] = lastwarn;
  assert( strcmp( wa , 'cachedOwner:notWorthCaching' ) , 'test_cache:adviseWarn' , ...
          'una CP constante con ADVISE deberia acabar avisando y llego "%s".' , wa );
  %APAGADO POR DEFECTO: sin ADVISE no se mide ni se avisa (si no, el veredicto
  %dependeria de lo rapida que sea la maquina y los goldens serian una loteria)
  [ ~ , u5 ] = fileparts( tempname );
  bn2 = matlab.lang.makeValidName( [ 'adv' u5 ] );
  lastwarn( '' );
  for i = 1:6
    Bv = cart();  Bv.DEBUG = false;            % ADVISE se queda en false
    Bv = Bv.Define( bn2 , @(c) 1 );
    evalc( 'x = Bv.( bn2 ); x = Bv.( bn2 );' );
  end
  [ ~ , wb ] = lastwarn;
  assert( isempty( wb ) , 'test_cache:adviseOff' , ...
          'con ADVISE apagado (el default) no deberia avisar nada; llego "%s".' , wb );
  %y el memo de "caso cerrado" hace lo que dice
  cachedOwner.adviseSettled( 'zzz_test' , true );
  assert(  cachedOwner.adviseSettled( 'zzz_test' ) , 'test_cache:settledYes' );
  assert( ~cachedOwner.adviseSettled( 'zzz_otra' ) , 'test_cache:settledNo'  );
  fprintf( 'OK  consejero o.ADVISE avisa de la CP que no se sostiene; apagado por defecto\n' );

  %% -------- 29) InfoCP: el informe de UNA CP, con historial
  %Reune lo que estaba repartido (tabla, trazas, ProfileCPs) y anade lo que no
  %habia: cuantas veces se ha recalculado, servido, replayado y caido.
  I1 = collatz();  I1 = I1.Reserve( 200 );
  v = I1.lens;  v = I1.best;  v = I1.lens;                                %#ok<NASGU>
  I1 = I1.Reserve( 5000 );                 % lens queda pendiente, best CAE
  v = I1.best;  v = I1.best;                                              %#ok<NASGU>
  SI = I1.InfoCP( 'lens' );
  assert( strcmp( SI.state , 'fresh' ) && SI.bytes > 0 , 'test_cache:infoState' , ...
          'lens deberia estar fresh y ocupar algo; dijo "%s" %d bytes.' , SI.state , SI.bytes );
  assert( isequal( fieldnames( SI.events ).' , { 'grow' } ) && ...
          strcmp( SI.events.grow , 'incremental' ) , 'test_cache:infoEvents' , ...
          'InfoCP deberia decir que grow es incremental.' );
  assert( isempty( SI.deps ) && isequal( SI.usedBy , { 'best' } ) , 'test_cache:infoDeps' , ...
          'lens no sale de nadie y la usa best; dijo sale de {%s} la usan {%s}.' , ...
          strjoin( SI.deps , ',' ) , strjoin( SI.usedBy , ',' ) );
  assert( SI.counts(1) >= 1 && SI.counts(2) >= 1 && SI.counts(3) >= 1 , ...
          'test_cache:infoCounts' , ...
          'lens deberia llevar >=1 recalculo, >=1 HIT y >=1 replay; lleva %s.' , mat2str(SI.counts) );
  SB = I1.InfoCP( 'best' );
  assert( isequal( SB.deps , { 'lens' } ) && SB.counts(4) >= 1 , 'test_cache:infoBest' , ...
          'best sale de lens y ha caido al menos una vez; dijo sale de {%s}, caidas %d.' , ...
          strjoin( SB.deps , ',' ) , SB.counts(4) );
  %MEDIR NO PUEDE FALSEAR EL HISTORIAL: cronometrar lee la CP decenas de veces
  %(timeit), y sin restaurar los contadores el informe mentia sobre si mismo
  %(visto: el contador de HIT se iba a 42 de una sola llamada)
  c1 = I1.InfoCP( 'lens' );  c2 = I1.InfoCP( 'lens' );
  assert( isequal( c1.counts , c2.counts ) , 'test_cache:infoPure' , ...
          'dos InfoCP seguidos dan historiales distintos (%s vs %s): medir esta contaminando.' , ...
          mat2str( c1.counts ) , mat2str( c2.counts ) );
  p1 = I1.InfoCP( 'lens' );  I1.ProfileCPs();  p2 = I1.InfoCP( 'lens' );
  assert( isequal( p1.counts , p2.counts ) , 'test_cache:profilePure' , ...
          'ProfileCPs deja el historial tocado (%s vs %s).' , ...
          mat2str( p1.counts ) , mat2str( p2.counts ) );
  %y la puerta del proxy imprime
  ti = evalc( 'I1.CP.lens.info' );
  assert( contains( ti , 'historial' ) && contains( ti , 'la usan' ) , 'test_cache:infoProxy' , ...
          'o.CP.<n>.info deberia imprimir el informe; salio: %s' , ti );
  expectError( @() I1.InfoCP( 'noExiste' ) , 'cachedOwner:noCP' );
  fprintf( 'OK  InfoCP    estado, eventos, deps, quien la usa, historial y coste -- sin falsear el historial\n' );

  fprintf( '\nTODOS LOS TESTS DEL PATRON DE CACHE: OK\n' );
end

function A = growIt( o )
  A = o;  A(3) = o;          % la via mas comun por la que nace un array
end

function o = setCP( o )
  o.total = 5;               % una CP no se asigna
end

function expectError( f , id )
  try
    f();
  catch ME
    assert( strcmp( ME.identifier , id ) , 'test_cache:wrongError' , ...
            'se esperaba %s y llego %s (%s)' , id , ME.identifier , ME.message );
    return;
  end
  error( 'test_cache:noError' , 'se esperaba el error %s y no hubo ninguno.' , id );
end
