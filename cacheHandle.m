classdef cacheHandle < handle
%cacheHandle  Almacen (por REFERENCIA) de valores derivados cacheados.
%
%   Es el "puntero" del patron: un handle que un objeto DUENO guarda como
%   propiedad privada y que sus copias comparten. Cada clave guarda una entrada
%   con ESTADO:
%     'fresh'    valor valido, se sirve tal cual
%     'pending'  valor viejo + log de eventos sin aplicar (replay perezoso: el
%                dueno lo resuelve en el proximo acceso, aplicando sobre el
%                valor viejo los handlers de evento de la definicion)
%   Las DEFINICIONES (que computa cada clave, y como reacciona a cada evento) NO
%   viven aqui: son asunto del dueno. Esta clase guarda VALORES y nada mas, y por
%   eso no sabe nada del dominio en el que se usa.
%
%   COPY-ON-WRITE. El uso previsto es que el dueno sea una clase de VALOR: al
%   editar sus datos NO muta este objeto, sino que construye un handle NUEVO con
%   las entradas que sobreviven o quedan pendientes, de modo que las copias
%   hermanas conservan su cache intacta. Si en cambio se MUTA este handle (borrar
%   un valor, resolver un pendiente, depositar un MISS) el efecto es compartido y
%   benigno: es estado de RENDIMIENTO, nunca dato.
%
%   EL ALMACEN ES UN STRUCT, no un containers.Map: se da por supuesto que las
%   claves ya son nombres MATLAB validos (quien registra las definiciones lo
%   exige) y el acceso por campo dinamico sale ~5x mas barato que un lookup de
%   Map -- importa porque el HIT esta en el camino caliente de toda lectura. De
%   paso hace `clone` UNA asignacion en vez de un bucle: el struct es un valor y
%   MATLAB ya lo copia perezosamente. Sirve la entrada COMPLETA (tryGet/entry)
%   para que quien lee no pague un lookup por campo.
%
%   LA PILA (`stk`) anota las claves en curso de resolucion, y sirve para dos
%   cosas: si el compute (o un handler) de una clave vuelve a pedir esa misma
%   clave, `enter` lanza un error claro en vez de recursar hasta agotar la pila;
%   y mientras hay alguien calculando, cada lectura de otra clave se apunta como
%   DEPENDENCIA suya (`noteDep`). Cada entrada guarda asi de que otras claves
%   salio su valor -- transitivamente --, que es lo que permite al dueno tirar
%   las CPs COMPUESTAS cuando cae aquella de la que salieron. Las dependencias
%   se DESCUBREN al calcular (lo que se leyo de verdad), no se declaran.
%
%   CICLO DE VIDA: no hay recolector, ni hace falta. Al derivar de `handle`,
%   MATLAB lo destruye por conteo de referencias en cuanto no queda ninguna, y el
%   unico puntero es la propiedad privada del dueno: cuando ese valor muere
%   (sobreescrito, fuera de scope, borrado) se libera el almacen entero.
%
% See also cacheProxy, cachedOwner.

  properties (Access = private)
    store = struct()   % clave -> struct('state','fresh'|'pending','value',v,'log',{eventos},'deps',{cellstr})
    stk   = {}         % PILA de claves en curso de resolucion (guarda de ciclos)
    acc   = {}         % lo que ha leido cada nivel de la pila (paralelo a stk)
    cnt   = struct()   % clave -> [miss hit replay caida]. VIVE FUERA DEL STORE a
                       % proposito: si viviera en la entrada, invalidar la clave
                       % borraria su historial, que es justo lo que se quiere
                       % contar. Lo copia `clone`, asi que el historial sigue a
                       % la estirpe del objeto (no al valor).
  end

  properties (SetAccess = private)
    %los dos se leen desde el camino caliente del dueno, asi que son PROPIEDADES
    %y no metodos: en un HIT hay que poder preguntar "¿me esta leyendo alguien?"
    %sin pagar una llamada
    nIn     = 0        % profundidad de la pila (0 = nadie esta calculando)
    hasDeps = false    % ¿alguna entrada registro dependencias? Es el atajo que
                       % permite que fire no pague NADA en las clases que no
                       % componen CPs. Puede quedarse en true de mas (borrar la
                       % ultima entrada con deps no lo baja): eso solo cuesta
                       % una pasada de mas, nunca correccion
  end

  methods
    function tf = has( obj , key )
      tf = isfield( obj.store , key );
    end

    function [ tf , e ] = tryGet( obj , key )
      %"esta? y dame la entrada" de una vez: el HIT no paga has+state+value
      tf = isfield( obj.store , key );
      if tf, e = obj.store.( key ); else, e = []; end
    end

    function e = entry( obj , key ), e = obj.store.( key ); end

    function setFresh( obj , key , v , d )
      %el {v} es obligatorio: struct('value',C) con C cell crearia un struct
      %ARRAY en vez de un escalar con el cell dentro
      if nargin < 4, d = {}; end                 %sin declarar = sin dependencias
      obj.store.( key ) = struct( 'state','fresh' , 'value',{v} , 'log',{{}} , 'deps',{d} );
      if ~isempty( d ), obj.hasDeps = true; end
    end

    function setPending( obj , key , v , L , d )
      %el valor pendiente sigue siendo EL VIEJO, asi que arrastra las
      %dependencias con las que se calculo
      if nargin < 5, d = {}; end
      obj.store.( key ) = struct( 'state','pending' , 'value',{v} , 'log',{L} , 'deps',{d} );
      if ~isempty( d ), obj.hasDeps = true; end
    end

    function d = deps( obj , key )
      %de que otras claves depende el valor guardado en `key`
      if isfield( obj.store , key ), d = obj.store.( key ).deps; else, d = {}; end
    end

    function ks = dropDependents( obj , key )
      %fuera todas las entradas calculadas A PARTIR de `key`. Lo usan el
      %recalculo forzado y el sembrado a mano: los dos cambian un valor sin
      %pasar por un evento, asi que quien lo hubiera usado queda rancio.
      ks = {};
      if ~obj.hasDeps, return; end
      for k = fieldnames( obj.store ).' , n = k{1};
        if any( strcmp( key , obj.store.( n ).deps ) ), ks{end+1} = n; end   %#ok<AGROW>
      end
      for i = 1:numel( ks ), obj.store = rmfield( obj.store , ks{i} ); end
    end

    function setEntry( obj , key , e )
      %copia una entrada TAL CUAL. Ya no lo usa el motor de la base: desde que
      %fire clona el store, la clave que sobrevive no se toca (era justo el
      %coste que hacia crecer una edicion con el numero de CPs cacheadas). Se
      %conserva porque SI lo usa el motor propio de msh, que todavia reconstruye
      %el store campo a campo. Cuando msh migre, esto se va.
      obj.store.( key ) = e;
    end

    function remove( obj , key )
      if isfield( obj.store , key ), obj.store = rmfield( obj.store , key ); end
    end

    function ks = keys( obj )
      ks = fieldnames( obj.store ).';   %FILA (como devolvia containers.Map)
    end

    function bump( obj , key , i )
      %+1 al contador i de `key` (1=miss 2=hit 3=replay 4=caida). Se llama desde
      %el camino caliente, asi que es lo mas corto que se pudo escribir.
      if isfield( obj.cnt , key )
        obj.cnt.( key )( i ) = obj.cnt.( key )( i ) + 1;
      else
        v = zeros( 1 , 4 );  v( i ) = 1;  obj.cnt.( key ) = v;
      end
    end

    function s = counts( obj , key )
      if isfield( obj.cnt , key ), s = obj.cnt.( key ); else, s = zeros( 1 , 4 ); end
    end

    function setCounts( obj , key , v )
      %devolver los contadores a un estado anterior. Existe por una razon muy
      %concreta: CRONOMETRAR UNA CP LA LEE, y `timeit` la lee decenas de veces,
      %asi que InfoCP y ProfileCPs inflaban justo el historial que estan
      %enseniando (medido: un InfoCP dejaba el contador de HIT en 42). Miden,
      %y devuelven la cuenta a donde estaba.
      obj.cnt.( key ) = v;
    end

    function ks = usedBy( obj , key )
      %quien SALE de `key`: el reverso de deps. Lo pregunta el informe por CP;
      %no esta en ningun camino caliente, asi que se recorre y ya.
      ks = {};
      if ~obj.hasDeps, return; end
      for k = fieldnames( obj.store ).' , n = k{1};
        if any( strcmp( key , obj.store.( n ).deps ) ), ks{end+1} = n; end   %#ok<AGROW>
      end
    end

    function c2 = clone( obj )
      %copia superficial: el store es un VALOR, basta asignarlo. La pila NO se
      %copia: es estado de UNA resolucion en curso, no del almacen. Los
      %CONTADORES si: son historial de esta estirpe de objetos.
      c2 = cacheHandle();
      c2.store   = obj.store;
      c2.hasDeps = obj.hasDeps;
      c2.cnt     = obj.cnt;
    end

    function c2 = cloneWithout( obj , key )
      c2 = obj.clone();
      c2.remove( key );
    end
  end

  methods   % ------------------------- guarda de ciclos Y captura de dependencias
    %La pila de claves en curso sirve para dos cosas a la vez: cortar los ciclos
    %(una clave que se pide a si misma) y saber QUIEN esta calculando cuando se
    %lee otra clave -- que es toda la maquinaria de dependencias. Por eso viven
    %juntas: la segunda sale casi gratis de la primera.
    function enter( obj , key )
      %marca `key` en curso de resolucion. Reentrar = el compute (o un handler)
      %de esa clave la vuelve a pedir: se corta con un error legible en vez de
      %recursar hasta agotar la pila.
      if any( strcmp( key , obj.stk ) )
        error('cacheHandle:cycle', ...
              'ciclo en la clave ''%s'': su compute (o un handler de evento) la vuelve a pedir.', key );
      end
      obj.stk{ end+1 } = key;
      obj.acc{ end+1 } = {};
      obj.nIn = numel( obj.stk );
    end

    function leave( obj , key )
      i = find( strcmp( key , obj.stk ) , 1 , 'last' );
      if isempty( i ), return; end
      obj.stk(i) = [];  obj.acc(i) = [];
      obj.nIn = numel( obj.stk );
    end

    function noteDep( obj , key )
      %se ha LEIDO `key` mientras se calculaba otra cosa: apuntala como
      %dependencia de todos los niveles en curso, JUNTO CON las dependencias que
      %`key` ya tenga registradas. Ese arrastre es lo que deja la lista
      %TRANSITIVAMENTE cerrada, y por eso a fire le basta UNA pasada: si A leyo
      %a B y B se calculo leyendo a C, A queda dependiendo de B y de C -- que es
      %lo que hace falta cuando el que cambia es C y B estaba cacheado.
      d = { key };
      if isfield( obj.store , key ), d = [ d , obj.store.( key ).deps ]; end
      for i = 1:obj.nIn
        if ~strcmp( obj.stk{i} , key ), obj.acc{i} = [ obj.acc{i} , d ]; end
      end
    end

    function d = takeDeps( obj )
      %lo leido por el nivel que esta a punto de terminar (se pide ANTES de
      %leave, con la clave todavia en la pila)
      if obj.nIn == 0, d = {}; return; end
      d = unique( obj.acc{ obj.nIn } , 'stable' );
    end
  end
end
