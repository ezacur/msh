classdef cachedOwner
%cachedOwner  Clase BASE para duenos de propiedades cacheadas: hereda y listo.
%
%   El patron completo son tres clases:
%     cacheHandle   plano de DATOS    -- el store compartido por las copias
%     cacheProxy    plano de CONTROL  -- reenvia o.CP.<nombre>[...] al dueno
%     cachedOwner   TODO LO DEMAS     -- los ganchos Y la politica (esta clase)
%   Una clase nueva hereda y aporta solo su dominio:
%
%     classdef miClase < cachedOwner
%       properties (Constant, Hidden)
%         CP_EVENTS = { 'miEvento' , 'general' }   % los eventos que la clase dispara
%       end
%       ...datos, editores que llaman a fire(), y sus Define...
%     end
%
%   LAS TRES REGLAS DE USO:
%     1. fire() en CADA metodo que edite datos, eventos de ESPECIFICO a GENERAL:
%            obj = obj.fire( { 'esp' , 'gen' } , argsDelEvento... )
%        La base no puede saber cuando cambia tu dominio; este es el paso que,
%        si se olvida, deja cache rancia. Los nombres se validan contra
%        CP_EVENTS: una errata es un error, no una cache rancia en silencio.
%     2. En Define, LA FIRMA DEL HANDLER ES LA POLITICA:
%            @(v,o,args...)  3+ argumentos = INCREMENTAL (uno por edit, en orden;
%                            el evento debe dispararse CON args). UNO POR CP:
%                            declarar dos es un error (el replay solo puede
%                            seguir una cadena incremental). Y NO PUEDE LEER `o`:
%                            se declara @(v,~,args...) -- ver el parrafo de abajo.
%            @(v,o)          2 argumentos  = SYNC ABSOLUTO (contra el estado
%                            actual; uno subsume N edits)
%            []              invalidar (cae y se recomputa al leer)
%            no declararlo   insensible (sobrevive intacto)
%     3. En Define, EL ORDEN DE DECLARACION ES LA PRIORIDAD del replay: declara
%        de especifico a general. Y DE UNA EDICION DECIDE EL EVENTO MAS
%        ESPECIFICO que la CP declare: su politica es la que se aplica, sea []
%        (cae), un sync absoluto o el atajo incremental. Un incremental
%        declarado en un evento GENERAL no adelanta al absoluto del ESPECIFICO.
%
%   NO CACHEES LO BARATO. Una lectura de CP cuesta ~110-135 us solo en llegar
%   (medido; una propiedad normal de una clase que sobrecarga subsref ya cuesta
%   ~105 us, y un campo de struct 1.6 us): el peaje es del despacho de MATLAB,
%   no del patron, pero se paga igual. Cachear un derivado que se computa en
%   menos que eso es PERDIDA NETA. Regla gruesa: de 0.5 ms para arriba. Los
%   ejemplos cart/route cachean sumas de tres elementos A PROPOSITO -- ensenan
%   el mecanismo con valores que caben en la cabeza --, no son una recomendacion.
%
%   POR QUE TODO VIVE AQUI (comprobado empiricamente, ver el tutorial): tanto
%   `builtin('subsref',obj,s)` como escribir propiedades privadas exigen estar
%   en un metodo de la clase del objeto -- y un metodo HEREDADO cuenta. Una
%   clase ajena no puede con ninguna de las dos cosas: una version anterior del
%   patron separaba la politica en un `cacheEngine` estatico y pagaba el peaje
%   de pasarle el estado empaquetado en cada llamada.
%
%   SAVE/LOAD sigue la asimetria del patron: DEFS viaja en el .mat (las
%   definiciones son DATO), CACHE es Transient y loadobj la revive vacia (los
%   valores son RENDIMIENTO: se recomputan al leer). TRES CONSECUENCIAS:
%     - los compute/handlers se SERIALIZAN. Un handle anonimo se lleva su
%       workspace capturado, asi que escribelos leyendo del objeto -- @(o)f(o.X)
%       -- y no capturando datos: @(o)f(X) mete una copia de X en cada .mat.
%     - si la subclase define su PROPIO loadobj, tiene que llamar al de aqui
%       (obj = loadobj@cachedOwner(obj)) o la cache llegara muerta.
%     - LA LOGICA GUARDADA GANA A LA NUEVA: un objeto cargado usa los compute y
%       handlers DEL FICHERO, no los del codigo actual (comprobado: mejorar un
%       compute no afecta a los .mat viejos). Y si RENOMBRAS un evento, las CPs
%       guardadas declaran el nombre viejo y pareceran insensibles: dato rancio
%       en silencio -- loadobj lo detecta y avisa (cachedOwner:staleEvents).
%       Si tu clase va a EVOLUCIONAR con .mat vivos, sobrescribe el gancho
%       PROTEGIDO refreshDefs y re-registra ahi (con Define) tus definiciones:
%       loadobj lo llama al cargar y el CODIGO -- no el fichero -- vuelve a ser
%       la fuente de la logica. Lo que no re-registres (las CPs ad hoc del
%       usuario) conserva lo del .mat.
%
%   DECLARAR DE MENOS ES EL DESCUIDO MAS FACIL, y no falla: la CP sobrevive al
%   evento que no declaro y sirve un valor rancio en silencio. Solo se puede
%   juzgar automaticamente el olvido TOTAL -- una CP que LEE datos del objeto y
%   no declara ni un evento --, y de ese avisa Define (cachedOwner:noEvents,
%   silenciable, y marcado tambien en la tabla). El olvido PARCIAL no se puede
%   ver desde aqui: haria falta saber que propiedad toca cada evento, y eso solo
%   lo sabe quien escribe los editores. La herramienta para ese caso es comparar
%   la lectura con el recalculo forzado: o.<nombre> frente a o.<nombre>_.
%
%   UNA CP PUEDE SALIR DE OTRA, Y EL MOTOR SE ENTERA SOLO. Si al calcular una CP
%   se lee otra, esa lectura queda apuntada como DEPENDENCIA (se descubre, no se
%   declara: cuenta lo que de verdad se leyo), y cuando un evento se lleva por
%   delante a la segunda, la primera CAE con ella. Sin eso, la compuesta -- que
%   normalmente no declara ningun evento -- SOBREVIVIA sirviendo un valor sacado
%   de datos que ya no valen, en silencio. Tres cosas que conviene saber:
%     - una dependencia rota GANA a cualquier declaracion: el handler de un
%       evento puede poner al dia por ESE evento, pero no puede saber que ha
%       cambiado la CP de la que sale su valor. Por eso la compuesta cae aunque
%       tuviera handler. Y la dependencia SOBREVIVE a los replays: poner el
%       valor al dia no cambia DE QUE CP salio, asi que las deps viejas se
%       conservan al refrescar (perderlas dejaba a la compuesta con handler
%       sobrevivir a la caida de su base a partir del primer replay).
%     - lo mismo vale para los dos gestos que cambian un valor SIN evento: el
%       recalculo forzado (o.<n>_) y el sembrado a mano (o.CP.<n>.set) tiran lo
%       que se hubiera calculado a partir de ellos.
%     - DENTRO DE LA CLASE NO SE PUEDE LEER CON UN PUNTO. MATLAB no pasa por
%       subsref en un metodo de la clase -- y una anonima escrita ahi hereda ese
%       contexto --, asi que @(o) o.otra*2 muere diciendo que "otra" no existe
%       (miente: existe). Desde dentro se lee o.access('otra'), que es la misma
%       puerta y tambien registra la dependencia; desde fuera, o.otra vale.
%       El error lo traduce runCompute (cachedOwner:cpFromInside).
%   El coste vive donde debe: si NINGUNA CP compone, el almacen lo sabe
%   (hasDeps) y ni fire ni la tabla hacen trabajo de mas.
%
%   UN HANDLER INCREMENTAL NO PUEDE LEER EL DUENO. El log guarda EDICIONES, no
%   estados: cuando el replay recorre N edits pendientes, el `o` que recibe el
%   handler es SIEMPRE el objeto de AHORA -- el de despues de las N -- y no el
%   que habia cuando ocurrio cada una. Un handler que lo lea da un valor MAL, y
%   en silencio (comprobado con la forma natural de una media incremental,
%   @(v,o,x) (v*(numel(o.VAL)-1)+x)/numel(o.VAL): con [2 4 6] devolvia 3.78 en
%   vez de 4, y con UNA sola edicion pendiente acertaba -- o sea, pasa los tests
%   de juguete y falla en produccion). Por eso:
%     - Define RECHAZA (cachedOwner:incrementalUsesOwner) el handler anonimo
%       cuyo cuerpo nombra su segundo argumento. La forma correcta es declararlo
%       no usado: @(v,~,args...). Si de verdad necesitas mirar el objeto, lo que
%       tienes no es un incremental sino un SYNC ABSOLUTO: @(v,o).
%     - lo que Define NO puede ver -- un handler OPACO (@funcion, @Clase.metodo:
%       func2str no trae cuerpo) o uno que llega de un .mat anterior a esta
%       guarda -- lo cierra el replay: en cuanto hay MAS DE UNA edicion
%       pendiente renuncia al atajo y degrada (dato correcto, atajo perdido) y
%       avisa una vez (cachedOwner:incrementalOpaque / :incrementalUsesOwner).
%       Con UNA sola edicion el dueno si es el correcto y el atajo se usa.
%
%   UN VALOR QUE ES HANDLE ROMPE EL AISLAMIENTO: todo lo demas de este patron
%   descansa en que MATLAB copia los valores, pero un objeto handle cacheado se
%   comparte por REFERENCIA entre las copias, asi que mutarlo in place (tipico en
%   un handler incremental) corrompe la cache de las hermanas -- y eso ya es
%   DATO malo, no rendimiento. Cachea valores; si tiene que ser un handle,
%   copialo en el compute y devuelve la copia desde el handler. Se avisa al
%   detectarlo (cachedOwner:handleValue).
%
%   LA CLASE ES ESCALAR, Y SE OBLIGA. El motor indexa por OBJETO, no por
%   elemento: en un array no hay a que clave cachear, y antes de cerrarlo un
%   array daba errores cripticos de MATLAB ("Too many input arguments") en vez de
%   decir que pasaba. Se bloquean las tres puertas por las que nace un array
%   (horzcat, vertcat, cat) y el () en subsref y subsasgn, todo con el mismo
%   error cachedOwner:notScalar. Se deja pasar o(1), que es un no-op inofensivo
%   del codigo generico. Si necesitas varios, mete los objetos en un CELL.
%
%   Si la subclase necesita su PROPIO subsref (arrays, alias legados, reentrada
%   en cadenas... el caso de msh), lo sobreescribe y llama a los PROTEGIDOS de
%   aqui: access (lectura perezosa), recompute (forzada), cachedAccess (el
%   plano de control) y, si le hace falta, liveCache y dbg.
%
%   HERENCIA EN CADENA: CP_EVENTS es Constant, asi que la clase HOJA es la que
%   la declara -- y debe listar tambien los eventos que disparen sus ancestros.
%
% See also cacheHandle, cacheProxy, cart, route, serie, facto, ivp.

  properties (Abstract, Constant, Hidden)
    CP_EVENTS   % cell con TODOS los eventos que la clase dispara. Es el
                % vocabulario valido de Define: caza las erratas de evento
  end

  properties
    DEBUG  = false  % true -> traza MISS/HIT/RPLAY/EVENT por consola
    ADVISE = false  % true -> el CONSEJERO: avisa (una vez) de la CP que no se
                    % sostiene, midiendo el compute en los MISS. No cambia el
                    % comportamiento, solo avisa. Ver noteCost.
                    % OPT-IN a proposito: cronometrar tiene coste (~1 us por
                    % MISS) y, sobre todo, un aviso que aparece segun lo rapida
                    % que sea la maquina rompe cualquier salida comparada byte a
                    % byte -- los goldens de este repo, sin ir mas lejos.
                    % Enciendelo cuando estes decidiendo, no de continuo.
  end

  properties (Access = private)
    DEFS = struct()              % nombre -> struct('compute',fh,'events',struct)
  end
  properties (Access = private, Transient)
    CACHE = cacheHandle.empty    % el "puntero" a los VALORES (compartido por copias)
  end

  %% ==================================================== REGISTRO
  methods
    function obj = cachedOwner()
      obj.CACHE = cacheHandle();          % handle vivo desde el inicio
    end

    function obj = Define( obj , name , computeFcn , varargin )
      %obj = obj.Define( nombre , @(o)... [, evento , handler|[] , ...] )
      %
      %   Registra (o REDEFINE, descartando el valor previo por COW) una CP.
      %   Ver las reglas 2 y 3 de la cabecera: firma = politica, orden de
      %   declaracion = prioridad.
      if ~isvarname( name )
        error('cachedOwner:name','"%s" no es un identificador MATLAB valido.', name );
      end
      %guardas de SOMBRA. subsref consulta las definiciones ANTES de caer al
      %builtin, asi que una CP taparia -- en silencio -- cualquier cosa que viva
      %en ese espacio de nombres: una propiedad, un metodo, la puerta CP, o el
      %recalculo forzado de otra CP.
      if name(end) == '_'
        %id PROPIO y no el generico de nombre invalido: son dos arreglos
        %distintos (uno es "eso no es un identificador", el otro "ese sufijo es
        %del recalculo forzado") y quien lo capture querra distinguirlos
        error('cachedOwner:reservedSuffix', ...
              '"%s": los nombres de CP no pueden terminar en ''_'' (sufijo reservado al recalculo forzado).', name );
      end
      %properties()/methods() OMITEN los miembros Hidden y los privados, asi que
      %la guarda dejaba pasar -- comprobado -- CP_EVENTS, DEFS, CACHE, loadobj,
      %cachedAccess y displayCachedView: precisamente el espacio de nombres que
      %dice cubrir. El metaclass los lista TODOS.
      %Hacen falta LAS DOS listas: la de la subclase no incluye las privadas de
      %la base (son inaccesibles alli), asi que sin ?cachedOwner se colaban DEFS
      %y CACHE -- comprobado.
      mc    = metaclass( obj );
      mb    = ?cachedOwner;
      taken = [ { 'CP' } ; { mc.PropertyList.Name }.' ; { mc.MethodList.Name }.' ; ...
                           { mb.PropertyList.Name }.' ; { mb.MethodList.Name }.' ];
      if ismember( name , taken )
        error('cachedOwner:shadow', ...
              '"%s" ya es propiedad, metodo o puerta de %s: la CP lo sombrearia.', name , class( obj ) );
      end
      if ~isa( computeFcn , 'function_handle' )
        error('cachedOwner:compute','computeFcn debe ser un function handle @(o)...');
      end
      if mod( numel( varargin ) , 2 )
        error('cachedOwner:events','los eventos van en pares nombre,handler.');
      end
      %GUARDA DEL OLVIDO TOTAL. Una CP que lee DATOS del objeto y no declara ni
      %un evento es casi siempre un descuido: sobrevive a todas las ediciones y
      %sirve un valor rancio en silencio. Solo se puede juzgar el olvido TOTAL --
      %con eventos declarados, saber si faltan exigiria un mapa propiedad->evento
      %que solo conoce quien escribe los editores --, asi que es un AVISO y no un
      %error, y con `warning('off','cachedOwner:noEvents')` se calla. Para el
      %olvido PARCIAL la herramienta sigue siendo comparar con o.<nombre>_.
      %Se miran solo las propiedades NO Constant y que no sean de la base: una
      %constante no cambia nunca, y DEBUG/DEFS/CACHE no son dominio.
      %definicion nueva, historial limpio: lo que se hubiera anotado sobre esta
      %CP (handler roto, valor handle, olvido de eventos) era de la definicion
      %ANTERIOR y ya no dice nada de esta
      cachedOwner.warnOnce( '' , [ class( obj ) '_' name ] );
      if isempty( varargin )
        pl   = mc.PropertyList;
        cand = setdiff( { pl( ~[ pl.Constant ] ).Name } , { mb.PropertyList.Name } );
        rd   = cachedOwner.propsRead( computeFcn , cand );
        if ~isempty( rd )
          cachedOwner.warnOnce( 'cachedOwner:noEvents' , [ class( obj ) '_' name ] , ...
            [ 'la CP "%s" de %s lee datos del objeto (%s) y NO declara ningun evento: sobrevivira ' ...
              'intacta a TODAS las ediciones y servira un valor rancio en cuanto eso cambie. ' ...
              'Declara en Define los eventos que le afectan ([] para invalidarla, o un handler ' ...
              'para ponerla al dia). Si de verdad no le afecta ninguno, silencia el aviso con ' ...
              'warning(''off'',''cachedOwner:noEvents'').' ] , ...
            name , class( obj ) , strjoin( rd , ', ' ) );
        end
      end
      ev = struct();                       %un struct CONSERVA el orden de insercion
      for a = 1:2:numel( varargin )
        en = varargin{a};
        if ~ismember( en , obj.CP_EVENTS )
          error('cachedOwner:unknownEvent','evento desconocido "%s" (validos: %s).', ...
                en , strjoin( obj.CP_EVENTS , ', ' ) );
        end
        h = varargin{a+1};
        if ~( isempty( h ) || isa( h , 'function_handle' ) )
          error('cachedOwner:handler','el handler de %s debe ser un handle o [] (invalidar).', en );
        end
        %FIRMA AMBIGUA. La politica se lee del nargin, y con un varargin que
        %empieza antes del tercer argumento no hay forma de saber la intencion:
        %@(varargin) da -1 y @(v,varargin) da -2, y los dos caian del lado
        %ABSOLUTO en silencio -- el handler se invocaba como h(v,o), asi que si
        %el autor lo queria incremental el replay fallaba y degradaba a
        %recompute, con el atajo ahi de adorno (comprobado). Se exige declararlo.
        if ~isempty( h )
          na = nargin( h );
          if na == -1 || na == -2
            error('cachedOwner:ambiguousHandler', ...
                  [ 'el handler de %s tiene firma ambigua (nargin=%d): con varargin ahi no se ' ...
                    'puede saber si es INCREMENTAL o SYNC ABSOLUTO. Declaralo: @(v,o) para ' ...
                    'absoluto, @(v,o,varargin) para incremental.' ] , en , na );
          end
        end
        ev.( en ) = h;
      end
      %UN SOLO handler incremental por CP: el replay sigue UNA cadena, asi que
      %un segundo quedaria muerto sin avisar (auditado: el log colapsaba a
      %marcador y la lectura acababa recomputando, con el handler ahi de adorno).
      inc = fieldnames( ev ).';
      inc = inc( cellfun( @(f) ~isempty( ev.(f) ) && cachedOwner.isInc( ev.(f) ) , inc ) );
      if numel( inc ) > 1
        error('cachedOwner:twoIncremental', ...
              '"%s" declara %d handlers incrementales (%s): solo se admite uno (los demas serian @(v,o) absolutos).', ...
              name , numel( inc ) , strjoin( inc , ', ' ) );
      end
      %EL INCREMENTAL SE RESUELVE AQUI, UNA VEZ. Antes lo buscaba `incOf` en cada
      %edicion (una vez por clave pendiente) y en cada replay, y buscarlo cuesta
      %un nargin() por handler: medido, ~180 us frente a ~5 us por leer un campo.
      %Registrar es frio; editar y leer no.
      %Y SE COMPRUEBA QUE NO LEA EL DUENO (ver el parrafo de la cabecera): en el
      %replay de VARIAS ediciones ese argumento es el estado FINAL, no el de cada
      %edicion, asi que un handler que lo lea da dato malo en silencio. Aqui se
      %caza el caso visible -- el handler ANONIMO --; del opaco se encarga el
      %replay. La respuesta es 'no' | 'yes' | 'unknown'.
      if isempty( inc ), ie = '';  ih = [];  io = 'no';
      else
        ie = inc{1};  ih = ev.( ie );
        [ io , onm ] = cachedOwner.ownerUse( ih );
        if strcmp( io , 'yes' )
          error('cachedOwner:incrementalUsesOwner', ...
                [ 'el handler incremental de "%s" (evento %s) lee su segundo argumento "%s", que ' ...
                  'es el DUENO: en el replay de varias ediciones pendientes ese argumento es el ' ...
                  'estado FINAL, no el de cada edicion, y el valor saldria mal EN SILENCIO. Si no ' ...
                  'lo necesitas, declaralo no usado: @(v,~,args...). Si lo necesitas, lo que tienes ' ...
                  'es un SYNC ABSOLUTO: @(v,o). (Si "%s" solo aparece dentro de un texto, renombra ' ...
                  'el argumento.)' ] , name , ie , onm , onm );
        end
      end
      existed = isfield( obj.DEFS , name );
      obj.DEFS.( name ) = struct( 'compute' , computeFcn , 'events' , ev , ...
                                  'inc' , ie , 'incH' , ih , 'incOwn' , io );
      c = obj.liveCache();
      if existed && ~isempty( c )
        %REDEFINIR ES CAMBIAR LA VERDAD, y arrastra a los que salieron de ella.
        %Tirar el valor viejo ya se hacia (era de OTRA definicion); lo que
        %faltaba era tirar tambien las CPs COMPUESTAS que se calcularon con el:
        %su valor salio de un compute que ya no existe, asi que es rancio.
        %Comprobado antes del arreglo: con base=10*N y viaAccess=base+2,
        %redefinir base a @(o)999 dejaba viaAccess sirviendo 12 con la verdad en
        %1001. Es el mismo agujero que ya tapaban recompute y .set, por la
        %puerta de al lado.
        c2 = c.clone();                       % COW: la hermana no se entera
        c2.remove( name );
        c2.dropDependents( name );
        obj.CACHE = c2;
      end
      obj.dbg( 'DEF   "%s" (eventos: %s)' , name , cachedOwner.lst( fieldnames(ev).' ) );
    end

    function varargout = VerifyCPs( obj )
      %[malas, R] = obj.VerifyCPs()   compara lo que SIRVE la cache con un
      %calculo desde cero, CP a CP.
      %VerifyCPs( obj )               lo mismo, pero IMPRIME el informe
      %
      %   Las dos formas no son capricho: con sintaxis de PUNTO, subsref pide
      %   siempre una salida -- es la misma regla que hace que `o.total` sin
      %   punto y coma muestre `ans` --, asi que `o.VerifyCPs` nunca llega aqui
      %   con nargout==0 y devuelve la lista. Para ver el informe, la llamada
      %   con sintaxis de FUNCION: VerifyCPs( o ), que no pasa por el despacho.
      %
      %   ES LA HERRAMIENTA PARA LOS DOS FALLOS QUE NINGUNA GUARDA PUEDE VER:
      %     - un handler que MIENTE: no falla, devuelve un valor equivocado (el
      %       que revienta ya se caza solo: cachedOwner:handlerFailed).
      %     - el olvido PARCIAL de eventos: declarar de menos, y la CP sobrevive
      %       a una edicion que si le afectaba (del olvido TOTAL avisa Define).
      %   Los dos se ven igual desde fuera -- la cache sirve una cosa y la
      %   verdad es otra -- y los dos son invisibles hasta que alguien compara.
      %
      %   NO TOCA LA CACHE: el replay de un pendiente se ejecuta SIN guardar el
      %   resultado, asi que el objeto queda igual y esto se puede repetir. Lo
      %   unico que si muta: calcular la verdad LEE, y leer resuelve pendientes
      %   de las CPs de las que dependa la que se esta verificando.
      %
      %   Solo mira las CPs CON VALOR: una que no se ha calculado no puede estar
      %   mal. Cuesta un compute por CP cacheada -- es una herramienta de
      %   depuracion y de test, no un guardia en caliente. En un test:
      %       assert( isempty( obj.VerifyCPs() ) , 'hay CPs rancias' )
      %
      %   OJO CON LOS VALORES HANDLE: se comparan por IDENTIDAD, asi que un
      %   compute que devuelva un handle nuevo saldra siempre como discrepancia.
      %   El informe lo marca; es falsa.
      %
      %   Y OJO CON LO COMPUESTO: la verdad de una CP que sale de otra se
      %   calcula con el valor CACHEADO de esa otra. Si la rancia es la de
      %   abajo, la que aparece en el informe es ella, que es la que hay que
      %   arreglar.
      ns = sort( fieldnames( obj.DEFS ).' );
      c  = obj.liveCache();
      R  = struct( 'name' , {} , 'state' , {} , 'ok' , {} , 'served' , {} , 'truth' , {} , 'err' , {} );
      if ~isempty( c )
        for i = 1:numel( ns ), n = ns{i};
          [ hit , e ] = c.tryGet( n );
          if ~hit, continue; end
          k = numel( R ) + 1;
          R(k).name = n;  R(k).state = e.state;  R(k).err = '';
          R(k).served = [];  R(k).truth = [];  R(k).ok = false;
          try
            if strcmp( e.state , 'fresh' ), sv = e.value;
            else, sv = obj.replay( n , obj.DEFS.( n ) , e.value , e.log );   % SIN guardar
            end
            tv = obj.truthOf( obj.DEFS.( n ) , n , c );
            R(k).served = sv;  R(k).truth = tv;  R(k).ok = isequaln( sv , tv );
          catch ME
            R(k).err = sprintf( '%s: %s' , ME.identifier , ME.message );
          end
        end
      end
      bad = { R( ~[ R.ok ] ).name };
      if nargout == 0
        cachedOwner.showVerify( class( obj ) , R );
      else
        out = { bad , R };
        varargout = out( 1:max( min( nargout , 2 ) , 1 ) );
      end
    end

    function varargout = ProfileCPs( obj )
      %R = obj.ProfileCPs()   mide, CP a CP, si cachear compensa
      %ProfileCPs( obj )      lo mismo, pero IMPRIME el informe
      %
      %   Las dos formas, por la misma razon que VerifyCPs: con sintaxis de
      %   PUNTO subsref siempre pide una salida, asi que o.ProfileCPs devuelve
      %   la tabla; para ver el informe, la llamada de FUNCION ProfileCPs(o).
      %
      %   QUE MIDE Y QUE NO. Mide dos cosas por CP, las dos en esta maquina y
      %   por el camino REAL (el HIT se cronometra llamando al subsref, no al
      %   motor por dentro, para que incluya el peaje del despacho):
      %       compute   lo que cuesta calcularla desde cero
      %       HIT       lo que cuesta servirla ya cacheada
      %   Y de ahi sale lo unico que importa de verdad: cada lectura EXTRA
      %   entre dos ediciones te ahorra (compute - HIT). Si eso es negativo,
      %   cachear esa CP es perdida neta -- recomputar sale mas barato que
      %   servirla --, y da igual el resto.
      %
      %   LO QUE NO PUEDE SABER, y por eso no da un si/no absoluto: CUANTAS
      %   VECES la lees entre ediciones. Ese numero lo pones tu, y es el que
      %   manda: si lees una sola vez por edicion, cachear no compensa NUNCA
      %   (recalculas igual y encima cada clave cacheada engorda el fire de
      %   cada edicion); si lees muchas, compensa en cuanto el ahorro sea
      %   positivo. El informe da el ahorro por lectura extra; multiplicalo por
      %   tus lecturas.
      %
      %   EL COSTE QUE SE MIDE ES EL MARGINAL, con las dependencias ya
      %   calientes: de una CP COMPUESTA se cronometra lo que cuesta ELLA, no
      %   lo que cuesta ella mas rehacer aquellas de las que sale. Es lo que
      %   hace falta para decidir -- la pregunta es "si NO la cacheo, que me
      %   cuesta cada lectura", y ahi el resto de la cache sigue funcionando.
      %   (El primer calculo de todos, con la cache entera fria, cuesta mas;
      %   ese numero no sirve para esta decision.)
      %
      %   OJO, Y ES LA DIFERENCIA CON VerifyCPs: esto SI toca la cache -- la
      %   deja CALIENTE, porque para cronometrar un HIT hace falta que haya
      %   valor. No cambia ningun resultado (leer es idempotente), pero si el
      %   estado: una CP que estaba sin calcular acaba calculada.
      %   Y CUESTA: ejecuta el compute de cada CP al menos una vez. Es una
      %   herramienta de banco, no de produccion.
      ns = sort( fieldnames( obj.DEFS ).' );
      R  = struct( 'name' , {} , 'work' , {} , 'compute' , {} , 'hit' , {} , ...
                   'bytes' , {} , 'gain' , {} , 'verdict' , {} );
      %cronometrar LEE las CPs (timeit, decenas de veces por clave), asi que se
      %guarda el historial y se devuelve al final: medir no puede falsear los
      %contadores que luego enseña InfoCP
      cc = obj.liveCache();  snap = struct();
      if ~isempty( cc )
        for i = 1:numel( ns ), snap.( ns{i} ) = cc.counts( ns{i} ); end
      end
      for i = 1:numel( ns )
        n = ns{i};  k = numel( R ) + 1;
        R(k).name = n;
        %POBLAR PRIMERO, Y NO ES UN DETALLE DE ORDEN: si se cronometra el
        %compute con la cache fria, una CP COMPUESTA arrastra el calculo de
        %aquellas de las que sale y sale un numero que no significa nada para
        %decidir (medido: `best` de collatz daba 68 ms -- 40 de ellos eran
        %calcular `lens` -- cuando su coste propio son microsegundos). Leerla
        %antes deja calientes sus dependencias, y entonces el compute mide su
        %coste MARGINAL: lo que costaria de verdad NO tenerla cacheada,
        %suponiendo que el resto de la cache se comporta con normalidad. Ese es
        %el numero con el que se decide.
        v = obj.access( n );
        w = whos( 'v' );  R(k).bytes = w.bytes;
        %COMPUTE: sin guardar (misma puerta que VerifyCPs), para medir el
        %camino frio y no un replay
        R(k).compute = cachedOwner.timeAdaptive( ...
                         @() obj.truthOf( obj.DEFS.( n ) , n , obj.liveCache() ) );
        %POR EL CAMINO REAL: llamar a subsref explicitamente SI pasa por el
        %despacho sobrecargado (dentro de un metodo, o.(n) no pasaria), asi que
        %el numero incluye el peaje que paga el usuario
        s = substruct( '.' , n );
        R(k).hit  = cachedOwner.timeAdaptive( @() subsref( obj , s ) );
        R(k).gain = R(k).compute - R(k).hit;
        %EL TRABAJO DE VERDAD, sin la maquinaria: se llama al compute PELADO.
        %Hace falta porque hay DOS preguntas distintas y confundirlas da
        %consejos absurdos. `compute` incluye el camino MISS entero (buscar la
        %definicion, pila, onCleanup, deposito), que en esta maquina son un par
        %de cientos de microsegundos PASE LO QUE PASE: por eso una CP
        %constante como @(c)'EUR' parecia "PAGA 4x". Comparando el TRABAJO con
        %el HIT se ve la otra pregunta, que es la que suele importar: si
        %calcularlo cuesta menos que servirlo cacheado, esto no deberia ser una
        %CP -- deberia ser un metodo, y te ahorras hasta el despacho.
        cf = obj.DEFS.( n ).compute;
        R(k).work = cachedOwner.timeAdaptive( @() cf( obj ) );
        if R(k).work < R(k).hit
          R(k).verdict = sprintf( ...
            'NO DEBERIA SER CP: calcularla (%s) cuesta menos que servirla (%s). Hazla un metodo.' , ...
            cachedOwner.fmtT( R(k).work ) , cachedOwner.fmtT( R(k).hit ) );
        elseif R(k).gain <= 0
          R(k).verdict = 'NO PAGA: recomputar cuesta menos que servirlo';
        elseif R(k).work < 2 * R(k).hit
          R(k).verdict = 'MARGINAL: el trabajo real apenas supera al peaje de servirlo';
        else
          R(k).verdict = sprintf( 'PAGA (%.0fx): cada lectura extra ahorra %s' , ...
                                  R(k).compute / R(k).hit , cachedOwner.fmtT( R(k).gain ) );
        end
      end
      if ~isempty( cc )
        for i = 1:numel( ns ), cc.setCounts( ns{i} , snap.( ns{i} ) ); end
      end
      if nargout == 0
        cachedOwner.showProfile( class( obj ) , R );
      else
        varargout = { R };
      end
    end

    function varargout = InfoCP( obj , name )
      %S = obj.InfoCP( 'x' )   TODO lo que el motor sabe de una CP, en un struct
      %InfoCP( obj , 'x' )     lo mismo, impreso
      %                        (tambien: obj.CP.x.info, que siempre imprime)
      %
      %   Reune en un sitio lo que estaba repartido entre la tabla, las trazas y
      %   ProfileCPs, y anade lo que no habia: el HISTORIAL. Campos:
      %     state     'fresh' | 'pending' (con cuantas ediciones) | 'sin calcular'
      %     bytes     lo que ocupa el valor ahora
      %     events    struct evento -> 'incremental' | 'sync' | 'invalida'
      %     deps      de que CPs SALE su valor (descubierto, no declarado)
      %     usedBy    que CPs salen DE ELLA (el reverso: quien cae si esta cae)
      %     counts    [miss hit replay caida] desde que nacio esta estirpe
      %     cost      coste del compute y de un HIT, medidos AHORA (o NaN si se
      %               pide sin medir: ver el argumento de abajo)
      %
      %   LOS CONTADORES SIGUEN A LA ESTIRPE, no al valor: viven fuera del store
      %   (si vivieran en la entrada, invalidar la clave borraria justo el
      %   historial que se quiere ver) y `clone` los copia, asi que una copia
      %   hermana hereda el historial y sigue contando por su cuenta.
      %
      %   MEDIR CUESTA (ejecuta el compute una vez), asi que solo se hace con la
      %   forma que imprime o pidiendo el struct de una CP concreta. Es la misma
      %   herramienta que ProfileCPs pero de una en una y con historial.
      if ~isfield( obj.DEFS , name )
        error('cachedOwner:noCP','no hay CP "%s" en %s.', name , class( obj ) );
      end
      d = obj.DEFS.( name );
      c = obj.liveCache();
      S = struct( 'name' , name , 'class' , class( obj ) , 'state' , 'sin calcular' , ...
                  'pending' , 0 , 'bytes' , 0 , 'events' , struct() , 'deps' , {{}} , ...
                  'usedBy' , {{}} , 'counts' , zeros(1,4) , 'work' , NaN , 'hit' , NaN );
      for f = fieldnames( d.events ).'
        S.events.( f{1} ) = cachedOwner.polDesc( d.events.( f{1} ) );
      end
      if ~isempty( c )
        S.counts = c.counts( name );
        S.deps   = c.deps( name );
        S.usedBy = c.usedBy( name );
        [ hit , e ] = c.tryGet( name );
        if hit
          S.state = e.state;
          if strcmp( e.state , 'pending' ), S.pending = numel( e.log ); end
          v = e.value;  w = whos( 'v' );  S.bytes = w.bytes;
        end
      end
      %los tiempos, con el mismo criterio que ProfileCPs (y su misma cautela:
      %poblar antes, para medir el coste MARGINAL de una CP compuesta)
      try
        v  = obj.access( name );                                          %#ok<NASGU>
        cf = d.compute;
        S.work = cachedOwner.timeAdaptive( @() cf( obj ) );
        s2 = substruct( '.' , name );
        S.hit  = cachedOwner.timeAdaptive( @() subsref( obj , s2 ) );
      catch
        %informar no puede romper nada: si el compute revienta, se deja NaN
      end
      %MEDIR HA LEIDO LA CP (timeit, decenas de veces): devolver el historial a
      %donde estaba, o el informe mentiria sobre si mismo
      if ~isempty( c ), c.setCounts( name , S.counts ); end
      if nargout == 0
        cachedOwner.showInfo( S );
      else
        varargout = { S };
      end
    end

    function obj = RemoveCP( obj , name )
      %obj = obj.RemoveCP( nombre )   borra definicion Y valor
      if ~isfield( obj.DEFS , name )
        error('cachedOwner:noCP','no hay CP "%s".', name );
      end
      c = obj.liveCache();
      if isempty( c )
        c = cacheHandle();                    %handle muerto/vacio: REVIVIRLO, como
                                              %hace fire. Si no, asignar la vacia
                                              %dejaba la cache apagada para siempre
                                              %(access computa sin cachear: nunca
                                              %mas un HIT, comprobado)
      else
        %quitar la definicion se lleva por delante lo que salio de ella: si no,
        %una CP compuesta seguia sirviendo su valor viejo mientras que forzar su
        %recalculo REVENTABA (MATLAB:nonExistentField, porque el compute busca
        %una CP que ya no existe) -- las dos cosas a la vez, que es lo peor
        c = c.clone();
        c.remove( name );
        c.dropDependents( name );
      end
      obj.DEFS  = rmfield( obj.DEFS , name );
      obj.CACHE = c;
    end
  end

  %% ==================================================== DESPACHO
  methods
    function varargout = subsref( obj , s )
      %ESCALAR SIEMPRE (ver la nota de la cabecera). Va lo PRIMERO porque todo lo
      %que sigue lee obj.DEFS, y sobre un array eso se expande a cs-list y muere
      %con un "Too many input arguments" que no le dice nada a nadie.
      if ~isscalar( obj ), cachedOwner.errNotScalar( class( obj ) , 'leer' ); end
      while strcmp( s(1).type , '()' )
        %o(1) es un no-op inofensivo que aparece en codigo generico: se QUITA el
        %eslabon y se sigue con el despacho normal de abajo -- no se recursa,
        %porque recursar exigia una salida y un STATEMENT del plano de control
        %(o(1).CP.x.delete devuelve cero) moria en unassignedOutputs (auditado).
        %Cualquier otro indice ya es pensar en arrays, y aqui no los hay.
        if ~all( cellfun( @(k) isequal( k , 1 ) , s(1).subs ) )
          cachedOwner.errNotScalar( class( obj ) , 'indexar con ()' );
        end
        if isscalar( s ), varargout = { obj }; return; end
        s = s(2:end);
      end
      %las cuatro puertas, de mas caliente a mas fria, y el fallback:
      %   o.<nombre>          lectura perezosa (el camino que se hace mil veces)
      %   o.CP[...]           plano de control (proxy)
      %   o.<nombre>_         RECALCULO forzado
      %   o.<nombre>.algo     leer e indexar DENTRO del valor
      %   todo lo demas       propiedades y metodos, via builtin
      %Este metodo tiene que ser de la clase (heredado cuenta): el builtin de la
      %ultima linea solo esquiva la sobrecarga desde aqui; fuera, recursa.
      nm = '';
      if strcmp( s(1).type , '.' ) && ( ischar( s(1).subs ) || isstring( s(1).subs ) )
        nm = char( s(1).subs );
      end
      if ~isempty( nm )
        if isscalar( s ) && isfield( obj.DEFS , nm )     % ---- CAMINO CALIENTE
          varargout = { obj.access( nm ) };
          return;
        end
        if strcmp( nm , 'CP' )                           % ---- plano de control
          if isscalar( s ), varargout{1} = cacheProxy( obj ); return; end
          out = obj.cachedAccess( s(2:end) );
          if isempty( out ), varargout = {};             % statement (.delete)
          else, varargout = out( 1:max( min( nargout , numel(out) ) , 1 ) );
          end
          return;
        end
        if numel(nm) > 1 && nm(end) == '_' && isfield( obj.DEFS , nm(1:end-1) )
          v = obj.recompute( nm(1:end-1) );              % ---- sufijo '_'
          if numel( s ) > 1
            [ varargout{ 1:max(nargout,1) } ] = subsref( v , s(2:end) );
          else
            varargout = { v };
          end
          return;
        end
        if isfield( obj.DEFS , nm )                      % ---- indexar DENTRO
          %subsref y NO builtin: el resto de la cadena es asunto del VALOR, y si
          %el valor es un objeto con su propio subsref hay que respetarlo (con
          %builtin, una CP que devuelve un cart moria en noSuchMethodOrField).
          %No hay riesgo de recursion: la cadena siempre se acorta.
          %Y se piden TODAS las salidas que MATLAB pida: pedir siempre una
          %truncaba las listas separadas por comas en silencio (ver
          %numArgumentsFromSubscript).
          [ varargout{ 1:max(nargout,1) } ] = subsref( obj.access( nm ) , s(2:end) );
          return;
        end
      end
      [ varargout{ 1:max(nargout,1) } ] = builtin( 'subsref' , obj , s );
    end

    function n = numArgumentsFromSubscript( obj , s , ctx )
      %cuantas salidas pide MATLAB a subsref. El default pide numel(obj(1))=1
      %tambien para un STATEMENT que empieza con (), asi que o(1).CP.x.delete
      %-- que devuelve CERO salidas -- moria en unassignedOutputs aunque el
      %despacho lo tratara bien (auditado; la cadena con '.' inicial no pasaba
      %por aqui y funcionaba). Statement -> 0; el resto -> 1. Devolver una
      %salida de mas cuando piden 0 es legal y es lo que setea `ans`, asi que
      %los statements normales (o.total sin ';') siguen mostrando su valor.
      %Un [a,b] = ... multi-salida NO pasa por aqui: el numero de LHS manda.
      if ctx == matlab.mixin.util.IndexingContext.Statement, n = 0; else, n = 1; end
      %...PERO una cadena que indexa DENTRO del valor puede dar una LISTA
      %SEPARADA POR COMAS, y decir 1 la truncaba EN SILENCIO: {o.cp.campo} con
      %un struct array 1x3 devolvia un elemento en vez de tres, y la misma
      %expresion pasando el valor por una variable devolvia tres (comprobado).
      %El numero no lo sabemos nosotros: lo sabe EL VALOR, y se le pregunta --
      %vale para struct arrays, para cells y para objetos con su propio
      %subsref. MATLAB solo llama aqui cuando la lista importa (medido: una
      %lectura simple ni pasa por aqui), asi que el camino caliente no se entera.
      if numel( s ) < 2 || ~isscalar( obj ) || ~strcmp( s(1).type , '.' ) || ...
         ~( ischar( s(1).subs ) || isstring( s(1).subs ) )
        return;
      end
      nm = char( s(1).subs );
      try
        if strcmp( nm , 'CP' )                          % el plano de control
          n = obj.cachedNumArgs( s(2:end) , ctx );  return;
        end
        forced = numel(nm) > 1 && nm(end) == '_' && isfield( obj.DEFS , nm(1:end-1) );
        if forced, base = nm(1:end-1); else, base = nm; end
        if ~isfield( obj.DEFS , base ), return; end
        if forced
          %con '_' NO se computa solo para contar: la cadena va a recalcular de
          %todas formas, y forzar aqui el calculo lo haria DOS veces (un
          %statement sin ';' construiria el bvh dos veces). Se cuenta si hay
          %valor a mano; si no, se deja el default y a lo sumo se trunca.
          c = obj.liveCache();
          if isempty( c ), return; end
          [ hit , e ] = c.tryGet( base );
          if ~hit || ~strcmp( e.state , 'fresh' ), return; end
          v = e.value;
        else
          v = obj.access( base );          % lo va a leer igual: queda cacheado
        end
        n = cachedOwner.csCount( v , s(2:end) , ctx , n );
      catch
        %contar no puede romper nada: si algo falla, el default -- y que el
        %error salga (con su mensaje de verdad) cuando corra el subsref
      end
    end

    function obj = subsasgn( obj , s , v )
      %Dos guardas y al builtin. EL COSTE NO ESTA DONDE PARECE: MATLAB NO llama a
      %este metodo para las escrituras hechas DENTRO de un metodo de la clase
      %(comprobado con un contador: 0 llamadas tras dos escrituras internas), asi
      %que los editores del dominio -- que son los que escriben en bucle -- no
      %pagan nada. Solo lo paga la escritura EXTERNA (o.DEBUG = true), que pasa
      %de ~24 us a ~0.5 ms. Es un precio que se paga una vez y a mano.
      if strcmp( s(1).type , '()' )
        cachedOwner.errNotScalar( class( obj ) , 'asignar por indice' );
      end
      if isscalar( obj ) && strcmp( s(1).type , '.' ) && ...
         ( ischar( s(1).subs ) || isstring( s(1).subs ) )
        nm = char( s(1).subs );
        if isfield( obj.DEFS , nm )
          %sin esto MATLAB suelta "Unrecognized property" sobre algo que SI se
          %puede leer, que es de las cosas que mas despistan del patron
          error('cachedOwner:readOnlyCP', ...
                [ '"%s" es una CP (propiedad DERIVADA): no se le asigna. Para sembrar un ' ...
                  'valor a mano, o = o.CP.%s.set( valor ); para cambiar lo que calcula, ' ...
                  'o = o.Define( ''%s'' , @(o)... ).' ] , nm , nm , nm );
        end
      end
      obj = builtin( 'subsasgn' , obj , s , v );
    end

    %LAS TRES PUERTAS POR LAS QUE NACE UN ARRAY. Con estas mas el guard de () en
    %subsref/subsasgn se cerraron las 16 vias que probe (concat, cat(3,..),
    %repmat, A(3)=o, B(2)=o con B nueva, o(ones), o([1 1]), o(:), o(2),
    %horzcat, [s.campo] desde un struct array...). SOBREVIVE UNA: `empty`, que
    %es estatico y sigue dando un 0x0 -- se deja a proposito, porque MATLAB lo
    %usa por dentro y el 0x0 es INERTE (comprobado: .total, .CP, .Demo() y (1)
    %sobre el dan todos cachedOwner:notScalar).
    %Lo limpio seria heredar de matlab.mixin.Scalar, que existe en R2022a y
    %bloquea todo esto de oficio -- pero es INCOMPATIBLE con el patron: MATLAB
    %rechaza la clase al cargarla ("Class 'cachedOwner' is not allowed to define
    %a 'subsref' method"), y el patron ES un subsref. Se probo: no arranca.
    function varargout = horzcat( varargin ) , cachedOwner.errCat( varargin ); end  %#ok<STOUT>
    function varargout = vertcat( varargin ) , cachedOwner.errCat( varargin ); end  %#ok<STOUT>
    function varargout = cat( ~ , varargin ) , cachedOwner.errCat( varargin ); end  %#ok<STOUT>
  end

  methods (Hidden)   % --------------------- el contrato de cacheProxy
    function out = cachedAccess( obj , s )
      %despacho de o.CP.<nombre>[...]. Devuelve un CELL ({} para los statements)
      if ~strcmp( s(1).type , '.' )
        error('cachedOwner:proxy','use .CP.<nombre>.');
      end
      name = char( s(1).subs );
      if ~isfield( obj.DEFS , name )
        error('cachedOwner:noCP','no hay CP "%s" definida.', name );
      end

      if isscalar( s ), out = { obj.access( name ) }; return; end

      if strcmp( s(2).type , '.' )
        opn = char( s(2).subs );
        switch opn
          case 'delete'
            %borra el VALOR. Muta el handle COMPARTIDO (las hermanas tambien lo
            %pierden: benigno, como mucho alguien recomputa) y es un STATEMENT.
            if numel( s ) > 2, error('cachedOwner:usage','.delete no admite mas indexacion.'); end
            c = obj.liveCache();
            if ~isempty( c ), c.remove( name ); end
            obj.dbg( 'CACHE "%s" valor borrado (definicion intacta)' , name );
            out = {};  return;

          case 'info'
            %o.CP.<n>.info -> imprime TODO lo que el motor sabe de esa CP.
            %Siempre imprime (como .delete es siempre statement); para el struct
            %esta S = o.InfoCP('<n>')
            if numel( s ) > 2, error('cachedOwner:usage','.info no admite mas indexacion.'); end
            InfoCP( obj , name );
            out = {};  return;

          case 'removeCP'
            %borra definicion Y valor -> dueno NUEVO (hay que reasignar)
            if numel( s ) > 2, error('cachedOwner:usage','.removeCP no admite mas indexacion.'); end
            out = { obj.RemoveCP( name ) };  return;

          case 'set'
            %siembra un valor a mano. CONSERVADOR: clona (COW) para que la
            %mentira quede aislada en el dueno que se devuelve. NADIE la verifica.
            if numel( s ) ~= 3 || ~strcmp( s(3).type , '()' ) || numel( s(3).subs ) ~= 1
              error('cachedOwner:usage','uso: obj = obj.CP.%s.set( valor ).', name );
            end
            c = obj.liveCache();
            if isempty( c ), c = cacheHandle(); else, c = c.clone(); end
            cachedOwner.chkHandleValue( class( obj ) , name , s(3).subs{1} );
            c.dropDependents( name );      %el valor cambia sin pasar por un
                                           %evento: lo que salio de el es rancio
            c.setFresh( name , s(3).subs{1} );
            obj.CACHE = c;
            obj.dbg( 'CACHE "%s" valor sembrado a mano' , name );
            out = { obj };  return;

          otherwise
            %o.CP.<nombre>.<evento> -> el handler declarado, tal cual ([] si la
            %politica era invalidar). Se INVOCA si la cadena trae un '()'.
            ev = obj.DEFS.( name ).events;
            if isfield( ev , opn )
              h = ev.( opn );
              if ~( numel( s ) >= 3 && strcmp( s(3).type , '()' ) )
                out = { h };  return;
              end
              if isempty( h )
                error('cachedOwner:notCallable','el evento %s de "%s" es [] (invalidar): no es invocable.', opn , name );
              end
              v = h( s(3).subs{:} );
              if numel( s ) > 3, v = subsref( v , s(4:end) ); end
              out = { v };  return;
            end
        end
      end
      %cualquier otra cosa: indexar DENTRO del valor (subsref, no builtin: ver
      %la nota en el subsref de arriba). El cell se llena con TODAS las salidas
      %de la cadena -- si es una lista separada por comas, el proxy ya se queda
      %con las que le pidan (out(1:min(nargout,...)))
      v = obj.access( name );
      k = cachedOwner.csCount( v , s(2:end) , matlab.mixin.util.IndexingContext.Expression , 1 );
      [ out{ 1:max(k,1) } ] = subsref( v , s(2:end) );
    end

    function n = cachedNumArgs( obj , s , ctx )
      %cuantas salidas da una cadena del PLANO DE CONTROL (o.CP.<...>). Lo
      %pregunta el proxy -- que no conoce el vocabulario de operaciones -- y
      %tambien el numArgumentsFromSubscript de arriba para las cadenas que
      %empiezan por .CP. Las OPERACIONES (delete/removeCP/set y los handlers de
      %evento) devuelven lo de siempre; solo se le pregunta al valor cuando la
      %cadena es indexacion de verdad.
      if ctx == matlab.mixin.util.IndexingContext.Statement, n = 0; else, n = 1; end
      if numel( s ) < 2 || ~strcmp( s(1).type , '.' ), return; end
      name = char( s(1).subs );
      if ~isfield( obj.DEFS , name ), return; end
      if strcmp( s(2).type , '.' )
        op = char( s(2).subs );
        if any( strcmp( op , { 'delete' , 'removeCP' , 'set' , 'info' } ) ) || ...
           isfield( obj.DEFS.( name ).events , op )
          return;
        end
      end
      try
        n = cachedOwner.csCount( obj.access( name ) , s(2:end) , ctx , n );
      catch
      end
    end

    function displayCachedView( obj )
      %la tabla que imprime el disp del proxy
      ns = sort( fieldnames( obj.DEFS ).' );
      if isempty( ns ), fprintf('  (sin CPs definidas)\n\n'); return; end
      w  = max( cellfun( @numel , ns ) ) + 1;
      c  = obj.liveCache();
      vn = upper( class( obj ) );  vn = vn(1);   %cart -> C, msh -> M
      %LOS AVISOS, TAMBIEN AQUI. Los fallos que puede haber -- un evento huerfano
      %tras renombrar, un valor handle que el COW no aisla, una CP que lee datos
      %sin declarar eventos, y un handler roto (este ultimo no da dato malo: solo
      %deja la clase lenta en silencio) --
      %se avisan con `warning`, y un warning se pierde con facilidad: en -batch,
      %enterrado en un script largo, o con el warning('off') de alguien mas
      %arriba. Esta tabla es el sitio PERSISTENTE donde alguien mira el estado,
      %asi que el aviso vive tambien aqui. La columna del marcador solo aparece
      %si hay algo que marcar: una tabla limpia no cambia ni un byte.
      pl   = metaclass( obj ).PropertyList;
      mb   = ?cachedOwner;
      cand = setdiff( { pl( ~[ pl.Constant ] ).Name } , { mb.PropertyList.Name } );
      cls  = class( obj );
      note = repmat( { '' } , 1 , numel( ns ) );
      for i = 1:numel( ns )
        hasV = false;  val = [];  d = obj.DEFS.( ns{i} );
        if ~isempty( c ), [ hasV , e ] = c.tryGet( ns{i} );  if hasV, val = e.value; end, end
        rd = {};
        if isempty( fieldnames( d.events ) ), rd = cachedOwner.propsRead( d.compute , cand ); end
        %ATAJOS ROTOS: la cuenta la lleva el memo de avisos, y se pregunta aqui
        %porque un warning salta una vez y se lo lleva el viento; la tabla no
        bad = {};
        for f = fieldnames( d.events ).'
          nf = cachedOwner.warnOnce( 'cachedOwner:handlerFailed' , [ cls '_' ns{i} '_' f{1} ] );
          if nf > 0, bad{end+1} = sprintf( '%s (x%d)' , f{1} , nf ); end            %#ok<AGROW>
        end
        note{i} = cachedOwner.cpNote( d.events , obj.CP_EVENTS , hasV , val , rd , bad );
      end
      flag = ~cellfun( @isempty , note );
      pre  = repmat( { '  ' } , 1 , numel( ns ) );   pre( flag ) = { '! ' };
      fprintf( '  CPs (%d) -- leer %s.<n> | recalcular %s.<n>_ | control %s.CP.<n> :\n' , ...
               numel(ns) , vn , vn , vn );
      %la leyenda solo si hay algo que leyendar: lo unico que no se adivina es
      %que el ORDEN de los eventos es la PRIORIDAD del replay
      if any( cellfun( @(x) ~isempty( fieldnames( obj.DEFS.(x).events ) ) , ns ) )
        fprintf( '  (eventos en orden de PRIORIDAD; incremental = atajo por edicion, sync = absoluto, invalida = cae)\n' );
      end
      for i = 1:numel( ns ), n = ns{i};
        if isempty( c ),                      st = '(sin calcular)';
        else
          [ hit , e ] = c.tryGet( n );
          if ~hit,                            st = '(sin calcular)';
          elseif strcmp( e.state , 'fresh' ), st = [ 'fresh: ' cachedOwner.fmt( e.value ) ];
          else, st = sprintf( '(pendiente: %d edit(s) en el log)' , numel( e.log ) );
          end
        end
        %LA POLITICA, NO SOLO EL NOMBRE. La firma del handler ES la politica, y
        %una firma no se ve: la tabla decia "eventos: additem, changeprice" sin
        %distinguir el atajo incremental del sync absoluto ni del [] que
        %invalida -- justo lo que la aridad esconde. Van EN ORDEN DE PRIORIDAD,
        %que es el de declaracion (fieldnames conserva el de insercion) y es el
        %que sigue el replay.
        ee  = obj.DEFS.(n).events;
        evs = fieldnames( ee ).';
        if isempty( evs ), evs = { '(ninguno: solo MISS/HIT)' };
        else
          evs = cellfun( @(e) sprintf( '%s->%s' , e , cachedOwner.polDesc( ee.(e) ) ) , ...
                         evs , 'UniformOutput' , false );
        end
        %las dependencias solo se imprimen si LAS HAY: una clase que no compone
        %CPs saca exactamente la misma tabla que antes, byte a byte
        dp = '';
        if ~isempty( c ), d = c.deps( n );
          if ~isempty( d ), dp = sprintf( ' | sale de: %s' , strjoin( d , ', ' ) ); end
        end
        fprintf( '  %s%-*s %-34s eventos: %s%s\n' , pre{i} , w , n , st , strjoin( evs , ', ' ) , dp );
      end
      if any( flag )
        fprintf( '  ! %d CP(s) CON AVISO:\n' , sum( flag ) );
        for i = find( flag ), fprintf( '      %s %s\n' , ns{i} , note{i} ); end
      end
      fprintf( '\n' );
    end
  end

  methods (Static, Hidden)
    function obj = loadobj( obj )
      %tras load, CACHE (Transient) llega muerta: revivirla vacia. Las
      %definiciones si viajaron en el .mat; los valores se recomputan al leer
      if isa( obj , 'cachedOwner' )
        obj.CACHE = cacheHandle();
        %MIGRACION de .mat viejos: DEFS viaja, asi que un fichero guardado antes
        %de precomputar el incremental trae entradas SIN los campos inc/incH y
        %fire/replay moririan buscandolos. Se rellenan aqui, una sola vez, con el
        %mismo incOf que hacia el trabajo antes (por eso sigue existiendo).
        ns = fieldnames( obj.DEFS ).';
        for i = 1:numel( ns )
          d = obj.DEFS.( ns{i} );  dirty = false;
          if ~isfield( d , 'inc' )
            [ d.inc , d.incH ] = cachedOwner.incOf( d.events );
            dirty = true;
          end
          %incOwn: los .mat anteriores a la guarda del incremental-que-lee-el-
          %dueno no lo traen, y ademas pueden traer justo el handler que hoy
          %Define rechaza. Se clasifica al cargar; si sale 'yes', el replay lo
          %cazara en caliente (y avisara) en vez de servir dato malo.
          if ~isfield( d , 'incOwn' )
            if isempty( d.incH ), d.incOwn = 'no';
            else,                 d.incOwn = cachedOwner.ownerUse( d.incH );
            end
            dirty = true;
          end
          if dirty, obj.DEFS.( ns{i} ) = d; end
        end
        %GANCHO DE EVOLUCION: la subclase con .mat vivos sobrescribe refreshDefs
        %y re-registra ahi sus definiciones -- el CODIGO vuelve a ser la fuente
        %de la logica; lo que no re-registre (CPs ad hoc del usuario) conserva
        %lo del fichero. Va ANTES del aviso de eventos huerfanos: lo que el
        %refresh arregla ya no es rancio, y del resto se avisa igual.
        obj = refreshDefs( obj );
        %EVOLUCION: si la clase RENOMBRO un evento despues de guardarse este
        %objeto, sus CPs declaran el nombre viejo, el codigo nuevo dispara el
        %nuevo, y la CP parece INSENSIBLE: sobrevive y sirve dato rancio en
        %silencio (comprobado con dos versiones de una clase). Es la misma
        %caceria de erratas de fire, por la puerta de atras -- asi que se avisa.
        cachedOwner.warnStaleEvents( class( obj ) , obj.DEFS , obj.CP_EVENTS );
      end
    end

    function out = costMemo( key , t )
      %memo de costes del CONSEJERO, dos modos:
      %   n = costMemo( key , t )   anota una medida y devuelve CUANTAS van
      %   m = costMemo( key )       devuelve la MEDIANA de lo anotado
      %OJO CON ESAS MEDIDAS: vienen de un tic/toc AISLADO dentro del MISS, y
      %estan sesgadas hacia arriba (arranque en frio; medido, 270x en un compute
      %barato). Sirven para CONTAR misses -- que es lo que dispara al consejero
      %-- y como registro aproximado, NUNCA para dictar el veredicto: ese lo da
      %una medida en condiciones (ver noteCost). Hidden y no private para que la
      %bateria pueda comprobarlo, como ownerUse o propsRead.
      %Mediana y no media: la PRIMERA toma carga el calentamiento y puede ser
      %60x la siguiente (medido en linsys); una media quedaria secuestrada por
      %ese valor.
      persistent M
      if isempty( M ), M = struct(); end
      f = matlab.lang.makeValidName( key );
      if nargin < 2
        if isfield( M , f ), out = median( M.( f ) ); else, out = NaN; end
        return;
      end
      if isfield( M , f ), M.( f )( end+1 ) = t; else, M.( f ) = t; end
      out = numel( M.( f ) );
    end

    function varargout = advicePending( key , name )
      %cola (de uno) de veredictos por emitir. Dos modos:
      %   advicePending( key , name )   apuntar
      %   [k,n] = advicePending()       retirar el apuntado ('' si no hay)
      persistent K N
      if nargin > 0, K = key;  N = name;  return; end
      varargout = { K , N };
      K = '';  N = '';
    end

    function tf = adviseSettled( key , set )
      %¿esta CP ya se gano el puesto de calle? Dos modos:
      %   tf = adviseSettled( key )        preguntar
      %        adviseSettled( key , true ) cerrar el caso
      %Existe para no re-medir lo caro: evaluar cuesta un compute entero, y las
      %CPs que claramente compensan son justo las que mas duele recalcular.
      persistent S
      if isempty( S ), S = struct(); end
      f = matlab.lang.makeValidName( key );
      if nargin > 1, S.( f ) = set;  tf = set;  return; end
      tf = isfield( S , f ) && S.( f );
    end

    function pj = hitToll( obj , name )
      %EL PEAJE: lo que cuesta SERVIR una CP ya cacheada, por el camino real
      %(subsref incluido). Es practicamente el mismo para cualquier clave -- es
      %despacho y un lookup --, asi que se mide UNA vez por clase y por sesion y
      %se reusa. Se mide con la clave que se acaba de calcular, que por
      %definicion esta caliente.
      persistent P
      if isempty( P ), P = struct(); end
      f = matlab.lang.makeValidName( class( obj ) );
      if isfield( P , f ), pj = P.( f ); return; end
      s  = substruct( '.' , name );
      pj = cachedOwner.timeAdaptive( @() subsref( obj , s ) );
      P.( f ) = pj;
    end

    function [ u , nm ] = ownerUse( h )
      %'no' | 'yes' | 'unknown': ¿el handler INCREMENTAL lee su segundo
      %argumento (el DUENO)? Leerlo da dato malo en el replay de varias
      %ediciones -- ver la cabecera --, asi que Define lo rechaza y el replay
      %renuncia al atajo cuando no se puede saber. `nm` es el nombre de ese
      %argumento, para poder decirlo en el mensaje de error.
      %  Hidden y no private: es una REGLA del contrato y la suite la prueba
      %  caso a caso (como warnStaleEvents/cpNote).
      %  Solo se puede mirar dentro de una ANONIMA: func2str de @funcion o de
      %  @Clase.metodo devuelve el nombre, sin cuerpo -> 'unknown'.
      u = 'unknown';  nm = '';
      s = func2str( h );
      if numel( s ) < 3 || ~strcmp( s(1:2) , '@(' ), return; end
      k = find( s == ')' , 1 );          %la lista de argumentos no lleva parentesis
      if isempty( k ), return; end
      p = strtrim( strsplit( s(3:k-1) , ',' ) );
      if numel( p ) < 2, return; end     %no puede pasar en un incremental (nargin>=3)
      nm = p{2};
      if strcmp( nm , '~' ), u = 'no'; return; end          %declarado NO USADO
      %el nombre cuenta solo como VARIABLE: ni detras de un punto (s.o es un
      %campo, no el dueno) ni pegado a otras letras (otro, xo), ni dentro de un
      %literal de texto ('... o ...' es un mensaje, no un uso)
      b = cachedOwner.stripLits( s(k+1:end) );
      if isempty( regexp( b , [ '(?<![\w.])' nm '(?![\w])' ] , 'once' ) )
        u = 'no';
      else
        u = 'yes';
      end
    end

    function ps = propsRead( h , cand )
      %que propiedades de `cand` lee un compute @(o)... Mismo analisis que
      %ownerUse -- solo se ve dentro de una ANONIMA, y el nombre no cuenta si
      %esta dentro de un literal --, y devuelve {} cuando no se puede saber, que
      %es lo conservador: el aviso solo salta con pruebas.
      %Se buscan los `o.LOQUESEA` del primer argumento, y de esos solo los que
      %son PROPIEDADES: leer otra CP (o.otra) o llamar a un metodo
      %(o.access('otra')) no es leer datos -- de las CPs ya se encargan las
      %dependencias.
      ps = {};
      s = func2str( h );
      if numel( s ) < 3 || ~strcmp( s(1:2) , '@(' ), return; end
      k = find( s == ')' , 1 );
      if isempty( k ), return; end
      p  = strtrim( strsplit( s(3:k-1) , ',' ) );
      nm = p{1};
      if isempty( nm ) || strcmp( nm , '~' ), return; end
      b = cachedOwner.stripLits( s(k+1:end) );
      t = regexp( b , [ '(?<![\w.])' nm '\.(\w+)' ] , 'tokens' );
      if isempty( t ), return; end
      ps = unique( [ t{:} ] , 'stable' );
      ps = ps( ismember( ps , cand ) );
    end

    function warnIncOwner( cls , name , kind , n )
      %el aviso de que el replay ha renunciado al atajo incremental. UNA vez por
      %clase.clave: es un defecto de DISENO del handler, no de esta lectura.
      %(Hidden para poder probar las dos ramas sin fabricar un .mat viejo.)
      if strcmp( kind , 'yes' )
        cachedOwner.warnOnce( 'cachedOwner:incrementalUsesOwner' , [ cls '_' name ] , ...
          [ 'la CP "%s" de %s trae un handler incremental que LEE el dueno (de un .mat anterior a ' ...
            'la guarda de Define). Con %d ediciones pendientes ese handler daria un valor MAL, asi ' ...
            'que se descarta el atajo y se recalcula: el DATO es correcto, el atajo se pierde. ' ...
            'Re-registrala con Define: @(v,~,args...) si no necesita el dueno, o @(v,o) si si.' ] , ...
          name , cls , n );
      else
        cachedOwner.warnOnce( 'cachedOwner:incrementalOpaque' , [ cls '_' name ] , ...
          [ 'el handler incremental de la CP "%s" de %s no es anonimo (@funcion o @Clase.metodo), ' ...
            'asi que no se puede comprobar si lee el dueno -- y leerlo daria un valor MAL con ' ...
            'varias ediciones pendientes. Con %d se descarta el atajo y se recalcula: el DATO es ' ...
            'correcto, el atajo se pierde. Escribelo como anonimo: @(v,~,args...).' ] , ...
          name , cls , n );
      end
    end

    function s = cpNote( ev , valid , hasVal , v , rd , bad )
      %el aviso de UNA CP para la tabla ('' si no hay ninguno). Hidden y con los
      %inputs explicitos -- como warnStaleEvents -- para poder testear la rama
      %del evento huerfano, que de otro modo exigiria tener dos versiones de una
      %misma clase vivas en la sesion.
      %`rd` (opcional) son las propiedades que lee un compute que no declaro
      %ningun evento: el mismo descuido que avisa Define, aqui de forma
      %PERSISTENTE (un warning se pierde; la tabla no).
      %`bad` (opcional) son los handlers que han fallado y cuantas veces: no dan
      %dato malo -- la cascada recalcula -- pero dejan la clase silenciosamente
      %lenta, asi que van los ULTIMOS en la prioridad de esta nota.
      if nargin < 5, rd  = {}; end
      if nargin < 6, bad = {}; end
      s = '';
      orph = setdiff( fieldnames( ev ).' , valid );
      if ~isempty( orph )
        s = sprintf( [ 'declara evento(s) que YA NO ESTAN en CP_EVENTS (%s): no reacciona a ' ...
                       'las ediciones y sirve valores rancios. Re-registrala con Define.' ] , ...
                     strjoin( orph , ', ' ) );
      elseif hasVal && isa( v , 'handle' )
        s = sprintf( [ 'cachea un objeto handle (%s): el COW no aisla ahi, mutarlo in place ' ...
                       'corrompe la cache de las copias hermanas.' ] , class( v ) );
      elseif ~isempty( rd )
        s = sprintf( [ 'lee datos del objeto (%s) y no declara ningun evento: sobrevive a todas ' ...
                       'las ediciones y servira un valor rancio. Declaralos en Define.' ] , ...
                     strjoin( rd , ', ' ) );
      elseif ~isempty( bad )
        s = sprintf( [ 'ATAJO ROTO (no es dato malo): el handler de %s ha fallado y se recalcula ' ...
                       'por el camino caro. Arreglalo o quitalo.' ] , strjoin( bad , ', ' ) );
      end
    end

    function stale = warnStaleEvents( cls , defs , valid )   %(Hidden: lo llama
      %loadobj; con inputs explicitos para poder testearlo sin fabricar un .mat
      %de una version vieja de una clase)
      stale = {};
      for n = fieldnames( defs ).', n = n{1};
        bad = setdiff( fieldnames( defs.(n).events ).' , valid );
        if ~isempty( bad ), stale{end+1} = sprintf( '%s (%s)' , n , strjoin( bad , ',' ) ); end  %#ok<AGROW>
      end
      if ~isempty( stale )
        warning('cachedOwner:staleEvents', ...
                [ '%s: definiciones cargadas de un .mat declaran eventos que ya no estan en ' ...
                  'CP_EVENTS -- %s. Esas CPs NO reaccionaran a las ediciones (pareceran ' ...
                  'insensibles) y serviran valores rancios. Re-registralas tras cargar ' ...
                  '(obj.Define) o refresca las definiciones en el loadobj de tu clase.' ] , ...
                cls , strjoin( stale , '; ' ) );
      end
    end
  end

  %% ==================================================== EVENTOS
  methods (Access = protected)
    function obj = fire( obj , fired , varargin )
      %disparar UN evento de edicion: fired = nombres ESPECIFICO -> GENERAL,
      %varargin = argumentos del evento semantico (los consume el handler
      %incremental, uno por edit). Llamar en CADA metodo que edite datos.
      %COW: handle nuevo; cada entrada sobrevive / queda pendiente / cae.
      if ischar( fired ), fired = { fired }; end
      %VALIDAR LOS NOMBRES, y no es paranoia: una errata aqui no falla nunca --
      %simplemente ninguna CP reconoce el evento, las que solo declaraban ese
      %SOBREVIVEN, y se sirve un valor rancio en silencio (auditado: dbl=0 con
      %N=5). Es el unico fallo del patron que el usuario no puede ver, asi que
      %se paga un ismember por EDICION (no por lectura) para cerrarlo.
      bad = fired( ~ismember( fired , obj.CP_EVENTS ) );
      if ~isempty( bad )
        error('cachedOwner:unknownEvent', ...
              'fire: evento desconocido "%s" (CP_EVENTS de %s: %s).', ...
              bad{1} , class( obj ) , strjoin( obj.CP_EVENTS , ', ' ) );
      end
      args = varargin;
      cOld = obj.liveCache();
      if isempty( cOld ), obj.CACHE = cacheHandle(); return; end
      edit = struct( 'fired' , { fired } , 'args' , { args } );
      %COW SIN REHACER EL STORE: se CLONA y solo se tocan las claves AFECTADAS.
      %Antes se reconstruia campo a campo (dos llamadas a metodo de handle por
      %clave, tambien para las INSENSIBLES), asi que el coste de UNA edicion
      %crecia con el numero de CPs cacheadas -- medido: Add con 63 CPs = 12 ms
      %frente a 0.5 ms con 4, y el rebuild aislado es x10..x23 mas lento que el
      %clone. Ahora la clave que sobrevive cuesta CERO. El orden de los campos
      %se conserva (clone lo hereda; rmfield y sobreescribir no lo alteran).
      cNew  = cOld.clone();
      dbgOn = obj.DEBUG;                             %las listas solo alimentan la traza
      D     = obj.DEFS;                              %HOISTED, y no es cosmetica: leer
                                                     %la propiedad DOS veces por clave
                                                     %dominaba la edicion (medido: 180
                                                     %us/clave -> el bucle entero pasa
                                                     %de 10.9 ms a decimas con 60 CPs)
      %DEPENDENCIAS: `chg` reune las claves que NO sobreviven, para poder tirar
      %despues a las que se calcularon A PARTIR de ellas. Solo se lleva la cuenta
      %si alguna CP compone -- en una clase que no lo haga (cart, route, msh hoy)
      %esto es un booleano falso y ni el bucle ni la pasada de abajo existen.
      trackDeps = cOld.hasDeps;  chg = {};
      pend = {};  drop = {};  surv = {};
      for k = cOld.keys(), key = k{1};
        if ~isfield( D , key )                       % huerfana: no hay definicion
          cNew.remove( key );
          if dbgOn, drop{end+1} = key; end                               %#ok<AGROW>
          continue;
        end
        d   = D.( key );
        ev  = d.events;
        hit = fired( isfield( ev , fired ) );        % declarados, EN EL ORDEN de fired
        if isempty( hit )                            % insensible: SOBREVIVE tal cual
          if dbgOn, surv{end+1} = key; end                               %#ok<AGROW>
          continue;                                  % ya esta en cNew: NO se toca
        elseif isempty( ev.( hit{1} ) )              % decide el MAS ESPECIFICO; [] = CAE
          cNew.remove( key );
          cNew.bump( key , 4 );
          if dbgOn, drop{end+1} = key; end                               %#ok<AGROW>
          if trackDeps, chg{end+1} = key; end                            %#ok<AGROW>
          continue;
        else                                         % con handler: queda PENDIENTE
          e = cOld.entry( key );
          if strcmp( e.state , 'pending' ), L = e.log; else, L = {}; end
          %EL LOG SOLO GUARDA FIDELIDAD SI EL REPLAY INCREMENTAL SIGUE VIVO:
          %esta clave declara un handler incremental, este edit dispara su
          %evento CON args, y el log era todo-incremental. En cuanto algo de
          %eso falla, el incremental muere PARA SIEMPRE (exige que TODOS los
          %edits lo sean) y al replay solo le importa LA UNION de eventos ->
          %el log entero se funde en UN marcador absoluto (sin el evento
          %incremental en su `fired`, para que el replay ni lo intente).
          %INVARIANTE: un log es o todo edits incrementales (con args) o UN
          %marcador sin args -- `purelog` se decide mirando la ultima entrada.
          %TECHO DEL LOG: un log incremental crece una entrada por edicion, y
          %cada evento lo copia entero al handle nuevo -> O(n^2) en una racha de
          %ediciones sin lecturas (auditado: 200 AddPoint = log de 200). Al
          %llegar al techo se colapsa a marcador, que SIEMPRE es correcto (el
          %replay degrada a sync absoluto o recompute): se pierde el atajo, a
          %cambio de memoria O(1) por clave y copia O(n) en total.
          %EL ATAJO ES DEL EVENTO QUE GANA, NO DEL LOTE. Antes bastaba con que el
          %evento incremental estuviera EN EL LOTE, asi que un incremental
          %declarado en el evento GENERAL se ponia por delante del handler
          %absoluto del ESPECIFICO, que no llegaba a correr NUNCA (comprobado:
          %especifico @(v,o) 999 y general @(v,o,x) v+x daba 105, no 999). La
          %especificidad ya mandaba tres lineas mas arriba para el [] -- un []
          %especifico no lo pisa el handler de un general --; aqui no, y esa
          %incoherencia era el fallo. Decide hit{1} y punto.
          ie      = d.inc;                    %precomputado en Define
          incEdit = ~isempty( ie ) && ~isempty( args ) && strcmp( hit{1} , ie );
          purelog = isempty( L ) || ~isempty( L{end}.args );
          if incEdit && purelog && numel( L ) >= cachedOwner.MAXLOG
            incEdit = false;                   %techo alcanzado: al marcador
          end
          if incEdit && purelog
            L{end+1} = edit;
          else
            u = {};
            for i = 1:numel( L ), u = [ u , L{i}.fired ]; end            %#ok<AGROW>
            u = unique( [ u , fired ] , 'stable' );
            if ~isempty( ie ), u = u( ~strcmp( u , ie ) ); end
            L = { struct( 'fired' , { u } , 'args' , {{}} ) };
          end
          cNew.setPending( key , e.value , L , e.deps );
          if dbgOn, pend{end+1} = key; end                               %#ok<AGROW>
          if trackDeps, chg{end+1} = key; end                            %#ok<AGROW>
        end
      end
      %SEGUNDA PASADA: las CPs COMPUESTAS. Una CP que se calculo leyendo otra no
      %puede sobrevivir a un evento que se lleva a esa otra por delante -- su
      %valor salio de datos que ya no valen, y sin esto sobrevivia en silencio
      %(era el agujero: el motor solo miraba los eventos que cada CP DECLARA, y
      %la compuesta normalmente no declara ninguno). La lista de dependencias es
      %transitivamente cerrada, asi que UNA pasada basta. Una dependencia rota
      %GANA a cualquier declaracion: un handler puede poner al dia por el evento,
      %pero no puede saber que ha cambiado la CP de la que sale su valor.
      if trackDeps && ~isempty( chg )
        for k = cNew.keys(), key = k{1};
          d = cNew.deps( key );
          if ~isempty( d ) && any( ismember( d , chg ) )
            cNew.remove( key );
            cNew.bump( key , 4 );
            if dbgOn
              %OJO AL surv: una compuesta llega aqui habiendo pasado ya por la
              %primera pasada, donde se la clasifico de INSENSIBLE (no declara
              %eventos) y se la anoto como superviviente. Si no se la quita, la
              %traza la saca en las DOS listas a la vez -- "caen {best} |
              %sobreviven {best}" --, que es justo lo contrario de lo que pasa.
              %Lo destapo collatz, el primer ejemplo con CPs compuestas.
              surv( strcmp( surv , key ) ) = [];
              drop{end+1} = [ key '(compuesta)' ];                        %#ok<AGROW>
            end
          end
        end
      end
      obj.CACHE = cNew;
      if dbgOn                                 %el guard evita los lst() en vano
        obj.dbg( 'EVENT {%s}: pendientes {%s} | caen {%s} | sobreviven {%s}' , ...
                 strjoin( fired ,'+') , cachedOwner.lst(pend) , ...
                 cachedOwner.lst(drop) , cachedOwner.lst(surv) );
      end
    end
  end

  %% ==================================================== MOTOR
  %PROTECTED y no private: una subclase con subsref propio (msh) tiene que poder
  %llamarlos desde el suyo. Auditado: con private, `obj.access(nm)` desde la
  %subclase moria en noSuchMethodOrField -- es decir, el camino de extension que
  %documenta la cabecera no existia.
  methods (Access = protected)
    function c = liveCache( obj )
      %el handle si esta VIVO; cacheHandle.empty si no (objeto default-
      %inicializado, a medio deserializar, o handle borrado a mano)
      c = obj.CACHE;
      if ~isempty(c) && ~isvalid(c), c = cacheHandle.empty; end
    end

    function v = access( obj , name )
      %LECTURA perezosa: HIT | replay del pendiente | MISS
      %ORDEN DELIBERADO: primero la cache, y la DEFINICION solo si hace falta
      %(extraer la entrada con sus handlers anonimos es caro y en un HIT no se
      %usa: medido, era el coste dominante del camino caliente).
      tCP = [];                                      % lo rellena solo el MISS
      c = obj.liveCache();
      %EL CONSEJERO, si esta encendido y no hay nada en vuelo: aqui es donde se
      %puede medir de verdad (ver flushAdvice). Apagado -- que es el default --
      %esto es leer una propiedad y seguir.
      if obj.ADVISE && ~isempty( c ) && c.nIn == 0, obj.flushAdvice(); end
      if ~isempty( c )
        %DEPENDENCIA ENTRE CPs: si esta lectura ocurre DENTRO del calculo de otra
        %CP, queda apuntada como dependencia suya. La pregunta es una propiedad
        %escalar del handle justamente para que el HIT no pague una llamada: en
        %el uso normal nadie esta calculando y esto es un `if` y ya.
        if c.nIn > 0, c.noteDep( name ); end
        [ hit , e ] = c.tryGet( name );              % UN lookup: el HIT es el camino caliente
        if hit && strcmp( e.state , 'fresh' )
          if obj.DEBUG, obj.dbg( 'HIT   "%s"' , name ); end
          c.bump( name , 2 );
          v = e.value;  return;
        end
      end

      %SOLO SE LLEGA AQUI SIN DEFINICION POR UNA PUERTA: un compute que hace
      %o.access('otra') sobre una CP que ya no existe -- subsref comprueba DEFS
      %antes de llamar aqui. Pasa tras un RemoveCP: la compuesta cae (bien) y al
      %recalcularse busca lo que se llevaron. Sin esto el usuario recibia un
      %"Unrecognized field name" de las tripas de MATLAB. El isfield solo se
      %paga en el camino FRIO, nunca en un HIT.
      if ~isfield( obj.DEFS , name )
        error('cachedOwner:noCP', ...
              [ 'no hay ninguna CP "%s" en %s, y algo la ha pedido con access(): casi seguro es ' ...
                'el compute de otra CP que salia de ella y a la que le han hecho RemoveCP. ' ...
                'Vuelve a definirla, o cambia el compute que la usa.' ] , name , class( obj ) );
      end
      r = obj.DEFS.( name );                         % a partir de aqui SI se usa
      if isempty( c ), v = obj.runCompute( r , name ); return; end  % sin handle: computa sin cachear

      c.enter( name );                               % guarda de ciclos + captura
      done = onCleanup( @() c.leave( name ) );  %#ok<NASGU>   (libera aunque falle)
      if hit
        c.bump( name , 3 );
        v = obj.replay( name , r , e.value , e.log );
        %EL REPLAY CONSERVA LAS DEPS VIEJAS. El valor puesto al dia SIGUE
        %saliendo de lo que leyo su compute: el handler incremental no puede
        %leer nada, y el sync puede no releer la CP de la que salio el valor,
        %asi que takeDeps() a secas las perdia -- y una compuesta CON handler
        %sobrevivia a la caida de su base tras el PRIMER replay, sirviendo
        %rancio en silencio (comprobado: servia 305 con la verdad en 905).
        %Conservar de mas es el lado seguro: como mucho, una caida de sobra.
        d = unique( [ e.deps , c.takeDeps() ] , 'stable' );
      else
        %EL CONSEJERO: cronometrar el compute sale GRATIS -- ya se esta
        %ejecutando -- y es la unica forma de saber si esta CP se sostiene CON
        %LOS DATOS DE AHORA. El tic/toc cuesta ~1 us y solo se paga en el
        %camino FRIO, nunca en un HIT. Ver noteCost.
        c.bump( name , 1 );
        tCP = tic;
        v   = obj.runCompute( r , name );
        tCP = toc( tCP );
        obj.dbg( 'MISS  "%s" -> calculado' , name );
        d = c.takeDeps();
      end
      cachedOwner.chkHandleValue( class( obj ) , name , v );   % solo en el camino FRIO
      c.setFresh( name , v , d );                    % mutacion compartida BENIGNA
      %EL CONSEJO VA AQUI, DESPUES DE DEPOSITAR, y no es cosmetica: para saber
      %cuanto cuesta SERVIR la clave hay que leerla, y leerla antes de que este
      %fresca la haria entrar por el camino frio -- o sea, reentrar en la clave
      %que se esta resolviendo, que es exactamente lo que la guarda de ciclos
      %corta (comprobado: cacheHandle:cycle). Ya fresca, esa lectura es un HIT
      %limpio, que ni siquiera llama a enter.
      if ~isempty( tCP ), obj.noteCost( name , tCP ); end
    end

    function v = recompute( obj , name )
      %o.<nombre>_ : recalculo forzado, descarta valor y log
      r = obj.DEFS.( name );
      c = obj.liveCache();
      if ~isempty( c )
        c.enter( name );
        done = onCleanup( @() c.leave( name ) );  %#ok<NASGU>
      end
      if ~isempty( c ), c.bump( name , 1 ); end   % un recalculo es un recalculo
      v = obj.runCompute( r , name );
      obj.dbg( 'RECMP "%s" -> recalculado a la fuerza' , name );
      if ~isempty( c )
        cachedOwner.chkHandleValue( class( obj ) , name , v );
        dp = c.takeDeps();
        %el valor cambia SIN pasar por un evento, asi que quien lo hubiera usado
        %para calcularse queda rancio: se tira. Solo cuesta si alguien compone.
        gone = c.dropDependents( name );
        if ~isempty( gone )
          obj.dbg( 'RECMP "%s" -> tiradas las CPs que salian de el: %s' , name , cachedOwner.lst( gone ) );
        end
        c.setFresh( name , v , dp );
      end
    end

    function noteCost( obj , name , t )
      %EL CONSEJERO (guarda barata, NO cambia el comportamiento). Cada MISS
      %aporta una medida del coste del compute; cuando una CP acumula
      %suficientes y TODAS quedan por debajo del peaje de servirla, se avisa una
      %vez: esa CP no se sostiene y deberia ser un metodo.
      %
      %POR QUE ESTO NO PUEDE SER UNA DECISION AUTOMATICA, y si un aviso: cuando
      %el motor llega aqui el peaje YA esta pagado (el despacho, buscar la
      %definicion, montar la pila). No guardar el valor no lo devuelve; solo
      %hace que la siguiente lectura lo pague otra vez. Lo unico que arregla una
      %CP barata es que NO SEA una CP, y eso solo lo puede decidir quien escribe
      %la clase. Misma filosofia que noEvents y handlerFailed: avisar, nunca
      %decidir en silencio.
      %
      %SE ADAPTA AL TAMANO SOLO: se remide en cada MISS, y cada edicion produce
      %un MISS nuevo. La misma CP callara en la malla de 10 millones de caras y
      %avisara en la de 100 (medido con collatz: el veredicto cambia de signo
      %entre N=100 y N=1000).
      if ~obj.ADVISE, return; end
      k = [ class( obj ) '_' name ];
      n = cachedOwner.costMemo( k , t );
      %CADA N MISSES, no solo a los N: los datos crecen y menguan, y con una
      %sola evaluacion por sesion el veredicto se quedaria congelado en el
      %tamano que hubiera la primera vez.
      if mod( n , cachedOwner.ADVISE_N ) ~= 0, return; end
      %...salvo las que YA se ganaron el puesto de calle: a esas no se las
      %vuelve a medir nunca, porque medir cuesta un compute entero y son
      %justamente las caras. Las baratas se pueden remedir sin remordimiento.
      if cachedOwner.adviseSettled( k ), return; end
      %AQUI SOLO SE APUNTA. Medir DENTRO del camino vivo no sale: se esta dentro
      %del enter/onCleanup de esta misma clave y las tomas salen contaminadas --
      %comprobado, el consejero se tragaba casos que ProfileCPs caza con margen
      %7x. El veredicto se emite en la proxima lectura, cuando no haya nada en
      %vuelo y las condiciones sean las mismas que las de ProfileCPs.
      cachedOwner.advicePending( k , name );
    end

    function flushAdvice( obj )
      %emite los veredictos apuntados, YA FUERA de cualquier resolucion en
      %curso: es el sitio donde medir significa algo (ver noteCost). Lo llama
      %access al empezar, y solo cuando ADVISE esta encendido y no hay nada en
      %vuelo. El guard de reentrada es obligatorio: medir LEE la CP, y leerla
      %volveria a entrar aqui.
      persistent BUSY
      if isempty( BUSY ), BUSY = false; end
      if BUSY, return; end
      [ k , name ] = cachedOwner.advicePending();
      if isempty( k ), return; end
      %el try/catch de abajo ya garantiza que se llega a apagar BUSY: un
      %onCleanup aqui sobra Y ADEMAS revienta (su destructor asignaria sobre el
      %workspace de una funcion que ya se fue: MATLAB:class:DestructorError)
      BUSY = true;
      try
        cf   = obj.DEFS.( name ).compute;
        work = cachedOwner.timeAdaptive( @() cf( obj ) );  % el compute PELADO
        pj   = cachedOwner.hitToll( obj , name );          % lo que cuesta SERVIRLA
        %CON MARGEN, y no es timidez: un veredicto que sale a cara o cruz es
        %peor que no darlo (medido: cart.total llego a dar 157 us contra 186,
        %un 15%, y avisaba unas veces si y otras no). Solo se opina cuando la CP
        %pierde de CALLE; para el resto esta ProfileCPs, que mide todo.
        if work * cachedOwner.ADVISE_MARGIN < pj
          cachedOwner.warnOnce( 'cachedOwner:notWorthCaching' , k , ...
            [ 'la CP "%s" de %s no se sostiene con estos datos: calcularla cuesta %s y SERVIRLA ' ...
              'cacheada cuesta %s. Y no vale con dejar de guardarla: cuando el motor puede medir ' ...
              'esto el peaje del despacho YA esta pagado, asi que lo que sobra es la CP entera. ' ...
              'Hazla un metodo normal. Con datos mas grandes el veredicto se da la vuelta y este ' ...
              'aviso no saldra; ProfileCPs(o) lo mide todo. Silenciable: ' ...
              'warning(''off'',''cachedOwner:notWorthCaching''), o o.ADVISE = false.' ] , ...
            name , class( obj ) , cachedOwner.fmtT( work ) , cachedOwner.fmtT( pj ) );
        elseif pj * cachedOwner.ADVISE_MARGIN < work
          cachedOwner.adviseSettled( k , true );   % se gana el puesto: caso cerrado
        end
      catch
        %opinar no puede romper una lectura, nunca
      end
      BUSY = false;
    end

    function v = truthOf( obj , r , name , c )
      %el valor calculado DESDE CERO y SIN guardarlo, con la misma guarda de
      %ciclos que una lectura normal (lo usa VerifyCPs: verificar no puede
      %cambiar lo que verifica)
      if ~isempty( c )
        c.enter( name );
        done = onCleanup( @() c.leave( name ) );  %#ok<NASGU>
      end
      v = obj.runCompute( r , name );
    end

    function v = runCompute( obj , r , name )
      %llama al compute Y TRADUCE el fallo mas desconcertante del patron. Dentro
      %de un metodo de la clase MATLAB NO pasa por subsref -- y una funcion
      %anonima escrita ahi hereda ese contexto --, asi que un compute definido
      %DENTRO de la clase que lea otra CP (@(o) o.otra*2) muere en un
      %"Unrecognized method, property, or field 'otra'" que MIENTE: la CP existe.
      %Desde FUERA de la clase la misma linea funciona. Aqui se detecta que el
      %nombre que falta es una CP y se dice que hacer.
      try
        v = r.compute( obj );
      catch ME
        if strcmp( ME.identifier , 'MATLAB:noSuchMethodOrField' )
          bad = regexp( ME.message , '''(\w+)''' , 'tokens' , 'once' );
          if ~isempty( bad ) && isfield( obj.DEFS , bad{1} ) && ~strcmp( bad{1} , name )
            error('cachedOwner:cpFromInside', ...
                  [ 'el compute de "%s" ha intentado leer la CP "%s" con un punto, y no se puede: ' ...
                    'MATLAB NO pasa por subsref dentro de un metodo de la clase (y una funcion ' ...
                    'anonima escrita ahi hereda ese contexto), asi que "%s" no se resuelve como CP. ' ...
                    'Desde DENTRO de la clase leela con o.access(''%s''); desde fuera, o.%s vale. ' ...
                    'Se registrara como dependencia: si "%s" cae, "%s" cae con ella.' ] , ...
                  name , bad{1} , bad{1} , bad{1} , bad{1} , bad{1} , name );
          end
        end
        rethrow( ME );
      end
    end

    function v = replay( obj , name , r , v0 , L )
      %REPLAY PEREZOSO con CASCADA: incremental -> sync absoluto -> recompute.
      %El incremental es el handler de 3+ args (si TODOS los edits del log
      %disparan su evento, con args); los sync absolutos se prueban EN EL ORDEN
      %DE DECLARACION de Define. Cada escalon con try/catch: un handler que
      %falla DEGRADA al siguiente, nunca rompe la lectura -- pero NO EN SILENCIO:
      %avisa una vez, lleva la cuenta y sale marcado en la tabla como ATAJO ROTO
      %(cachedOwner:handlerFailed). Degradar es correcto; que no se note, no: la
      %clase se quedaba dando el valor bueno por el camino caro para siempre.
      ie = r.inc;  ih = r.incH;                     %precomputados en Define
      useInc = ~isempty( ie ) && ...
               all( cellfun( @(ed) ~isempty( ed.args ) && any( strcmp( ie , ed.fired ) ) , L ) );
      %EL HANDLER VE EL DUENO DE AHORA, NO EL DE SU EDICION. Con UNA edicion
      %pendiente coinciden y no hay problema; con VARIAS, un handler que lea el
      %dueno da un valor MAL en silencio (ver la cabecera). Define ya rechaza el
      %handler anonimo que lo hace, asi que aqui solo llegan los que NO se pueden
      %inspeccionar -- opacos o venidos de un .mat viejo --: se renuncia al
      %atajo, que degrada a sync absoluto o recompute, siempre correctos.
      if useInc && numel( L ) > 1 && ~strcmp( r.incOwn , 'no' )
        useInc = false;
        obj.dbg( 'RPLAY "%s" -> atajo incremental descartado: %d ediciones y el handler puede leer el dueno' , ...
                 name , numel( L ) );
        cachedOwner.warnIncOwner( class( obj ) , name , r.incOwn , numel( L ) );
      end
      if useInc
        try
          v = v0;
          for i = 1:numel( L ), v = ih( v , obj , L{i}.args{:} ); end
          obj.dbg( 'RPLAY "%s" -> %d incremental(es) sobre el valor viejo' , name , numel(L) );
          return;
        catch ME
          obj.dbg( 'RPLAY "%s" -> incremental fallo (%s), probando sync absoluto' , name , ME.identifier );
          cachedOwner.warnHandlerFailed( class( obj ) , name , ie , ME );
        end
      end
      fired = {};
      for i = 1:numel( L ), fired = [ fired , L{i}.fired ]; end   %#ok<AGROW>
      for o = fieldnames( r.events ).'                            % ORDEN DE DECLARACION
        h = r.events.( o{1} );
        if isempty( h ) || cachedOwner.isInc( h ), continue; end  % solo los ABSOLUTOS
        if any( strcmp( o{1} , fired ) )
          try
            v = h( v0 , obj );
            obj.dbg( 'RPLAY "%s" -> sync absoluto via %s' , name , o{1} );
            return;
          catch ME
            obj.dbg( 'RPLAY "%s" -> handler de %s fallo (%s)' , name , o{1} , ME.identifier );
            cachedOwner.warnHandlerFailed( class( obj ) , name , o{1} , ME );
          end
        end
      end
      v = obj.runCompute( r , name );
      obj.dbg( 'RPLAY "%s" -> sin handler viable, recalculado' , name );
    end

    function obj = refreshDefs( obj )
      %gancho de EVOLUCION para las clases con .mat vivos: loadobj lo llama al
      %cargar, DESPUES de revivir la cache y ANTES del aviso de eventos
      %huerfanos. Sobrescribelo y re-registra ahi (con Define) las definiciones
      %de tu clase: el CODIGO, y no el fichero, vuelve a ser la fuente de la
      %logica -- sin esto, un objeto cargado usa los compute y handlers DEL
      %.mat para siempre. Lo que no re-registres (las CPs ad hoc que el usuario
      %definiera sobre el objeto) conserva lo del fichero. El default no hace
      %nada: el fichero manda, que es lo que siempre hizo.
    end

    function dbg( obj , fmt , varargin )
      if obj.DEBUG, fprintf( '  [%s] %s\n' , class( obj ) , sprintf( fmt , varargin{:} ) ); end
    end
  end

  properties (Constant, Access = private)
    MAXLOG    = 32   % techo del log incremental (ver fire)
    ADVISE_N  = 5    % misses que hay que acumular antes de opinar. Con menos, el
                     % ruido de una sola toma (y el calentamiento de la primera)
                     % daria consejos a cara o cruz
    ADVISE_MARGIN = 2  % cuanto tiene que perder la CP para que se opine: el
                     % compute tiene que costar menos de la MITAD que servirla.
                     % Medir dentro del camino vivo es ruidoso y un veredicto
                     % al 15% no es un veredicto (ver noteCost)
  end

  methods (Static, Access = private)
    function tf = isInc( h )
      %la FIRMA es la politica: 3+ argumentos (o varargin equivalente) = handler
      %INCREMENTAL @(v,o,args...); 2 = sync absoluto @(v,o)
      n  = nargin( h );
      tf = ( n >= 3 ) || ( n <= -3 );
    end

    function [ en , h ] = incOf( ev )
      %el PRIMER evento declarado con handler incremental ('' si no hay).
      %YA NO esta en el camino caliente: Define lo precomputa en DEFS.inc/.incH.
      %Sobrevive porque loadobj lo necesita para migrar los .mat guardados antes
      %de ese cambio, que traen DEFS sin esos campos.
      en = '';  h = [];
      for f = fieldnames( ev ).', f = f{1};
        g = ev.( f );
        if ~isempty( g ) && cachedOwner.isInc( g )
          en = f;  h = g;  return;
        end
      end
    end

    function n = csCount( v , s , ctx , dflt )
      %cuantas salidas da subsref(v,s). NO se calcula: SE LE PREGUNTA AL VALOR,
      %que es quien conoce sus reglas -- y funciona igual con un struct array
      %(numel de la lista), con un cell {:} y con un objeto que tenga su propio
      %numArgumentsFromSubscript. Si el valor no sabe responder, el default.
      try
        n = numArgumentsFromSubscript( v , s , ctx );
        if ~( isnumeric(n) && isscalar(n) && n >= 0 ), n = dflt; end
      catch
        n = dflt;
      end
    end

    function b = stripLits( b )
      %vacia los literales de texto ('...' y "...") para que un nombre que solo
      %aparece dentro de un mensaje no cuente como uso. La comilla simple es
      %ambigua (transpuesta): abre literal solo si NO viene detras de algo que
      %se pueda transponer -- identificador, cierre, punto u otra comilla.
      keep = true( 1 , numel( b ) );  q = 0;  i = 1;
      while i <= numel( b )
        c = b(i);
        if q == 0
          if c == '"'
            q = 2;  keep(i) = false;
          elseif c == '''' && ~( i > 1 && ( isletter( b(i-1) ) || ...
                                            any( b(i-1) == '0123456789_)]}.''' ) ) )
            q = 1;  keep(i) = false;
          end
        else
          keep(i) = false;
          if q == 1 && c == ''''
            if i < numel( b ) && b(i+1) == ''''      %'' escapada: sigue dentro
              keep(i+1) = false;  i = i + 2;  continue;
            end
            q = 0;
          elseif q == 2 && c == '"'
            q = 0;
          end
        end
        i = i + 1;
      end
      b = b( keep );
    end

    function n = warnOnce( id , key , fmt , varargin )
      %MEMO de avisos, en tres modos:
      %   warnOnce( id , key , fmt , ... )  avisa la PRIMERA vez y CUENTA todas
      %   n = warnOnce( id , key )          cuantas van (0 = ninguna): lo
      %                                     pregunta la tabla, que es donde el
      %                                     aviso no se pierde
      %   warnOnce( '' , key )              olvida lo anotado sobre esa clave
      %                                     (Define, al registrar: definicion
      %                                     nueva, historial limpio)
      %Los avisos de este patron denuncian un DISENO (un valor handle, un
      %handler que lee el dueno, un atajo roto), no una llamada concreta:
      %repetirlos en cada lectura seria ruido y acabarian silenciados. Pero la
      %CUENTA si importa -- fallar una vez y fallar siempre no es lo mismo --,
      %asi que se lleva aunque no se avise.
      %OJO: la clave es CLASE.cp, no el objeto. Dos objetos de la misma clase
      %comparten aviso; es la misma aproximacion que ya tenia handleValue.
      persistent seen
      if isempty( seen ), seen = struct(); end
      if isempty( id )                               % --- olvidar
        f  = fieldnames( seen );
        f  = f( startsWith( f , matlab.lang.makeValidName( [ key '_' ] ) ) );
        if ~isempty( f ), seen = rmfield( seen , f ); end
        n = 0;  return;
      end
      k = matlab.lang.makeValidName( [ key '_' id ] );
      if isfield( seen , k ), n = seen.( k ); else, n = 0; end
      if nargin < 3, return; end                     % --- solo preguntaba
      seen.( k ) = n + 1;
      if n == 0, warning( id , fmt , varargin{:} ); end
      n = n + 1;
    end

    function warnHandlerFailed( cls , name , ev , ME )
      %LA CASCADA SE TRAGA EL ERROR A PROPOSITO (que un handler falle no puede
      %romper una lectura), y ese era el precio: un atajo roto no se notaba
      %NUNCA -- la clase seguia dando el valor correcto, por el camino caro,
      %para siempre. Este aviso es la unica senal, y la cuenta que lleva el memo
      %alimenta la marca de la tabla.
      m = regexprep( strtrim( ME.message ) , '\s+' , ' ' );
      if numel( m ) > 140, m = [ m(1:140) '...' ]; end
      cachedOwner.warnOnce( 'cachedOwner:handlerFailed' , [ cls '_' name '_' ev ] , ...
        [ 'el handler de %s de la CP "%s" (%s) ha FALLADO y se ha degradado: el valor sigue siendo ' ...
          'correcto -- se recalcula por el camino caro -- pero el atajo esta roto y ya no se usa. ' ...
          'Causa: %s: %s. La cascada se traga el error para no romper la lectura, asi que este ' ...
          'aviso (y la marca de la tabla) es la unica senal. Silenciable: ' ...
          'warning(''off'',''cachedOwner:handlerFailed'').' ] , ...
        ev , name , cls , ME.identifier , m );
    end

    function chkHandleValue( cls , name , v )
      %EL UNICO CASO EN QUE LA CACHE DEJA DE SER SOLO RENDIMIENTO. Si el valor
      %es un objeto HANDLE, el COW no aisla: fire copia la REFERENCIA al handle
      %nuevo, asi que un handler (o cualquiera) que lo mute IN PLACE corrompe la
      %cache de las copias hermanas -- comprobado: la hermana, con un item,
      %devolvia el valor de la que tenia dos. No se puede impedir desde aqui (el
      %valor lo produce el dominio), asi que se AVISA, una vez por clase.clave y
      %solo en el camino frio. Silenciable: warning('off','cachedOwner:handleValue').
      %Deteccion best-effort: ve un handle desnudo, no uno metido en un cell o
      %en un campo de struct.
      if ~isa( v , 'handle' ), return; end
      cachedOwner.warnOnce( 'cachedOwner:handleValue' , [ cls '_' name ] , ...
              [ 'la CP "%s" de %s cachea un objeto handle (%s): las copias comparten esa ' ...
                'referencia, asi que mutarlo in place corrompe la cache de las hermanas. ' ...
                'Devuelve un valor, o copia el handle en el compute y en el handler.' ] , ...
              name , cls , class( v ) );
    end

    function errNotScalar( cls , what )
      error('cachedOwner:notScalar', ...
            [ '%s es una clase ESCALAR y no admite arrays (intentaste %s). El motor de ' ...
              'cache trabaja objeto a objeto: en un array no hay a que clave indexar. ' ...
              'Guarda los objetos en un CELL, o si ya tienes el array por otra via, lee ' ...
              'con arrayfun( @(x) x.<cp> , A ).' ] , cls , what );
    end

    function errCat( args )
      cls = 'cachedOwner';
      for i = 1:numel( args )
        if isa( args{i} , 'cachedOwner' ), cls = class( args{i} ); break; end
      end
      cachedOwner.errNotScalar( cls , 'concatenar' );
    end

    function s = lst( c )
      if isempty( c ), s = '-'; else, s = strjoin( c , ', ' ); end
    end

    function showVerify( cls , R )
      %el informe de VerifyCPs. Mismo criterio que la tabla: el marcador solo
      %aparece si hay algo que marcar
      if isempty( R )
        fprintf( '  VerifyCPs de %s: ninguna CP tiene valor cacheado (nada que verificar).\n\n' , cls );
        return;
      end
      bad = ~[ R.ok ];
      w   = max( cellfun( @numel , { R.name } ) ) + 1;
      fprintf( '  VerifyCPs de %s: %d CP(s) con valor, %d discrepa(n)\n' , cls , numel(R) , sum(bad) );
      for i = 1:numel( R )
        if R(i).ok
          fprintf( '    %-*s %-9s OK\n' , w , R(i).name , R(i).state );
        elseif ~isempty( R(i).err )
          fprintf( '  ! %-*s %-9s ERROR al verificar -- %s\n' , w , R(i).name , R(i).state , R(i).err );
        else
          ex = '';
          if isa( R(i).served , 'handle' )
            ex = '  (valor handle: se compara por IDENTIDAD, la discrepancia puede ser falsa)';
          end
          fprintf( '  ! %-*s %-9s sirve %s | verdad %s%s\n' , w , R(i).name , R(i).state , ...
                   cachedOwner.fmt( R(i).served ) , cachedOwner.fmt( R(i).truth ) , ex );
        end
      end
      if any( bad )
        fprintf( [ '  Una discrepancia = la cache sirve algo distinto de lo que da su compute:\n' ...
                   '  o un handler MIENTE, o a esa CP le falta declarar el evento que la afecta.\n' ] );
      end
      fprintf( '\n' );
    end

    function t = timeAdaptive( f )
      %cronometra `f` gastando lo justo. timeit es lo preciso, pero repite la
      %llamada un monton de veces: sobre un compute de 5 segundos (un bvh, sin
      %ir mas lejos) eso es inaceptable. Asi que se prueba UNA vez y solo se
      %afina con timeit lo que sale barato.
      %LA PRIMERA LLAMADA MIENTE, y no poco: incluye el CALENTAMIENTO de MATLAB
      %(JIT del camino nuevo, carga de la clase del valor...). Medido: la
      %factorizacion de linsys daba 45 ms en la primera toma y 0.7 ms en la
      %segunda -- un factor 60 --, y con eso ProfileCPs llegaba a soltar el
      %veredicto CONTRARIO al correcto. Por eso, cuando la primera toma sale
      %"cara" (y por tanto no se va a afinar con timeit) se repite UNA vez mas
      %ya en caliente y se toma la mejor. timeit no necesita este cuidado: ya
      %calienta por dentro antes de medir.
      t0 = tic;  f();  t = toc( t0 );
      if t < 0.02
        ok = false;
        try, t = timeit( f );  ok = true; catch, end
        if ~ok
          %timeit no siempre puede (segun el handle), y quedarse con la toma
          %CRUDA seria quedarse justo con la contaminada por el arranque en
          %frio: se repite en caliente y se toma la mejor
          t0 = tic;  f();  t = min( t , toc( t0 ) );
        end
      else
        t0 = tic;  f();  t = min( t , toc( t0 ) );
      end
    end

    function s = fmtT( t )
      %un tiempo, en la unidad que se lee bien
      if     t >= 1,     s = sprintf( '%.2f s'  , t       );
      elseif t >= 1e-3,  s = sprintf( '%.2f ms' , t*1e3   );
      else,              s = sprintf( '%.0f us' , t*1e6   );
      end
    end

    function s = fmtB( b )
      if     b >= 1e6, s = sprintf( '%.1f MB' , b/1e6 );
      elseif b >= 1e3, s = sprintf( '%.0f KB' , b/1e3 );
      else,            s = sprintf( '%d B'    , b     );
      end
    end

    function showInfo( S )
      %el informe de UNA CP (lo imprimen InfoCP y o.CP.<n>.info)
      st = S.state;
      if strcmp( st , 'pending' ), st = sprintf( 'pending (%d edicion(es) sin aplicar)' , S.pending ); end
      if S.bytes > 0, st = sprintf( '%s, %s' , st , cachedOwner.fmtB( S.bytes ) ); end
      fprintf( '\n  CP "%s" de %s\n' , S.name , S.class );
      fprintf( '    estado     %s\n' , st );
      ev = fieldnames( S.events ).';
      if isempty( ev )
        fprintf( '    eventos    ninguno: no reacciona a nada (solo MISS/HIT)\n' );
      else
        p = cellfun( @(f) sprintf( '%s -> %s' , f , S.events.( f ) ) , ev , 'UniformOutput' , false );
        fprintf( '    eventos    %s\n' , strjoin( p , ' | ' ) );
        fprintf( '               (en orden de PRIORIDAD; de una edicion decide el MAS ESPECIFICO)\n' );
      end
      fprintf( '    sale de    %s\n' , cachedOwner.lst( S.deps ) );
      if isempty( S.usedBy )
        fprintf( '    la usan    -\n' );
      else
        fprintf( '    la usan    %s   <- estas caen si esta cae\n' , cachedOwner.lst( S.usedBy ) );
      end
      n = S.counts;
      fprintf( '    historial  %d recalculo(s) | %d HIT | %d replay | %d caida(s) por evento\n' , ...
               n(1) , n(2) , n(3) , n(4) );
      if ~isnan( S.work )
        fprintf( '    coste      calcularla %s | servirla %s' , ...
                 cachedOwner.fmtT( S.work ) , cachedOwner.fmtT( S.hit ) );
        if S.work < S.hit
          fprintf( '   -> NO DEBERIA SER CP: hazla un metodo\n' );
        else
          fprintf( '   -> cachear PAGA (%.0fx)\n' , S.work / S.hit );
        end
      end
      if n(1) > 0 && n(2) == 0
        fprintf( '    OJO        se ha recalculado %d vez(ces) y NO se ha servido ni una: o la\n' , n(1) );
        fprintf( '               invalidas en cada edicion, o no te esta sirviendo de nada\n' );
      end
      fprintf( '\n' );
    end

    function showProfile( cls , R )
      %el informe de ProfileCPs
      if isempty( R )
        fprintf( '  ProfileCPs de %s: no hay CPs definidas.\n\n' , cls );
        return;
      end
      w   = max( [ cellfun( @numel , { R.name } ) , 7 ] ) + 1;
      bad = [ R.work ] < [ R.hit ] | [ R.gain ] <= 0;
      fprintf( [ '  ProfileCPs de %s: %d CP(s), medido en esta maquina.\n' ...
                 '    trabajo  el compute PELADO: lo que cuesta el calculo y nada mas\n' ...
                 '    MISS     el camino frio ENTERO por el patron (trabajo + maquinaria)\n' ...
                 '    HIT      servirla ya cacheada, por el camino real (con despacho)\n' ] , ...
               cls , numel( R ) );
      fprintf( '  %-*s %9s %10s %9s %12s %9s  %s\n' , w , 'CP' , 'trabajo' , 'MISS' , 'HIT' , ...
               'ahorro/lect' , 'memoria' , 'veredicto' );
      for i = 1:numel( R )
        if bad(i), mk = '! '; else, mk = '  '; end
        fprintf( '%s%-*s %9s %10s %9s %12s %9s  %s\n' , mk , w , R(i).name , ...
                 cachedOwner.fmtT( R(i).work ) , cachedOwner.fmtT( R(i).compute ) , ...
                 cachedOwner.fmtT( R(i).hit ) , cachedOwner.fmtT( R(i).gain ) , ...
                 cachedOwner.fmtB( R(i).bytes ) , R(i).verdict );
      end
      fprintf( [ '  DOS PREGUNTAS DISTINTAS, y conviene no mezclarlas:\n' ...
                 '   1. DEBERIA SER UNA CP?  compara TRABAJO con HIT. Si calcularla cuesta menos\n' ...
                 '      que servirla, el patron solo puede restar: hazla un metodo normal.\n' ...
                 '   2. CONVIENE GUARDARLA?  el ahorro es POR LECTURA EXTRA entre dos ediciones,\n' ...
                 '      asi que multiplicalo por las que hagas de verdad. Si lees una sola vez por\n' ...
                 '      edicion no compensa nunca: recalculas igual y cada clave cacheada engorda\n' ...
                 '      el fire de cada edicion.\n' ] );
      if any( bad )
        fprintf( '  Las marcadas con ! no se sostienen con estos datos (mira el veredicto).\n' );
      end
      fprintf( [ '  Y OJO CON EL TAMANO: esto vale para los datos que tiene el objeto AHORA. Una\n' ...
                 '  malla de 100 caras y otra de 10 millones dan veredictos opuestos -- vuelve a\n' ...
                 '  medir con datos representativos.\n\n' ] );
    end

    function s = polDesc( h )
      %la politica de UN evento, para la tabla. Son las tres que puede declarar
      %Define; la cuarta -- insensible -- es no declararlo, y por eso no aparece
      if isempty( h ),                s = 'invalida';
      elseif cachedOwner.isInc( h ),  s = 'incremental';
      else,                           s = 'sync';
      end
    end

    function s = fmt( v )
      if isnumeric(v) && isscalar(v), s = sprintf('%g',v);
      elseif ischar(v),               s = [ '''' v '''' ];
      else,                           s = sprintf('%s %s',mat2str(size(v)),class(v));
      end
    end
  end
end
