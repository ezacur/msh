classdef msh < matlab.mixin.CustomDisplay
%MSH  Contenedor de mallas (value class) con derivados cacheados perezosamente,
%     auto-invalidados/auto-actualizados mediante EVENTOS de edicion.
%
%   M = msh( V , F )              desde coordenadas + conectividad
%   M = msh( S )                  desde un struct legado (.xyz/.tri/.xyzF*/.triA*,
%                                 tambien .vertices/.faces)
%   M = msh( M0 )                 copia (value semantics; comparte cache mientras
%                                 los datos no cambien)
%   M = msh()                     malla vacia
%
%   DATOS PRINCIPALES
%     M.V        coordenadas nV x 2|3 (2 columnas = malla GENUINAMENTE 2D; la
%                dimension nsd es el numero de columnas, nunca se infiere de
%                los valores)
%     M.F        conectividad nF x k, k in {1,2,3,4}: puntos, segmentos,
%                triangulos, TETRAEDROS; mixtas 0-padded (ceros al final)
%     TIPOS OBLIGATORIOS: V es SIEMPRE double y F SIEMPRE int32 (los setters
%     validan y convierten cualquier entrada numerica); los ATTS (campos por
%     nodo/cara) conservan el tipo que les des.
%     escalares: M.nsd (2|3), M.nV, M.nF, M.ct (codigo VTK; alias M.celltype)
%     PLANAR / FLAT: una malla con nsd==3 cuyos vertices caen todos en UN
%     plano es "planar" (M.IsPlanar); si ese plano es exactamente z == 0 es
%     "flat" (M.IsFlat). Es descripcion, no tipo: nsd sigue siendo 3.
%     campos por nodo/cara: M.AddField('xyzFOO',v) / M.GetField('xyzFOO') ...
%
%   CPs (cached props, el corazon de la clase): derivados definidos por un
%   REGISTRO nombre -> { computeFcn , eventos }. Se calculan al pedirlos, se
%   guardan, y ante cada EVENTO de edicion o bien caen, o bien quedan
%   PENDIENTES de una actualizacion barata (replay perezoso en el proximo
%   acceso), o bien sobreviven intactos si el evento no les afecta.
%
%     M.bvh                     LEE (perezoso: HIT / replay / computa+guarda)
%     M.bvh_                    RECALCULA a la fuerza (descarta valor y log,
%                               computa fresco, guarda y devuelve)
%     M.bvh.frame / M.bvh_.frame   indexar dentro del valor
%     M.CP                      tabla de definiciones y estados (plano de control)
%     M.CP.bvh                  el valor (lectura perezosa, como M.bvh)
%     M.CP.bvh.delete           borra el valor (la definicion queda)
%     M = M.CP.bvh.removeProp   borra definicion y valor
%     M = M.CP.bvh.set( x )     siembra un valor a mano (aislado, COW)
%     M.CP.bvh.changeCoords     el handler del evento (invocable)
%
%     M = M.DefineCP( 'foo' , @(m) ... , evento , handler|[] , ... )
%     M = M.RemoveCP( 'foo' )
%
%   EVENTOS (de especifico a general): las ediciones disparan
%     M.V = ...    -> [changeNodeCount] [changeDim] changeCoords
%     M.F = ...    -> [changeFaceCount] changeConnectivity
%     Transform()  -> transform(T) + changeCoords
%   Para cada CP: evento no declarado = sobrevive; declarado con [] = invalida;
%   declarado con handler = queda pendiente y el handler lo actualiza en el
%   proximo acceso (transform es INCREMENTAL @(v,m,T); el resto son ABSOLUTOS
%   @(v,m), sincronizan contra la malla actual). Si un handler falla se degrada
%   en cascada: transform -> sync absoluto -> recompute.
%
%   CPs DE FABRICA (sobreescribibles con DefineCP):
%     bvh        blob BVH; transform->plegado O(1), changeCoords->refit O(n)
%     boundary   facetas del borde (MeshBoundary)
%     edges, esup, psup, esue, bbox, surfCent ([area,centroide])
%     triNormals normales por cara; transform->rotacion O(n) de las filas
%
%   CONVENCION DE NOMBRES POR CASO (tres espacios disjuntos):
%     MAYUSCULAS  = propiedades (V, F, VIZ, INFO, DEBUG y el proxy CP)
%     Capitalized = metodos publicos (Plot, Transform, Tidy, DefineCP, ...)
%     minuscula   = CPs (bvh, esup, boundary, ...); DefineCP lo exige
%   Excepciones documentadas: los contadores nsd/nV/nF/ct (M.ct = codigo VTK,
%   trivial, directo, alias M.celltype), los metodos del protocolo MATLAB
%   (subsref, disp, loadobj...) y los alias transicionales xyz/tri/celltype.
%   surfCent devuelve [area, centroide]; triNormals es la normal canonica por
%   cara (el campo NORMALS del usuario viaja aparte y Transform() lo rota);
%   normales por nodo: meshNormals(ToStruct(M),'angle'|...) o define la tuya.
%
%   QUERIES (usan/rellenan el bvh cacheado):
%     [e,cp,d,bc,F] = M.ClosestElement( P , Dmax )
%     [xyz,cell,t,rid] = M.IntersectRay( rays , MODE )
%
%   M.DEBUG = true  narra por consola la cadena de procesos: HIT/MISS/RPLAY/
%   RECMP de cada CP (con tiempos), eventos (que cae, que queda pendiente, que
%   sobrevive), Transform y queries. Tambien msh(V,F,'DEBUG',true).
%
%   OTROS: M.VIZ (preferencias de plot: defaults < VIZ < args explicitos),
%   M.INFO (metadatos libres; la textura vive en M.INFO.texture), M.Plot(...),
%   M.PlotBVH(), y delegaciones al toolbox legado: M.Tidy(...), M.Append(...),
%   M.RemoveFaces/Nodes(...). Puente: S = M.ToStruct() (habla .xyz/.tri).
%   OJO: las formas funcionales ya no despachan a la clase: plot(M) o
%   transform(M,T) caen a las funciones del path (transform "medio funciona"
%   via los alias xyz/tri pero DEVUELVE STRUCT, no msh) -- usa M.Plot() y
%   M.Transform(T).
%
% See also cacheHandle, cacheView, BVH, bvhClosestElement, bvhIntersectRay,
%          msh_CLASS_TUTORIAL.md.

  %% ------------------------------------------------------------------ DATOS
  % convencion por caso: MAYUSCULAS = propiedades (publicas y privadas),
  % Capitalized = metodos publicos, minuscula = CPs. El sufijo '_' sobre un
  % nombre de CP significa RECALCULAR (M.bvh_); el nombre desnudo LEE (M.bvh).
  properties (Access = private)
    VERTICES   = zeros(0,3)          % almacenamiento real de .V (SIEMPRE double)
    FACES      = zeros(0,3,'int32')  % almacenamiento real de .F (SIEMPRE int32)
    VATTS      = struct()            % atributos por nodo (legado .xyzNAME, sin prefijo)
    FATTS      = struct()            % atributos por cara (legado .triNAME, sin prefijo)
    cachePROPS = struct()            % REGISTRO de CPs: nombre -> struct
                                     %   .compute  @(m) valor
                                     %   .events   struct evento -> handler | []
  end

  properties   % --- datos publicos sin efecto en la cache de derivados -------
    VIZ   = struct()   % preferencias de visualizacion para Plot
    INFO  = struct()   % metadatos arbitrarios del usuario (textura: INFO.texture)
    DEBUG = false      % true -> narra cada paso: cache (HIT/MISS/RPLAY/RECMP),
                       % eventos (que cae/pende/sobrevive), Transform, queries
  end

  properties (Access = private, Transient)   % NO se serializa (save/load)
    CACHE = cacheHandle.empty    % handle a los VALORES cacheados (el "puntero")
  end

  properties (Dependent)
    V          % nV x 2|3 coordenadas (SIEMPRE double)
    F          % nF x k   conectividad 0-padded (SIEMPRE int32)
    nsd        % 2 o 3 (numero de columnas de V)
    nV         % numero de vertices
    nF         % numero de caras
    ct         % codigo VTK (directo, sin cache: trivial); contador-excepcion
               % como nsd/nV/nF (alias: M.celltype, via subsref); el resto de
               % derivados son CPs (nombre desnudo lee, sufijo '_' recalcula)
  end

  %% ====================================================== CONSTRUCCION
  methods
    function obj = msh( varargin )
      obj.cachePROPS   = msh.defaultRegistry();      % cachedProps de fabrica
      obj.CACHE = cacheHandle();                 % handle vivo desde el inicio
      if nargin == 0, return; end
      if nargin == 1 && isa( varargin{1} , 'msh' )
        obj = varargin{1};                     % copia por valor (cache compartida)
        return;
      end
      S = msh.parseInputs( varargin{:} );
      obj = obj.setFromStruct( S );
    end
  end

  methods (Static, Access = private)
    function S = parseInputs( varargin )
      if numel( varargin ) == 1
        S = varargin{1};
        if isnumeric( S ), S = struct( 'xyz' , S , 'tri' , zeros(0,3) ); end
        if ~isstruct( S )
          error('msh:input','unsupported input: use msh(nodes,faces), msh(struct) or msh(msh).');
        end
        if isfield( S , 'vertices' ), S.xyz = S.vertices;  S = rmfield( S ,'vertices'); end
        if isfield( S , 'faces' )
          if ~isempty( S.faces ) && min( S.faces(:) ) == 0, S.faces = S.faces + 1; end
          S.tri = S.faces;  S = rmfield( S ,'faces');
        end
        if ~isfield( S , 'xyz' ), error('msh:input','the struct needs .xyz (or .vertices).'); end
        if ~isfield( S , 'tri' ), S.tri = zeros(0,3); end
      else
        S = struct( 'xyz' , varargin{1} , 'tri' , varargin{2} );
        for a = 3:2:numel( varargin )           % pares nombre/valor extra
          S.( varargin{a} ) = varargin{a+1};
        end
      end
    end
  end

  %% =============================================== DOT-DISPATCH (subsref)
  % orden de despacho:  M.CP...  ->  M.<cp>_ (RECALCULA)  ->  M.<cp> (LEE)
  % -> alias legados -> builtin (props reales, metodos, encadenados).
  % NB: dentro de los metodos de la clase este subsref NO corre (regla de
  % MATLAB) -- el codigo interno usa accessCached/recomputeCached directamente.
  methods
    function varargout = subsref( M , s )
      if strcmp( s(1).type , '.' ) && ( ischar( s(1).subs ) || isstring( s(1).subs ) )
        nm = char( s(1).subs );
        if strcmp( nm , 'CP' )                 %plano de control (proxy)
          if isscalar( s ), varargout{1} = cacheView( M ); return; end
          out = M.cachedAccess( s(2:end) );
          if isempty( out ), varargout = {};
          else, varargout = out( 1:max( min( nargout , numel(out) ) , 1 ) );
          end
          return;
        end
        if numel( nm ) > 1 && nm(end) == '_' && isfield( M.cachePROPS , nm(1:end-1) )
          v = M.recomputeCached( nm(1:end-1) );          %sufijo '_': RECALCULA
          if numel( s ) > 1
            [ varargout{ 1:max( nargout , 1 ) } ] = builtin( 'subsref' , v , s(2:end) );
          else
            varargout{1} = v;
          end
          return;
        end
        if isfield( M.cachePROPS , nm )        %nombre desnudo: LEE (perezoso)
          v = M.accessCached( nm );            %(ops .delete/.set/...: solo M.CP)
          if numel( s ) > 1
            [ varargout{ 1:max( nargout , 1 ) } ] = builtin( 'subsref' , v , s(2:end) );
          else
            varargout{1} = v;
          end
          return;
        end
        %alias legados (solo LECTURA, via subsref): permiten que codigo viejo
        %y las funciones en transicion (@msh\private) lean un msh como si
        %fuera el struct legado; el saneamiento futuro los retirara
        switch nm
          case 'celltype', s(1).subs = 'ct';
          case 'xyz',      s(1).subs = 'V';
          case 'tri',      s(1).subs = 'F';
        end
      end
      [ varargout{ 1:max( nargout , 1 ) } ] = builtin( 'subsref' , M , s );
    end
  end

  %% ================================================= GET/SET + EVENTOS
  methods
    function v = get.V(obj), v = obj.VERTICES; end
    function v = get.F(obj), v = obj.FACES;    end

    function obj = set.V( obj , val )
      if ~isnumeric( val ) || ~ismatrix( val ) || ...
         ( ~isempty( val ) && size( val ,2) ~= 2 && size( val ,2) ~= 3 )
        error('msh:nodes','V must be nV x 2 or nV x 3 numeric.');
      end
      val = double( val );
      if ~isempty( obj.FACES ) && size( val ,1) < max( [ 0 ; obj.FACES(:) ] )
        error('msh:nodes','F would reference removed vertices (edit F first, or use RemoveNodes).');
      end
      fired = {};
      if size( val ,1) ~= size( obj.VERTICES ,1), fired{end+1} = 'changeNodeCount'; end
      if size( val ,2) ~= size( obj.VERTICES ,2), fired{end+1} = 'changeDim';       end
      fired{end+1} = 'changeCoords';
      obj.dbg( 'SET   V %dx%d -> %dx%d' , ...
               size(obj.VERTICES,1) , size(obj.VERTICES,2) , size(val,1) , size(val,2) );
      obj.VERTICES = val;
      obj = obj.reconcileNodeFields();
      obj = obj.fireEvents( fired , {} );
    end

    function obj = set.F( obj , val )
      val = msh.castFaces( val , size( obj.VERTICES ,1) );
      fired = {};
      if size( val ,1) ~= size( obj.FACES ,1), fired{end+1} = 'changeFaceCount'; end
      fired{end+1} = 'changeConnectivity';
      obj.dbg( 'SET   F %dx%d -> %dx%d' , ...
               size(obj.FACES,1) , size(obj.FACES,2) , size(val,1) , size(val,2) );
      obj.FACES = val;
      obj = obj.reconcileFaceFields();
      obj = obj.fireEvents( fired , {} );
    end
  end

  %% ================================================= ESCALARES + AZUCAR
  methods
    function d = get.nsd(obj), d = size( obj.VERTICES ,2); end
    function n = get.nV(obj),  n = size( obj.VERTICES ,1); end
    function n = get.nF(obj),  n = size( obj.FACES ,1);    end
    function c = get.ct(obj),  c = meshCelltype( obj.ToStruct() ); end

    function tf = IsFlat( obj )
      %malla "flat": nsd==3 y TODOS los vertices con z == 0 EXACTO (el caso
      %"casi 2D"; sigue siendo 3D: puede deformarse fuera del plano)
      X = obj.VERTICES;
      tf = size( X ,2) == 3 && ~isempty( X ) && all( X(:,3) == 0 );
    end
    function tf = IsPlanar( obj )
      %malla "planar": nsd==3 pero todos los vertices caen en UN plano
      %(cualquiera, tolerancia relativa; las flat son planares)
      tf = size( obj.VERTICES ,2) == 3 && msh.planarInfo( obj.VERTICES ) > 0;
    end
  end

  %% ========================================================== QUERIES
  methods
    function varargout = ClosestElement( M , P , varargin )
      %[e,cp,d,bc,F] = M.ClosestElement( P [, Dmax] )  -- usa el bvh cacheado
      M.dbg( 'QUERY ClosestElement: %d puntos' , size( P ,1) );
      t0 = tic;
      [ varargout{ 1:max(nargout,1) } ] = ...
          bvhClosestElement( { M.ToStruct() , M.accessCached( 'bvh' ) } , P , varargin{:} );
      M.dbg( 'QUERY ClosestElement resuelta en %.2f ms' , 1e3*toc(t0) );
    end
    function varargout = IntersectRay( M , ray , varargin )
      %[xyz,cell,t,rid] = M.IntersectRay( rays [, MODE] )
      if M.DEBUG
        mode = 'first';
        if ~isempty( varargin ) && ~isempty( varargin{end} ) && ischar( varargin{end} )
          mode = varargin{end};
        end
        if size( ray ,2) == 6, nr = size( ray ,1); else, nr = size( ray ,3); end
        M.dbg( 'QUERY IntersectRay (%s): %d rayos' , mode , nr );
      end
      t0 = tic;
      [ varargout{ 1:max(nargout,1) } ] = ...
          bvhIntersectRay( { M.ToStruct() , M.accessCached( 'bvh' ) } , ray , varargin{:} );
      M.dbg( 'QUERY IntersectRay resuelta en %.2f ms' , 1e3*toc(t0) );
    end
  end

  %% ==================================================== TRANSFORMACION
  methods
    function M = Transform( M , T )
      %Aplica T (4x4/3x4/3x3 homogenea; 2D: 3x3 homogenea 2D) a las coordenadas.
      %Dispara el evento semantico 'transform' (ademas de changeCoords): las
      %CPs con handler de transform quedan pendientes de un update INCREMENTAL
      %barato (bvh: plegado O(1); triNormals: rotacion) en vez de invalidarse.
      %Los campos NORMALS del usuario se rotan aqui mismo (via tools\transform,
      %que aplica R/det(R)^(1/3)). El metodo Capitalized ya NO ensombrece a la
      %funcion del path: llamada directa.
      M.dbg( 'TRANS Transform() sobre %d nodos' , size( M.VERTICES ,1) );
      S2 = transform( M.ToStruct() , T );
      M.VERTICES = S2.xyz;                       % mismo tamano: el evento va aparte
      if isfield( S2 , 'xyzNORMALS' ), M.VATTS.NORMALS = S2.xyzNORMALS; end
      if isfield( S2 , 'triNORMALS' ), M.FATTS.NORMALS = S2.triNORMALS; end
      M = M.fireEvents( { 'transform' , 'changeCoords' } , { T } );
    end
  end

  %% ======================================= CACHEDPROPS: REGISTRO + MOTOR
  methods
    function M = DefineCP( M , name , computeFcn , varargin )
      %M = M.DefineCP( nombre , @(m)... [, evento , handler|[] , ...] )
      %
      %   Registra (o REDEFINE, descartando el valor previo) una CP: lectura
      %   M.<nombre> (perezosa), recalculo M.<nombre>_ , control M.CP.<nombre>.
      %   Eventos no declarados no la afectan; declarados con [] la invalidan;
      %   con handler queda pendiente y se actualiza perezosamente
      %   ('transform' es incremental @(v,m,T); el resto absolutos @(v,m)).
      %
      %   El nombre debe empezar en MINUSCULA (convencion por caso: mayusculas
      %   = propiedades, Capitalized = metodos) y no pisar los contadores ni
      %   los alias legados.
      if ~isvarname( name )
        error('msh:cached','''%s'' no es un identificador MATLAB valido.', name );
      end
      if ~( name(1) >= 'a' && name(1) <= 'z' )
        error('msh:cached', ...
              'las CPs empiezan en minuscula (''%s'' no): MAYUSCULAS = propiedades, Capitalized = metodos.', name );
      end
      if ismember( name , { 'nsd','nV','nF','ct' , 'xyz','tri','celltype' , 'cached' } )
        error('msh:cached','''%s'' esta reservado (contador, alias legado o nombre historico).', name );
      end
      if ~isa( computeFcn , 'function_handle' )
        error('msh:cached','computeFcn debe ser un function handle @(m)...');
      end
      if mod( numel( varargin ) , 2 )
        error('msh:cached','los eventos van en pares nombre,handler.');
      end
      ev = struct();
      for a = 1:2:numel( varargin )
        en = msh.eventName( varargin{a} );
        h  = varargin{a+1};
        if ~( isempty( h ) || isa( h , 'function_handle' ) )
          error('msh:cached','el handler de %s debe ser un function handle o [] (invalidar).', en );
        end
        ev.( en ) = h;
      end
      existed = isfield( M.cachePROPS , name );
      M.cachePROPS.( name ) = struct( 'compute' , computeFcn , 'events' , ev );
      if existed && ~isempty( M.CACHE ) && isvalid( M.CACHE ) && M.CACHE.has( name )
        M.CACHE = M.CACHE.cloneWithout( name );   %el valor era de OTRA definicion
      end
      if existed, w = 'redefinida'; else, w = 'definida'; end
      M.dbg( 'DEF   CP ''%s'' %s (eventos: %s)' , name , w , ...
             msh.lst( fieldnames( ev ).' ) );
    end

    function M = RemoveCP( M , name )
      %M = M.RemoveCP( nombre )   borra definicion Y valor
      if ~isfield( M.cachePROPS , name )
        error('msh:cached','no hay CP ''%s''.', name );
      end
      M.cachePROPS = rmfield( M.cachePROPS , name );
      if ~isempty( M.CACHE ) && isvalid( M.CACHE ) && M.CACHE.has( name )
        M.CACHE = M.CACHE.cloneWithout( name );
      end
      M.dbg( 'DEF   CP ''%s'' eliminada (definicion y valor)' , name );
    end
  end

  methods (Hidden)
    function out = cachedAccess( M , s )
      %despachador de M.CP.<nombre>[...] (la vista cacheView). Devuelve un
      %cell de outputs (vacio para .delete).
      if ~strcmp( s(1).type , '.' )
        error('msh:cached','use M.CP.<nombre> (o lee directo con M.<nombre>).');
      end
      name = char( s(1).subs );
      if ~isfield( M.cachePROPS , name )
        error('msh:cached','no hay CP ''%s'' definida (ver DefineCP).', name );
      end
      r = M.cachePROPS.( name );
      if isscalar( s )                             % M.CP.bvh -> el valor
        out = { M.accessCached( name ) };
        return;
      end
      if strcmp( s(2).type , '.' )
        opn = char( s(2).subs );
        switch opn
          case 'delete'          %borra el VALOR (handle compartido); statement
            if numel( s ) > 2, error('msh:cached','.delete no admite mas indexacion.'); end
            c = M.CACHE;
            if ~isempty( c ) && isvalid( c ), c.remove( name ); end
            M.dbg( 'CACHE ''%s'' valor borrado (definicion intacta)' , name );
            out = {};
            return;
          case 'removeProp'      %borra definicion + valor; devuelve el msh nuevo
            if numel( s ) > 2, error('msh:cached','.removeProp no admite mas indexacion.'); end
            out = { M.RemoveCP( name ) };
            return;
          case 'set'             %siembra un valor a mano (conservador: COW, aislado)
            if numel( s ) ~= 3 || ~strcmp( s(3).type , '()' ) || numel( s(3).subs ) ~= 1
              error('msh:cached','uso: M = M.CP.%s.set( valor ).', name );
            end
            M2 = M;
            if isempty( M2.CACHE ) || ~isvalid( M2.CACHE ), M2.CACHE = cacheHandle();
            else,                                             M2.CACHE = M2.CACHE.clone();
            end
            M2.CACHE.setFresh( name , s(3).subs{1} );
            M2.dbg( 'CACHE ''%s'' valor sembrado a mano (set)' , name );
            out = { M2 };
            return;
          otherwise
            if ismember( opn , msh.eventNames() )   %handler de un evento
              if ~isfield( r.events , opn )
                error('msh:cached','''%s'' no declara el evento %s.', name , opn );
              end
              h = r.events.( opn );
              if numel( s ) >= 3 && strcmp( s(3).type , '()' )      %invocacion
                if isempty( h )
                  error('msh:cached','el evento %s de ''%s'' es [] (invalidar): no es invocable.', opn , name );
                end
                v = h( s(3).subs{:} );
                if numel( s ) > 3, v = builtin( 'subsref' , v , s(4:end) ); end
                out = { v };
              else
                out = { h };
              end
              return;
            end
        end
      end
      %cualquier otra cosa: indexar DENTRO del valor
      v = M.accessCached( name );
      out = { builtin( 'subsref' , v , s(2:end) ) };
    end

    function displayCachedView( M )
      %tabla que imprime la vista M.CP
      names = sort( fieldnames( M.cachePROPS ).' );
      if isempty( names )
        fprintf( '  (sin CPs definidas)\n\n' );  return;
      end
      w = max( cellfun( @numel , names ) );
      c = M.CACHE;
      fprintf( '  CPs (%d) -- leer M.<nombre> | recalcular M.<nombre>_ | control M.CP.<nombre> :\n' , numel( names ) );
      for n = names, n = n{1};
        r = M.cachePROPS.( n );
        if isempty( c ) || ~isvalid( c ) || ~c.has( n ), st = '(sin calcular)';
        elseif strcmp( c.state( n ) , 'fresh' ),         st = msh.fmtVal( c.value( n ) );
        else,                                            st = '(pendiente de replay)';
        end
        fprintf( '    %-*s  %s\n' , w , n , st );
        ev = fieldnames( r.events ).';
        if isempty( ev )
          fprintf( '    %-*s    eventos: (ninguno: solo MISS/HIT)\n' , w , '' );
        else
          p = cellfun( @(e) sprintf( '%s->%s' , e , msh.hDesc( r.events.(e) ) ) , ev ,'uni',0);
          fprintf( '    %-*s    eventos: %s\n' , w , '' , strjoin( p , ', ' ) );
        end
      end
      fprintf( '\n' );
    end
  end

  methods (Access = private)
    function v = accessCached( obj , name )
      r = obj.cachePROPS.( name );
      c = obj.CACHE;
      if isempty( c ) || ~isvalid( c )            % p.ej. array default-inicializado
        v = r.compute( obj );  return;
      end
      if c.has( name )
        if strcmp( c.state( name ) , 'fresh' )
          obj.dbg( 'HIT   ''%s''' , name );
          v = c.value( name );
          return;
        end
        v = obj.replayEntry( name , r , c.value( name ) , c.log( name ) );
        c.setFresh( name , v );      %resolver el pendiente: mutacion compartida benigna
        return;
      end
      t0 = tic;
      v = r.compute( obj );
      obj.dbg( 'MISS  ''%s'' -> calculado en %.2f ms' , name , 1e3*toc(t0) );
      c.setFresh( name , v );
    end

    function v = recomputeCached( obj , name )
      %M.<nombre>_ : RECALCULO forzado -- descarta valor y log de replay,
      %computa fresco desde la malla actual, guarda y devuelve. Escribe al
      %handle COMPARTIDO (mutacion benigna: los que comparten CACHE no han
      %divergido -- el COW separa al editar -- asi que el recalculo les vale).
      r = obj.cachePROPS.( name );
      t0 = tic;
      v = r.compute( obj );
      obj.dbg( 'RECMP ''%s'' -> recalculado a la fuerza en %.2f ms' , name , 1e3*toc(t0) );
      c = obj.CACHE;
      if ~isempty( c ) && isvalid( c ), c.setFresh( name , v ); end
    end

    function v = replayEntry( obj , name , r , v0 , L )
      %REPLAY PEREZOSO: aplica el log de eventos pendientes sobre el valor viejo.
      %'transform' es INCREMENTAL (recibe cada T); el resto son ABSOLUTOS
      %(sincronizan contra la malla ACTUAL, que ya lo contiene todo) -> si el
      %log mezcla, un solo sync absoluto lo subsume. Cascada ante fallos:
      %transform -> sync absoluto -> recompute.
      t0 = tic;
      allT = isfield( r.events , 'transform' ) && ~isempty( r.events.transform );
      if allT
        for i = 1:numel( L )
          if ~ismember( 'transform' , L{i}.fired ), allT = false; break; end
        end
      end
      if allT
        try
          v = v0;
          for i = 1:numel( L )
            v = r.events.transform( v , obj , L{i}.args{:} );
          end
          obj.dbg( 'RPLAY ''%s'' -> %d transform(s) incremental(es) en %.2f ms' , ...
                   name , numel( L ) , 1e3*toc(t0) );
          return;
        catch
          obj.dbg( 'RPLAY ''%s'' -> handler de transform fallo, probando sync absoluto' , name );
        end
      end
      fired = {};
      for i = 1:numel( L ), fired = [ fired , L{i}.fired ]; end        %#ok<AGROW>
      order = { 'changeNodeCount' , 'changeDim' , 'changeCoords' , ...
                'changeFaceCount' , 'changeConnectivity' };
      for o = order
        if any( strcmp( o{1} , fired ) ) && isfield( r.events , o{1} ) && ~isempty( r.events.( o{1} ) )
          try
            t0 = tic;
            v = r.events.( o{1} )( v0 , obj );
            obj.dbg( 'RPLAY ''%s'' -> sync absoluto via %s en %.2f ms' , name , o{1} , 1e3*toc(t0) );
            return;
          catch
            obj.dbg( 'RPLAY ''%s'' -> handler de %s fallo' , name , o{1} );
          end
        end
      end
      t0 = tic;
      v = r.compute( obj );
      obj.dbg( 'RPLAY ''%s'' -> sin handler viable, recalculado en %.2f ms' , name , 1e3*toc(t0) );
    end

    function obj = fireEvents( obj , fired , args )
      %UN evento de edicion: fired = nombres disparados (especifico->general),
      %args = argumentos del evento semantico (transform: {T}).
      %COW: handle nuevo; cada entrada sobrevive / queda pendiente / cae.
      cOld = obj.CACHE;
      if isempty( cOld ) || ~isvalid( cOld )
        obj.CACHE = cacheHandle();
        return;
      end
      edit = struct( 'fired' , { fired } , 'args' , { args } );
      cNew = cacheHandle();
      pend = {};  drop = {};  surv = {};
      for k = cOld.keys(), key = k{1};
        if ~isfield( obj.cachePROPS , key ), drop{end+1} = key; continue; end   %#ok<AGROW>
        ev  = obj.cachePROPS.( key ).events;
        hit = fired( isfield( ev , fired ) );
        if isempty( hit )                                %insensible: sobrevive
          cNew.setEntry( key , cOld.entry( key ) );
          surv{end+1} = key;                             %#ok<AGROW>
        elseif any( cellfun( @(e) ~isempty( ev.(e) ) , hit ) )   %handler: pendiente
          if strcmp( cOld.state( key ) , 'pending' ), L = cOld.log( key ); else, L = {}; end
          if isempty( L ) || ~isequal( L{end} , edit ), L{end+1} = edit; end
          cNew.setPending( key , cOld.value( key ) , L );
          pend{end+1} = key;                             %#ok<AGROW>
        else                                             %solo [] declarados: cae
          drop{end+1} = key;                             %#ok<AGROW>
        end
      end
      obj.CACHE = cNew;
      obj.dbg( 'EVENT {%s}: pendientes {%s} | caen {%s} | sobreviven {%s}' , ...
               strjoin( fired ,'+') , msh.lst( pend ) , msh.lst( drop ) , msh.lst( surv ) );
    end

    function dbg( obj , fmt , varargin )
      if obj.DEBUG, fprintf( '[msh] %s\n' , sprintf( fmt , varargin{:} ) ); end
    end
  end

  %% ======================================================== CAMPOS
  methods
    function obj = AddField( obj , name , val )
      %M = M.AddField('xyzFOO',v) / ('triBAR',v) / ('FOO',v) inferido por filas
      [ where , bare ] = obj.resolveField( name , size( val ,1) );
      if strcmp( where , 'node' )
        if size( val ,1) ~= obj.nV
          error('msh:field','node field "%s": %d rows but nV = %d.', name , size(val,1) , obj.nV );
        end
        obj.VATTS.( bare ) = val;
      else
        if size( val ,1) ~= obj.nF
          error('msh:field','face field "%s": %d rows but nF = %d.', name , size(val,1) , obj.nF );
        end
        obj.FATTS.( bare ) = val;
      end
    end
    function v = GetField( obj , name )
      [ where , bare ] = obj.resolveField( name , NaN );
      if strcmp( where , 'node' ), v = obj.VATTS.( bare );
      else,                        v = obj.FATTS.( bare );
      end
    end
    function obj = RmField( obj , name )
      [ where , bare ] = obj.resolveField( name , NaN );
      if strcmp( where , 'node' ), obj.VATTS = rmfield( obj.VATTS , bare );
      else,                        obj.FATTS = rmfield( obj.FATTS , bare );
      end
    end
    function tf = HasField( obj , name )
      try, obj.resolveField( name , NaN );  tf = true;
      catch, tf = false;
      end
    end
    function L = FieldNames( obj )
      %struct con .node y .face: nombres (con prefijo legado) de los campos
      L.node = strcat( 'xyz' , fieldnames( obj.VATTS ).' );
      L.face = strcat( 'tri' , fieldnames( obj.FATTS ).' );
    end
  end

  methods (Access = private)
    function [ where , bare ] = resolveField( obj , name , nrows )
      if strncmp( name , 'xyz' , 3 ) && numel( name ) > 3
        where = 'node';  bare = name(4:end);
        if ~isnan( nrows ) || isfield( obj.VATTS , bare ), return; end
      elseif strncmp( name , 'tri' , 3 ) && numel( name ) > 3
        where = 'face';  bare = name(4:end);
        if ~isnan( nrows ) || isfield( obj.FATTS , bare ), return; end
      else
        bare = name;
        if isfield( obj.VATTS , bare ), where = 'node'; return; end
        if isfield( obj.FATTS , bare ), where = 'face'; return; end
        if ~isnan( nrows )                       %inferir por numero de filas
          if nrows == obj.nV && nrows ~= obj.nF, where = 'node'; return; end
          if nrows == obj.nF && nrows ~= obj.nV, where = 'face'; return; end
          error('msh:field','cannot infer node/face for "%s": use the xyz/tri prefix.', name );
        end
      end
      error('msh:field','field "%s" not found.', name );
    end
  end

  %% ============================================ DELEGACIONES AL TOOLBOX
  methods
    function M = Tidy( M , varargin )
      M = msh( MeshTidy( M.ToStruct() , varargin{:} ) );
    end
    function M = RemoveFaces( M , idx )
      M = msh( MeshRemoveFaces( M.ToStruct() , idx ) );
    end
    function M = RemoveNodes( M , idx )
      M = msh( MeshRemoveNodes( M.ToStruct() , idx ) );
    end
    function M = Append( M , varargin )
      others = cellfun( @(x) toS(x) , varargin , 'uni' , 0 );
      M = msh( MeshAppend( M.ToStruct() , others{:} ) );
      function s = toS( x ), if isa( x ,'msh'), s = x.ToStruct(); else, s = x; end, end
    end
    function h = Plot( M , varargin )
      %precedencia: defaults de plotMESH < M.VIZ < args explicitos
      %(VIZ va DELANTE como pares nombre/valor: lo explicito lo pisa)
      fn = fieldnames( M.VIZ ).';
      vp = [ fn ; struct2cell( M.VIZ ).' ];
      h  = plotMESH( M.ToStruct() , vp{:} , varargin{:} );
    end
    function h = PlotBVH( M , varargin )
      %el metodo Capitalized ya no ensombrece a BVH\plotBVH.m: llamada directa
      h = plotBVH( M.accessCached( 'bvh' ) , M.ToStruct() , varargin{:} );
    end
  end

  %% =============================================== PUENTE STRUCT LEGADO
  methods
    function S = ToStruct( obj )
      %struct legado (.xyz/.tri/.xyzF*/.triA*/.texture); viz/info NO van
      S = struct( 'xyz' , obj.VERTICES , 'tri' , obj.FACES );
      for f = fieldnames( obj.VATTS ).', f = f{1};
        S.( ['xyz' f] ) = obj.VATTS.( f );
      end
      for f = fieldnames( obj.FATTS ).', f = f{1};
        S.( ['tri' f] ) = obj.FATTS.( f );
      end
      if isfield( obj.INFO , 'texture' ) && ~isempty( obj.INFO.texture )
        S.texture = obj.INFO.texture;
      end
    end
  end

  methods (Static)
    function obj = loadobj( obj )
      %tras load la cache (Transient) esta vacia: reinicializarla. El registro
      %cachePROPS SI se serializa (las definiciones viajan con el valor).
      if isempty( obj.CACHE ) || ~isvalid( obj.CACHE )
        obj.CACHE = cacheHandle();
      end
    end
  end

  %% ============================================ CONSTRUCCION INTERNA
  methods (Access = private)
    function obj = setFromStruct( obj , S )
      X = double( S.xyz );
      T = S.tri;
      if isempty( T ), T = zeros(0,3,'int32'); end
      obj.VERTICES = X;
      obj.FACES = msh.castFaces( T , size( X ,1) );
      nV = size( X ,1);  nF = size( T ,1);
      for f = fieldnames( S ).', f = f{1};
        if any( strcmp( f , {'xyz','tri','celltype'} ) ), continue; end
        if strcmp( f , 'texture' ), obj.INFO.texture = S.texture; continue; end
        if strcmp( f , 'DEBUG' ),   obj.DEBUG = logical( S.DEBUG ); continue; end
        v = S.(f);
        if strncmp( f , 'xyz' , 3 ) && numel( f ) > 3
          if size( v ,1) ~= nV
            warning('msh:field','node field "%s": %d rows vs nV = %d (cropped/padded).', f , size(v,1) , nV );
            v = cropPad( v , nV );
          end
          obj.VATTS.( f(4:end) ) = v;
        elseif strncmp( f , 'tri' , 3 ) && numel( f ) > 3
          if size( v ,1) ~= nF
            warning('msh:field','face field "%s": %d rows vs nF = %d (cropped/padded).', f , size(v,1) , nF );
            v = cropPad( v , nF );
          end
          obj.FATTS.( f(4:end) ) = v;
        else
          obj.INFO.( f ) = v;                   %campos desconocidos -> info
        end
      end
      function v = cropPad( v , n )
        if size( v ,1) > n, v = v( 1:n ,:,:,:,:); else, v( end+1:n ,:,:,:,:) = NaN; end
      end
    end

    function obj = reconcileNodeFields( obj )
      nV = size( obj.VERTICES ,1);
      for f = fieldnames( obj.VATTS ).', f = f{1};
        v = obj.VATTS.( f );
        if size( v ,1) == nV, continue; end
        warning('msh:field','node field "%s" resized %d -> %d rows (crop/NaN-pad; use RemoveNodes/Tidy to remap).', f , size(v,1) , nV );
        if size( v ,1) > nV, v = v( 1:nV ,:,:,:,:); else, v( end+1:nV ,:,:,:,:) = NaN; end
        obj.VATTS.( f ) = v;
      end
    end
    function obj = reconcileFaceFields( obj )
      nF = size( obj.FACES ,1);
      for f = fieldnames( obj.FATTS ).', f = f{1};
        v = obj.FATTS.( f );
        if size( v ,1) == nF, continue; end
        warning('msh:field','face field "%s" resized %d -> %d rows (crop/NaN-pad; use RemoveFaces/Tidy to remap).', f , size(v,1) , nF );
        if size( v ,1) > nF, v = v( 1:nF ,:,:,:,:); else, v( end+1:nF ,:,:,:,:) = NaN; end
        obj.FATTS.( f ) = v;
      end
    end
  end

  %% ============================================ REGISTRO DE FABRICA
  methods (Static, Access = private)
    function val = castFaces( val , nV )
      %contrato de tipos: FACES es SIEMPRE int32 (valida antes de convertir)
      if ~isnumeric( val ) || ~ismatrix( val ) || ...
         ( ~isempty( val ) && ( size( val ,2) < 1 || size( val ,2) > 4 ) )
        error('msh:faces','faces must be nF x k numeric, k in 1..4.');
      end
      if ~isempty( val )
        if isfloat( val )
          if any( ~isfinite( val(:) ) ) || any( val(:) ~= round( val(:) ) )
            error('msh:faces','faces must be integers (0 = padding).');
          end
          if max( val(:) ) > double( intmax('int32') )
            error('msh:faces','face indices exceed the int32 range.');
          end
        end
        if any( val(:) < 0 )
          error('msh:faces','faces must be nonnegative integers (0 = padding).');
        end
        if double( max( val(:) ) ) > nV
          error('msh:faces','faces reference missing nodes.');
        end
      end
      val = int32( val );
    end

    function R = defaultRegistry()
      %CPs con las que nace toda malla (sobreescribibles):
      R = struct();
      R.bvh = struct( ...
        'compute' , @(m) BVH( ToStruct( m ) ) , ...
        'events'  , struct( ...
           'transform'          , @(B,m,T) BVH( B , T ) , ...           %plegado O(1)
           'changeCoords'       , @(B,m) BVH( B , ToStruct( m ) ) , ... %refit O(n)
           'changeConnectivity' , [] ) );
      R.boundary = struct( ...
        'compute' , @(m) MeshBoundary( m.F ) , ...
        'events'  , struct( 'changeConnectivity' , [] ) );
      R.edges = struct( ...
        'compute' , @(m) meshEdges( struct( 'tri' , m.F ) ) , ...
        'events'  , struct( 'changeConnectivity' , [] ) );
      R.esup = struct( ...
        'compute' , @(m) meshEsuP( ToStruct( m ) , 'sparse' ) , ...
        'events'  , struct( 'changeConnectivity' , [] , 'changeNodeCount' , [] ) );
      R.psup = struct( ...
        'compute' , @(m) meshPsuP( ToStruct( m ) , 'sparse' ) , ...
        'events'  , struct( 'changeConnectivity' , [] , 'changeNodeCount' , [] ) );
      R.esue = struct( ...
        'compute' , @(m) meshEsuE( ToStruct( m ) , 'sparse' ) , ...
        'events'  , struct( 'changeConnectivity' , [] ) );
      R.bbox = struct( ...
        'compute' , @(m) meshBB( ToStruct( m ) ) , ...
        'events'  , struct( 'changeCoords' , [] ) );
      R.surfCent = struct( ...
        'compute' , @(m) msh.surfCentOf( m ) , ...
        'events'  , struct( 'changeCoords' , [] , 'changeConnectivity' , [] ) );
      R.triNormals = struct( ...
        'compute' , @(m) meshNormals( ToStruct( m ) ) , ...
        'events'  , struct( ...
           'changeCoords'       , [] , ...
           'changeConnectivity' , [] , ...
           'transform'          , @(N,m,T) N * msh.rotOf( T , size( N ,2) ).' ) );
    end

    function v = surfCentOf( m )
      [ s , c ] = meshSurface( ToStruct( m ) );
      v = [ s , c ];
    end

    function R = rotOf( T , d )
      %matriz (ortonormal) que TRANSPORTA NORMALES bajo una SEMEJANZA, para
      %vectores de d columnas: cross(Ra,Rb) = det(R)*R*cross(a,b), asi que se
      %incluye el factor det(R) (reflexiones exactas). ERROR si T no es
      %semejanza (la cascada del replay degrada sola).
      T = double( T );
      if d == 2
        if size( T ,1) < 2 || size( T ,2) < 2, error('msh:notSimilarity','T invalida.'); end
        A = T( 1:2 , 1:2 );
      else
        if size( T ,1) < 3 || size( T ,2) < 3, error('msh:notSimilarity','T sin parte lineal 3x3.'); end
        A = T( 1:3 , 1:3 );
        d = 3;
      end
      s = abs( det( A ) ) ^ ( 1/d );
      if s == 0, error('msh:notSimilarity','T es singular.'); end
      R = A / s;
      if norm( R.' * R - eye( d ) , 'fro' ) > 1e-9
        error('msh:notSimilarity','T no es una semejanza.');
      end
      R = R * sign( det( R ) );
    end

    function n = eventNames()
      n = { 'transform' , 'changeCoords' , 'changeNodeCount' , 'changeDim' , ...
            'changeConnectivity' , 'changeFaceCount' };
    end

    function en = eventName( a )
      a0 = a;
      a  = lower( char( a ) );
      switch a
        case 'transform',                                          en = 'transform';
        case {'changecoords','coords','coordinates','coordenadas'},en = 'changeCoords';
        case {'changenodecount','nodecount'},                      en = 'changeNodeCount';
        case {'changedim','dim'},                                  en = 'changeDim';
        case {'changeconnectivity','connectivity','conectividad'}, en = 'changeConnectivity';
        case {'changefacecount','facecount'},                      en = 'changeFaceCount';
        otherwise
          error('msh:cached','evento desconocido ''%s'' (validos: %s).', ...
                char( a0 ) , strjoin( msh.eventNames() , ', ' ) );
      end
    end

    function s = lst( c )
      if isempty( c ), s = '-'; else, s = strjoin( c , ', ' ); end
    end

    function s = hDesc( h )
      if isempty( h ), s = 'invalida'; else, s = 'handler'; end
    end

    function [ p , ax , c ] = planarInfo( X )
      %p = 0 no planar | 1 planar generica | 2 planar con plano axis-aligned
      %(ax en 1..3, coordenada c). Tolerancia RELATIVA sobre la extension.
      p = 0;  ax = 0;  c = NaN;
      if size( X ,2) ~= 3 || size( X ,1) < 3 || any( ~isfinite( X(:) ) ), return; end
      Xc = X - mean( X ,1);
      C  = Xc.' * Xc;
      [ Ve , D ] = eig( ( C + C.' )/2 );
      [ d , o ] = sort( diag( D ) );                 %ascendente
      if d(3) <= 0 || d(1) > 1e-14 * d(3), return; end
      p = 1;
      n = Ve( : , o(1) );                            %normal del plano
      [ mx , axk ] = max( abs( n ) );
      if mx >= 1 - 1e-9
        p = 2;  ax = axk;  c = mean( X(:,axk) );
      end
    end

    %-------------------------------------------- formateadores para display
    function s = thousands( n )
      %entero con separadores de miles: 123456712 -> '123,456,712'
      s = regexprep( sprintf( '%d' , round( double( n ) ) ) , ...
                     '(\d)(?=(\d{3})+$)' , '$1,' );
    end

    function s = dimStr( X )
      %'2D' | '3D' | '3D, flat (z = 0)' | '3D, planar[: eje = c]'
      %flat   = 3 columnas y z == 0 EXACTO en todos los vertices
      %planar = 3 columnas pero coplanar (cualquier plano)
      %pura descripcion: la DIMENSION nsd sigue siendo el numero de columnas
      if size( X ,2) == 2, s = '2D'; return; end
      s = '3D';
      if ~isempty( X ) && all( X(:,3) == 0 )
        s = '3D, flat (z = 0)';  return;
      end
      [ p , ax , c ] = msh.planarInfo( X );
      if p == 2
        axes_ = 'xyz';
        s = sprintf( '3D, planar: %s = %.6g' , axes_(ax) , c );
      elseif p == 1
        s = '3D, planar';
      end
    end

    function s = cellStr( T )
      if isempty( T ), s = '(none)'; return; end
      nz  = sum( T ~= 0 , 2 );
      nms = { 'points' , 'segments' , 'triangles' , 'tets' };
      cnt = [ sum(nz==1) , sum(nz==2) , sum(nz==3) , sum(nz>=4) ];
      u   = find( cnt > 0 );
      if isempty( u ),  s = '(empty cells)';           return; end
      if isscalar( u ), s = sprintf( '(%s)' , nms{u} ); return; end
      p = arrayfun( @(k) sprintf( '%s %s' , msh.thousands( cnt(k) ) , nms{k} ) , u ,'uni',0);
      s = sprintf( '(mixed: %s)' , strjoin( p , ' + ' ) );
    end

    function s = fmtSize( v )
      sz = strjoin( arrayfun( @(n) sprintf('%d',n) , size(v) ,'uni',0) , 'x' );
      if issparse( v ), c = [ 'sparse ' class(v) ]; else, c = class( v ); end
      if iscell( v ), s = sprintf( '{%s cell}' , sz );
      else,           s = sprintf( '[%s %s]' , sz , c );
      end
    end

    function s = fmtVal( v )
      %escalares y vectores "cortos" (<=3 numeros) se imprimen; el resto se
      %resume como [tamano clase]
      if isstruct( v ) && isscalar( v ) && isfield( v ,'nE') && isfield( v ,'volume')
        s = sprintf( '[blob BVH: %s elems, %s]' , msh.thousands( v.nE ) , v.volume );
      elseif isstring( v ) && isscalar( v )
        s = msh.fmtVal( char( v ) );
      elseif isnumeric( v ) || islogical( v )
        if issparse( v )
          s = sprintf( '%s, nnz %d' , msh.fmtSize( v ) , nnz( v ) );
        elseif isscalar( v )
          if islogical( v )
            if v, s = 'true'; else, s = 'false'; end
          else
            s = sprintf( '%.6g' , double( v ) );
          end
        elseif ~isempty( v ) && numel( v ) <= 3 && isnumeric( v )
          s = [ '[' strtrim( sprintf( '%.6g ' , double( v ) ) ) ']' ];
        else
          s = msh.fmtSize( v );
        end
      elseif ischar( v ) && size( v ,1) <= 1
        if numel( v ) > 32, v = [ v(1:29) '...' ]; end
        s = [ '''' v '''' ];
      elseif isstruct( v ) && isscalar( v )
        fn = fieldnames( v );
        if numel( fn ) <= 4
          s = sprintf( '[struct: %s]' , strjoin( fn.' , ', ' ) );
        else
          s = sprintf( '[struct with %d fields]' , numel( fn ) );
        end
      else
        s = msh.fmtSize( v );
      end
    end

    function s = kvStr( S )
      fn = fieldnames( S ).';
      p  = cellfun( @(f) sprintf( '%s = %s' , f , msh.fmtVal( S.(f) ) ) , fn ,'uni',0);
      s  = strjoin( p , ', ' );
    end
  end

  %% ============================================================== DISPLAY
  % informativo pero sin evaluar NINGUN getter Dependent: solo lee lo que ya
  % existe (datos + entradas ya presentes en la cache). disp y display
  % comparten esta implementacion (matlab.mixin.CustomDisplay).
  methods (Access = protected)
    function displayScalarObject( obj )
      X = obj.VERTICES;  T = obj.FACES;
      if isempty( X ) && isempty( T )
        fprintf( '  msh: empty\n' );
      else
        fprintf( '  msh: %s vertices (%s), %s faces %s\n' , ...
                 msh.thousands( size(X,1) ) , msh.dimStr( X ) , ...
                 msh.thousands( size(T,1) ) , msh.cellStr( T ) );
      end
      fn = fieldnames( obj.VATTS ).';
      if ~isempty( fn )
        p = cellfun( @(f) sprintf( 'xyz%s %s' , f , msh.fmtSize( obj.VATTS.(f) ) ) , fn ,'uni',0);
        fprintf( '    node fields: %s\n' , strjoin( p , ', ' ) );
      end
      fn = fieldnames( obj.FATTS ).';
      if ~isempty( fn )
        p = cellfun( @(f) sprintf( 'tri%s %s' , f , msh.fmtSize( obj.FATTS.(f) ) ) , fn ,'uni',0);
        fprintf( '    face fields: %s\n' , strjoin( p , ', ' ) );
      end
      if ~isempty( fieldnames( obj.VIZ ) )
        fprintf( '    VIZ:  %s\n' , msh.kvStr( obj.VIZ ) );
      end
      if ~isempty( fieldnames( obj.INFO ) )
        fprintf( '    INFO: %s\n' , msh.kvStr( obj.INFO ) );
      end
      names = fieldnames( obj.cachePROPS ).';
      c = obj.CACHE;
      live = {};  none = {};
      for n = names, n = n{1};
        if ~isempty( c ) && isvalid( c ) && c.has( n ), live{end+1} = n;   %#ok<AGROW>
        else,                                           none{end+1} = n;   %#ok<AGROW>
        end
      end
      if ~isempty( live )
        live = sort( live );
        w = max( cellfun( @numel , live ) );
        fprintf( '    CPs:\n' );
        for n = live, n = n{1};
          if strcmp( c.state( n ) , 'fresh' )
            d = msh.fmtVal( c.value( n ) );
          else
            L = c.log( n );
            evs = cellfun( @(e) e.fired{1} , L ,'uni',0);
            d = sprintf( '(pendiente: %s)' , strjoin( evs , ', ' ) );
          end
          fprintf( '      %-*s  %s\n' , w , n , d );
        end
      end
      if ~isempty( none )
        fprintf( '    definidas sin calcular: %s\n' , strjoin( sort( none ) , ', ' ) );
      end
      if obj.DEBUG, fprintf( '    DEBUG: on\n' ); end
      fprintf( '\n' );
    end
  end
end
