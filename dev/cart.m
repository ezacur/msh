classdef cart < cachedOwner
%cart  Carrito de la compra: el ejemplo MINIMO del patron de cache.
%
%   PARA QUE SIRVE ESTE FICHERO: es la clase de juguete del tutorial
%   (cache_TUTORIAL.html, en la raiz). Si vas a escribir tu propia clase
%   cacheada, copia esta forma: todo el mecanismo vive en cachedOwner y aqui
%   solo hay DOMINIO -- los datos, los editores (que avisan con fire), las
%   definiciones y el vocabulario de eventos.
%
%   LA IDEA EN DOS LINEAS: los derivados caros (el total) se apuntan en un
%   "cuaderno" al calcularse por primera vez. Editar el carrito NO recalcula
%   nada: solo decide que hacer con cada apunte -- dejarlo (no le afecta),
%   tacharlo (se recalculara al leer), o anotarle el cambio al margen para
%   ponerse al dia barato cuando alguien lo lea.
%
%   COMO SE USA:
%     C = cart().Demo();        % registra las 4 CPs de ejemplo
%     C = C.Add('cafe',3,2);    % editar (clase de VALOR: hay que reasignar)
%     C.total                   % leer: calcula SOLO si hace falta
%     C.total_                  % recalcular a la fuerza (para verificar)
%     C.CP                      % la tabla: que CPs hay y en que estado
%
%   CPs de ejemplo (las registra Demo); cada una ensena UNA politica:
%     total     suma de precio*cantidad  -- se pone al dia SOLA al anadir (+d)
%     ticket    los nombres concatenados -- se invalida al anadir ([])
%     topPrice  el precio mayor          -- solo reacciona al evento GENERAL
%     currency  constante                -- insensible: nada la toca
%
% See also cachedOwner, cacheHandle, cacheProxy, route, serie, facto, ivp, cartDemo.

  properties (SetAccess = private)
    NAME  = {}      % 1xN nombres
    PRICE = []      % 1xN precios unitarios
    QTY   = []      % 1xN cantidades
  end

  properties (Constant, Hidden)
    CP_EVENTS = { 'additem' , 'changeprice' , 'changedata' }   % los eventos que dispara esta clase
  end

  methods
    function obj = cart()
      obj.DEBUG = true;                       % el ejemplo narra lo que hace
    end

    function obj = Add( obj , name , price , qty )
      %anade una linea. LA REGLA DE ORO: todo metodo que edita datos termina en
      %fire, eventos de ESPECIFICO a GENERAL. Aqui ademas se pasa un ARGUMENTO
      %(el subtotal de la linea nueva): es lo que permite que `total` se ponga
      %al dia sumando en vez de recorriendo el carrito
      obj.NAME{  end+1 } = name;
      obj.PRICE( end+1 ) = price;
      obj.QTY(   end+1 ) = qty;
      obj = obj.fire( { 'additem' , 'changedata' } , price*qty );
    end

    function obj = SetPrice( obj , name , price )
      %cambia un precio. SIN args: no hay delta util que pasar, asi que quien
      %quiera ponerse al dia tendra que mirar el estado final (sync absoluto)
      k = strcmp( obj.NAME , name );
      if ~any( k ), error('cart:noitem','no hay "%s" en el carrito.', name ); end
      obj.PRICE( k ) = price;
      obj = obj.fire( { 'changeprice' , 'changedata' } );
    end

    function obj = Demo( obj )
      %registra las 4 CPs de ejemplo. Recuerda: LA FIRMA DEL HANDLER ES LA
      %POLITICA -- @(v,~,args) incremental ("aplicame este cambio"), @(v,c)
      %sync absoluto ("ponme al dia mirando el objeto"), [] invalidar
      %("tirame"), no declarar = insensible ("ni me toques")
      %El incremental lleva ~ y no un nombre A PROPOSITO: en el replay ese
      %argumento es el carrito de AHORA, no el de su edicion, asi que leerlo da
      %dato malo en silencio (Define lo rechaza)
      obj = obj.Define( 'total' , @(c) sum( c.PRICE .* c.QTY ) , ...
                        'additem'     , @(v,~,d) v + d , ...                 % INCREMENTAL
                        'changeprice' , @(v,c)   sum( c.PRICE .* c.QTY ) );  % absoluto
      obj = obj.Define( 'ticket' , @(c) strjoin( c.NAME , ' + ' ) , ...
                        'additem' , [] );                                    % INVALIDAR
      obj = obj.Define( 'topPrice' , @(c) max( [ 0 , c.PRICE ] ) , ...
                        'changedata' , @(v,c) max( [ 0 , c.PRICE ] ) );       % solo el GENERAL
      obj = obj.Define( 'currency' , @(c) 'EUR' );                           % sin eventos
    end
  end
end
