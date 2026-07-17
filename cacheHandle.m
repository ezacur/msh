classdef cacheHandle < handle
%cacheHandle  Almacen (por referencia) de los valores de las cachedProps de `msh`.
%
%   Es el "puntero" del diseno: un handle privado compartido entre copias de
%   una malla. Cada clave guarda una entrada con ESTADO:
%     'fresh'    valor valido, se sirve tal cual
%     'pending'  valor viejo + log de eventos sin aplicar (replay perezoso:
%                msh lo resuelve en el proximo acceso con los handlers de
%                eventos de la definicion)
%   Las DEFINICIONES (que computa cada clave, y como reacciona a los eventos)
%   NO viven aqui sino en el registro `cachePROPS` del dueno (parte del valor).
%
%   Al EDITAR la malla, msh NO muta este objeto: construye un handle NUEVO con
%   las entradas que sobreviven o quedan pendientes (copy-on-write), de modo
%   que las copias hermanas conservan su cache intacta. Si mutan este handle
%   (delete de un valor, resolucion de un pendiente, deposito de un MISS) el
%   efecto es compartido y benigno: solo estado de rendimiento.
%
% See also msh, cacheView.

  properties (Access = private)
    store   % containers.Map: key -> struct('state','fresh'|'pending','value',v,'log',{editos})
  end

  methods
    function obj = cacheHandle()
      obj.store = containers.Map('KeyType','char','ValueType','any');
    end

    function tf = has( obj , key )
      tf = isKey( obj.store , key );
    end

    function st = state( obj , key )
      e = obj.store( key );  st = e.state;
    end

    function v = value( obj , key )
      e = obj.store( key );  v = e.value;
    end

    function L = log( obj , key )
      e = obj.store( key );  L = e.log;
    end

    function setFresh( obj , key , v )
      obj.store( key ) = struct( 'state','fresh' , 'value',{v} , 'log',{{}} );
    end

    function setPending( obj , key , v , L )
      obj.store( key ) = struct( 'state','pending' , 'value',{v} , 'log',{L} );
    end

    function e = entry( obj , key )
      e = obj.store( key );
    end

    function setEntry( obj , key , e )
      obj.store( key ) = e;
    end

    function remove( obj , key )
      if isKey( obj.store , key ), remove( obj.store , key ); end
    end

    function ks = keys( obj )
      ks = keys( obj.store );
    end

    function c2 = clone( obj )
      %copia superficial (las entradas son valores)
      c2 = cacheHandle();
      for k = keys( obj.store )
        c2.store( k{1} ) = obj.store( k{1} );
      end
    end

    function c2 = cloneWithout( obj , key )
      c2 = obj.clone();
      c2.remove( key );
    end

    function clear( obj )
      remove( obj.store , keys( obj.store ) );
    end
  end
end
