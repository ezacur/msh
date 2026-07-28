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
%   EL ALMACEN ES UN STRUCT, no un containers.Map: las claves ya son nombres
%   MATLAB validos (DefineCP lo exige) y el acceso por campo dinamico sale ~5x
%   mas barato que un lookup de Map -- importa porque el HIT esta en el camino
%   caliente de toda lectura de CP. De paso hace `clone` UNA asignacion en vez
%   de un bucle: el struct es un valor y MATLAB ya copia perezosamente.
%   Sirve la entrada COMPLETA (tryGet/entry) para que quien lee no pague un
%   lookup por campo.
%
%   `inFlight` anota las claves en curso de resolucion: si el compute (o un
%   handler) de una CP vuelve a pedir esa misma CP, `enter` lanza un error
%   claro en vez de recursar hasta agotar la pila.
%
% See also msh, cacheView.

  properties (Access = private)
    store    = struct()   % clave -> struct('state','fresh'|'pending','value',v,'log',{editos})
    inFlight = struct()   % clave -> true mientras se resuelve (guarda de ciclos)
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

    function e  = entry( obj , key ), e  = obj.store.( key );        end
    function st = state( obj , key ), st = obj.store.( key ).state;  end
    function v  = value( obj , key ), v  = obj.store.( key ).value;  end
    function L  = log(   obj , key ), L  = obj.store.( key ).log;    end

    function setFresh( obj , key , v )
      %el {v} es obligatorio: struct('value',C) con C cell crearia un struct
      %ARRAY en vez de un escalar con el cell dentro
      obj.store.( key ) = struct( 'state','fresh' , 'value',{v} , 'log',{{}} );
    end

    function setPending( obj , key , v , L )
      obj.store.( key ) = struct( 'state','pending' , 'value',{v} , 'log',{L} );
    end

    function setEntry( obj , key , e )
      obj.store.( key ) = e;
    end

    function remove( obj , key )
      if isfield( obj.store , key ), obj.store = rmfield( obj.store , key ); end
    end

    function ks = keys( obj )
      ks = fieldnames( obj.store ).';   %FILA (como devolvia containers.Map)
    end

    function c2 = clone( obj )
      %copia superficial: el store es un VALOR, basta asignarlo
      c2 = cacheHandle();
      c2.store = obj.store;
    end

    function c2 = cloneWithout( obj , key )
      c2 = obj.clone();
      c2.remove( key );
    end

    function clear( obj )
      obj.store = struct();
    end
  end

  methods   % --------------------------------------- guarda de ciclos
    function enter( obj , key )
      %marca `key` en curso de resolucion. Reentrar = el compute (o un handler)
      %de esa CP la vuelve a pedir: se corta con un error legible en vez de
      %recursar hasta agotar la pila.
      if isfield( obj.inFlight , key )
        error('cacheHandle:cycle', ...
              'ciclo en la cachedProp ''%s'': su compute (o un handler de evento) la vuelve a pedir.', key );
      end
      obj.inFlight.( key ) = true;
    end

    function leave( obj , key )
      if isfield( obj.inFlight , key ), obj.inFlight = rmfield( obj.inFlight , key ); end
    end
  end
end
