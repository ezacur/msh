classdef serie < cachedOwner
%serie  Una lista de numeros. El ejemplo de DISEÑO del patron de cache.
%
%   El dominio es deliberadamente minusculo -- una lista, y dos formas de
%   editarla -- para que todo lo que se discuta sea POLITICA DE CACHE y no
%   dominio:
%       o = o.Push( v )    anade un valor al final
%       o = o.Pop()        quita el ultimo
%
%   Tres derivados, elegidos porque cada uno pide una politica distinta:
%       avg     la media          -- barata de actualizar al anadir, no al quitar
%       top     el maximo         -- barata al anadir, IMPOSIBLE al quitar
%       sorted  la lista ordenada -- no hay atajo por ningun lado
%
%   Y DOS REGISTROS para los mismos tres derivados, que es la gracia del
%   fichero: la politica es separable del dominio.
%       o = serie().Safe()   CONSERVADOR: un solo evento, todo se invalida
%       o = serie().Fast()   con los atajos que cada derivado admite
%   Los dos son CORRECTOS. El primero es el que hay que escribir primero.
%
%   Ver serie_TUTORIAL.html (el porque de cada decision, paso a paso).
%
% See also cachedOwner, cart, route, facto, ivp, serieDemo.

  properties (SetAccess = private)
    X = []          % los numeros, 1xN
  end

  properties (Constant, Hidden)
    %ESPECIFICO -> GENERAL. 'changedata' existe para poder decir "me da igual
    %que haya pasado, tirame" sin tener que enumerar las ediciones: una CP que
    %lo declare con [] queda inmunizada contra editores FUTUROS.
    CP_EVENTS = { 'push' , 'pop' , 'changedata' }
  end

  methods
    %% ---------------------------------------------------------- EDITORES
    function o = Push( o , v )
      %anade un valor. Los args del evento son TODO lo que un handler
      %incremental va a poder saber de esta edicion -- no puede mirar el objeto
      %(veria el estado FINAL, no el de su edicion) --, asi que se le pasa el
      %valor anadido Y el numero de elementos QUE QUEDAN. Ese `n` es el que
      %hace que una cadena de varios Push se replique bien: cada edicion lleva
      %el suyo.
      o.X( end+1 ) = v;
      o = o.fire( { 'push' , 'changedata' } , v , numel( o.X ) );
    end

    function o = Pop( o )
      %quita el ultimo. Mismos args: el valor que SALE y el n que queda.
      if isempty( o.X ), error('serie:empty','la serie esta vacia.'); end
      x = o.X( end );  o.X( end ) = [];
      o = o.fire( { 'pop' , 'changedata' } , x , numel( o.X ) );
    end

    %% -------------------------------------------------------- REGISTRO 1
    function o = Safe( o )
      %CONSERVADOR: cada CP declara SOLO el evento general, con [] -- "pase lo
      %que pase, tirame". Cualquier edicion invalida todo; la unica ventaja que
      %queda (que no es poca) es la pereza: se calcula al pedirlo y se reusa
      %hasta la siguiente edicion. No hay handlers, asi que no hay nada que
      %pueda mentir.
      o = o.Define( 'avg'    , @(s) sum( s.X ) / max( numel( s.X ) , 1 ) , 'changedata' , [] );
      o = o.Define( 'top'    , @(s) max( [ -inf , s.X ] )                , 'changedata' , [] );
      o = o.Define( 'sorted' , @(s) sort( s.X )                          , 'changedata' , [] );
    end

    %% -------------------------------------------------------- REGISTRO 2
    function o = Fast( o )
      %CON ATAJOS, uno por derivado y solo donde existe de verdad.
      %OJO A LOS ~: el segundo argumento de un handler INCREMENTAL es el dueno
      %y NO se puede leer (en el replay seria el estado final, no el de su
      %edicion), asi que se declara no usado. Todo lo que el atajo necesita
      %llega por los args del evento.
      o = o.Define( 'avg' , @(s) sum( s.X ) / max( numel( s.X ) , 1 ) , ...
                    'push' , @(m,~,v,n) m + ( v - m ) / n , ...              % INCREMENTAL
                    'pop'  , @(~,s) sum( s.X ) / max( numel( s.X ) , 1 ) );  % SYNC absoluto
      o = o.Define( 'top' , @(s) max( [ -inf , s.X ] ) , ...
                    'push' , @(t,~,v,~) max( t , v ) , ...                   % INCREMENTAL
                    'pop'  , [] );                                           % INVALIDA
      o = o.Define( 'sorted' , @(s) sort( s.X ) , 'changedata' , [] );       % sin atajo
    end
  end
end
