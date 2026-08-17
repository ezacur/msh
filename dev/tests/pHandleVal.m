classdef pHandleVal < cachedOwner
%pHandleVal  sonda: una CP cuyo VALOR es un objeto handle.
%
%   Es el UNICO sitio donde el patron deja de ser solo rendimiento y puede
%   corromper DATO, y hasta ahora solo se probaba que el motor AVISA
%   (cachedOwner:handleValue) -- no que el dano sea real. Esta sonda lo hace
%   demostrable:
%     hmap    devuelve un containers.Map, que es HANDLE: las copias comparten
%             la MISMA referencia, asi que mutarlo por una se ve por la otra
%     svalue  el mismo dato como struct, que es VALOR: ahi el aislamiento vale
%   Las dos declaran su evento a proposito, para que la guarda del olvido total
%   (noEvents) no salte y ensucie la salida del test.
%
%   El motor NO PUEDE impedirlo: el valor lo produce el dominio. Por eso avisa.
%
% See also test_cache, cachedOwner.
  properties (SetAccess = private)
    N = 1
  end
  properties (Constant, Hidden)
    CP_EVENTS = { 'chg' }
  end
  methods
    function obj = pHandleVal()
      obj.DEBUG = false;
      obj = obj.Define( 'hmap'   , @(o) containers.Map( 'n' , o.N ) , 'chg' , [] );
      obj = obj.Define( 'svalue' , @(o) struct( 'n' , o.N )         , 'chg' , [] );
    end
    function obj = Bump( obj )
      obj.N = obj.N + 1;
      obj = obj.fire( 'chg' );
    end
  end
end
