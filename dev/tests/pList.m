classdef pList < cachedOwner
%pList  sonda: CPs cuyo VALOR da una lista separada por comas al indexar.
%
%   items  struct array 1xN -> items.n son N salidas
%   cellv  cell 1xN         -> cellv{:} son N salidas
%   plain  un escalar       -> una salida, como siempre
%   El grupo 22 de test_cache comprueba que la cadena NO se trunca: pedir una
%   sola salida devolvia el primer elemento y se perdia el resto EN SILENCIO,
%   con el agravante de que la misma expresion, pasando el valor por una
%   variable, si daba las N.
%
% See also test_cache, cachedOwner.
  properties (SetAccess = private)
    N = 3
  end
  properties (Constant, Hidden)
    CP_EVENTS = { 'chg' }
  end
  methods
    function o = pList()
      o.DEBUG = false;
      o = o.Define( 'items' , @(p) struct( 'n' , num2cell( 1:p.N ) ) , 'chg' , [] );
      o = o.Define( 'cellv' , @(p) num2cell( 10*(1:p.N) )           , 'chg' , [] );
      o = o.Define( 'plain' , @(p) 40 + p.N                         , 'chg' , [] );
    end
    function o = Grow( o )
      o.N = o.N + 1;
      o = o.fire( { 'chg' } );
    end
  end
end
