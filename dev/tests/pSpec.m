classdef pSpec < cachedOwner
%pSpec  sonda: quien decide cuando el ESPECIFICO y el GENERAL traen politicas
%       distintas. Las dos CPs son simetricas a proposito:
%
%     k   especifico 'add' -> ABSOLUTO   | general 'chg' -> INCREMENTAL
%     m   especifico 'add' -> INCREMENTAL| general 'chg' -> ABSOLUTO
%
%   Los handlers devuelven valores imposibles (999, 777) para que se vea CUAL
%   corrio. Manda siempre el mas especifico DECLARADO Y DISPARADO: con Add(x)
%   -- que dispara los dos -- k debe dar 999 (su absoluto) y m debe dar 100+x
%   (su incremental); con Bump -- que dispara solo el general -- k debe usar su
%   incremental (101) y m su absoluto (777).
%
% See also test_cache, cachedOwner.
  properties (SetAccess = private)
    N = 0
  end
  properties (Constant, Hidden)
    CP_EVENTS = { 'add' , 'chg' }
  end
  methods
    function o = pSpec()
      o.DEBUG = false;
      o = o.Define( 'k' , @(p) 100 , 'add' , @(v,p)   999   , 'chg' , @(v,~,x) v + x );
      o = o.Define( 'm' , @(p) 100 , 'add' , @(v,~,x) v + x , 'chg' , @(v,p)   777   );
    end
    function o = Add( o , x )          % ESPECIFICO -> GENERAL
      o.N = o.N + x;
      o = o.fire( { 'add' , 'chg' } , x );
    end
    function o = Bump( o )             % solo el GENERAL
      o.N = o.N + 1;
      o = o.fire( { 'chg' } , 1 );
    end
  end
end
