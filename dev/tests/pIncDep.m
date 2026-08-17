classdef pIncDep < cachedOwner
%pIncDep  sonda: CP COMPUESTA con handler INCREMENTAL -- las deps deben
%sobrevivir al replay.
%
%   `comp` sale de `base` (su compute la lee via access) Y ADEMAS declara un
%   incremental para 'bump', que solo toca B (a base ni la roza). El refresco
%   del replay guardaba takeDeps() a secas -- y un handler incremental no lee
%   nada --, asi que la dependencia de base se perdia en el PRIMER replay: a
%   la siguiente caida de base, comp sobrevivia sirviendo rancio en silencio
%   (comprobado: servia 305 con la verdad en 905). Grupo 19 de test_cache.
%
% See also test_cache, pComp, cachedOwner.
  properties (SetAccess = private)
    A = 3      % lo lee base
    B = 0      % lo lee comp (su termino propio)
  end
  properties (Constant, Hidden)
    CP_EVENTS = { 'bump' , 'chgA' }
  end
  methods
    function obj = pIncDep()
      obj.DEBUG = false;
      obj = obj.Define( 'base' , @(o) o.A * 100 , 'chgA' , [] );
      obj = obj.Define( 'comp' , @(o) o.access('base') + o.B , ...
                        'bump' , @(v,~,d) v + d );
    end
    function obj = Bump( obj , d )     % toca B: base NO se entera (correcto)
      obj.B = obj.B + d;
      obj = obj.fire( 'bump' , d );
    end
    function obj = ChgA( obj , a )     % toca A: base cae, y comp debe caer CON ella
      obj.A = a;
      obj = obj.fire( 'chgA' );
    end
  end
end
