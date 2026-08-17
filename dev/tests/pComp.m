classdef pComp < cachedOwner
%pComp  sonda: CPs COMPUESTAS -- una CP que se calcula leyendo otra.
%
%   Las tres formas, para el grupo 19 de test_cache:
%     base        CP normal, cae con el evento 'chg'
%     dotInside   compute definido DENTRO de la clase que lee la CP con un
%                 PUNTO: no puede funcionar (MATLAB salta el subsref en contexto
%                 de clase) y el motor lo traduce a cachedOwner:cpFromInside
%     viaAccess   la forma correcta desde dentro: o.access('base')
%     stubborn    ademas de salir de base, DECLARA un handler que se queda con
%                 el valor viejo: sirve para comprobar que una dependencia rota
%                 gana a la declaracion
%
% See also test_cache, cachedOwner.
  properties (SetAccess = private)
    N = 1
  end
  properties (Constant, Hidden)
    CP_EVENTS = { 'chg' }
  end
  methods
    function obj = pComp()
      obj.DEBUG = false;
      obj = obj.Define( 'base'      , @(o) 10*o.N , 'chg' , [] );
      obj = obj.Define( 'dotInside' , @(o) o.base + 1 );
      obj = obj.Define( 'viaAccess' , @(o) o.access('base') + 2 );
      obj = obj.Define( 'stubborn'  , @(o) o.access('base') + 3 , ...
                        'chg' , @(v,o) v );        % handler que MIENTE a proposito
    end
    function obj = Bump( obj )
      obj.N = obj.N + 1;
      obj = obj.fire( { 'chg' } );
    end
  end
end
