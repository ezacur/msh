classdef cacheView
%CACHEVIEW  Vista proxy de las CPs de un objeto: lo que devuelve M.CP.
%
%   No tiene estado propio: envuelve al dueno (p.ej. un msh) y reenvia toda la
%   cadena de subscripts a su despachador cachedAccess, de modo que
%       M.CP.bvh                 == valor (computa/replay si hace falta)
%       M.CP.bvh.frame           == indexar dentro del valor
%       M.CP.bvh.delete          == borrar el valor (no la definicion)
%       M = M.CP.bvh.removeProp  == borrar definicion + valor
%       M = M.CP.bvh.set( x )    == sembrar un valor a mano (aislado, COW)
%       M.CP.bvh.changeCoords    == handler del evento (o [] si invalida)
%   Mostrar `M.CP` sin mas imprime la tabla de definiciones y estados.
%   (La lectura corriente no necesita el proxy: M.bvh lee perezoso y
%   M.bvh_ recalcula a la fuerza.)
%
% See also msh, cacheHandle.

  properties (Access = private)
    OWNER = []   % el objeto dueno (valor; debe implementar cachedAccess y
                 % displayCachedView)
  end

  methods
    function obj = cacheView( M )
      if nargin, obj.OWNER = M; end
    end

    function varargout = subsref( obj , s )
      out = cachedAccess( obj.OWNER , s );
      if isempty( out )
        varargout = {};
      else
        varargout = out( 1:max( min( nargout , numel(out) ) , 1 ) );
      end
    end

    function disp( obj )
      displayCachedView( obj.OWNER );
    end
  end
end
