classdef cacheView
%CACHEVIEW  Vista proxy de las cachedProps de un objeto: lo que devuelve M.cached.
%
%   No tiene estado propio: envuelve al dueno (p.ej. un msh) y reenvia toda la
%   cadena de subscripts a su despachador cachedAccess, de modo que
%       M.cached.BVH                 == valor (computa/replay si hace falta)
%       M.cached.BVH.frame           == indexar dentro del valor
%       M.cached.BVH.delete          == borrar el valor (no la definicion)
%       M = M.cached.BVH.removeProp  == borrar definicion + valor
%       M = M.cached.BVH.set( x )    == sembrar un valor a mano (aislado, COW)
%       M.cached.BVH.changeCoords    == handler del evento (o [] si invalida)
%   Mostrar `M.cached` sin mas imprime la tabla de definiciones y estados.
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
