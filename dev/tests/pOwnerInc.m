function v = pOwnerInc( v , r , d )
%pOwnerInc  sonda: handler incremental OPACO que LEE el dueno.
%
%   Es la media de la longitud de los tramos de una `route`, escrita de la forma
%   natural -- y por eso mismo MAL: en el replay de varias ediciones pendientes
%   `r` es la ruta de AHORA, no la de cada edicion, asi que el numero de tramos
%   que lee es el final para todas.
%
%   Existe para probar LA RED de cachedOwner: al no ser anonimo, func2str no
%   trae cuerpo y Define no puede rechazarlo (lo mismo que pasaria con un
%   handler venido de un .mat viejo). Con UNA edicion pendiente el dueno si es
%   el correcto y el atajo se usa; con MAS, el replay renuncia al atajo.
%
% See also test_cache, cachedOwner.
  ns = size( r.XY , 1 ) - 1;            % tramos tras esta edicion
  v  = ( v*( ns - 1 ) + d ) / ns;
end
