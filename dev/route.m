classdef route < cachedOwner
%route  Segundo ejemplo del patron de cache: una polilinea de waypoints.
%
%   Existe para demostrar UNA afirmacion: que heredar de cachedOwner basta y
%   que la clase nueva es puro dominio. Se escribio de cero y funciono a la
%   primera. Ademas ensena un caso que cart no ejercita: una CP PENDIENTE que
%   SOBREVIVE a un evento que no le afecta, conservando su log.
%
%   Derivados (los registra Demo):
%     len    longitud total    -- incremental al anadir punto, sync al escalar
%     bbox   caja envolvente   -- solo sync absoluto
%     npts   numero de puntos  -- incremental al anadir; INSENSIBLE a escalar
%
%   Uso:  R = route().Demo();     % o mira routeDemo.m
%
% See also cachedOwner, cart, serie, facto, ivp, routeDemo.

  properties (SetAccess = private)
    XY = zeros(0,2)      % waypoints, una fila por punto
  end

  properties (Constant, Hidden)
    CP_EVENTS = { 'addpoint' , 'changegeom' }   % los eventos que dispara esta clase
  end

  methods
    function obj = AddPoint( obj , x , y )
      %anade un waypoint. Evento ESPECIFICO->GENERAL, con args: el tramo nuevo
      d = 0;
      if ~isempty( obj.XY ), d = norm( [ x y ] - obj.XY(end,:) ); end
      obj.XY( end+1 , : ) = [ x y ];
      obj = obj.fire( { 'addpoint' , 'changegeom' } , d );
    end

    function obj = Scale( obj , k )
      %escala la ruta. Evento absoluto: SIN args
      obj.XY = obj.XY * k;
      obj = obj.fire( { 'changegeom' } );
    end

    function obj = Demo( obj )
      %registra las 3 CPs de ejemplo
      %ojo al ~ del incremental: en el replay ese argumento seria la ruta de
      %AHORA, no la de su edicion, asi que no se puede leer (Define lo rechaza)
      obj = obj.Define( 'len' , @(r) sum( vecnorm( diff( r.XY ) , 2 , 2 ) ) , ...
                        'addpoint'   , @(v,~,d) v + d , ...                           % INCREMENTAL
                        'changegeom' , @(v,r) sum( vecnorm( diff( r.XY ) ,2,2) ) );   % absoluto
      obj = obj.Define( 'bbox' , @(r) [ min( r.XY ,[],1) , max( r.XY ,[],1) ] , ...
                        'changegeom' , @(v,r) [ min( r.XY ,[],1) , max( r.XY ,[],1) ] );
      obj = obj.Define( 'npts' , @(r) size( r.XY , 1 ) , ...
                        'addpoint' , @(v,~,d) v + 1 );
      %npts NO declara changegeom: escalar no cambia cuantos puntos hay -> sobrevive
    end
  end
end
