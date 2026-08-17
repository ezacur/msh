classdef ivp < cachedOwner
%ivp  Problema de valor inicial (dy/dt = F(t,y), y(0) = Y0) con solucion
%     cacheada INCREMENTALMENTE: el ejemplo que SI paga el precio del patron.
%
%   La CP `sol` guarda la tabla ordenada de pares [t,y] ya resueltos (fila 1 =
%   [0,Y0]). Pedir y(t) con At:
%     - si t ya esta -> HIT, cero integracion;
%     - si no, se integra (RK4 a paso fijo DT) SOLO desde el tiempo guardado
%       mas cercano por debajo, y el par nuevo se inserta EN SU SITIO. Con
%       soluciones en t = [0 2 4 6], pedir y(5) integra de 4 a 5 y nada mas.
%   El coste de integrar es proporcional al tramo, asi que aqui el atajo PUEDE
%   amortizar el peaje del patron -- al contrario que facto, que ensena la
%   misma idea con costes de juguete y sale perdiendo 163x.
%
%   CUANTO PAGA, MEDIDO (bench_incremental, 20 consultas crecientes): el peaje
%   del patron es FIJO (~1-3 ms por lectura entre despacho, fire y replay) y el
%   ahorro es proporcional al trabajo evitado, asi que todo depende del paso:
%       dt=1e-3 ( 400 pasos/consulta)  ->  1.1x   (apenas empata)
%       dt=3e-4 (1333 pasos/consulta)  ->  3.1x
%       dt=1e-4 (4000 pasos/consulta)  ->  6.5x   (techo teorico: 10x)
%   La pregunta que decide NO es "¿es caro mi computo?" sino "¿el trabajo que
%   me ahorro es bastante mayor que unos milisegundos?".
%
%   LA LECCION DE DISENO: el atajo NO es un handler incremental. Necesita leer
%   el objeto (F, DT y la lista de tiempos pedidos), asi que es un SYNC
%   ABSOLUTO @(v,o) -- pero un sync absoluto RECIBE el valor viejo y nadie le
%   obliga a tirarlo: extendSol REUSA los pares ya integrados y solo anade los
%   que falten. "Absoluto" dice contra que estado se pone al dia (el de ahora,
%   N ediciones incluidas), no cuanto trabajo repite.
%
%       P = ivp( @(t,y) cos(t).*y , 1 );    % exacta: y = exp(sin(t))
%       [P,y] = P.At( 4 )                   % integra 0 -> 4
%       [P,y] = P.At( 6 )                   % integra SOLO 4 -> 6
%       [P,y] = P.At( 5 )                   % integra SOLO 4 -> 5, inserta ordenado
%       P = P.SetY0( 2 )                    % cambia el problema: sol cae entera
%
%   LA PUREZA ES UNA PRECONDICION DEL PATRON, Y ESTA CLASE ESTUVO A PUNTO DE
%   INCUMPLIRLA. El compute de una CP tiene que ser funcion PURA de los datos
%   del objeto: VerifyCPs comprueba la cache RECALCULANDO desde cero, asi que
%   si el valor depende ademas de la HISTORIA -- de en que orden se pidieron
%   las cosas -- "la verdad" no esta bien definida y la herramienta canta
%   discrepancias que no lo son.
%   Aqui el peligro era real y sutil. La version obvia integra "el tramo
%   partido en n pasos", asi que la rejilla depende DE DONDE ARRANCAS: pedir
%   y(1) y luego y(1/3) daba 2.3197758575243266 y pedirlos al reves
%   2.3197762151071655 (comprobado con DT=0.1), y VerifyCPs marcaba `sol` como
%   rancia sin que nada estuviera mal. El demo se libraba por los pelos: 4, 5 y
%   6 son multiplos exactos de DT y las rejillas coincidian.
%   La cura es anclar la rejilla al ORIGEN (ver rk4) y AJUSTAR a ella los
%   tiempos pedidos (ver At): asi la y del punto k es siempre "k pasos de
%   tamano DT desde 0", venga de donde venga, y el orden deja de importar.
%   MORALEJA TRANSFERIBLE: si tu compute ACUMULA, preguntate si dos historias
%   distintas con los mismos datos finales dan el mismo valor. Si la respuesta
%   es no, no tienes una CP: tienes estado que se te ha colado de contrabando
%   dentro de la cache.
%
%   Ademas sobrescribe refreshDefs: al cargar un .mat, las definiciones
%   vuelven a ser las DEL CODIGO (sin eso, un fichero viejo usaria para
%   siempre el extendSol con el que se guardo).
%
%   Ver incremental_TUTORIAL.html (los numeros, y donde esta el cruce).
%
% See also cachedOwner, facto, cart, route, serie, ivpDemo, bench_incremental.

  properties (SetAccess = private)
    F  = @(t,y) cos(t) .* y    % la EDO: dy/dt = F(t,y)
    Y0 = 1                     % y en t = 0
    DT = 1e-3                  % paso EXACTO del RK4 (define la rejilla global)
    TS = zeros(1,0)            % tiempos pedidos, YA AJUSTADOS a la rejilla
                               % (> 0, unicos, ordenados)
  end

  properties (Constant, Hidden)
    CP_EVENTS = { 'newtime' , 'changeproblem' }
  end

  methods
    function obj = ivp( f , y0 , dt )
      if nargin >= 1, obj.F  = f;  end
      if nargin >= 2, obj.Y0 = y0; end
      if nargin >= 3, obj.DT = dt; end
      %el constructor NO narra (como route y serie): quien quiera trazas pone
      %o.DEBUG = true despues. Si narrara aqui, cualquier bucle que construya
      %objetos -- bench_incremental, sin ir mas lejos -- se llenaria de DEF.
      obj = obj.registerCPs();
    end

    function [ obj , y , t ] = At( obj , t )
      %[P,y,t] = P.At(t)   y(t), integrando SOLO desde el guardado mas cercano.
      %   El t que se DEVUELVE es el pedido AJUSTADO a la rejilla k*DT (ver el
      %   parrafo de la PUREZA en la cabecera): pedir 0.7 con DT=0.1 da y(0.7),
      %   pero pedir 0.75 da y(0.8). Sin ese ajuste el resultado dependeria del
      %   ORDEN de las peticiones y la CP no seria verificable.
      if ~isscalar( t ) || ~isfinite( t ) || t < 0
        error('ivp:t','t debe ser un escalar finito >= 0.');
      end
      t = round( t / obj.DT ) * obj.DT;      % A LA REJILLA, siempre
      if t > 0 && ~any( obj.TS == t )
        obj.TS = sort( [ obj.TS , t ] );    % ORDENADOS: los tempranos ayudan a los tardios
        obj = obj.fire( 'newtime' );
      end
      S = obj.access( 'sol' );              % dentro de la clase se lee con access
      y = S( S(:,1) == t , 2 );
    end

    function obj = SetY0( obj , y0 )
      %cambia el problema: TODO lo integrado deja de valer
      obj.Y0 = y0;
      obj = obj.fire( 'changeproblem' );
    end
  end

  methods (Hidden)
    function S = extendSol( obj , S )
      %pone la tabla S = [t,y] al dia con obj.TS: integra SOLO los tiempos que
      %falten, cada uno desde el guardado mas cercano por debajo -- que puede
      %ser uno recien insertado en esta misma pasada.
      %Es a la vez el compute (partiendo de [0,Y0]) y el handler de 'newtime'
      %(partiendo del valor viejo): la misma funcion, distinto punto de partida.
      for t = obj.TS
        if any( S(:,1) == t ), continue; end
        i = find( S(:,1) <= t , 1 , 'last' );
        [ y , n ] = obj.rk4( S(i,1) , S(i,2) , t );
        obj.dbg( 'ODE   y(%g): %d paso(s) de RK4 desde t=%g' , t , n , S(i,1) );
        S = [ S( 1:i , : ) ; t , y ; S( i+1:end , : ) ];
      end
    end
  end

  methods (Access = protected)
    function obj = refreshDefs( obj )
      %gancho de EVOLUCION (ver cachedOwner): al cargar, la logica es la del codigo
      obj = obj.registerCPs();
    end
  end

  methods (Access = private)
    function obj = registerCPs( obj )
      %el sync absoluto REUSA el valor viejo (ver la cabecera); [] en
      %changeproblem porque un Y0 nuevo invalida todo lo integrado
      obj = obj.Define( 'sol' , @(o) extendSol( o , [ 0 , o.Y0 ] ) , ...
                        'newtime'       , @(v,o) extendSol( o , v ) , ...
                        'changeproblem' , [] );
    end

    function [ y , n ] = rk4( obj , t0 , y0 , t1 )
      %RK4 clasico sobre la REJILLA GLOBAL k*DT, y esa es toda la diferencia con
      %la version obvia. El paso es SIEMPRE DT (no "el tramo partido en n"), y el
      %tiempo de cada paso se calcula desde el ORIGEN (t = k*DT) y no desde el
      %ancla (t0 + k*h). Asi, avanzar del indice i al j hace EXACTAMENTE las
      %mismas operaciones en el mismo orden vengas de donde vengas, y por
      %induccion y(j) sale bit a bit identico se haya integrado de una tirada o
      %en dos tramos. Es lo que hace que `sol` sea funcion PURA de (F,Y0,DT,TS)
      %-- ver el parrafo de la PUREZA en la cabecera.
      h  = obj.DT;
      k0 = round( t0 / h );  k1 = round( t1 / h );   % indices de rejilla
      n  = k1 - k0;
      f  = obj.F;  y = y0;
      for k = k0:k1-1
        t  = k * h;                                  % DESDE EL ORIGEN
        s1 = f( t       , y            );
        s2 = f( t + h/2 , y + h/2 * s1 );
        s3 = f( t + h/2 , y + h/2 * s2 );
        s4 = f( t + h   , y + h   * s3 );
        y  = y + h/6 * ( s1 + 2*s2 + 2*s3 + s4 );
      end
    end
  end
end
