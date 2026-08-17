classdef collatz < cachedOwner
%collatz  Cadenas de Collatz (3n+1) memoizadas: el ejemplo donde el atajo GANA.
%
%   El problema: la cadena de n es n -> n/2 (par) o 3n+1 (impar), hasta llegar
%   a 1. `lens(n)` es cuantos terminos tiene esa cadena, contando n y el 1.
%   Calcularla para 1..N se memoiza sola: al bajar de n te caes en algo que ya
%   sabes, asi que cada n cuesta solo su tramo propio.
%
%   POR QUE ESTE EJEMPLO EXISTE: es el que cierra el juego de los tres.
%     facto    el atajo PIERDE contra su compute (un cumprod vectorizado)
%     ivp      empata o gana poco, segun el trabajo que haya debajo
%     collatz  el atajo gana 100x, medido
%   La diferencia es que aqui el recompute NO PUEDE HACER TRAMPA: no existe
%   ningun builtin vectorizado que calcule esto, asi que el camino caro tambien
%   es un bucle y el atajo gana en la proporcion que le toca. Esa es la leccion:
%   un atajo gana cuando compite contra otro bucle, no contra una llamada
%   vectorizada de MATLAB.
%
%   ES EXACTO, y eso importa mas de lo que parece. Todo es aritmetica ENTERA
%   (n/2 y 3n+1), asi que el atajo y el recalculo forzado coinciden BIT A BIT
%   -- no "casi" --, y VerifyCPs puede comparar con isequaln sin reservas. En
%   ivp, en cambio, cualquier descuido con la rejilla se paga en los ultimos
%   bits (ver el parrafo de la PUREZA en ivp.m).
%   SIN DESBORDAMIENTO, con margen: los valores INTERMEDIOS de la cadena suben
%   mucho mas que n -- para n <= 1e6 el maximo es 56.991.483.520, en n=704511 --
%   pero eso deja un factor 158.000 hasta 2^53, donde un double deja de contar
%   enteros exactos. Techo recomendado: N <= 1e7. Por encima, comprueba.
%
%   DOS CPs, Y LA SEGUNDA SALE DE LA PRIMERA:
%     lens   la tabla [1..N] de longitudes -- atajo INCREMENTAL en 'grow'
%     best   [n,longitud] del campeon bajo N -- NO declara ningun evento y aun
%            asi cae cuando la tabla crece, porque el motor DESCUBRE que salio
%            de `lens` (mira la columna "sale de:" en la tabla de o.CP).
%   Es el unico ejemplo del repo con una CP COMPUESTA: la maquinaria de
%   dependencias descubiertas solo se veia en las sondas de test.
%
%       C = collatz();
%       [C,c]   = C.Steps( 27 )       % la cadena de 27 (famosa por larga)
%       [C,n,c] = C.Best( 1000 )      % el campeon bajo 1000
%       [C,n,c] = C.Best( 10000 )     % lens se EXTIENDE; best CAE y se rehace
%
%   Ver incremental_TUTORIAL.html (los tres ejemplos, y donde esta el cruce).
%
% See also cachedOwner, facto, ivp, cart, route, serie, collatzDemo,
%          bench_incremental.

  properties (SetAccess = private)
    N = 1        % hasta donde llega la tabla (solo CRECE)
  end

  properties (Constant, Hidden)
    CP_EVENTS = { 'grow' }
  end

  methods
    function obj = collatz()
      %el constructor NO narra (como route, serie, facto e ivp): las trazas las
      %enciende quien las quiera con o.DEBUG = true
      obj = obj.registerCPs();
    end

    function obj = Reserve( obj , n )
      %apunta que hara falta llegar a n SIN calcular nada
      if ~isscalar( n ) || n < 1 || n ~= round( n ) || ~isfinite( n )
        error('collatz:n','n debe ser un entero >= 1.');
      end
      if n > 1e7
        warning('collatz:big', ...
                [ 'n = %g: por encima de ~1e7 las excursiones de la cadena se acercan a 2^53 ' ...
                  'y el double deja de contar enteros exactos. Comprueba antes de fiarte.' ] , n );
      end
      if n > obj.N
        a = obj.N;  obj.N = n;
        obj = obj.fire( 'grow' , a , n );   % args: DESDE donde y HASTA donde
      end
    end

    function [ obj , c ] = Steps( obj , n )
      %[C,c] = C.Steps(n)   terminos de la cadena de n (contando n y el 1)
      obj = obj.Reserve( n );
      v   = obj.access( 'lens' );      % dentro de la clase se lee con access
      c   = v( n );
    end

    function [ obj , n , c ] = Best( obj , N )
      %[C,n,c] = C.Best(N)   el n <= N con la cadena mas larga, y su longitud
      obj = obj.Reserve( N );
      b   = obj.access( 'best' );
      n   = b(1);  c = b(2);
    end
  end

  methods (Access = protected)
    function obj = refreshDefs( obj )
      %gancho de EVOLUCION (ver cachedOwner): al cargar, manda el CODIGO
      obj = obj.registerCPs();
    end
  end

  methods (Access = private)
    function obj = registerCPs( obj )
      %`lens` con atajo INCREMENTAL de verdad: el handler recibe la tabla vieja
      %(que es a la vez el valor Y el memo que necesita) y el tramo que falta,
      %asi que NO necesita mirar el objeto -- por eso el ~ del segundo argumento
      obj = obj.Define( 'lens' , @(o) collatz.extend( [] , 0 , o.N ) , ...
                        'grow' , @(v,~,a,b) collatz.extend( v , a , b ) );
      %`best` NO declara eventos, y no le hace falta: sale de `lens`, el motor
      %lo descubre al calcularla, y cuando `lens` cambia esta cae con ella
      obj = obj.Define( 'best' , @(o) collatz.champion( o.access('lens') ) );
    end
  end

  methods (Static, Hidden)
    function v = extend( v , a , b )
      %extiende la tabla de longitudes de a+1 hasta b REUSANDO v (que ya cubre
      %1..a). Es a la vez el compute (con v vacia) y el handler de 'grow' (con
      %la tabla vieja): la misma funcion, distinto punto de partida.
      %EL INVARIANTE que lo hace barato: se recorre n ASCENDENTE, asi que en
      %cuanto la cadena baja por debajo de n ya sabemos su longitud -- de ahi
      %que el bucle sea `while x >= n` y no una caminata hasta el 1.
      if isempty( v ), v = 1;  a = max( a , 1 ); end   % v(1) = 1: la cadena de 1 es {1}
      if numel( v ) < b, v( b ) = 0; end               % crecer de UNA vez
      for n = a+1 : b
        x = n;  c = 0;
        while x >= n
          if mod( x , 2 ) == 0, x = x / 2; else, x = 3*x + 1; end
          c = c + 1;
        end
        v( n ) = c + v( x );        % lo andado + lo que ya sabiamos
      end
    end

    function b = champion( v )
      %[n,longitud] del maximo. max() devuelve la PRIMERA aparicion, asi que
      %ante un empate gana el n mas pequeno: deterministico, que es lo que
      %necesita un valor cacheado (y un golden)
      [ c , n ] = max( v );
      b = [ n , c ];
    end
  end
end
