classdef linsys < cachedOwner
%linsys  A*x = b con la FACTORIZACION cacheada: el ejemplo mas practico.
%
%   Resolver un sistema cuesta O(n^3) la primera vez y O(n^2) las siguientes,
%   PERO solo si te guardas la factorizacion. Eso es exactamente una CP:
%       dcmp   decomposition(A) -- cara, opaca, y se invalida si A cambia
%       invA   la inversa explicita, que SALE de dcmp (CP compuesta)
%   Medido con n=800: factorizar 13.5 ms, resolver con la factorizacion
%   cacheada 0.23 ms, y A\b desde cero 9.3 ms -> 40x por cada termino
%   independiente. Es el analogo directo del `bvh` de msh: una estructura cara
%   que se construye una vez y se consulta muchas.
%
%   COMO SE USA -- y fijate en que NO HAY QUE REASIGNAR:
%       L = linsys( A );
%       x = L \ b;              % MISS: factoriza y guarda
%       y = L \ b2;             % HIT: reusa la factorizacion
%       z = B / L;              % b/A, con la MISMA factorizacion (via d')
%       W = inv( L );           % CP compuesta invA (usala poco: ver abajo)
%       L = L.SetA( A2 );       % esto SI edita: dcmp cae, invA cae con ella
%
%   POR QUE `L \ b` NO DEVUELVE EL OBJETO, cuando facto/ivp/collatz si tenian
%   que hacerlo: aquellos EDITABAN datos al preguntar (subian N, anadian un
%   tiempo), y editar en una clase de valor obliga a reasignar. Aqui resolver
%   no toca los datos: solo LEE un derivado. Y el valor calculado se deposita
%   en el cacheHandle, que las copias COMPARTEN, asi que la factorizacion queda
%   caliente para el llamante sin devolver nada. Es la mejor demostracion de
%   por que el almacen es un handle y no una propiedad mas.
%
%   LA LECCION QUE NO ENSENA NINGUN OTRO EJEMPLO: hay valores que NO PUEDEN
%   viajar en un .mat. MATLAB se niega a serializar un `decomposition` -- avisa
%   "Saving a decomposition is not supported" y al cargar te devuelve, casi en
%   silencio, la factorizacion de una MATRIZ VACIA. Quien la guardara en una
%   propiedad normal tendria un objeto roto en el fichero sin enterarse. Como la
%   cache del patron es Transient, aqui el problema NO PUEDE DARSE: la
%   factorizacion nunca viaja y se rehace en la primera lectura tras el load.
%   La cache transitoria no es una optimizacion; es el UNICO sitio correcto
%   donde este valor puede vivir.
%
%   Y UNA POLITICA SIN ATAJO, que tambien hay que saber elegir: `dcmp` declara
%   'changematrix' -> [] y punto. No hay forma de "poner al dia" una
%   factorizacion cuando cambia la matriz; el atajo no existe, y forzar uno
%   seria inventarselo. No toda CP necesita un handler listo.
%
%   NO SOBRECARGUES `size` EN UNA CLASE DE ESTAS. Seria comodo que size(L)
%   diera el tamano de A, pero `isscalar` se calcula CON size, y el subsref de
%   cachedOwner empieza comprobando isscalar(obj): con size sobrecargado, toda
%   lectura de CP moriria en cachedOwner:notScalar. Si necesitas el tamano,
%   L.A da la matriz (y `inv` ya lo resuelve por dentro).
%
%   Ver incremental_TUTORIAL.html y la §9 de cache_TUTORIAL.html (ProfileCPs
%   sobre esta clase da el veredicto mas rotundo del repo).
%
% See also cachedOwner, collatz, facto, ivp, cart, route, serie, linsysDemo.

  properties (SetAccess = private)
    A = []       % la matriz del sistema
  end

  properties (Constant, Hidden)
    CP_EVENTS = { 'changematrix' }
  end

  methods
    function obj = linsys( A )
      obj = obj.registerCPs();
      if nargin >= 1, obj = obj.SetA( A ); end
    end

    function obj = SetA( obj , A )
      %cambia la matriz: TODO lo factorizado deja de valer
      if ~isnumeric( A ) || ~ismatrix( A ) || isempty( A )
        error('linsys:A','A debe ser una matriz numerica no vacia.');
      end
      obj.A = A;
      obj = obj.fire( 'changematrix' );
    end

    function x = Solve( obj , b )
      %x = L.Solve(b)   resuelve A*x = b reusando la factorizacion cacheada.
      %NO devuelve el objeto a proposito: no edita datos, solo lee un derivado,
      %y el valor se deposita en el handle COMPARTIDO (ver la cabecera).
      x = obj.access( 'dcmp' ) \ b;     % dentro de la clase se lee con access
    end

    %% ---------------------------------------------------- operadores
    function x = mldivide( a , b )
      %L \ b  ->  Solve(b). El objeto es el primer argumento.
      if ~isa( a , 'linsys' )
        error('linsys:mldivide','use L\\b con L a la izquierda.');
      end
      x = Solve( a , b );
    end

    function x = mrdivide( a , b )
      %b / L  ->  x*A = b, o sea A'*x' = b'. OJO: aqui el objeto es el SEGUNDO
      %argumento, que es la trampa clasica al sobrecargar mrdivide.
      %Y se reusa LA MISMA factorizacion: `d'\y` esta soportado (comprobado),
      %aunque `transpose(d)` a secas no lo este
      %(MATLAB:decomposition:TransposeNotSupported). Nada de factorizar A'.
      if ~isa( b , 'linsys' )
        error('linsys:mrdivide','use b/L con L a la derecha.');
      end
      x = ( b.access( 'dcmp' )' \ a' )';
    end

    function X = inv( obj )
      %inv(L)  ->  la CP compuesta invA (A \ eye(n)), cacheada.
      %NUMERICAMENTE, casi nunca es lo que quieres: para resolver un sistema usa
      %L\b, que es mas rapido y mas estable. Esta aqui porque a veces se
      %necesita la inversa de verdad -- y porque es una CP COMPUESTA de libro.
      X = obj.access( 'invA' );
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
      %SIN ATAJO POSIBLE: una factorizacion no se "pone al dia" cuando cambia
      %la matriz, asi que [] es la politica correcta y no una dejadez
      obj = obj.Define( 'dcmp' , @(o) decomposition( o.A ) , ...
                        'changematrix' , [] );
      %invA SALE de dcmp -- el motor lo descubre y la tira con ella --, pero
      %ademas DECLARA el evento. No es redundancia inutil: la dependencia es la
      %red de seguridad y la declaracion es la intencion escrita, y de paso
      %calla el aviso noEvents (su compute lee o.A para el tamano, y la
      %heuristica no puede saber que la dependencia ya la cubre).
      obj = obj.Define( 'invA' , @(o) o.access('dcmp') \ eye( size( o.A , 1 ) ) , ...
                        'changematrix' , [] );
    end
  end
end
