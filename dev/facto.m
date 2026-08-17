classdef facto < cachedOwner
%facto  Factoriales memoizados: el ejemplo TODO-INCREMENTAL del patron.
%
%   La CP `facts` guarda la tabla [0!, 1!, ..., N!] y N solo CRECE: pedir un
%   factorial que ya esta es un HIT, y pedir uno mas alla EXTIENDE la tabla
%   multiplicando solo los factores que faltan. El handler incremental de
%   'grow' recibe (a,b) -- "llegaba a a! y tiene que llegar a b!" -- y parte
%   de v(end), asi que F.Fact(7) despues de F.Fact(5) cuesta DOS productos,
%   no siete. Es la forma pura del atajo: los args del evento llevan TODO lo
%   que el handler necesita, y una cadena de ediciones sin lecturas entre
%   medias se repite en orden con los (a,b) de cada una.
%
%       F = facto();
%       [F,x] = F.Fact(5)     % MISS: computa 0!..5! y los guarda
%       [F,x] = F.Fact(7)     % replay incremental: SOLO 6 y 7
%       [F,x] = F.Fact(3)     % HIT: ya estaba
%       F = F.Reserve(12)     % apunta la edicion SIN leer (queda pendiente)
%
%   OJO, Y NO ES UNA ADVERTENCIA DE ADORNO: esta clase es un ANTI-EJEMPLO
%   medido. `bench_incremental` dice que resolver 20 factoriales por aqui
%   cuesta 163x lo que resolverlos a pelo (26 ms contra 0.2 ms), y que el
%   ATAJO en si cuesta MAS que el compute que sustituye -- porque el compute es
%   un cumprod VECTORIZADO y el atajo, que concatena la tabla entera, compite
%   contra el en desventaja. Sirve para ENSENAR el mecanismo con valores que
%   caben en la cabeza; no para copiarlo. El que si puede pagar es ivp (la
%   misma idea con una EDO), y su tutorial mide donde esta la raya.
%
%   Ver incremental_TUTORIAL.html (los numeros, y donde esta el cruce).
%
% See also cachedOwner, ivp, cart, route, serie, factoDemo, bench_incremental.

  properties (SetAccess = private)
    N = 0        % el factorial mas alto ya apuntado
  end

  properties (Constant, Hidden)
    CP_EVENTS = { 'grow' }
  end

  methods
    function obj = facto()
      %el constructor NO narra (como route y serie): quien quiera trazas pone
      %o.DEBUG = true despues. Si narrara aqui, cualquier bucle que construya
      %objetos -- bench_incremental, sin ir mas lejos -- se llenaria de DEF.
      %facts(k) = (k-1)!, y facts(1) = 0! = 1 SIEMPRE: asi v nunca esta vacia
      %y el handler puede partir de v(end) sin caso especial
      obj = obj.Define( 'facts' , @(o) cumprod( [ 1 , 1:o.N ] ) , ...
                        'grow'  , @(v,~,a,b) [ v , v(end) * cumprod( a+1 : b ) ] );
    end

    function obj = Reserve( obj , n )
      %apunta que hara falta llegar a n! SIN leer: la tabla queda pendiente y
      %el log acumula ediciones que el proximo Fact repetira en orden
      if ~isscalar( n ) || n < 0 || n ~= round( n )
        error('facto:n','n debe ser un entero >= 0.');
      end
      if n > obj.N
        a = obj.N;  obj.N = n;
        obj = obj.fire( 'grow' , a , n );
      end
    end

    function [ obj , f ] = Fact( obj , n )
      %[F,f] = F.Fact(n)   devuelve n!, calculando SOLO lo que falte
      obj = obj.Reserve( n );
      t = obj.access( 'facts' );            % dentro de la clase se lee con access
      f = t( n+1 );
    end
  end
end
