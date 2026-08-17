function R = bench_incremental( rounds )
%BENCH_INCREMENTAL  Cuanto ahorra de verdad el atajo incremental (facto e ivp).
%
%   bench_incremental()      imprime la tabla
%   R = bench_incremental()  ademas la devuelve (struct con los tiempos)
%
%   Existe porque el tutorial afirma que `ivp` SI paga el peaje del patron y
%   `facto` no, y esas dos frases tienen que ser comprobables en la maquina de
%   quien lea. Mide tres cosas por clase:
%     recompute  leer teniendo que calcularlo todo desde cero
%     atajo      leer teniendo ya la mitad hecha (replay incremental o sync)
%     HIT        leer algo que ya esta (solo el peaje del despacho)
%   ...y, en ivp, el mismo RK4 escrito como FUNCION SUELTA, para separar el
%   coste del patron del coste de escribir el integrador dentro de una clase.
%
%   RONDAS INTERCALADAS Y MEDIANA, igual que el test de rendimiento de la
%   bateria y por el mismo motivo: en una sola toma estos numeros se mueven
%   facil un 50% (medido; la primera timeit de la sesion paga el calentamiento).
%
%   Uso:  addpath('C:\repos\msh','C:\repos\msh\dev'); bench_incremental
%
% See also facto, ivp, incremental_TUTORIAL.html, test_cache.

  if nargin < 1, rounds = 3; end
  f  = @(t,y) cos(t) .* y;
  NQ = 20;                       % consultas crecientes del escenario realista

  NC = 1e5;                      % tamano de la tabla de collatz
  fa = zeros(rounds,5);  iv = zeros(rounds,6);  co = zeros(rounds,4);
  for r = 1:rounds
    %-------- collatz: tabla 1..NC, o extenderla los ultimos 2000
    C0 = collatz();                                        % vacia
    C1 = collatz();  C1 = C1.Reserve( NC-2000 );  v1 = C1.lens;  %#ok<NASGU>
    C2 = C1;         C2 = C2.Reserve( NC );       v2 = C2.lens;  %#ok<NASGU>
    co(r,1) = timeit( @() steps( C0 , NC ) );              % recompute
    co(r,2) = timeit( @() steps( C1 , NC ) );              % atajo
    co(r,3) = timeit( @() steps( C2 , NC ) );              % HIT
    %el mismo trabajo SIN el patron: la tabla desde cero, sin clase ni cache
    co(r,4) = timeit( @() collatz.extend( [] , 0 , NC ) );
    %-------- facto: tabla hasta 2000, o extenderla de 1000 a 2000
    %(facto/ivp nacen con DEBUG apagado: no hay que tocarlo, y menos DENTRO de
    %lo que se mide -- una escritura externa cuesta ~0.5 ms de subsasgn)
    F0 = facto();                                          % vacia
    F1 = facto();  F1 = F1.Fact( 1000 );                   % ya llega a 1000!
    F2 = F1;       F2 = F2.Fact( 2000 );                   % ya llega a 2000!
    fa(r,1) = timeit( @() fact( F0 , 2000 ) );             % recompute
    fa(r,2) = timeit( @() fact( F1 , 2000 ) );             % atajo
    fa(r,3) = timeit( @() fact( F2 , 2000 ) );             % HIT
    %el escenario REAL: NQ consultas crecientes reusando el objeto (la tabla se
    %extiende) frente a LO QUE ESCRIBIRIAS SIN EL PATRON (calcular cada
    %respuesta desde cero, sin clase ni cache). OJO: el "sin patron" NO puede
    %ser "construir un objeto nuevo cada vez" -- eso mide sobre todo el
    %`Define`, que hace introspeccion de metaclass y cuesta milisegundos; seria
    %cargarle al baseline un coste que no tiene nada que ver con cachear.
    fa(r,4) = timeit( @() seqFact(  NQ ) );
    fa(r,5) = timeit( @() seqFactRaw( NQ ) );

    %-------- ivp: integrar hasta t=8, o extender desde t=4
    P0 = ivp( f , 1 );                                     % nada integrado
    P1 = ivp( f , 1 );  P1 = P1.At( 4 );                   % 0..4 hecho
    P2 = P1;  P2 = P2.At( 8 );                             % 0..8 hecho
    iv(r,1) = timeit( @() at( P0 , 8 ) );                  % recompute
    iv(r,2) = timeit( @() at( P1 , 8 ) );                  % atajo
    iv(r,3) = timeit( @() at( P2 , 8 ) );                  % HIT
    iv(r,4) = timeit( @() rk4raw( f , 1 , 0 , 8 , 1e-3 ) );% el RK4 a pelo
    iv(r,5) = timeit( @() seqAt(    f , NQ ) );
    iv(r,6) = timeit( @() seqAtRaw( f , NQ ) );
  end
  fa = median( fa , 1 );  iv = median( iv , 1 );  co = median( co , 1 );

  R = struct( 'facto' , struct( 'recompute',fa(1) , 'atajo',fa(2) , 'hit',fa(3) , ...
                                'seqCache',fa(4) , 'seqFresh',fa(5) ) , ...
              'ivp'   , struct( 'recompute',iv(1) , 'atajo',iv(2) , 'hit',iv(3) , 'raw',iv(4) , ...
                                'seqCache',iv(5) , 'seqFresh',iv(6) ) );

  fprintf( '\n  bench_incremental  (mediana de %d rondas intercaladas, R%s)\n\n' , ...
           rounds , version('-release') );
  fprintf( '  UNA lectura\n' );
  fprintf( '  %-34s %10s %10s %10s\n' , '' , 'recompute' , 'atajo' , 'HIT' );
  row( 'facto: tabla hasta 2000!' , fa(1:3) );
  row( 'ivp:   y(8), dt=1e-3 (8000 pasos)' , iv(1:3) );
  row( sprintf('collatz: tabla 1..%g (+2000)',NC) , co(1:3) );
  fprintf( '\n  %d CONSULTAS CRECIENTES: el patron contra NO tenerlo\n' , NQ );
  fprintf( '  %-34s %10s %10s %10s\n' , '' , 'con CP' , 'sin patron' , 'ahorro' );
  fprintf( '  %-34s %9.1fms %9.1fms %9.1fx\n' , 'facto: 100!, 200!, ... 2000!' , ...
           fa(4)*1e3 , fa(5)*1e3 , fa(5)/fa(4) );
  fprintf( '  %-34s %9.1fms %9.1fms %9.1fx\n' , 'ivp:   y(0.4), y(0.8), ... y(8)' , ...
           iv(5)*1e3 , iv(6)*1e3 , iv(6)/iv(5) );
  %----- DONDE ESTA EL PUNTO DE CRUCE. Lo de arriba dice "gana/pierde" en UN
  %ajuste; esto dice POR QUE: el peaje del patron es FIJO por lectura (~1-3 ms
  %entre despacho, fire y replay) y el ahorro crece con el trabajo que hay
  %debajo. Barriendo DT se ve la asintota: con 20 consultas crecientes el
  %trabajo evitado es ~10x, y el ratio medido tiende a el segun el trabajo por
  %consulta supera al peaje. UNA ronda: timeit ya promedia por dentro y con
  %DT pequeno cada llamada cuesta segundos.
  DTS = [ 1e-3 , 3e-4 , 1e-4 ];
  sw  = zeros( numel(DTS) , 3 );
  for j = 1:numel( DTS )
    sw(j,1) = DTS(j);
    sw(j,2) = timeit( @() seqAt(    f , NQ , DTS(j) ) );
    sw(j,3) = timeit( @() seqAtRaw( f , NQ , DTS(j) ) );
  end
  R.sweep = sw;
  fprintf( '\n  ivp: %d consultas crecientes, barriendo el paso (donde esta el cruce)\n' , NQ );
  fprintf( '  %-14s %10s %10s %10s %10s\n' , 'dt' , 'pasos' , 'con CP' , 'sin patron' , 'ahorro' );
  for j = 1:size(sw,1)
    fprintf( '  %-14g %10d %9.1fms %9.1fms %9.1fx\n' , sw(j,1) , round(0.4/sw(j,1)) , ...
             sw(j,2)*1e3 , sw(j,3)*1e3 , sw(j,3)/sw(j,2) );
  end

  fprintf( '\n  facto  el ATAJO cuesta %.2fx lo que su propio recompute (el compute es un\n' , ...
           fa(2)/fa(1) );
  fprintf( '         cumprod VECTORIZADO y el atajo, que concatena la tabla entera, compite\n' );
  fprintf( '         contra el en desventaja) y el patron entero cuesta %.0fx lo que resolver\n' , ...
           fa(4)/fa(5) );
  fprintf( '         las %d consultas a pelo -> NO PAGA, ni el atajo ni cachear.\n' , NQ );
  fprintf( '  ivp    con dt=1e-3 apenas empata (%.1fx): el trabajo por consulta (~%.1f ms) no\n' , ...
           sw(1,3)/sw(1,2) , sw(1,2)*1e3/NQ );
  fprintf( '         llega a cubrir el peaje. Al bajar dt el ahorro tiende a ~%.0fx, que es el\n' , ...
           NQ*(NQ+1)/2/NQ );
  fprintf( '         trabajo que de verdad se evita -> PAGA cuando hay bastante debajo.\n' );
  fprintf( '  ivp    el mismo RK4 como funcion suelta: %.1f ms (integrar 0..8) frente a %.1f ms\n' , ...
           iv(4)*1e3 , iv(1)*1e3 );
  fprintf( '         por el camino completo: meter el integrador en un classdef no penaliza,\n' );
  fprintf( '         lo que se paga es el despacho, UNA vez por lectura y no por paso.\n' );
  fprintf( '  colltz el atajo ahorra %.0fx (%.1f ms -> %.2f ms) extendiendo los ultimos 2000 de\n' , ...
           co(1)/co(2) , co(1)*1e3 , co(2)*1e3 );
  fprintf( '         %g. Y aqui el recompute NO PUEDE HACER TRAMPA: no hay builtin vectorizado\n' , NC );
  fprintf( '         que calcule esto (la tabla a pelo cuesta %.1f ms, igual que por la clase),\n' , co(4)*1e3 );
  fprintf( '         asi que el camino caro tambien es un bucle -> EL ATAJO GANA DE VERDAD.\n\n' );
  if nargout == 0, clear R; end
end

function seqFact( nq )
  %nq consultas CRECIENTES reusando el objeto: la tabla se va extendiendo
  F = facto();
  for k = 1:nq, [ F , ~ ] = F.Fact( 100*k ); end
end

function seqFactRaw( nq )
  %las mismas nq respuestas SIN el patron: cada una desde cero
  for k = 1:nq
    t = cumprod( [ 1 , 1:100*k ] );  x = t(end);  %#ok<NASGU>
  end
end

function seqAt( f , nq , dt )
  if nargin < 3, dt = 1e-3; end
  P = ivp( f , 1 , dt );
  for k = 1:nq, [ P , ~ ] = P.At( 0.4*k ); end
end

function seqAtRaw( f , nq , dt )
  %las mismas nq respuestas SIN el patron: integrar 0 -> t cada vez
  if nargin < 3, dt = 1e-3; end
  for k = 1:nq, y = rk4raw( f , 1 , 0 , 0.4*k , dt ); end                  %#ok<NASGU>
end

function row( name , t )
  fprintf( '  %-34s %9.2fms %9.2fms %9.3fms\n' , name , t(1)*1e3 , t(2)*1e3 , t(3)*1e3 );
end

function x = fact(  F , n ), [ ~ , x ] = F.Fact(  n ); end
function x = steps( C , n ), [ ~ , x ] = C.Steps( n ); end
function y = at(   P , t ), [ ~ , y ] = P.At(   t ); end

function y = rk4raw( f , y0 , t0 , t1 , dt )
  %el mismo integrador de ivp.rk4, pero suelto: separa el coste del PATRON del
  %coste de escribir el integrador dentro de un classdef
  n = max( 1 , ceil( ( t1 - t0 ) / dt ) );
  h = ( t1 - t0 ) / n;  y = y0;
  for k = 0:n-1
    t  = t0 + k*h;
    k1 = f( t       , y            );
    k2 = f( t + h/2 , y + h/2 * k1 );
    k3 = f( t + h/2 , y + h/2 * k2 );
    k4 = f( t + h   , y + h   * k3 );
    y  = y + h/6 * ( k1 + 2*k2 + 2*k3 + k4 );
  end
end
