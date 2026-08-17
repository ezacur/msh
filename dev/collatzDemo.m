%collatzDemo  Sesion paso a paso del ejemplo donde el atajo GANA (collatz).
%
%   Dos cosas que ningun otro demo ensena a la vez:
%     - un atajo incremental que de verdad compensa (aqui se ve el mecanismo;
%       los numeros los da bench_incremental: 100x extendiendo una tabla de 1e6)
%     - una CP COMPUESTA: `best` sale de `lens`, no declara ningun evento, y
%       aun asi cae cuando la tabla crece -- el motor lo descubre solo.
%
%   OJO SI LO EDITAS: la salida es un GOLDEN -- test_cache la compara byte a
%   byte contra collatzDemo_expected.txt. test_cache('record') regenera el
%   golden, SOLO tras un cambio deliberado. Por eso los tamanos son pequenos:
%   un demo tiene que ser instantaneo y deterministico.
%
%   Uso:  addpath('C:\repos\msh','C:\repos\msh\dev'); collatzDemo
%
% See also collatz, incremental_TUTORIAL.html, factoDemo, ivpDemo.

function collatzDemo()

  sep( 0 , 'C = collatz()   -- dos CPs: `lens` (incremental) y `best` (SIN eventos)' );
  C = collatz();  C.DEBUG = true;        % las trazas [collatz] las pide el demo
  disp( C.CP );

  sep( 1 , '[C,c] = C.Steps(27)   -- MISS: memoiza 1..27 de una vez' );
  [ C , c ] = C.Steps( 27 );
  fprintf( '  la cadena de 27 tiene %d terminos\n' , c );

  sep( 2 , '[C,n,c] = C.Best(1000)   -- lens se EXTIENDE (28..1000); best MISS' );
  [ C , n , c ] = C.Best( 1000 );
  fprintf( '  campeon bajo 1000: n=%d con %d terminos\n' , n , c );

  sep( 3 , 'la tabla: fijate en "sale de: lens" -- nadie lo declaro' );
  disp( C.CP );

  sep( 4 , '[C,c] = C.Steps(500)   -- HIT: ya estaba, cero trabajo' );
  [ C , c ] = C.Steps( 500 );
  fprintf( '  la cadena de 500 tiene %d terminos\n' , c );

  sep( 5 , '[C,n,c] = C.Best(10000)   -- lens PENDIENTE; best CAE con ella' );
  [ C , n , c ] = C.Best( 10000 );
  fprintf( '  campeon bajo 10000: n=%d con %d terminos\n' , n , c );

  sep( 6 , 'C.Reserve(20000) y C.Reserve(30000)   -- dos ediciones sin leer' );
  C = C.Reserve( 20000 );
  C = C.Reserve( 30000 );
  disp( C.CP );

  sep( 7 , 'C.best   -- el replay repite las 2 ediciones sobre lens, y best se rehace' );
  b = C.best;
  fprintf( '  campeon bajo 30000: n=%d con %d terminos\n' , b(1) , b(2) );

  sep( 8 , 'ES EXACTO: el atajo coincide BIT A BIT con el recalculo forzado' );
  fprintf( '  lens identica al recompute: %d\n' , isequal( C.lens , C.lens_ ) );
  fprintf( '  VerifyCPs no marca nada:    %d\n' , isempty( C.VerifyCPs() ) );

  fprintf( '\n' );
end

function sep( n , t )
  fprintf( '\n%s\n== PASO %d. %s\n' , repmat('-',1,76) , n , t );
end
