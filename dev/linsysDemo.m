%linsysDemo  Sesion paso a paso del ejemplo mas practico (linsys).
%
%   Ensena tres cosas que ningun otro demo junta:
%     - una CP cuyo valor es un objeto CARO Y OPACO (una factorizacion), que es
%       el analogo directo del `bvh` de msh;
%     - que leer NO obliga a reasignar cuando no se editan datos: `x = L\b`
%       deja la factorizacion caliente porque el almacen es un handle
%       compartido;
%     - operadores sobrecargados (\, / e inv) apoyados en la misma CP.
%
%   OJO SI LO EDITAS: la salida es un GOLDEN -- test_cache la compara byte a
%   byte contra linsysDemo_expected.txt. Por eso NO se imprime ningun tiempo
%   (variarian en cada maquina) y los residuos salen como booleanos contra una
%   tolerancia, no como numeros.
%
%   Uso:  addpath('C:\repos\msh','C:\repos\msh\dev'); linsysDemo
%
% See also linsys, incremental_TUTORIAL.html, collatzDemo.

function linsysDemo()

  rng( 0 );                                  % deterministico: esto es un golden
  n = 400;
  A = randn( n ) + n*eye( n );               % bien condicionada
  b = randn( n , 1 );  b2 = randn( n , 1 );  B = randn( 3 , n );

  sep( 0 , 'L = linsys(A)   -- dos CPs: dcmp (la factorizacion) e invA (compuesta)' );
  L = linsys( A );  L.DEBUG = true;          % las trazas [linsys] las pide el demo
  disp( L.CP );

  sep( 1 , 'x = L\b   -- MISS: factoriza (O(n^3)) y guarda' );
  x = L \ b;
  fprintf( '  residuo |A*x-b|/|b| < 1e-10 : %d\n' , norm(A*x-b)/norm(b) < 1e-10 );

  sep( 2 , 'y = L\b2   -- HIT: NO vuelve a factorizar (esto es todo el asunto)' );
  y = L \ b2;
  fprintf( '  residuo |A*y-b2|/|b2| < 1e-10 : %d\n' , norm(A*y-b2)/norm(b2) < 1e-10 );

  sep( 3 , 'z = B/L   -- b/A con la MISMA factorizacion (via d''), otro HIT' );
  z = B / L;
  fprintf( '  residuo |z*A-B|/|B| < 1e-10 : %d\n' , norm(z*A-B,'fro')/norm(B,'fro') < 1e-10 );

  sep( 4 , 'la tabla: dcmp fresh, invA aun sin calcular' );
  disp( L.CP );

  sep( 5 , 'W = inv(L)   -- MISS de invA, que SALE de dcmp (mira "sale de:")' );
  W = inv( L );
  fprintf( '  |A*W - I| < 1e-10 : %d\n' , norm( A*W - eye(n) , 'fro' ) < 1e-10 );
  disp( L.CP );

  sep( 6 , 'inv(L) otra vez   -- HIT: la inversa ya estaba' );
  W2 = inv( L );
  fprintf( '  identica a la anterior: %d\n' , isequal( W , W2 ) );

  sep( 7 , 'M = L; L = L.SetA(A2)   -- dcmp cae, invA cae con ella; M no se entera' );
  M  = L;
  A2 = randn( n ) + n*eye( n );
  L  = L.SetA( A2 );
  disp( L.CP );

  sep( 8 , 'L\b usa la matriz NUEVA; M\b sigue con la vieja (COW)' );
  xn = L \ b;
  xm = M \ b;
  fprintf( '  L resuelve con A2 : %d\n' , norm(A2*xn-b)/norm(b) < 1e-10 );
  fprintf( '  M resuelve con A  : %d\n' , norm(A *xm-b)/norm(b) < 1e-10 );

  fprintf( '\n' );
end

function sep( n , t )
  fprintf( '\n%s\n== PASO %d. %s\n' , repmat('-',1,76) , n , t );
end
