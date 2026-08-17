%factoDemo  Sesion paso a paso del ejemplo TODO-INCREMENTAL (facto).
%
%   La CP `facts` solo CRECE: cada paso ensena que pedir mas alla extiende la
%   tabla con los productos que faltan y ni uno mas, y que pedir lo que ya
%   esta es un HIT. Comparar con ivpDemo, donde el mismo gesto se implementa
%   como sync absoluto que reusa el valor viejo.
%
%   OJO SI LO EDITAS: la salida es un GOLDEN -- test_cache la compara byte a
%   byte contra factoDemo_expected.txt. test_cache('record') regenera el
%   golden, SOLO tras un cambio deliberado.
%
%   Uso:  addpath('C:\repos\msh','C:\repos\msh\dev'); factoDemo
%
% See also facto, ivpDemo, cartDemo, cachedOwner.

function factoDemo()

  sep( 0 , 'F = facto()   -- registra `facts` (incremental en ''grow''); cache VACIA' );
  F = facto();  F.DEBUG = true;        % las trazas [facto] las pide el demo
  disp( F.CP );

  sep( 1 , '[F,x] = F.Fact(5)   -- MISS: computa 0!..5! de una vez' );
  [ F , x ] = F.Fact( 5 );
  fprintf( '  5! = %d\n' , x );

  sep( 2 , '[F,x] = F.Fact(7)   -- replay INCREMENTAL: solo los productos 6 y 7' );
  [ F , x ] = F.Fact( 7 );
  fprintf( '  7! = %d\n' , x );

  sep( 3 , '[F,x] = F.Fact(3)   -- HIT: ya estaba en la tabla' );
  [ F , x ] = F.Fact( 3 );
  fprintf( '  3! = %d\n' , x );

  sep( 4 , 'F.Reserve(9) y F.Reserve(12)   -- dos ediciones SIN leer: log de 2' );
  F = F.Reserve( 9 );
  F = F.Reserve( 12 );
  disp( F.CP );

  sep( 5 , '[F,x] = F.Fact(12)   -- el replay repite las 2 ediciones EN ORDEN' );
  [ F , x ] = F.Fact( 12 );
  fprintf( '  12! = %d\n' , x );

  sep( 6 , 'F.facts contra F.facts_ (recalculo forzado): el atajo no miente' );
  fprintf( '  coinciden: %d\n' , isequal( F.facts , F.facts_ ) );

  fprintf( '\n' );
end

function sep( n , t )
  fprintf( '\n%s\n== PASO %d. %s\n' , repmat('-',1,76) , n , t );
end
