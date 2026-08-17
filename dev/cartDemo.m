%cartDemo  Sesion paso a paso del patron de cache (17 pasos, con trazas).
%
%   Ejemplo DIDACTICO, independiente de msh: el dueno es cart (carrito de la
%   compra). Cada PASO imprime lo que hace y las trazas [cart] del motor;
%   la seccion 5 de cache_TUTORIAL.html (raiz del repo) es esta misma sesion
%   comentada paso a paso, con el estado del store dibujado en cada uno --
%   correr esto con el tutorial al lado es la mejor forma de aprender el patron.
%
%   OJO SI LO EDITAS: la salida de esta sesion es un GOLDEN -- test_cache la
%   compara byte a byte contra cartDemo_expected.txt. Cambiar cualquier texto
%   impreso rompe el test (a proposito: asi se detecta cualquier cambio de
%   comportamiento del motor). test_cache('record') regenera el golden, SOLO
%   tras un cambio deliberado.
%
%   Uso:  addpath('C:\repos\msh','C:\repos\msh\dev'); cartDemo
%
% See also cart, cachedOwner, cacheHandle, cacheProxy, routeDemo.

function cartDemo()

  sep( 0 , 'C = cart().Demo()   -- registra 4 definiciones; la cache queda VACIA' );
  C = cart().Demo();
  disp( C.CP );

  sep( 1 , 'C = C.Add(''cafe'',3,2)   -- evento con la cache vacia: nada que hacer' );
  C = C.Add( 'cafe' , 3 , 2 );

  sep( 2 , 'leer las 4 CPs   -- MISS: computa y deposita en el store' );
  fprintf( '  total=%g  ticket=''%s''  topPrice=%g  currency=''%s''\n' , ...
           C.total , C.ticket , C.topPrice , C.currency );

  sep( 3 , 'C.total otra vez   -- HIT (no computa nada)' );
  fprintf( '  total=%g\n' , C.total );
  disp( C.CP );

  sep( 4 , 'D = C   -- la COPIA comparte el MISMO cacheHandle' );
  D = C;

  sep( 5 , 'C = C.Add(''pan'',1,1)   -- COW: handle NUEVO para C; D no se entera' );
  C = C.Add( 'pan' , 1 , 1 );
  disp( C.CP );

  sep( 6 , 'C.total -> REPLAY incremental (6+1=7).  D.total -> HIT del viejo (6)' );
  fprintf( '  C.total=%g\n' , C.total );
  fprintf( '  D.total=%g\n' , D.total );

  sep( 7 , 'C.ticket   -- cayo con el additem ([] = invalidar), asi que MISS' );
  fprintf( '  ticket=''%s''\n' , C.ticket );

  sep( 8 , 'C = C.SetPrice(''cafe'',4)   -- ticket SOBREVIVE (no declara el evento)' );
  C = C.SetPrice( 'cafe' , 4 );
  disp( C.CP );

  sep( 9 , 'C.total   -- el log no es todo-incremental: cascada a SYNC ABSOLUTO' );
  fprintf( '  total=%g\n' , C.total );

  sep( 10 , 'C.topPrice   -- 2 edits fundidos en UN marcador: un solo sync absoluto' );
  fprintf( '  topPrice=%g\n' , C.topPrice );

  sep( 11 , 'C.total_   -- recalculo FORZADO (descarta valor y log)' );
  fprintf( '  total=%g\n' , C.total_ );

  sep( 12 , 'C.CP.total.delete   -- borra el VALOR (la definicion queda)' );
  C.CP.total.delete
  disp( C.CP );

  sep( 13 , 'C = C.CP.total.set( 999 )   -- siembra a mano: NADIE lo verifica' );
  C = C.CP.total.set( 999 );
  fprintf( '  total=%g   <-- mentira consentida\n' , C.total );

  sep( 14 , 'h = C.CP.total.additem   -- el handler declarado, tal cual' );
  h = C.CP.total.additem;
  fprintf( '  %s\n' , func2str( h ) );

  sep( 15 , 'P = C.CP   -- el proxy en una VARIABLE: la cadena rota sigue valiendo' );
  P = C.CP;
  fprintf( '  class(P)=%s   P.ticket=''%s''\n' , class(P) , P.ticket );

  sep( 16 , 'guarda de ciclos: una CP que se pide a si misma' );
  C = C.Define( 'bucle' , @(c) c.bucle );
  try
    disp( C.bucle );
  catch ME
    fprintf( '  ERROR ATRAPADO [%s]\n  %s\n' , ME.identifier , ME.message );
  end

  fprintf( '\n' );
end

function sep( n , t )
  fprintf( '\n%s\n== PASO %d. %s\n' , repmat('-',1,76) , n , t );
end
