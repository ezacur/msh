%routeDemo  Sesion paso a paso de route: el patron de cache como puro dominio.
%
%   Ver route_TUTORIAL.html en la raiz del repo para la explicacion linea a
%   linea del fichero y de cada paso de esta sesion.
%
%   Uso:  addpath('C:\repos\msh','C:\repos\msh\dev'); routeDemo
%
% See also route, cachedOwner, cartDemo.

function routeDemo()

  sep( 0 , 'R = route().Demo()   -- registra 3 definiciones; la cache queda VACIA' );
  R = route().Demo();
  R.DEBUG = true;
  disp( R.CP );

  sep( 1 , 'R = R.AddPoint(0,0); R = R.AddPoint(3,4)   -- eventos con la cache vacia' );
  R = R.AddPoint( 0 , 0 );
  R = R.AddPoint( 3 , 4 );

  sep( 2 , 'R.len y R.npts   -- MISS: computa y deposita. len = |(3,4)| = 5' );
  fprintf( '  len=%g   npts=%d\n' , R.len , R.npts );

  sep( 3 , 'R2 = R   -- la COPIA comparte el MISMO cacheHandle' );
  R2 = R;

  sep( 4 , 'R = R.AddPoint(3,8)   -- COW: len y npts quedan PENDIENTES con el edit (d=4)' );
  R = R.AddPoint( 3 , 8 );
  disp( R.CP );

  sep( 5 , 'R.len -> REPLAY incremental (5+4=9).  R2.len -> HIT del viejo (5)' );
  fprintf( '  R.len=%g\n' , R.len );
  fprintf( '  R2.len=%g\n' , R2.len );

  sep( 6 , 'R = R.Scale(2)   -- npts SOBREVIVE al evento... CON su log pendiente intacto' );
  R = R.Scale( 2 );
  disp( R.CP );

  sep( 7 , 'R.npts -> el edit del paso 4 se resuelve AHORA, incremental (2+1=3)' );
  fprintf( '  npts=%d\n' , R.npts );

  sep( 8 , 'R.len -> el log no es incremental (Scale no llevo args): sync absoluto (18)' );
  fprintf( '  len=%g\n' , R.len );

  sep( 9 , 'R.bbox -> nadie lo habia pedido: MISS recien ahora (pereza)' );
  fprintf( '  bbox=[%g %g %g %g]\n' , R.bbox );

  sep( 10 , 'R.len_ -> recalculo FORZADO: verifica que los handlers no mintieron' );
  fprintf( '  len=%g\n' , R.len_ );

  fprintf( '\n' );
end

function sep( n , t )
  fprintf( '\n%s\n== PASO %d. %s\n' , repmat('-',1,76) , n , t );
end
