%ivpDemo  Sesion paso a paso del ejemplo que PAGA el precio del patron (ivp).
%
%   dy/dt = cos(t)*y con y(0)=1: la exacta es y = exp(sin(t)), asi que cada
%   paso puede ensenar el error de verdad. El coste de integrar es
%   proporcional al tramo, y las trazas [ivp] ODE dicen cuantos pasos de RK4
%   costo cada peticion: ahi se ve el ahorro (At(6) tras At(4) integra 2000
%   pasos, no 6000).
%
%   OJO SI LO EDITAS: la salida es un GOLDEN -- test_cache la compara byte a
%   byte contra ivpDemo_expected.txt. test_cache('record') regenera el golden,
%   SOLO tras un cambio deliberado.
%
%   Uso:  addpath('C:\repos\msh','C:\repos\msh\dev'); ivpDemo
%
% See also ivp, factoDemo, cartDemo, cachedOwner.

function ivpDemo()

  exact = @(t,y0) y0 * exp( sin( t ) );    % la EDO es lineal: y = Y0*exp(sin t)

  sep( 0 , 'P = ivp( @(t,y) cos(t).*y , 1 )   -- registra `sol`; cache VACIA' );
  P = ivp( @(t,y) cos(t).*y , 1 );  P.DEBUG = true;   % las trazas las pide el demo
  disp( P.CP );

  sep( 1 , '[P,y] = P.At(4)   -- MISS: integra 0 -> 4 (4000 pasos)' );
  [ P , y ] = P.At( 4 );
  fprintf( '  y(4) = %.6f   (clava la exacta a 1e-9: %d)\n' , y , abs( y - exact(4,1) ) < 1e-9 );

  sep( 2 , '[P,y] = P.At(6)   -- sync que REUSA: integra SOLO 4 -> 6 (2000 pasos)' );
  [ P , y ] = P.At( 6 );
  fprintf( '  y(6) = %.6f   (clava la exacta a 1e-9: %d)\n' , y , abs( y - exact(6,1) ) < 1e-9 );

  sep( 3 , '[P,y] = P.At(5)   -- 4 -> 5 (1000 pasos), y se inserta EN SU SITIO' );
  [ P , y ] = P.At( 5 );
  fprintf( '  y(5) = %.6f   (clava la exacta a 1e-9: %d)\n' , y , abs( y - exact(5,1) ) < 1e-9 );
  fprintf( '  tiempos guardados: %s\n' , mat2str( P.sol(:,1).' ) );

  sep( 4 , 'P.At(5) otra vez   -- HIT: cero integracion' );
  [ P , y ] = P.At( 5 );
  fprintf( '  y(5) = %.6f\n' , y );

  sep( 5 , 'D = P; P = P.SetY0(2)   -- cambia el problema: sol CAE para P, no para D' );
  D = P;
  P = P.SetY0( 2 );
  disp( P.CP );

  sep( 6 , 'P.At(4) recomputa TODO con Y0=2; D conserva su solucion intacta (COW)' );
  [ P , y ] = P.At( 4 );
  fprintf( '  P: y(4) = %.6f   (clava la exacta a 1e-9: %d)\n' , y , abs( y - exact(4,2) ) < 1e-9 );
  [ D , y ] = D.At( 4 );
  fprintf( '  D: y(4) = %.6f   (HIT del valor viejo, con Y0=1)\n' , y );

  fprintf( '\n' );
end

function sep( n , t )
  fprintf( '\n%s\n== PASO %d. %s\n' , repmat('-',1,76) , n , t );
end
