function serieDemo()
%serieDemo  La sesion del tutorial de diseño (serie_TUTORIAL.html), paso a paso.
%
%   Cada paso imprime lo que hace y lo que el motor narra (DEBUG on). Compara
%   los DOS registros de la misma clase -- Safe (conservador) y Fast (con
%   atajos) -- sobre las mismas ediciones.
%
%   Uso:  addpath('C:\repos\msh','C:\repos\msh\dev');  serieDemo
%
%   OJO: la salida de este demo es un GOLDEN de test_cache (byte a byte).
%   Cualquier cambio aqui, en serie.m o en lo que imprime el motor exige
%   regrabar con test_cache('record') -- solo tras un cambio deliberado.
%
% See also serie, cachedOwner, cartDemo, routeDemo.

  sep( 0 , 'S = serie().Safe()   -- registro CONSERVADOR: todo se invalida' );
  S = serie().Safe();  S.DEBUG = true;
  disp( S.CP );

  sep( 1 , 'S = S.Push(4); S = S.Push(8)   -- editar NO calcula nada' );
  S = S.Push(4);  S = S.Push(8);

  sep( 2 , 'S.avg   -- primera lectura: MISS, se calcula y se guarda' );
  fprintf( '  avg = %g\n' , S.avg );

  sep( 3 , 'S.avg   -- segunda lectura: HIT, no se recalcula' );
  fprintf( '  avg = %g\n' , S.avg );

  sep( 4 , 'S = S.Push(6)   -- una edicion tira TODO lo cacheado' );
  S = S.Push(6);
  disp( S.CP );

  sep( 5 , 'F = serie().Fast()   -- mismos derivados, politicas distintas' );
  F = serie().Fast();  F.DEBUG = true;
  F = F.Push(4);  F = F.Push(8);
  fprintf( '  avg = %g , top = %g\n' , F.avg , F.top );
  disp( F.CP );

  sep( 6 , 'F = F.Push(6)   -- avg y top quedan PENDIENTES; sorted cae' );
  F = F.Push(6);
  disp( F.CP );

  sep( 7 , 'F.avg   -- replay INCREMENTAL: se pone al dia sin recorrer la lista' );
  fprintf( '  avg = %g   (verdad: %g)\n' , F.avg , mean([4 8 6]) );

  sep( 8 , 'F = F.Pop()   -- al quitar: avg se SINCRONIZA, top se INVALIDA' );
  F = F.Pop();
  disp( F.CP );
  fprintf( '  avg = %g , top = %g\n' , F.avg , F.top );

  sep( 9 , 'VerifyCPs( F )   -- comparar lo que sirve la cache con la verdad' );
  VerifyCPs( F );

  sep( 10 , 'un handler MAL escrito, y como se caza' );
  W = serie().Fast();  W.DEBUG = false;
  W = W.Define( 'malo' , @(s) sum(s.X)/max(numel(s.X),1) , ...
                'push' , @(m,~,v,n) (m+v)/2 );        % NO es la media incremental
  W = W.Push(2);  W = W.Push(4);  W = W.Push(6);
  m = W.malo;                                                             %#ok<NASGU>
  W = W.Push(8);  W = W.Push(10);
  fprintf( '  la cache sirve %g, la verdad es %g -- y no falla ni avisa nadie\n' , ...
           W.malo , mean([2 4 6 8 10]) );
  VerifyCPs( W );
end

function sep( n , t )
  fprintf( '\n%s\n%2d) %s\n' , repmat('=',1,74) , n , t );
end
