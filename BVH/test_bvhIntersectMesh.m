function test_bvhIntersectMesh
%TEST_BVHINTERSECTMESH  Curva de interseccion: topologia, exactitud y la regla de
%   paridad en las zonas coplanares.
%
%   El invariante que se comprueba en todos los casos es TOPOLOGICO, no metrico:
%   una curva de interseccion valida tiene todos los nodos de grado 2 (lazos
%   cerrados) o grado 2 salvo los extremos (curvas abiertas sobre superficies con
%   borde), y NUNCA nodos de grado > 2 -- un grado 3 significa que la curva se
%   bifurca, que es como se manifiestan los fallos de soldadura.
  root = fileparts( fileparts( mfilename('fullpath') ) );
  addpath( root , fullfile(root,'MESH') , fullfile(root,'tools') , fullfile(root,'BVH') );

  %% 1) dos esferas: la interseccion es un CIRCULO conocido analiticamente
  S1 = sphereMesh(3);
  S2 = sphereMesh(3);  S2.xyz = S2.xyz + [0.8 0 0];
  [ C , K ] = bvhIntersectMesh( S1 , S2 );
  chk( C , '1) dos esferas' , 0 , 1 );          %0 extremos, 1 componente
  x0 = 0.4;  r0 = sqrt(1-x0^2);
  assert( max(abs( C.xyz(:,1) - x0 )) < 0.01 , '1) fuera del plano del circulo' );
  assert( max(abs( sqrt(sum(C.xyz(:,2:3).^2,2)) - r0 )) < 0.01 , '1) radio' );
  assert( size(K,1) == size(C.xyz,1) , '1) K debe tener una fila por nodo' );
  assert( all( ismember( K(:,1) , [1 2 3] ) ) , '1) kind invalido' );

  %% 2) el orden de los argumentos no cambia la curva
  C2 = bvhIntersectMesh( S2 , S1 );
  assert( isequal( sortrows(round(C.xyz*1e10)) , sortrows(round(C2.xyz*1e10)) ) , ...
          '2) A x B y B x A deben dar la misma curva' );

  %% 3) mallas disjuntas -> curva vacia
  S3 = sphereMesh(2);  S3.xyz = S3.xyz + [10 0 0];
  C3 = bvhIntersectMesh( S1 , S3 );
  assert( isempty( C3.tri ) , '3) mallas disjuntas no deben cortar' );

  %% 4) superficie ABIERTA: la curva tiene 2 extremos
  [uu,vv] = meshgrid( linspace(-1.4,0.3,14) , linspace(-1.4,1.4,14) );
  Pt = struct('xyz',[ uu(:) , zeros(numel(uu),1) , vv(:) ] , 'tri',delaunay(uu(:),vv(:)) );
  C4 = bvhIntersectMesh( S1 , Pt );
  chk( C4 , '4) esfera x parche abierto' , 2 , 1 );

  %% 5) DOS lazos disjuntos: esfera atravesada por un cilindro
  Cy = cylMesh( 0.5 , 2 , 24 , 9 );
  C5 = bvhIntersectMesh( S1 , Cy );
  chk( C5 , '5) esfera x cilindro' , 0 , 2 );

  %% 6) resoluciones distintas
  S6 = sphereMesh(1);  S6.xyz = S6.xyz*0.9 + [0.5 0 0];
  C6 = bvhIntersectMesh( S1 , S6 );
  chk( C6 , '6) resoluciones distintas' , 0 , 1 );

  %% 7) LA REGLA DE PARIDAD: coplanar + cruce en las mismas mallas.
  %%    Sin cancelar por paridad, la zona coplanar aporta toda su red interna de
  %%    aristas (medido: >1000 segmentos de ruido y cientos de nodos de grado >2).
  [xg,yg] = meshgrid( linspace(0,2,17) , linspace(0,1,9) );
  Tg = delaunay( xg(:) , yg(:) );
  A7 = struct('xyz',[ xg(:) , yg(:) , zeros(numel(xg),1) ] , 'tri',Tg );
  for shift = [ 0 , 0.031 ]                     %coincidente y desplazado
    x2 = xg(:) + shift;  y2 = yg(:) + shift/2;
    z2 = zeros(numel(x2),1);  w = x2 > 1;
    z2(w) = 0.4*sin( 2*pi*( x2(w)-1 ) );
    B7 = struct('xyz',[ x2 , y2 , z2 ] , 'tri',Tg );
    C7 = bvhIntersectMesh( A7 , B7 );
    d  = degrees( C7 );
    assert( ~isempty( C7.tri ) , '7) shift %g: deberia haber curva' , shift );
    assert( sum( d > 2 ) == 0 , ...
            '7) shift %g: %d nodo(s) de grado >2: la paridad no cancelo el interior' , ...
            shift , sum(d>2) );
  end

  %% 8) las claves son CANONICAS: ningun nodo duplicado, ni por clave ni por sitio
  assert( size(unique(K,'rows'),1) == size(K,1) , '8) claves duplicadas' );
  assert( size(unique(C.xyz,'rows'),1) == size(C.xyz,1) , '8) nodos coincidentes' );

  %% 9) EL RECORRIDO contra la FUERZA BRUTA.
  %%    OJO con la invariante correcta: los dos son filtros CONSERVADORES pero
  %%    NINGUNO contiene al otro. La fuerza bruta compara AABBs en MUNDO; el
  %%    recorrido compara la AABB del triangulo de A en el BUILD FRAME de B contra
  %%    las cajas de nodo, y la AABB de un triangulo no es invariante bajo la
  %%    rotacion del frame. Lo que si tiene que cumplirse: todo par que se corte
  %%    DE VERDAD esta en los dos, y la CURVA resultante es identica.
  if exist( 'bvhTriPairs_mx' ,'file') == 3
    for cs = { {S1,S2} , {S1,Pt} , {S1,Cy} , {S1,S6} }
      A9 = cs{1}{1};  B9 = cs{1}{2};
      Bb  = BVH( struct('xyz',B9.xyz,'tri',B9.tri) );
      Pt9 = bvhTriPairs_mx( intoFrame( A9.xyz , Bb ) , A9.tri , Bb );
      Pb9 = bruteRef( A9 , B9 );
      %los pares de la fuerza bruta que SE CORTAN tienen que estar en el recorrido
      [ ~,~,~,~, Cb ] = triTriPairs_mx( A9.xyz,A9.tri, B9.xyz,B9.tri, Pb9 );
      real9 = Pb9( Cb > 0 ,:);
      assert( all( ismember( real9 , Pt9 , 'rows' ) ) , ...
              '9) el recorrido pierde %d par(es) que SI se cortan' , ...
              sum( ~ismember( real9 , Pt9 , 'rows' ) ) );
      %y la curva tiene que salir identica por los dos caminos
      [ Ka , Xa , Sa ] = triTriPairs_mx( A9.xyz,A9.tri, B9.xyz,B9.tri, Pt9 );  %#ok<ASGLU>
      [ Kb , Xb , Sb ] = triTriPairs_mx( A9.xyz,A9.tri, B9.xyz,B9.tri, Pb9 );  %#ok<ASGLU>
      assert( isequal( sortrows(unique(Ka,'rows')) , sortrows(unique(Kb,'rows')) ) , ...
              '9) recorrido y fuerza bruta dan claves distintas' );
      assert( size(Sa,1) == size(Sb,1) , '9) numero de segmentos distinto' );
    end
    fprintf( '  recorrido validado contra fuerza bruta (misma curva)\n' );
  else
    fprintf( '  (bvhTriPairs_mx sin compilar: recorrido no validado)\n' );
  end

  %% 10) AUTOINTERSECCION. Una esfera cerrada NO se autocorta: curva vacia. Sin
  %%     excluir el par identico y la adyacencia devolvia su red de aristas entera
  %%     (134 de 162 nodos con grado >2).
  C10 = bvhIntersectMesh( S1 , S1 );
  assert( isempty( C10.tri ) , ...
          '10) una esfera cerrada no se autocorta: %d segmento(s) de mas' , size(C10.tri,1) );

  %%     y una malla que SI se autocorta: las dos esferas de arriba unidas en UNA.
  %%     Su autointerseccion tiene que ser exactamente el mismo circulo del caso 1.
  U = struct( 'xyz' , [ S1.xyz ; S2.xyz ] , ...
              'tri' , [ S1.tri ; S2.tri + size(S1.xyz,1) ] );
  C11 = bvhIntersectMesh( U , U );
  chk( C11 , '10) autointerseccion de la union' , 0 , 1 );
  assert( size(C11.tri,1) == size(C.tri,1) , ...
          '10) la autointerseccion deberia dar el mismo circulo: %d vs %d segmentos' , ...
          size(C11.tri,1) , size(C.tri,1) );
  assert( isequal( sortrows(round(C11.xyz*1e10)) , sortrows(round(C.xyz*1e10)) ) , ...
          '10) los nodos de la autointerseccion no coinciden con los del caso 1' );

  %%     B deformada NO cuenta como autointerseccion aunque comparta numeracion
  Sd = S1;  Sd.xyz = S1.xyz .* [1 1 0.6];
  C12 = bvhIntersectMesh( S1 , Sd );
  assert( ~isempty( C12.tri ) , ...
          '10) una malla deformada comparte ids pero es OTRA superficie: debe cortar' );

  %% 11) el bundle {M,B} REUTILIZA el blob (antes lo aceptaba y lo tiraba)
  if exist( 'bvhTriPairs_mx' ,'file') == 3
    Bb11 = BVH( struct('xyz',S2.xyz,'tri',S2.tri) );
    C13  = bvhIntersectMesh( S1 , {S2,Bb11} );
    assert( isequal( C13.xyz , C.xyz ) && isequal( C13.tri , C.tri ) , ...
            '11) con blob y sin blob deben dar la misma curva' );
    %un blob que NO sirve (otra malla) debe ignorarse y reconstruirse, no mentir
    Bbad = BVH( struct('xyz',S6.xyz,'tri',S6.tri) );
    C14  = bvhIntersectMesh( S1 , {S2,Bbad} );
    assert( isequal( C14.xyz , C.xyz ) && isequal( C14.tri , C.tri ) , ...
            '11) un blob ajeno debe descartarse, no usarse' );
  end

  fprintf( 'ALL bvhIntersectMesh tests passed.\n' );
end

function P = bruteRef( A , B )
%referencia: todos los pares cuyas AABB de cara se solapan
  a = boxesOf(A);  b = boxesOf(B);  I = cell(size(a,1),1);
  for i = 1:size(a,1)
    w = a(i,1)<=b(:,4) & a(i,4)>=b(:,1) & a(i,2)<=b(:,5) & a(i,5)>=b(:,2) & ...
        a(i,3)<=b(:,6) & a(i,6)>=b(:,3);
    if any(w), I{i} = [ repmat(i,sum(w),1) , find(w) ]; end
  end
  P = vertcat( I{:} );  if isempty(P), P = zeros(0,2); end
end
function Bx = boxesOf( M )
  X = double(M.xyz); T = double(M.tri); Bx = zeros(size(T,1),6);
  for c = 1:3
    V = [ X(T(:,1),c) , X(T(:,2),c) , X(T(:,3),c) ];
    Bx(:,c) = min(V,[],2);  Bx(:,c+3) = max(V,[],2);
  end
end
function X = intoFrame( X , B )
%al BUILD frame del blob (el mex espera las coordenadas alli)
  X = double(X);  X(:,end+1:3) = 0;
  if isequal( B.frame , eye(4) ), return; end
  Af = B.frame(1:3,1:3);  s2 = trace(Af.'*Af)/3;
  X  = ( X - B.frame(1:3,4).' ) * ( Af / s2 );
end

function d = degrees( C )
  if isempty(C.tri), d = zeros(size(C.xyz,1),1); return; end
  d = accumarray( C.tri(:) ,1,[size(C.xyz,1) 1]);
end

function chk( C , tag , nEnds , nComp )
  assert( ~isempty( C.tri ) , '%s: curva vacia' , tag );
  assert( all( C.tri(:) >= 1 & C.tri(:) <= size(C.xyz,1) ) , '%s: indice fuera de rango' , tag );
  assert( ~any( C.tri(:,1) == C.tri(:,2) ) , '%s: segmento degenerado' , tag );
  d = degrees( C );
  assert( ~any( d == 0 ) , '%s: nodo huerfano' , tag );
  assert( sum( d > 2 ) == 0 , '%s: %d nodo(s) de grado >2 (la curva se bifurca)' , tag , sum(d>2) );
  assert( sum( d == 1 ) == nEnds , '%s: %d extremos, esperaba %d' , tag , sum(d==1) , nEnds );
  nc = ncomp( C );
  assert( nc == nComp , '%s: %d componentes, esperaba %d' , tag , nc , nComp );
end

function n = ncomp( C )
  nV = size( C.xyz ,1);
  G = sparse( [C.tri(:,1);C.tri(:,2)] , [C.tri(:,2);C.tri(:,1)] , 1 , nV , nV );
  [ ~ , ~ , r ] = dmperm( G + speye(nV) );
  n = numel(r) - 1;
end

function M = cylMesh( r , h , nt , nz )
  th = linspace(0,2*pi,nt+1).';  th(end) = [];
  zz = linspace(-h,h,nz);
  [TH,ZZ] = meshgrid(th,zz);
  V = [ r*cos(TH(:)) , r*sin(TH(:)) , ZZ(:) ];
  T = [];
  for i = 1:nz-1
    for j = 1:nt
      j2 = mod(j,nt)+1;
      a = (j-1)*nz+i;  b = (j2-1)*nz+i;
      T = [ T ; a b b+1 ; a b+1 a+1 ];                                     %#ok<AGROW>
    end
  end
  M = struct('xyz',V,'tri',T);
end
