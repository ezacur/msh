function bench_vertexSeed
%BENCH_VERTEXSEED  ¿Acelera la busqueda sembrar una cota superior por punto?
%
%   Heuristica: el vertice mas cercano (BVH de PUNTOS sobre los vertices) da
%   una cota superior d_v de la distancia real; los triangulos de su 1-ring
%   (EsuP) la aprietan. La cota se siembra como Dmax POR PUNTO en
%   bvhClosestElement y poda la travesia desde la raiz.
%
%     H0: query normal (referencia)
%     H1: query de vertices  ->  siembra d_v
%     H2: H1 + distancia exacta al 1-ring  ->  siembra mas apretada
%
%   SEMANTICA (ojo, no es el Dmax de radio-de-busqueda): la cota es ALCANZABLE
%   (viene de un candidato real), asi que un miss (e=0/d=Inf/cp=NaN) bajo esta
%   siembra significa "el candidato ERA el ganador" y se RELLENA con el.
%   Garantia: si hubiera un elemento estrictamente mejor que la cota, la
%   travesia sembrada lo habria encontrado -> tras el relleno el resultado es
%   exacto. Ademas la cota se infla (1+1e-9) para que el propio elemento de la
%   cota se re-encuentre y los misses sean rarisimos. El bench VERIFICA
%   igualdad de d contra H0 en todos los casos.

  here = fileparts( mfilename('fullpath') );
  addpath( fullfile( here , '..' , 'MESH' ) );
  addpath( fullfile( here , '..' , 'tools' ) );
  rng(11);

  %malla: esfera triangulada y luego ABOLLADA (conectividad de la esfera,
  %geometria no convexa y de densidad no uniforme)
  nV = 26000;
  V0 = randn( nV ,3);  V0 = V0 ./ sqrt( sum( V0.^2 ,2) );
  T  = convhulln( V0 );
  V  = V0 .* ( 1 + 0.25*sin( 4*V0(:,1) ) .* cos( 3*V0(:,2) ) );
  M  = struct( 'xyz',V , 'tri',T );
  nF = size( T ,1);
  scale = norm( max(V,[],1) - min(V,[],1) );
  tiny  = 1e-12 * scale;

  fprintf( 'malla: %d vertices, %d triangulos (esfera abollada)\n' , nV , nF );

  %estructuras (one-time)
  tic;  B  = BVH( M );                                  tB  = toc;
  Mv = struct( 'xyz',V , 'tri',(1:nV).' );
  tic;  Bv = BVH( Mv );                                 tBv = toc;
  tic;  EsuP = meshEsuP( M ,'sparse');                  tE  = toc;
  [ tj0 , vj0 ] = find( EsuP );
  repTri = accumarray( vj0 , tj0 , [nV,1] , @min );     %un tri incidente por vertice
  fprintf( 'one-time: blob tri %.0f ms | blob vertices %.0f ms | EsuP %.0f ms\n\n' , ...
           1e3*tB , 1e3*tBv , 1e3*tE );

  %conjuntos de puntos
  nPn = 2e5;  nPm = 2e5;  nPf = 1e5;
  Pnear = surfPoints( V , T , nPn ) + 0.005*scale*randn( nPn ,3);
  Pmid  = ( rand( nPm ,3) - 0.5 ) .* ( 1.2*( max(V,[],1)-min(V,[],1) ) ) + ...
          ( max(V,[],1)+min(V,[],1) )/2;
  Pfar  = 10*scale*randn( nPf ,3);

  sets = { 'near-surface' , Pnear ; 'caja media' , Pmid ; 'far-field' , Pfar };

  fprintf( '%-13s %8s | %8s %8s %8s | %8s %8s | %9s %9s | %s\n' , ...
           'set' , 'H0' , 'q.vert' , 'ring' , 'seeded' , 'H1 tot' , 'H2 tot' , ...
           'H1/H0' , 'H2/H0' , 'exact%%/miss' );
  fprintf( '%s\n' , repmat( '-' , 1 , 110 ) );

  for s = 1:size( sets ,1)
    P  = sets{s,2};
    nP = size( P ,1);

    %---- H0: referencia
    t0 = Inf;
    for r = 1:3, tic; [ e0 , ~ , d0 ] = bvhClosestElement( {M,B} , P );  t0 = min( t0 , toc ); end

    %---- query de vertices (compartida por H1 y H2)
    tv = Inf;
    for r = 1:3, tic; [ ev , ~ , dv ] = bvhClosestElement( {Mv,Bv} , P );  tv = min( tv , toc ); end

    %---- H1: siembra con la distancia al vertice mas cercano
    seed1 = dv*(1+1e-9) + tiny;
    t1 = Inf;
    for r = 1:3, tic; [ e1 , ~ , d1 ] = bvhClosestElement( {M,B} , P , seed1 );  t1 = min( t1 , toc ); end
    m1 = e1 == 0;                                    %miss => el VERTICE era el ganador
    e1(m1) = repTri( ev(m1) );  d1(m1) = dv(m1);

    %---- H2: + cota del 1-ring (EsuP del vertice mas cercano)
    tr = Inf;
    for r = 1:3
      tic;
      cols = EsuP( : , ev );                         %nF x nP sparse
      [ ti , pi ] = find( cols );                    %pares (triangulo, punto)
      d2q = pTriD2( P(pi,:) , V(T(ti,1),:) , V(T(ti,2),:) , V(T(ti,3),:) );
      %min por punto SIN accumarray(@min): asignacion en orden descendente
      %(la ultima escritura por punto = su minimo)
      [ ds , o ] = sort( d2q , 'descend' );
      dub = inf( nP ,1);  ec = zeros( nP ,1);
      dub( pi(o) ) = ds;  ec( pi(o) ) = ti( o );
      dub = sqrt( dub );
      tr = min( tr , toc );
    end
    seed2 = dub*(1+1e-9) + tiny;
    t2 = Inf;
    for r = 1:3, tic; [ e2 , ~ , d2 ] = bvhClosestElement( {M,B} , P , seed2 );  t2 = min( t2 , toc ); end
    m2 = e2 == 0;                                    %miss => el candidato del ring gana
    e2(m2) = ec(m2);  d2(m2) = dub(m2);

    %---- EXACTITUD (contra H0, ambas rutas deben ser identicas)
    err1 = max( abs( d1 - d0 ) );
    err2 = max( abs( d2 - d0 ) );
    assert( err1 <= 1e-9*scale , 'H1 no exacta (err %g)' , err1 );
    assert( err2 <= 1e-9*scale , 'H2 no exacta (err %g)' , err2 );
    ringExact = mean( abs( dub - d0 ) <= 1e-12*scale + 1e-9*d0 );  %el 1-ring YA era el ganador

    us = @(t) 1e6*t/nP;
    fprintf( '%-13s %8.3f | %8.3f %8.3f %8.3f | %8.3f %8.3f | %8.2fx %8.2fx | %4.1f%% / %d+%d\n' , ...
             sets{s,1} , us(t0) , us(tv) , us(tr) , us(t2) , ...
             us(tv)+us(t1) , us(tv)+us(tr)+us(t2) , ...
             t0/(tv+t1) , t0/(tv+tr+t2) , 100*ringExact , nnz(m1) , nnz(m2) );
  end

  fprintf( '\n(µs/punto, best-of-3, single-thread; "seeded" = query sembrada de H2;\n' );
  fprintf( ' exact%% = puntos donde el 1-ring ya contenia al ganador; miss = rellenos H1+H2)\n' );
end

%---------------------------------------------------------------- helpers
function P = surfPoints( V , T , n )
%puntos uniformes-ish sobre la superficie: triangulo aleatorio + baricentricas
  r  = randi( size(T,1) , n ,1);
  w  = -log( rand( n ,3) );  w = w ./ sum( w ,2);
  P  = w(:,1).*V(T(r,1),:) + w(:,2).*V(T(r,2),:) + w(:,3).*V(T(r,3),:);
end

function d2 = pTriD2( Q , A , B , C )
%distancia^2 punto-triangulo, vectorizada (Ericson, todas las regiones)
  ab = B - A;  ac = C - A;  aq = Q - A;
  d1 = sum( ab.*aq ,2);  d2_ = sum( ac.*aq ,2);
  cp = A;  done = ( d1 <= 0 & d2_ <= 0 );                       %region A

  bq = Q - B;  d3 = sum( ab.*bq ,2);  d4 = sum( ac.*bq ,2);
  m = ~done & d3 >= 0 & d4 <= d3;                               %region B
  cp(m,:) = B(m,:);  done = done | m;

  vc = d1.*d4 - d3.*d2_;
  m = ~done & vc <= 0 & d1 >= 0 & d3 <= 0;                      %arista AB
  v = d1 ./ ( d1 - d3 );
  cp(m,:) = A(m,:) + v(m).*ab(m,:);  done = done | m;

  cq = Q - C;  d5 = sum( ab.*cq ,2);  d6 = sum( ac.*cq ,2);
  m = ~done & d6 >= 0 & d5 <= d6;                               %region C
  cp(m,:) = C(m,:);  done = done | m;

  vb = d5.*d2_ - d1.*d6;
  m = ~done & vb <= 0 & d2_ >= 0 & d6 <= 0;                     %arista AC
  w = d2_ ./ ( d2_ - d6 );
  cp(m,:) = A(m,:) + w(m).*ac(m,:);  done = done | m;

  va = d3.*d6 - d5.*d4;
  m = ~done & va <= 0 & ( d4 - d3 ) >= 0 & ( d5 - d6 ) >= 0;    %arista BC
  w = ( d4 - d3 ) ./ ( ( d4 - d3 ) + ( d5 - d6 ) );
  cp(m,:) = B(m,:) + w(m).*( C(m,:) - B(m,:) );  done = done | m;

  m = ~done;                                                    %interior
  den = 1 ./ ( va + vb + vc );
  v = vb .* den;  w = vc .* den;
  cp(m,:) = A(m,:) + v(m).*ab(m,:) + w(m).*ac(m,:);

  d2 = sum( ( Q - cp ).^2 ,2);
end
