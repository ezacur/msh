function test_fanClosestElement
%TEST_FANCLOSESTELEMENT  fanClosestElement: garantias, formas de llamada y guards.

  addpath( fullfile( fileparts( mfilename('fullpath') ) , '..' , 'MESH' ) );
  rng(19);

  %% nube de puntos (celltype 1): el fan ES exacto
  Vp = randn( 20000 ,3);
  Mp = struct( 'xyz',Vp , 'tri',(1:2e4).' );
  P  = randn( 5000 ,3)*2;
  [ ~ , ~ , d0 ] = bvhClosestElement( Mp , P );
  [ ea , cpa , da ] = fanClosestElement( Mp , P );
  assert( max( abs( da - d0 ) ) < 1e-12 , 'nube: debe ser exacto' );
  assert( max( abs( sqrt(sum((cpa-Vp(ea,:)).^2,2)) ) ) < 1e-12 , 'nube: cp = el vertice' );
  fprintf( 'nube-puntos   ok  (fan == exacto)\n' );

  %% triangulos: cota superior + hit alto + salidas completas
  V = randn( 8000 ,3);  V = V ./ sqrt( sum( V.^2 ,2) );
  M = struct( 'xyz',V , 'tri',convhulln( V ) );
  P = [ elemPts(M,20000) + 0.01*randn(20000,3) ; randn(10000,3)*3 ];
  [ e0 , ~ , d0 ] = bvhClosestElement( M , P );
  [ ea , cpa , da , bca , Fa ] = fanClosestElement( M , P );
  assert( all( da >= d0 - 1e-12 ) , 'tri: cota superior violada' );
  hit = abs( da - d0 ) <= 1e-12 + 1e-9*d0;
  assert( mean( hit ) > 0.85 , 'tri: hit rate sospechosamente bajo (%.1f%%)' , 100*mean(hit) );
  assert( all( ea >= 1 ) , 'tri: sin misses sin Dmax' );
  %bc REGION-EXACTAS del MEX: suman 1, en [0,1], reconstruyen cp
  assert( max( abs( sum( bca ,2) - 1 ) ) < 1e-13 , 'tri: bc no suman 1' );
  assert( min(bca(:)) >= 0 && max(bca(:)) <= 1 , 'tri: bc fuera de [0,1]' );
  recA = bca(:,1).*V(M.tri(ea,1),:) + bca(:,2).*V(M.tri(ea,2),:) + bca(:,3).*V(M.tri(ea,3),:);
  assert( max(max( abs( recA - cpa ) )) < 1e-12 , 'tri: bc no reconstruyen cp' );
  assert( isfield( Fa ,'type') && isfield( Fa ,'onBoundary') && all( Fa.type >= 1 ) , 'tri: F' );
  dcp = sqrt( sum( ( P - cpa ).^2 ,2) );
  assert( max( abs( dcp - da ) ) < 1e-9 , 'tri: |P-cp| debe igualar d' );
  fprintf( 'triangulos    ok  (cota + hit %.1f%% + bc region-exacta/F)\n' , 100*mean(hit) );

  %% las 4 formas: nodos dados (vector/escalar) y Bnodes reproducen la etapa 1
  Mn = struct( 'xyz',M.xyz , 'tri',(1:size(M.xyz,1)).' );
  nodes = bvhClosestElement( Mn , P );
  [ e1 , ~ , d1 ] = fanClosestElement( {M,nodes} , P );
  assert( isequal( e1 , ea ) && isequal( d1 , da ) , 'seed vector: no reproduce' );
  Bn = BVH( Mn );
  [ e2 , ~ , d2 ] = fanClosestElement( {M,Bn} , P );
  assert( isequal( e2 , ea ) && isequal( d2 , da ) , 'Bnodes: no reproduce' );
  %seed escalar: el fan de UN nodo -- cota valida y cp dentro de su abanico
  [ e3 , ~ , d3 ] = fanClosestElement( {M,nodes(1)} , P(1:100,:) );
  assert( all( d3 >= d0(1:100) - 1e-12 ) , 'seed escalar: viola la cota' );
  fan1 = find( any( M.tri == nodes(1) ,2) );
  assert( all( ismember( e3 , fan1 ) ) , 'seed escalar: gano un elemento fuera del abanico' );
  fprintf( 'formas        ok  (vector/Bnodes reproducen; escalar barre su abanico)\n' );

  %% Dmax escalar y vectorial (corta por la distancia AL NODO en etapa 1)
  [ eD , ~ , dD ] = fanClosestElement( M , P , 0.05 );
  assert( all( isinf( dD( eD == 0 ) ) ) && any( eD == 0 ) , 'Dmax: semantica de miss' );
  assert( all( abs( dD(eD>0) - da(eD>0) ) < 1e-12 ) , 'Dmax: los encontrados no cambian' );
  seed = inf( size(P,1) ,1);  seed(2:2:end) = 1e-12;
  [ eV , ~ , dV ] = fanClosestElement( M , P , seed );
  assert( all( eV(1:2:end) >= 1 ) && max( abs( dV(1:2:end) - da(1:2:end) ) ) < 1e-12 , ...
          'Dmax vectorial: la mitad Inf debe reproducir la query normal' );
  assert( all( eV(2:2:end) == 0 ) , 'Dmax vectorial: la mitad ~0 debe fallar' );
  %con nodos DADOS, Dmax filtra por la distancia al ELEMENTO
  [ eG , ~ , dG ] = fanClosestElement( {M,nodes} , P , 0.05 );
  assert( all( eG( da < 0.05 ) > 0 ) && all( eG( da >= 0.05 ) == 0 ) , ...
          'Dmax con nodos dados: debe filtrar por d del elemento' );
  assert( all( isinf( dG( eG == 0 ) ) ) , 'Dmax con nodos dados: miss -> Inf' );
  fprintf( 'Dmax          ok  (etapa 1 por d_nodo; con nodos dados por d_elemento)\n' );

  %% segmentos y tets: cota + funcionamiento
  s = linspace( 0 , 6*pi , 5001 ).';
  Ms = struct( 'xyz',[cos(s),sin(s),0.1*s] , 'tri',[(1:5e3).',(2:5e3+1).'] );
  P2 = randn( 8000 ,3)*2;
  [ ~ , ~ , d0 ] = bvhClosestElement( Ms , P2 );
  [ ~ , ~ , da ] = fanClosestElement( Ms , P2 );
  assert( all( da >= d0 - 1e-12 ) , 'seg: cota violada' );

  Xt = randn( 3000 ,3);
  Mt = struct( 'xyz',Xt , 'tri',delaunayn( Xt ) );
  P3 = [ randn( 5000 ,3)*0.8 ; randn( 3000 ,3)*4 ];          %interior + exterior
  [ ~ , ~ , d0 ] = bvhClosestElement( Mt , P3 );
  [ ~ , ~ , da ] = fanClosestElement( Mt , P3 );
  assert( all( da >= d0 - 1e-12 ) , 'tet: cota violada' );
  fprintf( 'seg/tets      ok  (cota superior en ambos)\n' );

  %% mixta 0-padded + vertice suelto como semilla -> miss
  Mm = struct( 'xyz',[ M.xyz ; Ms.xyz + 5 ] , ...
               'tri',[ M.tri ; [ Ms.tri + size(M.xyz,1) , zeros(5e3,1) ] ] );
  [ ~ , ~ , d0 ] = bvhClosestElement( Mm , P2 );
  [ ~ , ~ , da ] = fanClosestElement( Mm , P2 );
  assert( all( da >= d0 - 1e-12 ) , 'mixta: cota violada' );
  Miso = struct( 'xyz',[ M.xyz ; 99 99 99 ] , 'tri',M.tri );  %nodo sin elementos
  ni = size( Miso.xyz ,1);
  [ ei , ~ , di ] = fanClosestElement( {Miso,ni} , P(1:5,:) );
  assert( all( ei == 0 ) && all( isinf( di ) ) , 'vertice suelto: debe dar miss' );
  fprintf( 'mixta/suelto  ok\n' );

  %% guards: Bnodes ajeno/rancio, nodos malos
  try
    fanClosestElement( {M,BVH(M)} , P(1:10,:) );             %blob de ELEMENTOS
    error('test:foreign','debia errar');
  catch ME
    assert( strcmp( ME.identifier , 'fanClosestElement:B' ) , 'foreign: %s' , ME.identifier );
  end
  M2 = M;  M2.xyz = M.xyz * 2;                               %malla editada, Bnodes viejo
  try
    fanClosestElement( {M2,Bn} , P(1:10,:) );
    error('test:stale','debia errar');
  catch ME
    assert( strcmp( ME.identifier , 'bvhClosestElement:staleBVH' ) , 'stale: %s' , ME.identifier );
  end
  try
    fanClosestElement( {M,[1;2]} , P(1:10,:) );              %vector de nodos corto
    error('test:count','debia errar');
  catch ME
    assert( strcmp( ME.identifier , 'fanClosestElement:seed' ) , 'count: %s' , ME.identifier );
  end
  try
    fanClosestElement( {M,1.5} , P(1:10,:) );                %nodo no entero
    error('test:int','debia errar');
  catch ME
    assert( strcmp( ME.identifier , 'fanClosestElement:seed' ) , 'int: %s' , ME.identifier );
  end
  fprintf( 'guards        ok  (blob ajeno + Bnodes rancio + nodos malos)\n' );

  fprintf( 'ALL fanClosestElement tests passed.\n' );
end

function P = elemPts( M , n )
  T = M.tri;  r = randi( size(T,1) , n ,1);
  w = -log( rand( n ,3) );  w = w ./ sum( w ,2);
  P = w(:,1).*M.xyz(T(r,1),:) + w(:,2).*M.xyz(T(r,2),:) + w(:,3).*M.xyz(T(r,3),:);
end
