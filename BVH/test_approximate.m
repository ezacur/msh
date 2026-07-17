function test_approximate
%TEST_APPROXIMATE  approximateClosestElement: garantias, semantica y firma.

  addpath( fullfile( fileparts( mfilename('fullpath') ) , '..' , 'MESH' ) );
  rng(19);

  %% nube de puntos (celltype 1): el aproximado ES exacto
  Vp = randn( 20000 ,3);
  Mp = struct( 'xyz',Vp , 'tri',(1:2e4).' );
  P  = randn( 5000 ,3)*2;
  [ ~ , ~ , d0 ] = bvhClosestElement( Mp , P );
  [ ea , cpa , da ] = approximateClosestElement( Mp , P );
  assert( max( abs( da - d0 ) ) < 1e-12 , 'nube: debe ser exacto' );
  assert( max( abs( sqrt(sum((cpa-Vp(ea,:)).^2,2)) ) ) < 1e-12 , 'nube: cp = el vertice' );
  fprintf( 'nube-puntos   ok  (aproximado == exacto)\n' );

  %% triangulos: cota superior + hit alto + salidas completas
  V = randn( 8000 ,3);  V = V ./ sqrt( sum( V.^2 ,2) );
  M = struct( 'xyz',V , 'tri',convhulln( V ) );
  P = [ elemPts(M,20000) + 0.01*randn(20000,3) ; randn(10000,3)*3 ];
  Ba = approximateClosestElement( M );                       %forma constructor
  [ e0 , ~ , d0 ] = bvhClosestElement( M , P );
  [ ea , cpa , da , bca , Fa ] = approximateClosestElement( {M,Ba} , P );
  assert( all( da >= d0 - 1e-12 ) , 'tri: cota superior violada' );
  hit = abs( da - d0 ) <= 1e-12 + 1e-9*d0;
  assert( mean( hit ) > 0.85 , 'tri: hit rate sospechosamente bajo (%.1f%%)' , 100*mean(hit) );
  assert( all( ea >= 1 ) , 'tri: sin misses sin Dmax' );
  %bc REGION-EXACTAS del MEX (igual que bvhClosestElement): suma 1, en [0,1],
  %y reconstruyen cp a precision de maquina
  assert( max( abs( sum( bca ,2) - 1 ) ) < 1e-13 , 'tri: bc no suman 1' );
  assert( min(bca(:)) >= 0 && max(bca(:)) <= 1 , 'tri: bc fuera de [0,1]' );
  recA = bca(:,1).*V(M.tri(ea,1),:) + bca(:,2).*V(M.tri(ea,2),:) + bca(:,3).*V(M.tri(ea,3),:);
  assert( max(max( abs( recA - cpa ) )) < 1e-12 , 'tri: bc no reconstruyen cp' );
  assert( isfield( Fa ,'type') && isfield( Fa ,'onBoundary') && all( Fa.type >= 1 ) , 'tri: F' );
  dcp = sqrt( sum( ( P - cpa ).^2 ,2) );
  assert( max( abs( dcp - da ) ) < 1e-9 , 'tri: |P-cp| debe igualar d' );
  fprintf( 'triangulos    ok  (cota + hit %.1f%% + bc region-exacta/F)\n' , 100*mean(hit) );

  %% Dmax escalar y vectorial (corta por la distancia AL VERTICE)
  [ eD , ~ , dD ] = approximateClosestElement( {M,Ba} , P , 0.05 );
  assert( all( isinf( dD( eD == 0 ) ) ) && any( eD == 0 ) , 'Dmax: semantica de miss' );
  assert( all( abs( dD(eD>0) - da(eD>0) ) < 1e-12 ) , 'Dmax: los encontrados no cambian' );
  %vectorial: OJO, el corte es por la distancia AL VERTICE (etapa 1), no al
  %elemento -- una cota de elemento (da) NO es alcanzable para la etapa 1.
  %Probamos la mecanica por punto con una mezcla Inf / ~0:
  seed = inf( size(P,1) ,1);  seed(2:2:end) = 1e-12;
  [ eV , ~ , dV ] = approximateClosestElement( {M,Ba} , P , seed );
  assert( all( eV(1:2:end) >= 1 ) && max( abs( dV(1:2:end) - da(1:2:end) ) ) < 1e-12 , ...
          'Dmax vectorial: la mitad Inf debe reproducir la query normal' );
  assert( all( eV(2:2:end) == 0 ) , 'Dmax vectorial: la mitad ~0 debe fallar' );
  fprintf( 'Dmax          ok  (escalar + vector por punto; corta por d_vertice)\n' );

  %% segmentos y tets: cota + funcionamiento
  s = linspace( 0 , 6*pi , 5001 ).';
  Ms = struct( 'xyz',[cos(s),sin(s),0.1*s] , 'tri',[(1:5e3).',(2:5e3+1).'] );
  P2 = randn( 8000 ,3)*2;
  [ ~ , ~ , d0 ] = bvhClosestElement( Ms , P2 );
  [ ~ , ~ , da ] = approximateClosestElement( Ms , P2 );
  assert( all( da >= d0 - 1e-12 ) , 'seg: cota violada' );

  Xt = randn( 3000 ,3);
  Mt = struct( 'xyz',Xt , 'tri',delaunayn( Xt ) );
  P3 = [ randn( 5000 ,3)*0.8 ; randn( 3000 ,3)*4 ];          %interior + exterior
  [ ~ , ~ , d0 ] = bvhClosestElement( Mt , P3 );
  [ ~ , ~ , da ] = approximateClosestElement( Mt , P3 );
  assert( all( da >= d0 - 1e-12 ) , 'tet: cota violada' );
  fprintf( 'seg/tets      ok  (cota superior en ambos)\n' );

  %% mixta 0-padded
  Mm = struct( 'xyz',[ M.xyz ; Ms.xyz + 5 ] , ...
               'tri',[ M.tri ; [ Ms.tri + size(M.xyz,1) , zeros(5e3,1) ] ] );
  [ ~ , ~ , d0 ] = bvhClosestElement( Mm , P2 );
  [ ~ , ~ , da ] = approximateClosestElement( Mm , P2 );
  assert( all( da >= d0 - 1e-12 ) , 'mixta: cota violada' );
  fprintf( 'mixta         ok\n' );

  %% blob rancio -> ERROR; blob ajeno -> ERROR
  M2 = M;  M2.xyz = M.xyz * 2;
  try
    approximateClosestElement( {M2,Ba} , P(1:10,:) );
    error('test:stale','debia errar');
  catch ME
    assert( strcmp( ME.identifier , 'approximateClosestElement:staleBVH' ) , 'stale: %s' , ME.identifier );
  end
  try
    approximateClosestElement( M , P(1:10,:) , BVH( M ) );   %blob del motor exacto
    error('test:foreign','debia errar');
  catch ME
    assert( strcmp( ME.identifier , 'approximateClosestElement:B' ) , 'foreign: %s' , ME.identifier );
  end
  fprintf( 'guardias      ok  (stale + blob ajeno)\n' );

  fprintf( 'ALL approximateClosestElement tests passed.\n' );
end

function P = elemPts( M , n )
  T = M.tri;  r = randi( size(T,1) , n ,1);
  w = -log( rand( n ,3) );  w = w ./ sum( w ,2);
  P = w(:,1).*M.xyz(T(r,1),:) + w(:,2).*M.xyz(T(r,2),:) + w(:,3).*M.xyz(T(r,3),:);
end
