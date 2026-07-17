function bench_approximate
%BENCH_APPROXIMATE  Evaluacion exhaustiva de approximateClosestElement:
%   velocidad y precision vs bvhClosestElement (exacto), sobre 9 mallas
%   (todos los celltypes, uniformes/abolladas/no-uniformes/anisotropas/
%   grandes) x 4 regimenes de puntos (sobre la superficie, cerca, caja
%   media, lejos).
%
%   Metricas por combinacion (nP = 1e5, best-of-3, single-thread):
%     t_apx , t_ex   µs/punto
%     x              speedup = t_ex / t_apx
%     hit%%          puntos con el resultado EXACTO (|d_apx - d_ex| <= tol)
%     p99/h, max/h   percentil 99 y maximo del error (d_apx - d_ex),
%                    normalizados por h = mediana de la primera arista
%   El bench ASSERTA en cada combinacion la propiedad de cota superior:
%   d_apx >= d_ex (el abanico es un subconjunto real de la malla).

  here = fileparts( mfilename('fullpath') );
  addpath( fullfile( here , '..' , 'MESH' ) );
  rng(13);

  MESHES = {};

  %1) triangulos uniformes (esfera 52k)
  V0 = randn( 26000 ,3);  V0 = V0 ./ sqrt( sum( V0.^2 ,2) );
  T0 = convhulln( V0 );
  MESHES(end+1,:) = { 'tri-uniforme' , struct('xyz',V0,'tri',T0) };

  %2) triangulos abollados (no convexa, densidad variable)
  Vb = V0 .* ( 1 + 0.25*sin( 4*V0(:,1) ) .* cos( 3*V0(:,2) ) );
  MESHES(end+1,:) = { 'tri-abollada' , struct('xyz',Vb,'tri',T0) };

  %3) triangulos NO uniformes (muestreo concentrado en un polo: tamanos ~x30)
  Du = randn( 4000 ,3);   Du = Du ./ sqrt( sum( Du.^2 ,2) );
  Dc = randn( 22000 ,3)*0.12 + [ 0 0 1 ];  Dc = Dc ./ sqrt( sum( Dc.^2 ,2) );
  Vn = [ Du ; Dc ];
  MESHES(end+1,:) = { 'tri-NOuniforme' , struct('xyz',Vn,'tri',convhulln(Vn)) };

  %4) triangulos anisotropos (abollada estirada x6 en z: slivers)
  MESHES(end+1,:) = { 'tri-anisotropa' , struct('xyz',Vb.*[1 1 6],'tri',T0) };

  %5) segmentos (helice, 100k)
  s = linspace( 0 , 60*pi , 100001 ).';
  MESHES(end+1,:) = { 'seg-helice' , struct('xyz',[cos(s),sin(s),0.02*s] , ...
                                            'tri',[(1:1e5).',(2:1e5+1).']) };

  %6) tetraedros (delaunay de 15k puntos)
  Xt = randn( 15000 ,3);
  MESHES(end+1,:) = { 'tets-delaunay' , struct('xyz',Xt,'tri',delaunayn(Xt)) };

  %7) mixta (triangulos + segmentos, 0-padded)
  Vm1 = randn( 10000 ,3);  Vm1 = Vm1 ./ sqrt( sum( Vm1.^2 ,2) );
  Tm1 = convhulln( Vm1 );
  sm  = linspace( 0 , 20*pi , 30001 ).';
  Vm2 = [ cos(sm)+3 , sin(sm) , 0.05*sm ];
  Tm2 = [ (1:3e4).' , (2:3e4+1).' ] + size( Vm1 ,1);
  MESHES(end+1,:) = { 'mixta-tri+seg' , struct('xyz',[Vm1;Vm2] , ...
                                    'tri',[Tm1 ; [Tm2 , zeros(3e4,1)]]) };

  %8) nube de puntos (celltype 1): el aproximado debe ser EXACTO
  Vp = randn( 100000 ,3);
  MESHES(end+1,:) = { 'nube-puntos' , struct('xyz',Vp,'tri',(1:1e5).') };

  %9) triangulos grandes (200k tris)
  Vg = randn( 100000 ,3);  Vg = Vg ./ sqrt( sum( Vg.^2 ,2) );
  MESHES(end+1,:) = { 'tri-grande-200k' , struct('xyz',Vg,'tri',convhulln(Vg)) };

  nP = 1e5;
  fprintf( '%d mallas x 4 sets de %g puntos | µs/punto, best-of-3, 1 hilo\n' , size(MESHES,1) , nP );

  for im = 1:size( MESHES ,1)
    nm = MESHES{im,1};  M = MESHES{im,2};
    V  = M.xyz;  T = M.tri;
    scale = norm( max(V,[],1) - min(V,[],1) );

    %h local: mediana de la primera arista (elementos con >= 2 nodos)
    wk = size(T,2) >= 2 && any( T(:,min(2,end)) > 0 );
    if size( T ,2) >= 2
      m2 = T(:,2) > 0;
      h = median( sqrt( sum( ( V(T(m2,2),:) - V(T(m2,1),:) ).^2 ,2) ) );
    else
      h = NaN;                                  %nube de puntos: sin aristas
    end

    tic;  Bx = BVH( M );                          tBx = toc;
    tic;  Ba = approximateClosestElement( M );    tBa = toc;

    fprintf( '\n== %-16s  %7d elems, %6d verts | build: exacto %.0f ms, aprox %.0f ms | h = %.3g\n' , ...
             nm , size(T,1) , size(V,1) , 1e3*tBx , 1e3*tBa , h );
    fprintf( '   %-6s | %7s %7s %6s | %6s %9s %9s\n' , ...
             'set' , 't_apx' , 't_ex' , 'x' , 'hit%' , 'p99/h' , 'max/h' );

    SETS = { 'surf' , elemPoints( V , T , nP ) ; ...
             'near' , elemPoints( V , T , nP ) + 0.005*scale*randn( nP ,3) ; ...
             'mid'  , ( rand( nP ,3) - 0.5 ) .* ( 1.2*( max(V,[],1)-min(V,[],1) ) ) + ( max(V,[],1)+min(V,[],1) )/2 ; ...
             'far'  , 10*scale*randn( nP ,3) };

    for is = 1:size( SETS ,1)
      P = SETS{is,2};

      tx = Inf;  for r = 1:3, tic; [ ~ , ~ , d0 ] = bvhClosestElement( {M,Bx} , P );          tx = min( tx , toc ); end
      ta = Inf;  for r = 1:3, tic; [ ~ , ~ , da ] = approximateClosestElement( {M,Ba} , P );  ta = min( ta , toc ); end

      err = da - d0;
      assert( min( err ) >= -1e-9*scale , '%s/%s: COTA VIOLADA (%g)' , nm , SETS{is,1} , min(err) );
      hit = err <= 1e-12*scale + 1e-9*max( d0 , 0 );
      es  = sort( err );
      p99 = es( ceil( 0.99*numel(es) ) );
      us  = @(t) 1e6*t/nP;
      fprintf( '   %-6s | %7.3f %7.3f %6.2f | %5.1f%% %9.2e %9.2e\n' , ...
               SETS{is,1} , us(ta) , us(tx) , tx/ta , 100*mean(hit) , p99/h , max(err)/h );
    end
  end

  fprintf( ['\n(hit%% = resultado identico al exacto; p99/max del ERROR de distancia\n' ...
            ' normalizado por la arista mediana h; cota d_apx >= d_ex ASSERTADA en todo)\n'] );
end

function P = elemPoints( V , T , n )
%puntos sobre/dentro de elementos aleatorios (baricentricas exponenciales)
  r = randi( size(T,1) , n ,1);
  k = size( T ,2);
  w = -log( rand( n , k ) );
  w( T(r,:) == 0 ) = 0;                       %lanes 0-padded no cuentan
  w = w ./ sum( w ,2);
  P = zeros( n ,3);
  X = V;  X(:,end+1:3) = 0;
  for c = 1:k
    m = T(r,c) > 0;
    P(m,:) = P(m,:) + w(m,c) .* X( T(r(m),c) ,:);
  end
end
