function bench_seedCeiling
%BENCH_SEEDCEILING  Techo teorico de CUALQUIER heuristica de siembra de cotas.
%
%   Pregunta: si la heuristica costara CERO y diera la cota PERFECTA (d0
%   exacto, inflado), ¿cuanto acelera la query exacta? Y ¿que papel juega el
%   WARM START del kernel (el ganador del punto Morton-vecino anterior)?
%
%   Cuatro mediciones por set (mex directo, blobs 'noframe', 1 hilo):
%     plain          warm start ON , sin siembra      (el motor tal cual)
%     plain-nw       warm start OFF, sin siembra      (cuanto vale el warm start)
%     oracle         warm start ON , siembra d0*(1+1e-9)  (techo con heuristica gratis)
%     oracle-nw      warm start OFF, siembra perfecta (la siembra sola)
%
%   Si oracle ~ plain y oracle-nw ~ plain, el kernel ya se auto-siembra (via
%   warm start / descenso best-first) y NINGUNA heuristica de cota puede
%   ganar en modo exacto. Ademas: modo APROXIMADO (devolver el candidato del
%   1-ring sin verificar): coste y tasa de acierto.

  here = fileparts( mfilename('fullpath') );
  addpath( fullfile( here , '..' , 'MESH' ) );
  rng(11);

  nV = 26000;
  V0 = randn( nV ,3);  V0 = V0 ./ sqrt( sum( V0.^2 ,2) );
  T  = convhulln( V0 );
  V  = V0 .* ( 1 + 0.25*sin( 4*V0(:,1) ) .* cos( 3*V0(:,2) ) );
  M  = struct( 'xyz',V , 'tri',T );
  scale = norm( max(V,[],1) - min(V,[],1) );
  tiny  = 1e-12 * scale;

  B  = BVH( M , [] , 'noframe' );          %frame identidad: mex directo valido
  Mv = struct( 'xyz',V , 'tri',(1:nV).' );
  Bv = BVH( Mv , [] , 'noframe' );
  EsuP = meshEsuP( M ,'sparse');

  nPn = 2e5;  nPm = 2e5;  nPf = 1e5;
  Pnear = surfPoints( V , T , nPn ) + 0.005*scale*randn( nPn ,3);
  Pmid  = ( rand( nPm ,3) - 0.5 ) .* ( 1.2*( max(V,[],1)-min(V,[],1) ) ) + ...
          ( max(V,[],1)+min(V,[],1) )/2;
  Pfar  = 10*scale*randn( nPf ,3);
  sets  = { 'near-surface' , Pnear ; 'caja media' , Pmid ; 'far-field' , Pfar };

  fprintf( 'malla: %d vertices, %d triangulos | µs/punto, best-of-3, 1 hilo\n\n' , nV , size(T,1) );
  fprintf( '%-13s | %8s %9s | %8s %10s | %s\n' , 'set' , 'plain' , 'plain-nw' , ...
           'oracle' , 'oracle-nw' , 'aprox: vq+ring, acierto, err max' );
  fprintf( '%s\n' , repmat( '-' , 1 , 100 ) );

  for s = 1:size( sets ,1)
    P  = sets{s,2};
    nP = size( P ,1);
    us = @(t) 1e6*t/nP;

    [ ~ , ~ , d0 ] = bvhClosestElement_mx( P , B , 1 , Inf );
    oracleSeed = d0*(1+1e-9) + tiny;

    t = [ Inf Inf Inf Inf ];
    for r = 1:3
      tic; [ ~,~,dA ] = bvhClosestElement_mx( P , B , 1 , Inf );              t(1)=min(t(1),toc);
      tic; [ ~,~,dB ] = bvhClosestElement_mx( P , B , 1 , Inf        , 1 );   t(2)=min(t(2),toc);
      tic; [ eC,~,dC ] = bvhClosestElement_mx( P , B , 1 , oracleSeed );      t(3)=min(t(3),toc);
      tic; [ eD,~,dD ] = bvhClosestElement_mx( P , B , 1 , oracleSeed , 1 );  t(4)=min(t(4),toc);
    end
    %exactitud (la siembra perfecta inflada no debe perder a nadie)
    assert( max( abs( dA - d0 ) ) == 0 , 'plain no deterministico' );
    assert( max( abs( dB - d0 ) ) < 1e-12*scale , 'plain-nw difiere' );
    okC = eC > 0;  okD = eD > 0;
    assert( nnz(~okC) == 0 && nnz(~okD) == 0 , 'oracle: misses inesperados (%d,%d)' , nnz(~okC) , nnz(~okD) );
    assert( max( abs( dC - d0 ) ) < 1e-12*scale && max( abs( dD - d0 ) ) < 1e-12*scale , 'oracle difiere' );

    %modo APROXIMADO: candidato del 1-ring del vertice mas cercano, SIN verificar
    tq = Inf;
    for r = 1:3, tic; [ ev , ~ , ~ ] = bvhClosestElement_mx( P , Bv , 1 , Inf );  tq = min( tq , toc ); end
    trg = Inf;
    for r = 1:3
      tic;
      cols = EsuP( : , ev );
      [ ti , pi ] = find( cols );
      d2q = pTriD2( P(pi,:) , V(T(ti,1),:) , V(T(ti,2),:) , V(T(ti,3),:) );
      [ ds , o ] = sort( d2q , 'descend' );
      dub = inf( nP ,1);  dub( pi(o) ) = ds;  dub = sqrt( dub );
      trg = min( trg , toc );
    end
    hit = abs( dub - d0 ) <= 1e-12*scale + 1e-9*d0;
    errmax = max( dub - d0 );

    fprintf( '%-13s | %8.3f %9.3f | %8.3f %10.3f | %5.2f+%4.2f µs, %5.1f%%, %.2e\n' , ...
             sets{s,1} , us(t(1)) , us(t(2)) , us(t(3)) , us(t(4)) , ...
             us(tq) , us(trg) , 100*mean(hit) , errmax );
  end

  fprintf( ['\nlectura: oracle = techo de una heuristica GRATIS y PERFECTA;\n' ...
            'plain-nw vs oracle-nw = lo que la siembra externa aporta cuando NO hay warm start.\n'] );
end

function P = surfPoints( V , T , n )
  r  = randi( size(T,1) , n ,1);
  w  = -log( rand( n ,3) );  w = w ./ sum( w ,2);
  P  = w(:,1).*V(T(r,1),:) + w(:,2).*V(T(r,2),:) + w(:,3).*V(T(r,3),:);
end

function d2 = pTriD2( Q , A , B , C )
  ab = B - A;  ac = C - A;  aq = Q - A;
  d1 = sum( ab.*aq ,2);  d2_ = sum( ac.*aq ,2);
  cp = A;  done = ( d1 <= 0 & d2_ <= 0 );
  bq = Q - B;  d3 = sum( ab.*bq ,2);  d4 = sum( ac.*bq ,2);
  m = ~done & d3 >= 0 & d4 <= d3;  cp(m,:) = B(m,:);  done = done | m;
  vc = d1.*d4 - d3.*d2_;
  m = ~done & vc <= 0 & d1 >= 0 & d3 <= 0;
  v = d1 ./ ( d1 - d3 );  cp(m,:) = A(m,:) + v(m).*ab(m,:);  done = done | m;
  cq = Q - C;  d5 = sum( ab.*cq ,2);  d6 = sum( ac.*cq ,2);
  m = ~done & d6 >= 0 & d5 <= d6;  cp(m,:) = C(m,:);  done = done | m;
  vb = d5.*d2_ - d1.*d6;
  m = ~done & vb <= 0 & d2_ >= 0 & d6 <= 0;
  w = d2_ ./ ( d2_ - d6 );  cp(m,:) = A(m,:) + w(m).*ac(m,:);  done = done | m;
  va = d3.*d6 - d5.*d4;
  m = ~done & va <= 0 & ( d4 - d3 ) >= 0 & ( d5 - d6 ) >= 0;
  w = ( d4 - d3 ) ./ ( ( d4 - d3 ) + ( d5 - d6 ) );
  cp(m,:) = B(m,:) + w(m).*( C(m,:) - B(m,:) );  done = done | m;
  m = ~done;
  den = 1 ./ ( va + vb + vc );
  v = vb .* den;  w = vc .* den;
  cp(m,:) = A(m,:) + v(m).*ab(m,:) + w(m).*ac(m,:);
  d2 = sum( ( Q - cp ).^2 ,2);
end
