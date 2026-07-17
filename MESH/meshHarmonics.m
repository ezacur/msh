function [ Phi , lambda , L , M ] = meshHarmonics( C , K , metric )
% [ Phi , lambda , L , M ] = meshHarmonics( C , K , metric )
%
% Modos propios del operador de Laplace-Beltrami (FEM P1) de la malla C,
% resolviendo el problema generalizado    L * phi = lambda * M * phi .
% (los autovectores son los "armonicos" de la malla / manifold harmonics)
%
% Entradas:
%   C      : malla con campos .xyz y .tri. celltype 3 (segmentos),
%            5 (triangulos) o 10 (tetraedros).
%   K      : numero de modos a calcular. Inf (por defecto) -> todos.
%   metric : matriz de masa, se pasa tal cual a meshMetric:
%            'lumped' (def.) | 'consistent' | 'voronoi' | 'connectivity'.
%
% Salidas:
%   Phi    : [Nverts x K] autofunciones, M-ortonormales (Phi.'*M*Phi = I).
%   lambda : [K x 1]       autovalores en orden creciente (lambda(1) ~ 0).
%   L      : [Nverts x Nverts] matriz de rigidez (stiffness), dispersa.
%   M      : [Nverts x Nverts] matriz de masa ( = meshMetric(C,metric) ).

  if nargin < 2 || isempty( K )      , K = Inf;          end
  if nargin < 3 || isempty( metric ) , metric = 'lumped'; end

  celltype = meshCelltype( C );
  n = size( C.xyz ,1);

  % ---- matriz de masa (metrica) ----
  M = meshMetric( C , metric );
  M = ( M + M.' )/2;                         % simetriza por seguridad numerica
  if any( full(diag(M)) <= 0 )
    warning('meshHarmonics:singularMass',...
      'Hay vertices con masa <= 0 (no referenciados o celdas degeneradas); M es singular.');
  end

  % ---- matriz de rigidez (stiffness) FEM P1 ----
  L = stiffnessMatrix( C , celltype , n );
  L = ( L + L.' )/2;

  % ---- numero de modos ----
  if ~isfinite( K ), K = n; end
  K = max( 1 , min( round(K) , n ) );

  % ---- problema de autovalores generalizado   L phi = lambda M phi ----
  if K >= n-1
    % todos (o casi todos) los modos: denso
    [ Phi , D ] = eig( full(L) , full(M) );
    lambda = real( diag(D) );
  else
    % pocos modos: disperso con shift-invert.
    %   resolvemos (L + s*M) phi = (lambda + s) M phi  y restamos s.
    %   s = 1e-8 * (escala tipica de autovalor) regulariza la singularidad
    %   de L (cuyo nucleo es la funcion constante) sin sesgar el espectro.
    s = 1e-8 * ( full(sum(diag(L))) / full(sum(diag(M))) );
    [ Phi , D ] = eigs( L + s*M , M , K , 'smallestabs' );
    lambda = real( diag(D) ) - s;
  end

  % ---- ordenar ascendente y recortar a K ----
  [ lambda , o ] = sort( lambda , 'ascend' );
  Phi    = real( Phi(:,o) );
  lambda = lambda(1:K);
  Phi    = Phi(:,1:K);

  % ---- M-ortonormalizar  (Phi.'*M*Phi = I) ----
  s = sqrt( max( sum( Phi .* (M*Phi) ,1) , realmin ) );
  Phi = Phi ./ s;

  % ---- fijar signo (componente de mayor magnitud positiva) -> reproducible ----
  nc = size( Phi ,2);
  [ ~ , im ] = max( abs(Phi) ,[],1 );
  sg = sign( Phi( sub2ind( size(Phi) , im , 1:nc ) ) );
  sg( sg == 0 ) = 1;
  Phi = Phi .* sg;

end


function L = stiffnessMatrix( C , celltype , n )
% Matriz de rigidez FEM P1 (Laplaciano de masa-cero, filas suman 0, PSD).

  t = C.tri;
  if celltype == 3

    P1 = C.xyz( t(:,1) ,:);  P2 = C.xyz( t(:,2) ,:);
    len = sqrt( sum( (P2-P1).^2 ,2) );
    v   = 1 ./ len;                               % rigidez local 1D = 1/l
    I = [ t(:,1) ; t(:,2) ; t(:,1) ; t(:,2) ];
    J = [ t(:,1) ; t(:,2) ; t(:,2) ; t(:,1) ];
    S = [   v    ;   v    ;  -v    ;  -v    ];

  elseif celltype == 5

    % Laplaciano cotangente (FEM P1 sobre superficie)
    ang = meshQuality( C , 'angles' );            % [Nt x 3] grados en v1,v2,v3
    w12 = cotd( ang(:,3) )/2;    % arista (1,2) opuesta a v3
    w23 = cotd( ang(:,1) )/2;    % arista (2,3) opuesta a v1
    w31 = cotd( ang(:,2) )/2;    % arista (3,1) opuesta a v2
    I = [ t(:,1) ; t(:,2) ; t(:,3) ;  t(:,1) ; t(:,2) ;  t(:,2) ; t(:,3) ;  t(:,3) ; t(:,1) ];
    J = [ t(:,1) ; t(:,2) ; t(:,3) ;  t(:,2) ; t(:,1) ;  t(:,3) ; t(:,2) ;  t(:,1) ; t(:,3) ];
    S = [ w12+w31 ; w12+w23 ; w23+w31 ;  -w12 ; -w12 ;  -w23 ; -w23 ;  -w31 ; -w31 ];

  elseif celltype == 10

    % rigidez FEM P1 del tetraedro:  K_ij = V * (grad phi_i . grad phi_j)
    P = C.xyz;  P(:,end+1:3) = 0;
    P1 = P( t(:,1) ,:); P2 = P( t(:,2) ,:); P3 = P( t(:,3) ,:); P4 = P( t(:,4) ,:);
    e1 = P2-P1;  e2 = P3-P1;  e3 = P4-P1;
    cr1 = cross( e2 , e3 , 2 );
    cr2 = cross( e3 , e1 , 2 );
    cr3 = cross( e1 , e2 , 2 );
    vol6 = sum( e1 .* cr1 ,2);                    % = 6 * volumen (con signo)
    g2 = cr1 ./ vol6;  g3 = cr2 ./ vol6;  g4 = cr3 ./ vol6;  % grad phi_2..4
    g1 = -( g2 + g3 + g4 );                                  % grad phi_1
    V  = abs( vol6 )/6;                                      % volumen
    g  = { g1 , g2 , g3 , g4 };

    nt = size( t ,1);
    I = zeros( 16*nt ,1 );  J = I;  S = I;
    k = 0;
    for a = 1:4
      for b = 1:4
        rows = k*nt + (1:nt);
        I(rows) = t(:,a);
        J(rows) = t(:,b);
        S(rows) = V .* sum( g{a} .* g{b} ,2);
        k = k + 1;
      end
    end

  else
    error('meshHarmonics:celltype','celltype %d no soportado (solo 3, 5 o 10).',celltype);
  end

  L = sparse( I , J , S , n , n );

end
