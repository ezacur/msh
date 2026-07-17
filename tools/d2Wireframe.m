function [ds, sid, cp, bc] = d2Wireframe(xyz, M, prec)
% d2Wireframe  Distancia minima de puntos a un wireframe (segmentos).
%
%   [ds, sid, cp, bc] = d2Wireframe(xyz, M)
%   [ds, sid, cp, bc] = d2Wireframe(xyz, M, 'single')
%
% Entradas:
%   xyz : (N x nsd) puntos de consulta (nsd = 2 o 3)
%   M   : struct con  .xyz (V x nsd) vertices
%                     .tri (S x 2)   segmentos (indices en .xyz)
%   'single' (opcional): hace el barrido punto-segmento en simple
%       precision (~2x menos trafico de memoria). El ganador se
%       reconstruye siempre en doble, asi que ds/cp/bc son exactos para
%       el segmento elegido; solo el desempate entre segmentos casi
%       equidistantes (d^2 dentro de ~1e-7 relativo) puede diferir.
%
% Salidas:
%   ds  : (N x 1)   distancia euclidea al segmento mas cercano
%   sid : (N x 1)   fila de M.tri del segmento mas cercano
%   cp  : (N x nsd) punto mas cercano sobre el wireframe
%   bc  : (N x 2)   baricentricas [w1 w2] tal que cp = w1*P1 + w2*P2
%
% Drop-in optimizado (la referencia se conserva como d2Wireframe_old),
% reformulado para no materializar arrays (nsd x N x S). MAS preciso que la
% referencia lejos del origen (centra en el centroide: err 3e-14 vs 3e-8 de
% _old con offset 1e7). Con w = x - p1 y t = clamp( dot(w,v)/|v|^2 , [0,1] ):
%       d^2 = |w|^2 - 2*t*dot(w,v) + t^2*|v|^2
% Todo el barrido se reduce a dos GEMM sobre coordenadas aumentadas
% [X 1], que pliegan las constantes por segmento dentro del producto:
%       t_raw = [X 1] * [V'/|v|^2 ; -(p1.v)/|v|^2]
%       G     = [X 1] * [2*P1'    ; -|p1|^2      ]   ( = 2*X*P1' - |p1|^2 )
%       score = t.*|v|^2.*(t - 2*t_raw) - G          ( = d^2 - |x|^2 )
% |x|^2 es constante por fila y se omite del argmin. El punto cercano se
% reconstruye solo para el ganador y ds se recalcula exacto a partir de
% el, de modo que la cancelacion del desarrollo de Gram solo puede
% afectar al desempate entre segmentos casi equidistantes (y se acota
% centrando las coordenadas en el centroide del wireframe).

  % --- Validacion
  if ~isstruct(M) || ~isfield(M,'xyz') || ~isfield(M,'tri')
    error('d2Wireframe:input', 'M debe ser struct con campos .xyz y .tri');
  end
  if size(M.tri,2) ~= 2
    error('d2Wireframe:input', 'M.tri debe tener 2 columnas (segmentos).');
  end
  if isempty(M.tri)
    error('d2Wireframe:input', 'M.tri esta vacio.');
  end
  if size(xyz,2) ~= size(M.xyz,2)
    error('d2Wireframe:input', 'xyz y M.xyz deben tener igual dimension.');
  end
  use_single = nargin >= 3 && strcmpi(prec, 'single');

  N   = size(xyz, 1);
  nsd = size(xyz, 2);
  S   = size(M.tri, 1);

  % --- Precalculo de segmentos, centrado en el centroide del wireframe
  mu = mean(M.xyz, 1);
  P1 = M.xyz(M.tri(:,1), :) - mu;        % S x nsd
  V  = M.xyz(M.tri(:,2), :) - mu - P1;   % S x nsd
  X  = xyz - mu;                         % N x nsd

  vv = sum(V.^2, 2).';                   % 1 x S, |v|^2
  inv_vv = 1 ./ vv;
  inv_vv(vv <= eps(max(vv))) = 0;        % degenerados -> t = 0 (distancia a P1)
  cs = sum(P1 .* V, 2).' .* inv_vv;      % 1 x S, (p1.v)/|v|^2
  np = sum(P1.^2, 2).';                  % 1 x S, |p1|^2

  W1 = [V.' .* inv_vv; -cs];             % (nsd+1) x S : t_raw = [X 1]*W1
  W2 = [2*P1.'; -np];                    % (nsd+1) x S : G = 2*X*P1' - |p1|^2

  Xs = X;
  if use_single
    Xs = single(X);  W1 = single(W1);  W2 = single(W2);  vv = single(vv);
  end

  % --- Salidas
  ds  = zeros(N, 1);
  sid = zeros(N, 1);
  cp  = zeros(N, nsd);
  bc  = zeros(N, 2);

  % --- Lotes: ~4 temporales (n x S) vivos -> pico ~ 4*8*MAX_ELEMS bytes
  MAX_ELEMS = 1e7;
  batch = max(1, min(N, floor(MAX_ELEMS / S)));

  for i = 1:batch:N
    j = min(i + batch - 1, N);
    n = j - i + 1;

    Xa = [Xs(i:j,:), ones(n, 1, 'like', Xs)];

    A = Xa * W1;                         % n x S, t sin clampear
    T = min(max(A, 0), 1);               % n x S, t clampeado

    % score = d^2 - |x|^2 = t.*vv.*(t - 2*t_raw) - (2*X*P1' - np)
    A = T - 2*A;
    A = A .* T;
    A = A .* vv;
    A = A - Xa * W2;

    [~, bi] = min(A, [], 2);

    tb  = double( T( (bi-1)*n + (1:n).' ) );   % t del segmento ganador
    cpb = P1(bi,:) + tb .* V(bi,:);            % punto cercano (frame centrado)

    D = X(i:j,:) - cpb;                  % ds exacto: sin error de Gram
    ds(i:j)   = sqrt(sum(D.^2, 2));
    sid(i:j)  = bi;
    cp(i:j,:) = cpb + mu;
    bc(i:j,:) = [1-tb, tb];
  end
end
