function M = MeshImproveAcutes( M , L , maxIts )

  try
    [~,~,A] = unique( M.triPIECE );
  catch
    A = zeros( size( M.tri ,1) , 1 );
  end
  for it = 1:5
    [ M.tri , map ] = limpiar_malla_atributos( M.xyz , M.tri , A , L , maxIts );
    for f = fieldnames( M ,'^tri.+').', f = f{1};
      M.(f) = M.(f)( map ,:,:,:,:);
    end
    A = A( map ,:);
    M = MeshTidy( M ,0,true);
  end


  return

  if nargin < 3, maxIts = 1e4; end
  it = 0;
  while 1, it = it+1; if it > maxIts, break; end
    A = min( meshQuality( M ,'angles') ,[],2);
    [l,t] = min( A );
    if l > L, break; end
    disp(l)
    M0 = M;

    [ M.tri , map ] = collapse_triangle_with_map( M.xyz , M.tri , t );
    for f = fieldnames( M ,'^tri.+').', f = f{1};
      M.(f) = M.(f)( map ,:,:,:,:);
    end

    figure(1);
    subplot(211); plotMESH( M0 ,'-','td','PIECE','EdgeAlpha',0.3,'light'); colormap parula; set(gca,'CameraPosition',[1116.5,23.7,80.6],'CameraTarget',[-2.4,-41.1,-9.4],'CameraViewAngle',0.5);
                 hplotMESH( MeshRemoveFaces( M0 ,{t},true) ,'4o-','EdgeColor','r');
    subplot(212); plotMESH( M  ,'-','td','PIECE','EdgeAlpha',0.3,'light'); colormap parula; set(gca,'CameraPosition',[1116.6,23.7,80.6],'CameraTarget',[-2.4,-41.1,-9.4],'CameraViewAngle',0.5); %,'CameraUpVector',[-0.08,-3.e-3,0.997]
    linkprop(gaa(gcf),{'View','XLim','YLIm','ZLim','CameraPosition','CameraTarget','CameraUpVector','CameraViewAngle','clim','Colormap'});
    drawnow


    X = MeshFlatten( M0 ); X0 = X;
    X.tri = M.tri;
    for f = fieldnames( X ,'^tri.+').', f = f{1};
      X.(f) = X.(f)( map ,:,:,:,:);
    end

    figure(2);
     silhouette( plotMESH( X0 ,'-','td','PIECE','EdgeColor','m','EdgeAlpha',0.3) ,'EdgeColor','g','LineWidth',3); colormap parula; drawnow; set(gca,'CameraPosition',[mean( X0.xyz( X0.tri(t,:) ,:) ,1),0]+[0,0,10],'CameraTarget',[mean( X0.xyz( X0.tri(t,:) ,:) ,1),0],'CameraViewAngle',10);
    hplotMESH( MeshRemoveFaces( X0 ,{t},true) ,'2','EdgeColor','r','FaceAlpha',0.4);
    figure(3);
     silhouette( plotMESH( X  ,'-','td','PIECE','EdgeColor','m','EdgeAlpha',0.3) ,'EdgeColor','g','LineWidth',3); colormap parula; drawnow; set(gca,'CameraPosition',[mean( X0.xyz( X0.tri(t,:) ,:) ,1),0]+[0,0,10],'CameraTarget',[mean( X0.xyz( X0.tri(t,:) ,:) ,1),0],'CameraViewAngle',10);
    hplotMESH( MeshRemoveFaces( X0 ,{t},true) ,'2','EdgeColor','r','FaceAlpha',0.4);
    linkprop(gaa([2,3]),{'View','XLim','YLIm','ZLim','CameraPosition','CameraTarget','CameraUpVector','CameraViewAngle','clim','Colormap'});
    drawnow
       drawnow


  end

  M = MeshTidy( M ,0,true);

end



function [F_new, map_indices] = collapse_triangle_with_map(V, F, t_idx)
% COLLAPSE_OPTIMAL_GEOMETRIC
% Colapsa el triangulo t_idx eligiendo la arista que minimiza la distorsion
% geometrica (la mas corta), respetando la topologia y los bordes.
%
% NO MODIFICA V. Solo reestructura F.
%
% Salidas:
%   F_new: Nueva topologia.
%   map_indices: Vector que mapea las caras de F_new a las de F original.

    if t_idx > size(F, 1) || t_idx < 1, error('Indice fuera de rango'); end

    % --- 1. PREPARACION ---
    % Objeto de triangulacion para consultas topologicas eficientes
    oldW = warning( 'off' , 'MATLAB:triangulation:PtsNotInTriWarnId' ); oldW = onCleanup( @()warning(oldW) );
    TR = triangulation(F, V);

%     % Detectar vertices de borde (Boundary Nodes)
%     % freeBoundary devuelve aristas, extraemos los vertices unicos.
%     b_edges = freeBoundary(TR);
%     if ~isempty(b_edges)
%         boundary_nodes = unique(b_edges(:));
%     else
%         boundary_nodes = [];
%     end

    boundary_nodes = unique( MeshBoundary( Mesh(V,F) ).tri );
    
    % Mapa logico para consulta rapida O(1)
    is_border = false(size(V,1), 1);
    is_border(boundary_nodes) = true;

    % Vertices del triangulo objetivo
    t_verts = F(t_idx, :); % [v1, v2, v3]
    
    % Las 3 posibles aristas (pares de vertices)
    candidates = [t_verts(1), t_verts(2);
                  t_verts(2), t_verts(3);
                  t_verts(3), t_verts(1)];
              
    best_keep = -1;
    best_remove = -1;
    min_distortion = inf; % Buscamos minimizar la longitud de la arista colapsada
    found_valid = false;
    
    % --- 2. EVALUACION DE CANDIDATOS ---
    for i = 1:3
        u = candidates(i, 1);
        v = candidates(i, 2);
        
        % A. Verificacion Topologica (Link Condition)
        if ~check_link_condition(TR, u, v)
            continue; % Si rompe la malla, saltamos esta arista
        end
        
        % B. Logica de Borde y Direccion del Colapso
        u_bd = is_border(u);
        v_bd = is_border(v);
        
        keep = -1; remove = -1;
        
        if u_bd && ~v_bd
            % U es borde, V es interno -> V debe morir, U sobrevive.
            keep = u; remove = v;
        elseif ~u_bd && v_bd
            % V es borde, U es interno -> U debe morir, V sobrevive.
            keep = v; remove = u;
        else
            % Ambos borde o ambos internos -> Podemos elegir cualquiera.
            % Por defecto elegimos U como superviviente para medir distancia.
            keep = u; remove = v;
        end
        
        % C. Evaluacion Geometrica (Distancia)
        % Como no movemos vertices (snap), la "distorsion" es proporcional
        % a la distancia entre los puntos. Cuanto mas cerca esten, mejor.
        dist_sq = sum((V(keep,:) - V(remove,:)).^2);
        
        % D. Seleccion del Minimo
        if dist_sq < min_distortion
            min_distortion = dist_sq;
            best_keep = keep;
            best_remove = remove;
            found_valid = true;
        end
    end
    
    % Vector de indices originales
    original_ids = (1:size(F, 1))';
    
    if ~found_valid
        warning('Triangulo %d: No se puede colapsar ninguna arista sin violar topologia.', t_idx);
        F_new = F;
        map_indices = original_ids;
        return;
    end
    
    % --- 3. EJECUCION DEL COLAPSO ---
%     fprintf('Colapsando T%d: Vertice %d absorbe al %d (Distancia^2: %.4f)\n', ...
%             t_idx, best_keep, best_remove, min_distortion);
        
    F_new = F;
    % Reemplazo global: el vertice eliminado se convierte en el superviviente
    F_new(F_new == best_remove) = best_keep;
    
    % --- 4. LIMPIEZA Y MAPEO ---
    % Identificar caras degeneradas (triangulos con vertices repetidos)
    deg_mask = (F_new(:,1) == F_new(:,2)) | ...
               (F_new(:,2) == F_new(:,3)) | ...
               (F_new(:,1) == F_new(:,3));
           
    % Eliminar filas degeneradas
    F_new(deg_mask, :) = [];
    
    % Generar el mapa de indices (para propiedades externas)
    map_indices = original_ids(~deg_mask);
    
end

% --- FUNCION AUXILIAR (Link Condition) ---
function is_safe = check_link_condition_old(TR, u, v)
    % 1. Vecinos de u y v
    att_u = vertexAttachments(TR, u); n_u = unique(TR.ConnectivityList(att_u{1},:));
    att_v = vertexAttachments(TR, v); n_v = unique(TR.ConnectivityList(att_v{1},:));
    
    % Excluirse mutuamente
    n_u(n_u==u | n_u==v) = [];
    n_v(n_v==v | n_v==u) = [];
    
    % 2. Interseccion de vecinos (Link geometrico)
    common_neighbors = intersect(n_u, n_v);
    
    % 3. Caras compartidas por la arista (Link topologico)
    edge_faces_ids = edgeAttachments(TR, u, v);
    
    if isempty(edge_faces_ids{1})
        % Arista flotante o invalida
        is_safe = false; return;
    end
    
    face_verts = TR.ConnectivityList(edge_faces_ids{1}, :);
    link_edge = unique(face_verts(:));
    link_edge(link_edge==u | link_edge==v) = [];
    
    % 4. Validacion: Los vecinos comunes deben ser EXACTAMENTE los opuestos de las caras compartidas
    is_safe = isequal(sort(common_neighbors(:)), sort(link_edge(:)));
end
























































function [F_new, map_indices] = limpiar_malla_atributos(V, F, A, min_angle_deg,max_iter)
% LIMPIAR_MALLA_ATRIBUTOS
% Elimina triangulos con angulos < min_angle_deg, respetando estrictamente
% las fronteras geometricas y las interfaces entre regiones definidas por A.
%
% Entradas:
%   V: Vertices (Nx3) - NO SE MODIFICAN
%   F: Caras (Mx3)
%   A: Vector de atributos (Mx1) - Usado para detectar interfaces
%   min_angle_deg: (Opcional) Default 2.0
%
% Salidas:
%   F_new: Nueva topologia
%   map_indices: Vector tal que A(map_indices) reconstruye los atributos de F_new

    if nargin < 4, min_angle_deg = 2.0; end
    
    % Copias de trabajo internas
    F_curr = F;
    A_curr = A; % Mantenemos A_curr internamente solo para detectar bordes rapido
    
    original_ids = (1:size(F, 1))'; % Este es tu 'map'
    min_angle_rad = deg2rad(min_angle_deg);
    
    has_changed = true;
    iter = 0;
    if nargin < 5, max_iter = size(F,1); end
    
    % Desactivar warning de triangulacion una sola vez
    oldW = warning('off', 'MATLAB:triangulation:PtsNotInTriWarnId');
    cleanW = onCleanup(@() warning(oldW));

    while has_changed && iter < max_iter
        iter = iter + 1;
        has_changed = false; % Asumimos false hasta encontrar candidato
        
        % 1. Topologia actual
        TR = triangulation(F_curr, V);
        
        % 2. Mapeo de Fronteras (Feature Edges)
        % Una arista es "feature" si es borde exterior O separa dos atributos distintos
        edges_list = edges(TR);
        att = edgeAttachments(TR, edges_list);
        
        % Matriz dispersa para consulta rapida O(1)
        num_v = size(V,1);
        is_feature_edge = sparse(edges_list(:,1), edges_list(:,2), false, num_v, num_v);
        is_feature_node = false(num_v, 1);
        
        for i = 1:size(edges_list, 1)
            faces_idx = att{i};
            u = edges_list(i,1);
            v = edges_list(i,2);
            
            es_borde = (length(faces_idx) == 1);
            es_interfaz = false;
            
            if length(faces_idx) == 2
                % Si los triangulos vecinos tienen distinto atributo, es una frontera interna
                if A_curr(faces_idx(1)) ~= A_curr(faces_idx(2))
                    es_interfaz = true;
                end
            end
            
            if es_borde || es_interfaz
                is_feature_edge(u,v) = true;
                is_feature_edge(v,u) = true;
                is_feature_node(u) = true;
                is_feature_node(v) = true;
            end
        end
        
        % 3. Buscar triangulos malos
        min_angles = get_min_angles(V, F_curr); disp( min( min_angles ) )
        bad_tris = find(min_angles < min_angle_rad);
        
        if isempty(bad_tris), break; end
        
        % Ordenar: intentar colapsar primero los triangulos mas agudos
        [~, sort_idx] = sort(min_angles(bad_tris));
        sorted_bad_tris = bad_tris(sort_idx);
        
        % 4. Intentar colapsar
        for k = 1:length(sorted_bad_tris)
            t_idx = sorted_bad_tris(k);
            verts = F_curr(t_idx, :);
            
            % Encontrar vertice con el angulo mas agudo
            angs_tri = get_triangle_angles(V(verts(1),:), V(verts(2),:), V(verts(3),:));
            [~, min_local_idx] = min(angs_tri);
            
            % Definir arista opuesta al angulo agudo
            % min_idx 1 (V1) -> opuesto es arista V2-V3
            if min_local_idx == 1, u = verts(2); v = verts(3);
            elseif min_local_idx == 2, u = verts(3); v = verts(1);
            else, u = verts(1); v = verts(2); end
            
            % --- CHEQUEOS DE SEGURIDAD ---
            
            % A. Topologia (Link Condition)
            if ~check_link_condition(TR, u, v), continue; end
            
            % B. Geometria y Atributos (Feature Preservation)
            u_feat = is_feature_node(u);
            v_feat = is_feature_node(v);
            uv_is_feat = is_feature_edge(u,v);
            
            keep = -1; remove = -1;
            
            % Logica de supervivencia de nodos
            if uv_is_feat
                % La arista es frontera: colapso permitido solo sobre la misma frontera
                keep = u; remove = v;
            elseif ~u_feat && ~v_feat
                % Ambos internos: colapso libre
                keep = u; remove = v;
            elseif u_feat && ~v_feat
                % U frontera, V interno: V muere
                keep = u; remove = v;
            elseif ~u_feat && v_feat
                % U interno, V frontera: U muere
                keep = v; remove = u;
            else
                % Ambos son nodos de frontera, pero la arista NO es frontera.
                % Esto cruzaria una region interior conectando dos bordes -> PROHIBIDO.
                continue; 
            end
            
            % --- EJECUTAR COLAPSO ---
            F_curr(F_curr == remove) = keep;
            
            % Limpieza de degenerados
            deg_mask = (F_curr(:,1) == F_curr(:,2)) | ...
                       (F_curr(:,2) == F_curr(:,3)) | ...
                       (F_curr(:,1) == F_curr(:,3));
            
            F_curr(deg_mask, :) = [];
            A_curr(deg_mask) = [];       % Sincronizar atributos internos
            original_ids(deg_mask) = []; % Sincronizar mapa
            
            has_changed = true;
            break; % Salir y recalcular topologia
        end
    end
    
    F_new = F_curr;
    map_indices = original_ids;
end

% --- AUXILIARES (Sin cambios en logica) ---

function angles = get_min_angles(V, F)
    v1 = V(F(:,1), :); v2 = V(F(:,2), :); v3 = V(F(:,3), :);
    L1 = sqrt(sum((v2-v1).^2, 2)); L2 = sqrt(sum((v3-v2).^2, 2)); L3 = sqrt(sum((v1-v3).^2, 2));
    arg1 = (L1.^2 + L3.^2 - L2.^2) ./ (2 .* L1 .* L3);
    arg2 = (L1.^2 + L2.^2 - L3.^2) ./ (2 .* L1 .* L2);
    arg3 = (L2.^2 + L3.^2 - L1.^2) ./ (2 .* L2 .* L3);
    angles = min(acos(max(min([arg1, arg2, arg3], 1), -1)), [], 2);
end

function angs = get_triangle_angles(p1, p2, p3)
    angs(1) = vec_ang(p2-p1, p3-p1);
    angs(2) = vec_ang(p3-p2, p1-p2);
    angs(3) = vec_ang(p1-p3, p2-p3);
end

function a = vec_ang(v1, v2)
    n = norm(v1)*norm(v2);
    if n==0, a=0; return; end
    a = acos(max(min(dot(v1,v2)/n, 1), -1));
end

function is_safe = check_link_condition(TR, u, v)
    att_u = vertexAttachments(TR, u); att_v = vertexAttachments(TR, v);
    n_u = unique(TR.ConnectivityList(att_u{1}, :)); 
    n_v = unique(TR.ConnectivityList(att_v{1}, :));
    n_u(n_u==u | n_u==v) = []; n_v(n_v==v | n_v==u) = [];
    common = intersect(n_u, n_v);
    edge_faces = edgeAttachments(TR, u, v);
    if isempty(edge_faces{1}), is_safe = false; return; end
    link = unique(TR.ConnectivityList(edge_faces{1}, :));
    link(link==u | link==v) = [];
    is_safe = isequal(sort(common(:)), sort(link(:)));
end
