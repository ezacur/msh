function M_out = MeshMorton( M )
    % applyMortonOrdering: Reordena vértices y caras de una malla usando la
    % curva de Morton 3D (Z-order) para maximizar la localidad en memoria Caché.
    % Detecta automáticamente los campos con prefijos 'xyz' y 'tri'.

    % Validaciones básicas
    if ~isfield(M, 'xyz') || ~isfield(M, 'tri')
        error('La estructura debe contener los campos .xyz y .tri');
    end

    num_V = size(M.xyz, 1);
    num_F = size(M.tri, 1);

    % --- 1. Reordenamiento de Vértices ---
    codes_V = computeMorton3D(M.xyz);
    [~, perm_V] = sort(codes_V);
    
    % Matriz de permutación inversa para actualizar las referencias en las caras
    % Tiene que ser double para que MATLAB la use como índice
    inv_perm_V = zeros(num_V, 1, 'double'); 
    inv_perm_V(perm_V) = 1:num_V;

    % --- 2. Reordenamiento de Caras ---
    % Calculamos el centroide de cada cara antes de alterar los vértices.
    % Usamos la topología y geometría original.
    A = M.xyz(M.tri(:,1), :);
    B = M.xyz(M.tri(:,2), :);
    C = M.xyz(M.tri(:,3), :);
    centroids = (A + B + C) / 3.0;
    
    codes_F = computeMorton3D(centroids);
    [~, perm_F] = sort(codes_F);

    % --- 3. Actualización de Topología (.tri) ---
    % Actualizamos los índices que conforman los triángulos usando la inversa
    M_out = M; % Trabajamos sobre una copia
    M_out.tri = cast(inv_perm_V(M.tri), class(M.tri));

    % --- 4. Permutación Dinámica de Campos ---
    fields = fieldnames(M_out);
    
    for i = 1:length(fields)
        fn = fields{i};
        val = M_out.(fn);
        
        % Omitimos el procesamiento de cosas que no sean matrices por seguridad
        if ~isnumeric(val) && ~islogical(val)
            continue;
        end

        if startsWith(fn, 'xyz') && size(val, 1) == num_V
            % Es un campo de vértices (nV x C)
            if ismatrix(val)
                M_out.(fn) = val(perm_V, :);
            elseif ndims(val) == 3
                M_out.(fn) = val(perm_V, :, :); % Soporte para campos 3D
            end
            
        elseif startsWith(fn, 'tri') && size(val, 1) == num_F
            % Es un campo de caras (nF x C)
            if ismatrix(val)
                M_out.(fn) = val(perm_F, :);
            elseif ndims(val) == 3
                M_out.(fn) = val(perm_F, :, :);
            end
        end
    end
end

%% Funciones Locales (Optimizadas para Vectorización)
function codes = computeMorton3D(pts)
    % Cuantiza puntos 3D a enteros de 21 bits y los intercala en un uint64
    if isempty(pts)
        codes = zeros(0, 1, 'uint64');
        return;
    end

    % Normalización dentro del Bounding Box
    p_min = min(pts, [], 1);
    p_max = max(pts, [], 1);
    extent = max(p_max - p_min, 1e-12); % Evita dividir por cero en mallas degeneradas
    
    % Escalamos al rango de 21 bits (0 a 2.097.151)
    scale = (2^21) - 1;
    q = uint64((pts - p_min) ./ extent .* scale);

    % Extraemos columnas
    x = q(:, 1);
    y = q(:, 2);
    z = q(:, 3);

    % Esparcir los bits (Vectorizado)
    x = spreadBits21(x);
    y = spreadBits21(y);
    z = spreadBits21(z);

    % Intercalar: z | (y << 1) | (x << 2)
    codes = bitor(z, bitor(bitshift(y, 1), bitshift(x, 2)));
end

function v = spreadBits21(v)
    % Separa los 21 bits bajos de cada número dejando 2 espacios vacíos
    % entre ellos. Basado en "Magic Numbers / Dilated Integers".
    v = bitand(v, uint64(0x00000000001FFFFF));
    v = bitand(bitor(v, bitshift(v, 32)), uint64(0x1F00000000FFFF));
    v = bitand(bitor(v, bitshift(v, 16)), uint64(0x1F0000FF0000FF));
    v = bitand(bitor(v, bitshift(v, 8)),  uint64(0x100F00F00F00F00F));
    v = bitand(bitor(v, bitshift(v, 4)),  uint64(0x10C30C30C30C30C3));
    v = bitand(bitor(v, bitshift(v, 2)),  uint64(0x1249249249249249));
end