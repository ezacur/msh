function G = meshMetric( M , type )
% G = meshMetric( M , type )
%
% Matriz de "metrica" (producto interno L2) de la malla M.
%   type = 'consistent'   (def.) matriz de masa consistente (FEM P1)
%        = 'lumped'              diagonal, suma de filas de la consistente
%        = 'voronoi'             diagonal de areas de Voronoi mixtas (Meyer 2003)
%        = 'connectivity'        combinatoria: masa consistente con peso
%                                unitario por celda (ignora la geometria;
%                                su patron de no-ceros es el grafo de adyacencia)
%
% Soporta celltype 3 (lineas), 5 (triangulos) y 10 (tetraedros).
% 'voronoi' solo esta implementada para celltype 3 y 5.

  if nargin < 2 || isempty( type ), type = 'consistent'; end
  type = lower( type );

  celltype = meshCelltype( M );

  switch type
    case 'consistent'
      G = consistentMass( M , celltype , [] );

    case 'lumped'
      % agrupamiento por suma de filas de la masa consistente
      % (diagonal, conserva la masa total)
      G = consistentMass( M , celltype , [] );
      n = size( G ,1);
      G = spdiags( full( sum( G ,2) ) , 0 , n , n );

    case 'voronoi'
      G = voronoiMass( M , celltype );

    case 'connectivity'
      % masa consistente con peso unitario por celda -> independiente de
      % la geometria. sum(G(:)) == numero de celdas.
      w = ones( size(M.tri,1) ,1 );
      G = consistentMass( M , celltype , w );

    otherwise
      error('meshMetric:type','tipo "%s" no soportado (consistent, lumped, voronoi o connectivity).',type);
  end

end


function G = consistentMass( M , celltype , w )
% w : peso por celda. Si esta vacio se usa la medida geometrica
%     (longitud/area/volumen); con w = ones(...) se obtiene la
%     version combinatoria, independiente de la geometria.

  if nargin < 3, w = []; end
  sz = [1,1] * size(M.xyz,1);
  if 0
  elseif celltype == 3

    if isempty( w ), w = meshQuality( M ,'length'); end

    G = accumarray( M.tri(:,[1,1]) , w , sz ,[],[],true) +...
        accumarray( M.tri(:,[2,2]) , w , sz ,[],[],true) +...
        accumarray( M.tri(:,[1,2]) , w , sz ,[],[],true);
    G = G / 6;
    G = G + G.';

  elseif celltype == 5

    if isempty( w ), w = meshQuality( M , 'area' ); end

    G = accumarray( M.tri(:,[1,1]) , w , sz ,[],[],true) +...
        accumarray( M.tri(:,[2,2]) , w , sz ,[],[],true) +...
        accumarray( M.tri(:,[3,3]) , w , sz ,[],[],true) +...
        accumarray( M.tri(:,[1,2]) , w , sz ,[],[],true) +...
        accumarray( M.tri(:,[1,3]) , w , sz ,[],[],true) +...
        accumarray( M.tri(:,[2,3]) , w , sz ,[],[],true);
    G = G / 12;
    G = G + G.';

  elseif celltype == 10

    if isempty( w ), w = meshQuality( M , 'volume' ); end

    G = accumarray( M.tri(:,[1,1]) , w , sz ,[],[],true) +...
        accumarray( M.tri(:,[2,2]) , w , sz ,[],[],true) +...
        accumarray( M.tri(:,[3,3]) , w , sz ,[],[],true) +...
        accumarray( M.tri(:,[4,4]) , w , sz ,[],[],true) +...
        accumarray( M.tri(:,[1,2]) , w , sz ,[],[],true) +...
        accumarray( M.tri(:,[1,3]) , w , sz ,[],[],true) +...
        accumarray( M.tri(:,[1,4]) , w , sz ,[],[],true) +...
        accumarray( M.tri(:,[2,3]) , w , sz ,[],[],true) +...
        accumarray( M.tri(:,[2,4]) , w , sz ,[],[],true) +...
        accumarray( M.tri(:,[3,4]) , w , sz ,[],[],true);
    G = G / 20;
    G = G + G.';

  else
    error('meshMetric:celltype','celltype %d no soportado (solo 3, 5 o 10).',celltype);
  end

end


function G = voronoiMass( M , celltype )

  n = size( M.xyz ,1);
  if celltype == 3

    % en una polilinea el dual de Voronoi son medios segmentos ( = lumped )
    L = meshQuality( M , 'length' );
    d = accumarray( M.tri(:) , [ L ; L ]/2 , [n,1] );

  elseif celltype == 5

    [ A , ang , L ] = meshQuality( M , 'area' , 'angles' , 'lengths' );
    %  A   : [Nt x 1] area por triangulo
    %  ang : [Nt x 3] angulos (grados) en los vertices v1,v2,v3
    %  L   : [Nt x 3] longitudes de arista [ l12 , l23 , l31 ]

    l12 = L(:,1); l23 = L(:,2); l31 = L(:,3);
    c1  = cotd( ang(:,1) );   % cot del angulo en v1
    c2  = cotd( ang(:,2) );   % cot del angulo en v2
    c3  = cotd( ang(:,3) );   % cot del angulo en v3

    % areas de Voronoi (circuncentricas) por vertice:
    %   cada arista contribuye (l^2 * cot(angulo_opuesto))/8 a sus 2 extremos
    AV = [ ( l12.^2.*c3 + l31.^2.*c2 ) , ...   % v1: aristas e12 (op. v3) y e31 (op. v2)
           ( l12.^2.*c3 + l23.^2.*c1 ) , ...   % v2: aristas e12 (op. v3) y e23 (op. v1)
           ( l23.^2.*c1 + l31.^2.*c2 ) ] / 8;  % v3: aristas e23 (op. v1) y e31 (op. v2)

    % correccion "mixta" (Meyer et al. 2003) para triangulos obtusos:
    %   vertice obtuso -> A/2 ; los otros dos -> A/4
    obtV = ang > 90;            % vertice obtuso (a lo sumo uno por fila)
    obtT = any( obtV ,2);       % triangulo obtuso
    AV( obtT ,:) = repmat( A(obtT)/4 , 1 , 3 );
    for c = 1:3
      w = obtV(:,c);
      AV( w , c ) = A(w)/2;
    end

    d = accumarray( M.tri(:) , AV(:) , [n,1] );

  else
    error('meshMetric:voronoi','metrica "voronoi" solo implementada para celltype 3 y 5 (no %d).',celltype);
  end

  G = spdiags( d , 0 , n , n );

end
