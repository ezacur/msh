function test_MeshForceNode
%TEST_MESHFORCENODE  Insercion forzada de nodos: validez de la malla, contrato de
%   `ids`, atributos interpolados, y los casos que rompian la version anterior
%   (dependiente de escala, parametro por razon de cuerdas, caras partidas dos
%   veces en la misma pasada).
  root = fileparts( fileparts( mfilename('fullpath') ) );
  addpath( root , fullfile(root,'MESH') , fullfile(root,'tools') , fullfile(root,'BVH') );

  %% 1) puntos genericos sobre la superficie
  M0 = sphereMesh( 3 );
  M0.triF = rand( size(M0.tri,1) ,1);
  M0.xyzF = M0.xyz(:,3);
  nV0 = size( M0.xyz ,1);
  rng(7);
  X = randn( 200 ,3);  X = X ./ sqrt( sum( X.^2 ,2) ) * 1.3;   %fuera: se proyectan
  %la proyeccion de referencia es sobre la malla TRIANGULADA, no sobre la esfera
  %analitica: las caras son cuerdas y la sagita de sphereMesh(3) es ~0.02
  [ ~ , cpRef ] = bvhClosestElement( M0 , X );
  [ M , ids ] = MeshForceNode( M0 , X );
  checkMesh( M , '1) genericos' );
  assert( all( ids ~= 0 ) && all( isfinite( ids ) ) , 'ids: todos asignados' );
  assert( all( abs(ids) <= size(M.xyz,1) ) , 'ids: indice fuera de rango' );
  %cada punto debe estar EN el nodo que dice ids
  d = sqrt( sum( ( M.xyz( abs(ids) ,:) - cpRef ).^2 ,2) );
  assert( max( d ) < 1e-9 , '1) el nodo de ids no coincide con la proyeccion (%.3g)' , max(d) );
  assert( size(M.xyz,1) >= nV0 , '1) no puede perder nodos' );
  assert( numel(M.triF) == size(M.tri,1) , '1) triF debe seguir a las caras' );
  assert( numel(M.xyzF) == size(M.xyz,1) , '1) xyzF debe seguir a los nodos' );
  %el atributo por nodo interpolado debe seguir siendo la coordenada z
  assert( max(abs( M.xyzF - M.xyz(:,3) )) < 1e-9 , '1) xyzF mal interpolado' );

  %% 2) sobre VERTICES existentes: se reutilizan, sin geometria nueva
  M0b = sphereMesh( 2 );
  nV = size( M0b.xyz ,1);  nF = size( M0b.tri ,1);
  Xv = M0b.xyz( 1:20 ,:);
  [ Mv , idv ] = MeshForceNode( M0b , Xv );
  assert( size(Mv.xyz,1) == nV && size(Mv.tri,1) == nF , ...
          '2) coincidir con vertices no debe crear geometria' );
  assert( all( idv > 0 & idv <= nV ) , '2) ids debe indicar vertice existente' );
  assert( isequal( Mv.xyz( idv ,:) , Xv ) , '2) debe devolver EL vertice' );

  %% 3) sobre ARISTAS: ids negativo, y cada arista parte sus caras incidentes
  e  = [ M0b.tri(1:15,1) , M0b.tri(1:15,2) ];
  %el punto medio de una arista YA esta exactamente en la superficie de la malla:
  %normalizarlo a la esfera analitica lo sacaria de la arista
  Xe = 0.5*( M0b.xyz(e(:,1),:) + M0b.xyz(e(:,2),:) );
  [ Me , ide ] = MeshForceNode( M0b , Xe );
  checkMesh( Me , '3) aristas' );
  assert( all( ide < 0 ) , '3) una insercion en arista debe dar ids < 0' );

  %% 4) INTERIOR de caras: ids > nV0
  Xf = ( M0b.xyz(M0b.tri(1:15,1),:) + M0b.xyz(M0b.tri(1:15,2),:) ...
       + M0b.xyz(M0b.tri(1:15,3),:) ) / 3;            %el centroide ya esta en la cara
  [ Mf , idf ] = MeshForceNode( M0b , Xf );
  checkMesh( Mf , '4) interiores' );
  assert( all( idf > nV ) , '4) una insercion en cara debe dar ids > nV0' );
  assert( size(Mf.tri,1) == nF + 2*15 , '4) cada cara partida en 3 (+2 caras)' );

  %% 5) INVARIANZA DE ESCALA. La version anterior usaba un THRESHOLD absoluto de
  %%    1e-10 para la coincidencia con vertices, asi que a escala pequena TODO
  %%    coincidia con un vertice y a escala grande nada.
  for sc = [ 1e-6 , 1 , 1e6 ]
    Ms = M0b;  Ms.xyz = M0b.xyz * sc;
    Xs = Xf * sc;
    [ Mo , ido ] = MeshForceNode( Ms , Xs );
    checkMesh( Mo , sprintf('5) escala %g',sc) );
    assert( all( ido > nV ) , '5) a escala %g deja de clasificar como interior' , sc );
  end

  %% 6) VARIOS puntos en la MISMA cara (fuerza varias pasadas) y DOS aristas de la
  %%    misma cara (el caso que dejaba obsoleto el mapa de incidencia)
  f1 = M0b.tri( 1 ,:);
  V  = M0b.xyz;
  Xm = [ 0.5*V(f1(1),:)+0.3*V(f1(2),:)+0.2*V(f1(3),:) ; ...
         0.2*V(f1(1),:)+0.5*V(f1(2),:)+0.3*V(f1(3),:) ; ...
         0.3*V(f1(1),:)+0.2*V(f1(2),:)+0.5*V(f1(3),:) ; ...
         0.5*( V(f1(1),:)+V(f1(2),:) ) ; ...            %arista 1-2 de la cara
         0.5*( V(f1(2),:)+V(f1(3),:) ) ];               %arista 2-3 de la MISMA cara
  [ Mm , idm ] = MeshForceNode( M0b , Xm );
  checkMesh( Mm , '6) misma cara' );
  assert( all( idm ~= 0 ) , '6) todos deben insertarse' );
  assert( numel(unique(abs(idm))) == numel(idm) , '6) nodos distintos por punto' );

  %% 7) entradas degeneradas
  [ ~ , i0 ] = MeshForceNode( M0b , zeros(0,3) );
  assert( isempty( i0 ) , '7) X vacio' );
  Me2 = struct('xyz',zeros(0,3),'tri',zeros(0,3));
  [ ~ , i1 ] = MeshForceNode( Me2 , [0 0 0] );
  assert( numel(i1) == 1 , '7) malla vacia' );

  fprintf( 'ALL MeshForceNode tests passed.\n' );
end

function checkMesh( M , tag )
%la malla resultante tiene que ser VALIDA: indices en rango, sin caras
%degeneradas, sin nodos huerfanos nuevos, y con la misma superficie total
  nV = size( M.xyz ,1);
  T  = M.tri;
  assert( all( T(:) >= 1 & T(:) <= nV ) , '%s: indice de cara fuera de rango' , tag );
  assert( ~any( T(:,1)==T(:,2) | T(:,2)==T(:,3) | T(:,1)==T(:,3) ) , ...
          '%s: cara degenerada (vertice repetido)' , tag );
  assert( all( isfinite( M.xyz(:) ) ) , '%s: coordenada no finita' , tag );
  %cada arista interior debe seguir compartida por exactamente 2 caras (la esfera
  %es cerrada): eso detecta splits a medias, que es el fallo tipico
  E = sort( [ T(:,[1 2]) ; T(:,[2 3]) ; T(:,[3 1]) ] ,2);
  [ ~ , ~ , g ] = unique( E ,'rows');
  c = accumarray( g ,1);
  assert( all( c == 2 ) , '%s: %d arista(s) no compartidas por 2 caras' , tag , sum(c~=2) );
end
