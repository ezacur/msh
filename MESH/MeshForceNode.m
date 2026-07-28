function [ M , ids ] = MeshForceNode( M , X )
%MESHFORCENODE  Insert points into a triangular mesh, forcing them to be nodes.
%
%   [ M , ids ] = MeshForceNode( M , X )
%
%   Modifies the triangular surface mesh M so that every point of X becomes a
%   VERTEX of the mesh. Each point is first PROJECTED onto the surface of M, and
%   then inserted according to where its projection lands:
%     * ON an existing vertex  -> reused, no new geometry.
%     * INSIDE a face          -> the face is split into 3 (a fan around the
%                                 new vertex): 1 triangle -> 3 triangles.
%     * ON an edge             -> the (up to 2) faces sharing that edge are
%                                 each split in 2.
%   So X are snapped to the surface, NOT kept at their input coordinates; the
%   inserted node sits at the projected location  M.xyz( abs(ids) , : ).
%
%   INPUTS
%     M : struct with fields .xyz (nV x 3) and .tri (nF x 3). Any extra field
%         named  xyz*  (e.g. xyzF, xyzRGB) is treated as a PER-VERTEX attribute
%         and interpolated barycentrically at inserted vertices; any field
%         named  tri*  (e.g. triF) is a PER-FACE attribute and is copied from
%         the parent face onto the child faces. (The plain 'xyz'/'tri' — 3-char
%         names — are the geometry; other 3-char fields are ignored.)
%     X : P x 3 list of points to force into the mesh.
%
%   OUTPUTS
%     M   : the densified mesh (more vertices / faces than the input).
%     ids : P x 1. For each input point i, WHERE it ended up. The node itself
%           is always  M.xyz( abs(ids(i)) , : ). The value encodes the kind:
%
%             ids(i) > 0  and <= nV0  : coincided with an EXISTING vertex
%                                       (nV0 = size of the ORIGINAL M.xyz;
%                                        record it before calling if you need
%                                        to tell this case apart from the next)
%             ids(i) >  nV0           : inserted INSIDE a face (face split 1->3)
%             ids(i) <  0             : inserted ON an edge; its vertex index
%                                       is  -ids(i)  (edge faces split in 2)
%
%           i.e. the SIGN flags edge-insertions (negative) vs vertex/interior
%           (positive), and the MAGNITUDE is always the vertex index in the
%           returned M. On return every ids(i) is nonzero and finite.
%
%   CLASSIFICATION -- no geometric tolerance. Where a projection lands is decided
%   by the BVH engine's own region cascade (bvhClosestElement's 5th output,
%   F.type: 1 = vertex, 2 = edge, 3 = face interior), which is the SAME cascade
%   that produced the closest point, so the answer is consistent with it by
%   construction. That replaces the hand-tuned absolute THRESHOLD this function
%   used to carry, which conflated two different units -- a DISTANCE (vertex
%   coincidence) and a DIMENSIONLESS barycentric (interior-vs-edge) -- in one
%   constant, and was therefore scale-dependent: 1e-10 is 0.1 nm on a mesh in
%   metres and enormous on a mesh in micrometres. The region classification is
%   combinatorial and scale-invariant (verified on a mesh scaled by 1e-6).
%
%   The edge parameter is taken from the region-exact barycentrics, NOT
%   re-derived as a chord-length ratio |X-A|/|B-A| (which is only correct when
%   the point lies exactly on the segment) as this function used to do.
%
%   One point is inserted per face and per edge per pass, so several points on
%   the same face/edge are handled in successive passes. Each pass assigns at
%   least one point, so the loop terminates in at most numel(X)/3 passes; there
%   is no tolerance relaxation and nothing is ever silently snapped.
%
%   No VTK. Uses the BVH engine (bvhClosestElement); add BVH\ to the path.
%
% See also bvhClosestElement, BVH, MeshSubdivide, sphereMesh.

  if isempty( M.xyz ) || isempty( M.tri ) || isempty( X )
    ids = zeros( size(X,1) ,1);  return;
  end
  if size( M.tri ,2) ~= 3
    error('MeshForceNode:tri','only TRIANGLE surfaces (nF x 3) are supported.');
  end

  %atributos por nodo (xyz*) y por cara (tri*); los nombres de 3 letras son la
  %geometria y cualquier otro campo de 3 letras se ignora
  xyzAtts = {};  triAtts = {};
  for f = fieldnames( M ).', f = f{1};
    if numel( f ) <= 3, continue; end
    if strncmp( f , 'xyz' ,3 )
      xyzAtts{end+1,1} = f;                                              %#ok<AGROW>
      if isnumeric( M.(f) ), M.(f) = double( M.(f) ); end
    elseif strncmp( f , 'tri' ,3 )
      triAtts{end+1,1} = f;                                              %#ok<AGROW>
    end
  end
  M.xyz = double( M.xyz );
  M.tri = double( M.tri );

  %proyectar UNA vez: insertar densifica la conectividad pero NO mueve la
  %superficie, asi que estas coordenadas siguen siendo validas en cada pasada
  [ ~ , cp ] = bvhClosestElement( struct('xyz',M.xyz,'tri',M.tri) , double( X ) );
  X = cp;  X(:,4) = 0;                       % col 4 = nodo asignado (0 = pendiente)

  maxPass = size( X ,1) + 2;                 % cota: cada pasada asigna >= 1 punto
  for pass = 1:maxPass
    p = find( ~X(:,4) );
    if isempty( p ), break; end

    S = struct( 'xyz' , M.xyz , 'tri' , M.tri );
    [ eid , ~ , ~ , bc , F ] = bvhClosestElement( S , X(p,1:3) );
    if any( eid == 0 )
      error('MeshForceNode:miss','%d point(s) found no closest element.', sum(eid==0) );
    end

    %TODA la clasificacion se extrae AHORA, contra la conectividad que vio la
    %consulta y antes de mutar nada. Es imprescindible: la rama de caras modifica
    %M.tri (cambia la 3a columna de las caras partidas), asi que si la rama de
    %aristas leyera M.tri despues, la cara de un punto-de-arista podria ya no ser
    %la que la consulta clasifico -> vertices e interpolacion equivocados.
    if any( F.type == 0 | F.type > 3 )
      error('MeshForceNode:type','unexpected region code %d.', max(F.type) );
    end
    isV = F.type == 1;
    isF = F.type == 3;
    isE = F.type == 2;

    %(a) vertices: que nodo es
    vNode = zeros( numel(p) ,1);
    if any( isV )
      [ ~ , kv ] = max( bc(isV,:) ,[],2);
      vNode(isV) = M.tri( sub2ind( size(M.tri) , eid(isV) , kv ) );
    end
    %(b) aristas: arista CANONICA (n1 < n2) y parametro, sacados del bc
    eN1 = zeros(numel(p),1); eN2 = zeros(numel(p),1); eT = zeros(numel(p),1);
    if any( isE )
      wE = find( isE );
      [ ~ , kmin ] = min( bc(wE,:) ,[],2);
      KE = [ 2 3 ; 1 3 ; 1 2 ];              %los dos que NO son el minimo
      ke = KE( kmin ,:);
      T3 = M.tri( eid(wE) ,:);
      B3 = bc( wE ,:);
      lin = (1:numel(wE)).';
      a1 = T3( sub2ind(size(T3), lin , ke(:,1) ) );
      a2 = T3( sub2ind(size(T3), lin , ke(:,2) ) );
      c1 = B3( sub2ind(size(B3), lin , ke(:,1) ) );
      c2 = B3( sub2ind(size(B3), lin , ke(:,2) ) );
      sw = a2 < a1;                           %arista CANONICA: a1 < a2
      [ a1(sw) , a2(sw) ] = deal( a2(sw) , a1(sw) );
      [ c1(sw) , c2(sw) ] = deal( c2(sw) , c1(sw) );
      eN1(wE) = a1;  eN2(wE) = a2;  eT(wE) = c2 ./ ( c1 + c2 );
    end

    %--- (1) cayo sobre un VERTICE: se reutiliza, sin geometria nueva
    if any( isV )
      X( p(isV) ,4) = vNode(isV);
    end

    %--- (2) INTERIOR de una cara: fan 1 -> 3.  Una por cara y pasada.
    w = find( isF );
    if ~isempty( w )
      [ ~ , u ] = unique( eid(w) , 'first' );
      w  = w( sort(u) );
      pp = p(w);  FF = eid(w);  BB = bc(w,:);
      nV = size( M.xyz ,1);  nP = numel( pp );
      X( pp ,4) = nV + (1:nP).';

      M.xyz = [ M.xyz ; X(pp,1:3) ];
      for a = xyzAtts(:).', a = a{1};
        M.(a)( nV+(1:nP) ,:,:,:,:) = M.(a)( M.tri(FF,1) ,:,:,:,:) .* BB(:,1) + ...
                                     M.(a)( M.tri(FF,2) ,:,:,:,:) .* BB(:,2) + ...
                                     M.(a)( M.tri(FF,3) ,:,:,:,:) .* BB(:,3);
      end
      %la cara original pasa a (v1,v2,new) y se anaden (v2,v3,new) y (v3,v1,new)
      newV  = X(pp,4);
      M.tri = [ M.tri ; M.tri(FF,[2 3]) , newV ; M.tri(FF,[3 1]) , newV ];
      M.tri( FF ,3) = newV;
      for a = triAtts(:).', a = a{1};
        M.(a) = [ M.(a) ; M.(a)( FF ,:,:,:) ; M.(a)( FF ,:,:,:) ];
      end
    end

    %--- (3) sobre una ARISTA: cada cara incidente se parte en 2.
    %  n1/n2/tt ya venian calculados ARRIBA, contra la conectividad original.
    w = find( isE );
    if isempty( w ), continue; end
    n1 = eN1(w);  n2 = eN2(w);  tt = eT(w);

    %una por arista y pasada
    [ ~ , u ] = unique( [n1 n2] , 'rows' , 'first' );
    u  = sort(u);
    cp_ = p( w(u) );  cn1 = n1(u);  cn2 = n2(u);  ct = tt(u);

    %incidencia arista -> caras, construida UNA vez por pasada. Antes esto era un
    %ismember sobre TODO M.tri por cada punto (O(P*nF)) y ademas sobre un M.tri
    %que crecia dentro del bucle.
    nF   = size( M.tri ,1);
    Eall = sort( [ M.tri(:,[1 2]) ; M.tri(:,[2 3]) ; M.tri(:,[3 1]) ] ,2);
    Fall = repmat( (1:nF).' ,3,1);
    [ Eu , ~ , g ] = unique( Eall , 'rows' );
    [ gs , o ] = sort( g );
    Fs   = Fall( o );
    st   = find( [ true ; diff(gs) ~= 0 ] );
    cnt  = diff( [ st ; numel(gs)+1 ] );
    [ okE , gq ] = ismember( [cn1 cn2] , Eu , 'rows' );
    if ~all( okE )
      error('MeshForceNode:edge','an edge belongs to no face (corrupt connectivity?).');
    end

    %SELECCION: ademas de una por arista, ninguna CARA puede partirse dos veces en
    %la misma pasada -- si no, la incidencia de arriba se queda obsoleta a mitad
    %del bucle (partir la cara (a,b,c) por [b,c] la deja como (a,m,c), y la arista
    %[a,b] pasa a vivir en la cara NUEVA, que el mapa no conoce). Los puntos no
    %seleccionados esperan a la pasada siguiente.
    claimed = false( nF ,1);
    sel = false( numel(cp_) ,1);
    for e = 1:numel( cp_ )
      fs = Fs( st(gq(e)) : st(gq(e)) + cnt(gq(e)) - 1 );
      if any( claimed(fs) ), continue; end
      claimed(fs) = true;  sel(e) = true;
    end
    pp = cp_(sel);  n1 = cn1(sel);  n2 = cn2(sel);  tt = ct(sel);  gq = gq(sel);

    nV = size( M.xyz ,1);
    newV = nV + (1:numel(pp)).';
    X( pp ,4) = -newV;
    M.xyz = [ M.xyz ; X(pp,1:3) ];
    for a = xyzAtts(:).', a = a{1};
      M.(a)( newV ,:,:,:,:) = M.(a)( n1 ,:,:,:,:) .* (1-tt) + ...
                              M.(a)( n2 ,:,:,:,:) .*    tt;
    end

    for e = 1:numel( pp )
      fs = Fs( st(gq(e)) : st(gq(e)) + cnt(gq(e)) - 1 );
      F0 = M.tri( fs ,:);  F1 = F0;
      F0( F0 == n1(e) ) = newV(e);    M.tri( fs ,:) = F0;    %mitad hacia n2
      F1( F1 == n2(e) ) = newV(e);    M.tri = [ M.tri ; F1 ]; %mitad hacia n1
      for a = triAtts(:).', a = a{1};
        M.(a) = [ M.(a) ; M.(a)( fs ,:,:,:) ];
      end
    end
  end

  if ~all( X(:,4) )
    error('MeshForceNode:stuck','%d point(s) could not be inserted.', sum(~X(:,4)) );
  end
  ids = X(:,4);

end
