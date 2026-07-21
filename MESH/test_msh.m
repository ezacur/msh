%TEST_MSH  The msh value class: lazy cache, COW invalidation, queries, fields.
%
%  API por caso: MAYUSCULAS = propiedades (V,F,VIZ,INFO,DEBUG,CP), Capitalized
%  = metodos (Plot,Transform,DefineCP,...), minuscula = CPs (bvh,esup,...).
%  Acceso a CPs: M.bvh LEE (perezoso), M.bvh_ RECALCULA, M.CP.bvh = control.

function test_msh
  rng(31);
  root = fileparts( fileparts( mfilename('fullpath') ) );
  addpath( root );                              %@msh + cacheHandle + cacheView
  addpath( fullfile( root , 'tools' ) );        %transform
  addpath( fullfile( root , 'BVH' ) );          %motor BVH/bvhClosestElement/...

  %% 1) construction forms
  V = randn( 700 ,3);  V = V ./ sqrt( sum( V.^2 ,2) );
  F = convhulln( V );
  M = msh( V , F );
  assert( isequal( M.V , V ) && isequal( M.F , F ) , 'ctor: V/F' );
  assert( M.nV == 700 && M.nF == size(F,1) && M.nsd == 3 , 'ctor: scalars' );
  assert( isa( M.V ,'double') && isa( M.F ,'int32') , 'tipos: contrato' );
  Mty = msh( single( V ) , int64( F ) );          %cualquier numerico entra...
  assert( isa( Mty.V ,'double') && isa( Mty.F ,'int32') , 'tipos: coercion' );
  Mty = msh( V , F , 'xyzColor' , single( rand(700,3) ) );
  assert( isa( Mty.GetField('xyzColor') ,'single') , 'tipos: los campos conservan el suyo' );

  S = struct( 'xyz',V , 'tri',F , 'xyzTemp',rand(700,1) , 'triLabel',rand(size(F,1),1) );
  M2 = msh( S );
  assert( M2.HasField('xyzTemp') && M2.HasField('triLabel') , 'ctor: struct fields' );
  S2 = M2.ToStruct();
  assert( isequal( S2.xyzTemp , S.xyzTemp ) && isequal( S2.triLabel , S.triLabel ) , ...
          'ToStruct: field roundtrip' );

  Mv = msh( struct( 'vertices',V , 'faces',F ) );        %patch-style aliases
  assert( isequal( Mv.F , F ) , 'ctor: vertices/faces aliases' );
  Mc = msh( M );                                          %copy
  assert( isequal( Mc.V , M.V ) , 'ctor: copy' );
  M0 = msh();
  assert( M0.nV == 0 , 'ctor: empty' );
  X2 = rand( 80 ,2);
  M2d = msh( X2 , delaunayn( X2 ) );
  assert( M2d.nsd == 2 , '2D: nsd from columns' );

  %% 2) CPs: valores correctos + cacheo real (lectura con nombre desnudo)
  B1 = M.bbox;
  assert( isequal( B1 , [ min(V,[],1) ; max(V,[],1) ] ) , 'bbox: wrong' );
  [ sRef , cRef ] = meshSurface( struct('xyz',V,'tri',F) );
  sc = M.surfCent;                                        %[area, centroide]
  assert( abs( sc(1) - sRef ) < 1e-12 && max(abs( sc(2:end) - cRef )) < 1e-12 , ...
          'surfCent: wrong' );
  assert( isequal( M.boundary , MeshBoundary( F ) ) , 'boundary: wrong' );
  assert( isempty( M.boundary ) , 'boundary: closed hull must be empty' );
  Mopen = msh( V , F( V(F(:,1),3) < 0.5 ,:) );            %open cap
  assert( ~isempty( Mopen.boundary ) , 'boundary: open mesh must have one' );
  assert( M.ct == 5 , 'ct' );
  assert( M.celltype == 5 , 'celltype: alias via subsref' );
  assert( isequal( M.xyz , M.V ) && isequal( M.tri , M.F ) , 'xyz/tri: alias de lectura' );
  assert( isequal( size( M.esup ) , [ size(F,1) , 700 ] ) , 'esup: size' );

  %caching is real: second bvh access must be a hit (>=20x faster)
  Vb = randn( 26000 ,3);  Vb = Vb ./ sqrt( sum( Vb.^2 ,2) );
  Mb = msh( Vb , convhulln( Vb ) );
  tic;  B = Mb.bvh;   t1 = toc;                           %build (~40 ms)
  tic;  B = Mb.bvh;   t2 = toc;                           %cache hit
  assert( t2 < t1/20 , 'cache: bvh second access should be a hit (%.1f vs %.1f ms)' , 1e3*t2 , 1e3*t1 );

  %% 3) the acceptance example: M1/M2/M3 cache semantics
  M1 = msh( V , F );
  B1 = M1.bvh;                                  %fills the SHARED handle
  Mcopy = M1;                                   %sibling copy
  T = [ 0.9*eye(3) , [1;2;3] ; 0 0 0 1 ];
  M3 = M1.Transform( T );                       %similarity: bvh FOLDS, not rebuilt
  B3 = M3.bvh;
  assert( isequal( B3.child4 , B1.child4 ) && isequal( B3.X , B1.X ) , ...
          'Transform: the bvh must be the folded one (same hierarchy & X)' );
  assert( ~isequal( B3.frame , B1.frame ) , 'Transform: the frame must have moved' );

  M1.F( 1:10 ,:) = [];                          %edit M1 connectivity
  Bc = Mcopy.bvh;                               %sibling keeps its cache intact
  assert( isequal( Bc.child4 , B1.child4 ) , 'COW: sibling cache was disturbed' );
  Bn = M1.bvh;                                  %edited instance rebuilds
  assert( Bn.nE == size(F,1) - 10 , 'invalidation: bvh must see the new faces' );

  %deformation: hierarchy survives (lazy REFIT, not rebuild)
  Md = msh( V , F );
  Bd1 = Md.bvh;
  Md.V = V + 0.02*sin( 5*V(:,[2 3 1]) );        %small deformation
  Bd2 = Md.bvh;                                 %lazy refit on access
  assert( isequal( Bd2.child4 , Bd1.child4 ) && isequal( Bd2.perm , Bd1.perm ) , ...
          'deform: hierarchy must survive (refit)' );
  assert( ~isequal( Bd2.bounds4 , Bd1.bounds4 ) , 'deform: bounds must refresh' );

  %% 4) queries through the class == direct engine calls
  P = randn( 500 ,3) * 1.5;
  [ e1 , cp1 , d1 , bc1 , F1 ] = M.ClosestElement( P );
  Bref = BVH( struct('xyz',V,'tri',F) );
  [ e2 , ~ , d2 ] = bvhClosestElement( { struct('xyz',V,'tri',F) , Bref } , P );
  assert( max( abs( d1 - d2 ) ) < 1e-12 , 'ClosestElement: differs from engine' );
  assert( isstruct( F1 ) && numel( F1.type ) == 500 , 'ClosestElement: F output' );
  [ ~ , ~ , dD ] = M.ClosestElement( P , 0.2 );           %Dmax through the class
  assert( all( isinf( dD( d1 >= 0.2 ) ) ) , 'ClosestElement: Dmax misses' );

  rays = [ randn(200,3)*3 , randn(200,3)*0.2 ];
  [ xyz1 , c1 , t1 ] = M.IntersectRay( rays , 'first' );
  [ ~ , c2 , t2 ] = bvhIntersectRay( struct('xyz',V,'tri',F) , rays , 'first' );
  w = c1 > 0;
  assert( isequal( w , c2 > 0 ) && max( abs( t1(w) - t2(w) ) ) < 1e-12 , ...
          'IntersectRay: differs from engine' );

  %% 5) fields API
  Mf = msh( V , F );
  Mf = Mf.AddField( 'xyzTemp' , rand(700,2) );
  Mf = Mf.AddField( 'Label' , rand( size(F,1) ,1) );      %inferred -> face
  assert( isequal( size( Mf.GetField('xyzTemp') ) , [700 2] ) , 'field: get' );
  assert( Mf.HasField('triLabel') , 'field: inferred face + prefixed access' );
  L = Mf.FieldNames();
  assert( any( strcmp( L.node ,'xyzTemp') ) && any( strcmp( L.face ,'triLabel') ) , 'field: names' );
  Mf = Mf.RmField( 'xyzTemp' );
  assert( ~Mf.HasField('xyzTemp') , 'field: rm' );

  %field reconciliation on node-count change (crop/pad with warning)
  ws = warning( 'off' , 'msh:field' );
  Mf = Mf.AddField( 'xyzA' , rand(700,1) );
  Mf.V = [ Mf.V ; 0 0 0 ];
  v = Mf.GetField( 'xyzA' );
  warning( ws );
  assert( size( v ,1) == 701 && isnan( v(end) ) , 'field: pad on node growth' );

  %% 6) VIZ / INFO / Plot / textura (en INFO.texture)
  Mp = msh( V , F );
  Mp.VIZ.FaceColor = 'r';
  Mp.INFO.description  = 'una malla';
  assert( strcmp( Mp.INFO.description , 'una malla' ) , 'INFO' );
  Mp.INFO.texture = uint8( zeros( 4 ,4 ,3) );
  assert( isequal( size( Mp.INFO.texture ) , [4 4 3] ) , 'texture en INFO' );
  Sp = Mp.ToStruct();
  assert( isfield( Sp ,'texture') , 'ToStruct exporta INFO.texture' );
  h = Mp.Plot();  close( ancestor( h ,'figure') );        %smoke test

  %% 7) save/load: cache is Transient, data survives, queries work after load
  Ms = msh( V , F );
  Ms = Ms.AddField( 'xyzT' , rand(700,1) );
  [ ~ , ~ , dA ] = Ms.ClosestElement( P(1:50,:) );
  fn = [ tempname , '.mat' ];
  save( fn , 'Ms' );  L2 = load( fn );  delete( fn );
  assert( isequal( L2.Ms.V , Ms.V ) && L2.Ms.HasField('xyzT') , 'save/load: data' );
  [ ~ , ~ , dB ] = L2.Ms.ClosestElement( P(1:50,:) );     %cache rebuilt lazily
  assert( max( abs( dA - dB ) ) < 1e-12 , 'save/load: queries after load' );

  %% 8) guards
  Mg = msh( V , F );
  try
    Mg.V = V( 1:10 ,:);                         %faces would dangle
    error('test:guard','shrinking V under F must error');
  catch ME
    assert( strcmp( ME.identifier , 'msh:nodes' ) , 'guard: wrong error %s' , ME.identifier );
  end
  try
    Mg.F = [ 1 2 , 9999 ];                      %missing node
    error('test:guard','F referencing missing vertices must error');
  catch ME
    assert( strcmp( ME.identifier , 'msh:faces' ) , 'guard: wrong error %s' , ME.identifier );
  end
  try
    Mg.F = [ 1 2 2.5 ];                         %non-integer face index
    error('test:guard','non-integer F must error');
  catch ME
    assert( strcmp( ME.identifier , 'msh:faces' ) , 'guard: wrong error %s' , ME.identifier );
  end

  %% 9) DEBUG: narracion de la cadena de procesos (eventos + replay)
  Mdbg = msh( V , F , 'DEBUG' , true );
  assert( Mdbg.DEBUG , 'DEBUG: por constructor' );
  out = evalc( 'b9 = Mdbg.boundary; b9 = Mdbg.boundary;' );
  assert( contains( out , 'MISS' ) && contains( out , 'HIT' ) , 'DEBUG: MISS/HIT' );
  out = evalc( 'Mdbg.V = Mdbg.V * 1.1;' );
  assert( contains( out , 'SET' ) && contains( out , 'EVENT' ) && ...
          contains( out , 'changeCoords' ) , 'DEBUG: SET/EVENT' );
  out = evalc( 'B9 = Mdbg.bvh; Mdbg.V = Mdbg.V + 0.01; B9 = Mdbg.bvh;' );
  assert( contains( out , 'MISS' ) && contains( out , 'RPLAY' ) && ...
          contains( out , 'sync absoluto' ) , 'DEBUG: build + replay refit' );
  out = evalc( 'M9 = Mdbg.Transform( [ eye(3) , [1;0;0] ; 0 0 0 1 ] ); B9 = M9.bvh;' );
  assert( contains( out , 'TRANS' ) && contains( out , 'incremental' ) , 'DEBUG: fold perezoso' );
  Mdbg.DEBUG = false;
  out = evalc( 'b9 = Mdbg.boundary;' );
  assert( isempty( out ) , 'DEBUG off: no debe imprimir' );

  %% 10) display informativo (disp y display comparten displayScalarObject)
  Mdisp = msh( V , F );
  Mdisp = Mdisp.AddField( 'xyzTemp' , rand(700,2) );
  Mdisp.VIZ.FaceColor = [0.8 0.2 0.2];
  Mdisp.INFO.id = 'demo';
  s10 = Mdisp.surfCent;   %#ok<NASGU> rellena la cache
  out = evalc( 'disp( Mdisp )' );
  assert( contains( out , '(triangles)' ) , 'disp: tipo de celdas' );
  assert( contains( out , 'xyzTemp [700x2 double]' ) , 'disp: campos con dimensiones' );
  assert( contains( out , 'FaceColor = [0.8 0.2 0.2]' ) , 'disp: vizProp con valor' );
  assert( contains( out , 'id = ''demo''' ) , 'disp: info' );
  assert( ~isempty( regexp( out , 'surfCent\s+\[1x4 double\]' ,'once') ) , 'disp: vector largo resumido' );
  assert( contains( out , 'definidas sin calcular' ) , 'disp: registradas sin valor' );

  Mflat = msh( [ X2 , zeros(80,1) ] , delaunayn( X2 ) );   %flat: z == 0 exacto
  out = evalc( 'disp( Mflat )' );
  assert( contains( out , 'flat (z = 0)' ) , 'disp: malla flat' );
  assert( Mflat.IsFlat() && Mflat.IsPlanar() , 'flat => planar' );
  Mpla = msh( [ X2 , 2.5 + zeros(80,1) ] , delaunayn( X2 ) );   %planar, no flat
  out = evalc( 'disp( Mpla )' );
  assert( contains( out , 'planar: z = 2.5' ) , 'disp: malla planar' );
  assert( ~Mpla.IsFlat() && Mpla.IsPlanar() , 'planar y no flat' );
  assert( ~M.IsPlanar() && ~M.IsFlat() , 'la esfera no es planar' );

  Mmix = msh( V , [ 1 2 0 ; 4 5 6 ] );                     %mixta 0-padded
  out = evalc( 'disp( Mmix )' );
  assert( contains( out , 'mixed: 1 segments + 1 triangles' ) , 'disp: mixta' );

  out = evalc( 'disp( Mb )' );                             %26k vertices (bloque 2)
  assert( contains( out , '26,000 vertices' ) , 'disp: separadores de miles' );

  %% 11) CPs genericas: define / eventos / proxy CP / override
  Mq = msh( V , F );
  Mq = Mq.DefineCP( 'halfN' , @(m) m.nV / 2 , 'changeNodeCount' , [] );
  assert( Mq.CP.halfN == 350 , 'CP: computa (proxy)' );
  assert( Mq.halfN == 350 , 'CP: lectura desnuda' );

  Mq.DEBUG = true;
  out = evalc( 'v11 = Mq.halfN;' );
  assert( contains( out , 'HIT' ) , 'CP: hit' );
  out = evalc( 'Mq.V = Mq.V + 0.5; v11 = Mq.halfN;' );
  assert( contains( out , 'HIT' ) , 'CP: insensible a changeCoords' );
  out = evalc( 'Mq.V = [ Mq.V ; 0 0 0 ]; v11 = Mq.halfN;' );
  assert( contains( out , 'MISS' ) && v11 == 350.5 , 'CP: invalida con changeNodeCount' );
  Mq.DEBUG = false;

  %proxy: delete (statement, borra solo el valor) / set (conservador) / removeProp
  Mq.CP.halfN.delete;
  Mq.DEBUG = true;
  out = evalc( 'v11 = Mq.halfN;' );
  assert( contains( out , 'MISS' ) , 'proxy: delete borra el valor' );
  Mq.DEBUG = false;
  Mq2 = Mq.CP.halfN.set( 999 );
  assert( Mq2.halfN == 999 , 'proxy: set siembra' );
  assert( Mq.halfN == 350.5 , 'proxy: set es aislado (COW)' );
  Mq3 = Mq.CP.halfN.removeProp;
  try
    Mq3.CP.halfN;
    error('test:proxy','removeProp debia dejarla indefinida');
  catch ME
    assert( strcmp( ME.identifier , 'msh:cached' ) , 'proxy: removeProp' );
  end
  assert( Mq.halfN == 350.5 , 'proxy: removeProp no toca al original' );

  %handler de evento: obtener e invocar; indexar dentro del valor
  h11 = Mq.CP.bvh.changeCoords;
  assert( isa( h11 , 'function_handle' ) , 'proxy: handler getter' );
  Bh = h11( Mq.bvh , Mq );                               %refit manual
  assert( isequal( Bh.child4 , Mq.bvh.child4 ) , 'proxy: handler invocable' );
  assert( isequal( size( Mq.CP.bvh.frame ) , [4 4] ) , 'proxy: drill al valor' );
  assert( Mq.bvh.nE == size( F ,1) , 'drill via lectura desnuda' );

  %override de una definicion de fabrica (aislado de la instancia original)
  Msp = Mq.DefineCP( 'bvh' , @(m) BVH( ToStruct( m ) , [] , 'sphere' ) , ...
                     'changeConnectivity' , [] );
  assert( strcmp( Msp.bvh.volume , 'sphere' ) , 'override: nueva definicion' );
  assert( strcmp( Mq.bvh.volume , 'aabb' )   , 'override: la original intacta' );

  %replay perezoso: N ediciones -> UN solo sync en el acceso
  Mz = msh( V , F );  B0 = Mz.bvh;
  Mz.DEBUG = true;
  out = evalc( 'Mz.V = Mz.V + 0.01; Mz.V = Mz.V + 0.01; Bz = Mz.bvh;' );
  assert( numel( strfind( out , 'EVENT' ) ) == 2 , 'replay: dos eventos anotados' );
  assert( numel( strfind( out , 'RPLAY' ) ) == 1 && contains( out , 'sync absoluto' ) , ...
          'replay: un solo sync para N ediciones' );
  assert( isequal( Bz.child4 , B0.child4 ) && ~isequal( Bz.bounds4 , B0.bounds4 ) , ...
          'replay: jerarquia intacta, bounds refrescados' );

  %% 12) semantica nueva: nombre desnudo LEE, sufijo '_' RECALCULA
  Mn2 = msh( V , F );
  Mn2 = Mn2.DefineCP( 'rndTag' , @(m) rand() );        %sin eventos: sobrevive a todo
  r1 = Mn2.rndTag;                                     %MISS: computa y guarda
  assert( Mn2.rndTag == r1 , 'bare: segunda lectura debe ser HIT (mismo valor)' );
  r2 = Mn2.rndTag_;                                    %RECMP: valor nuevo
  assert( r2 ~= r1 , 'sufijo: debe recalcular (rand nuevo)' );
  assert( Mn2.rndTag == r2 , 'sufijo: el recalculo queda guardado' );
  Mn2.DEBUG = true;
  out = evalc( 'v12 = Mn2.rndTag_;' );
  assert( contains( out , 'RECMP' ) , 'sufijo: narra RECMP' );
  out = evalc( 'v12 = Mn2.rndTag;' );
  assert( contains( out , 'HIT' ) , 'bare tras recalculo: HIT' );
  Mn2.DEBUG = false;
  %set sembrado + '_' lo pisa con el valor verdadero
  Mn3 = Mn2.CP.rndTag.set( -1 );
  assert( Mn3.rndTag == -1 , 'set: el valor sembrado se lee desnudo' );
  v12b = Mn3.rndTag_;
  assert( v12b ~= -1 && Mn3.rndTag == v12b , 'sufijo: pisa el valor sembrado' );
  %drill en ambos accesos
  assert( isequal( size( Mn2.bvh.frame ) , [4 4] ) , 'bare: drill al valor' );
  assert( Mn2.bvh_.nE == size( F ,1) , 'sufijo: drill tras recalcular' );
  %validacion DefineCP: caso y reservados
  try
    Mn2.DefineCP( 'Bad' , @(m) 1 );
    error('test:val','DefineCP debia rechazar Capitalized');
  catch ME
    assert( strcmp( ME.identifier , 'msh:cached' ) , 'val caso: %s' , ME.identifier );
  end
  try
    Mn2.DefineCP( 'xyz' , @(m) 1 );
    error('test:val','DefineCP debia rechazar un alias reservado');
  catch ME
    assert( strcmp( ME.identifier , 'msh:cached' ) , 'val reservado: %s' , ME.identifier );
  end
  %el nombre historico 'cached' ya no resuelve (renombrado a M.CP)
  ok12 = false;
  try, Mn2.cached; catch, ok12 = true; end                 %#ok<VUNUS>
  assert( ok12 , 'M.cached debe errar (el proxy es M.CP)' );

  fprintf( 'ALL msh class tests passed.\n' );
end
