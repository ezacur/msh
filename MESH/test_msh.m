%TEST_MSH  The msh value class: lazy cache, COW invalidation, queries, fields.

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
  assert( isa( Mty.getField('xyzColor') ,'single') , 'tipos: los campos conservan el suyo' );

  S = struct( 'xyz',V , 'tri',F , 'xyzTemp',rand(700,1) , 'triLabel',rand(size(F,1),1) );
  M2 = msh( S );
  assert( M2.hasField('xyzTemp') && M2.hasField('triLabel') , 'ctor: struct fields' );
  S2 = M2.toStruct();
  assert( isequal( S2.xyzTemp , S.xyzTemp ) && isequal( S2.triLabel , S.triLabel ) , ...
          'toStruct: field roundtrip' );

  Mv = msh( struct( 'vertices',V , 'faces',F ) );        %patch-style aliases
  assert( isequal( Mv.F , F ) , 'ctor: vertices/faces aliases' );
  Mc = msh( M );                                          %copy
  assert( isequal( Mc.V , M.V ) , 'ctor: copy' );
  M0 = msh();
  assert( M0.nV == 0 , 'ctor: empty' );
  X2 = rand( 80 ,2);
  M2d = msh( X2 , delaunayn( X2 ) );
  assert( M2d.nsd == 2 , '2D: nsd from columns' );

  %% 2) cachedProps: valores correctos + cacheo real (acceso con sufijo)
  B1 = M.bbox_;
  assert( isequal( B1 , [ min(V,[],1) ; max(V,[],1) ] ) , 'bbox: wrong' );
  [ sRef , cRef ] = meshSurface( struct('xyz',V,'tri',F) );
  sc = M.surfCent_;                                       %[area, centroide]
  assert( abs( sc(1) - sRef ) < 1e-12 && max(abs( sc(2:end) - cRef )) < 1e-12 , ...
          'surfCent: wrong' );
  assert( isequal( M.boundary_ , MeshBoundary( F ) ) , 'boundary: wrong' );
  assert( isempty( M.boundary_ ) , 'boundary: closed hull must be empty' );
  Mopen = msh( V , F( V(F(:,1),3) < 0.5 ,:) );            %open cap
  assert( ~isempty( Mopen.boundary_ ) , 'boundary: open mesh must have one' );
  assert( M.ct == 5 , 'ct' );
  assert( M.celltype == 5 , 'celltype: alias via subsref' );
  assert( isequal( M.xyz , M.V ) && isequal( M.tri , M.F ) , 'xyz/tri: alias de lectura' );
  assert( isequal( size( M.EsuP_ ) , [ size(F,1) , 700 ] ) , 'EsuP: size' );

  %caching is real: second BVH access must be a hit (>=20x faster)
  Vb = randn( 26000 ,3);  Vb = Vb ./ sqrt( sum( Vb.^2 ,2) );
  Mb = msh( Vb , convhulln( Vb ) );
  tic;  B = Mb.BVH_;   t1 = toc;                          %build (~40 ms)
  tic;  B = Mb.BVH_;   t2 = toc;                          %cache hit
  assert( t2 < t1/20 , 'cache: BVH second access should be a hit (%.1f vs %.1f ms)' , 1e3*t2 , 1e3*t1 );

  %% 3) the acceptance example: M1/M2/M3 cache semantics
  M1 = msh( V , F );
  B1 = M1.BVH_;                                 %fills the SHARED handle
  Mcopy = M1;                                   %sibling copy
  T = [ 0.9*eye(3) , [1;2;3] ; 0 0 0 1 ];
  M3 = M1.transform( T );                       %similarity: BVH FOLDS, not rebuilt
  B3 = M3.BVH_;
  assert( isequal( B3.child4 , B1.child4 ) && isequal( B3.X , B1.X ) , ...
          'transform: the BVH must be the folded one (same hierarchy & X)' );
  assert( ~isequal( B3.frame , B1.frame ) , 'transform: the frame must have moved' );

  M1.F( 1:10 ,:) = [];                          %edit M1 connectivity
  Bc = Mcopy.BVH_;                              %sibling keeps its cache intact
  assert( isequal( Bc.child4 , B1.child4 ) , 'COW: sibling cache was disturbed' );
  Bn = M1.BVH_;                                 %edited instance rebuilds
  assert( Bn.nE == size(F,1) - 10 , 'invalidation: BVH must see the new faces' );

  %deformation: hierarchy survives (lazy REFIT, not rebuild)
  Md = msh( V , F );
  Bd1 = Md.BVH_;
  Md.V = V + 0.02*sin( 5*V(:,[2 3 1]) );        %small deformation
  Bd2 = Md.BVH_;                                %lazy refit on access
  assert( isequal( Bd2.child4 , Bd1.child4 ) && isequal( Bd2.perm , Bd1.perm ) , ...
          'deform: hierarchy must survive (refit)' );
  assert( ~isequal( Bd2.bounds4 , Bd1.bounds4 ) , 'deform: bounds must refresh' );

  %% 4) queries through the class == direct engine calls
  P = randn( 500 ,3) * 1.5;
  [ e1 , cp1 , d1 , bc1 , F1 ] = M.closestElement( P );
  Bref = BVH( struct('xyz',V,'tri',F) );
  [ e2 , ~ , d2 ] = bvhClosestElement( { struct('xyz',V,'tri',F) , Bref } , P );
  assert( max( abs( d1 - d2 ) ) < 1e-12 , 'closestElement: differs from engine' );
  assert( isstruct( F1 ) && numel( F1.type ) == 500 , 'closestElement: F output' );
  [ ~ , ~ , dD ] = M.closestElement( P , 0.2 );           %Dmax through the class
  assert( all( isinf( dD( d1 >= 0.2 ) ) ) , 'closestElement: Dmax misses' );

  rays = [ randn(200,3)*3 , randn(200,3)*0.2 ];
  [ xyz1 , c1 , t1 ] = M.intersectRay( rays , 'first' );
  [ ~ , c2 , t2 ] = bvhIntersectRay( struct('xyz',V,'tri',F) , rays , 'first' );
  w = c1 > 0;
  assert( isequal( w , c2 > 0 ) && max( abs( t1(w) - t2(w) ) ) < 1e-12 , ...
          'intersectRay: differs from engine' );

  %% 5) fields API
  Mf = msh( V , F );
  Mf = Mf.addField( 'xyzTemp' , rand(700,2) );
  Mf = Mf.addField( 'Label' , rand( size(F,1) ,1) );      %inferred -> face
  assert( isequal( size( Mf.getField('xyzTemp') ) , [700 2] ) , 'field: get' );
  assert( Mf.hasField('triLabel') , 'field: inferred face + prefixed access' );
  L = Mf.fieldNames();
  assert( any( strcmp( L.node ,'xyzTemp') ) && any( strcmp( L.face ,'triLabel') ) , 'field: names' );
  Mf = Mf.rmField( 'xyzTemp' );
  assert( ~Mf.hasField('xyzTemp') , 'field: rm' );

  %field reconciliation on node-count change (crop/pad with warning)
  ws = warning( 'off' , 'msh:field' );
  Mf = Mf.addField( 'xyzA' , rand(700,1) );
  Mf.V = [ Mf.V ; 0 0 0 ];
  v = Mf.getField( 'xyzA' );
  warning( ws );
  assert( size( v ,1) == 701 && isnan( v(end) ) , 'field: pad on node growth' );

  %% 6) viz / INFO / plot / textura (en INFO.texture)
  Mp = msh( V , F );
  Mp.viz.FaceColor = 'r';
  Mp.INFO.description  = 'una malla';
  assert( strcmp( Mp.INFO.description , 'una malla' ) , 'INFO' );
  Mp.INFO.texture = uint8( zeros( 4 ,4 ,3) );
  assert( isequal( size( Mp.INFO.texture ) , [4 4 3] ) , 'texture en INFO' );
  Sp = Mp.toStruct();
  assert( isfield( Sp ,'texture') , 'toStruct exporta INFO.texture' );
  h = Mp.plot();  close( ancestor( h ,'figure') );        %smoke test

  %% 7) save/load: cache is Transient, data survives, queries work after load
  Ms = msh( V , F );
  Ms = Ms.addField( 'xyzT' , rand(700,1) );
  [ ~ , ~ , dA ] = Ms.closestElement( P(1:50,:) );
  fn = [ tempname , '.mat' ];
  save( fn , 'Ms' );  L2 = load( fn );  delete( fn );
  assert( isequal( L2.Ms.V , Ms.V ) && L2.Ms.hasField('xyzT') , 'save/load: data' );
  [ ~ , ~ , dB ] = L2.Ms.closestElement( P(1:50,:) );     %cache rebuilt lazily
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
  out = evalc( 'b9 = Mdbg.boundary_; b9 = Mdbg.boundary_;' );
  assert( contains( out , 'MISS' ) && contains( out , 'HIT' ) , 'DEBUG: MISS/HIT' );
  out = evalc( 'Mdbg.V = Mdbg.V * 1.1;' );
  assert( contains( out , 'SET' ) && contains( out , 'EVENT' ) && ...
          contains( out , 'changeCoords' ) , 'DEBUG: SET/EVENT' );
  out = evalc( 'B9 = Mdbg.BVH_; Mdbg.V = Mdbg.V + 0.01; B9 = Mdbg.BVH_;' );
  assert( contains( out , 'MISS' ) && contains( out , 'RPLAY' ) && ...
          contains( out , 'sync absoluto' ) , 'DEBUG: build + replay refit' );
  out = evalc( 'M9 = Mdbg.transform( [ eye(3) , [1;0;0] ; 0 0 0 1 ] ); B9 = M9.BVH_;' );
  assert( contains( out , 'TRANS' ) && contains( out , 'incremental' ) , 'DEBUG: fold perezoso' );
  Mdbg.DEBUG = false;
  out = evalc( 'b9 = Mdbg.boundary_;' );
  assert( isempty( out ) , 'DEBUG off: no debe imprimir' );

  %% 10) display informativo (disp y display comparten displayScalarObject)
  Mdisp = msh( V , F );
  Mdisp = Mdisp.addField( 'xyzTemp' , rand(700,2) );
  Mdisp.viz.FaceColor = [0.8 0.2 0.2];
  Mdisp.INFO.id = 'demo';
  s10 = Mdisp.surfCent_;   %#ok<NASGU> rellena la cache
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
  assert( Mflat.isFlat() && Mflat.isPlanar() , 'flat => planar' );
  Mpla = msh( [ X2 , 2.5 + zeros(80,1) ] , delaunayn( X2 ) );   %planar, no flat
  out = evalc( 'disp( Mpla )' );
  assert( contains( out , 'planar: z = 2.5' ) , 'disp: malla planar' );
  assert( ~Mpla.isFlat() && Mpla.isPlanar() , 'planar y no flat' );
  assert( ~M.isPlanar() && ~M.isFlat() , 'la esfera no es planar' );

  Mmix = msh( V , [ 1 2 0 ; 4 5 6 ] );                     %mixta 0-padded
  out = evalc( 'disp( Mmix )' );
  assert( contains( out , 'mixed: 1 segments + 1 triangles' ) , 'disp: mixta' );

  out = evalc( 'disp( Mb )' );                             %26k vertices (bloque 2)
  assert( contains( out , '26,000 vertices' ) , 'disp: separadores de miles' );

  %% 11) cachedProps genericas: define / alias / eventos / proxy / override
  Mq = msh( V , F );
  Mq = Mq.defineCachedProp( 'halfN' , @(m) m.nV / 2 , 'changeNodeCount' , [] );
  assert( Mq.cached.halfN == 350 , 'cachedProp: computa' );
  assert( Mq.halfN_ == 350 , 'cachedProp: alias sufijo' );

  Mq.DEBUG = true;
  out = evalc( 'v11 = Mq.halfN_;' );
  assert( contains( out , 'HIT' ) , 'cachedProp: hit' );
  out = evalc( 'Mq.V = Mq.V + 0.5; v11 = Mq.halfN_;' );
  assert( contains( out , 'HIT' ) , 'cachedProp: insensible a changeCoords' );
  out = evalc( 'Mq.V = [ Mq.V ; 0 0 0 ]; v11 = Mq.halfN_;' );
  assert( contains( out , 'MISS' ) && v11 == 350.5 , 'cachedProp: invalida con changeNodeCount' );
  Mq.DEBUG = false;

  %proxy: delete (statement, borra solo el valor) / set (conservador) / removeProp
  Mq.cached.halfN.delete;
  Mq.DEBUG = true;
  out = evalc( 'v11 = Mq.halfN_;' );
  assert( contains( out , 'MISS' ) , 'proxy: delete borra el valor' );
  Mq.DEBUG = false;
  Mq2 = Mq.cached.halfN.set( 999 );
  assert( Mq2.halfN_ == 999 , 'proxy: set siembra' );
  assert( Mq.halfN_ == 350.5 , 'proxy: set es aislado (COW)' );
  Mq3 = Mq.cached.halfN.removeProp;
  try
    Mq3.cached.halfN;
    error('test:proxy','removeProp debia dejarla indefinida');
  catch ME
    assert( strcmp( ME.identifier , 'msh:cached' ) , 'proxy: removeProp' );
  end
  assert( Mq.halfN_ == 350.5 , 'proxy: removeProp no toca al original' );

  %handler de evento: obtener e invocar; indexar dentro del valor
  h11 = Mq.cached.BVH.changeCoords;
  assert( isa( h11 , 'function_handle' ) , 'proxy: handler getter' );
  Bh = h11( Mq.BVH_ , Mq );                              %refit manual
  assert( isequal( Bh.child4 , Mq.BVH_.child4 ) , 'proxy: handler invocable' );
  assert( isequal( size( Mq.cached.BVH.frame ) , [4 4] ) , 'proxy: drill al valor' );
  assert( Mq.BVH_.nE == size( F ,1) , 'proxy: drill via alias' );

  %override de una definicion de fabrica (aislado de la instancia original)
  Msp = Mq.defineCachedProp( 'BVH' , @(m) BVH( toStruct( m ) , [] , 'sphere' ) , ...
                             'changeConnectivity' , [] );
  assert( strcmp( Msp.BVH_.volume , 'sphere' ) , 'override: nueva definicion' );
  assert( strcmp( Mq.BVH_.volume , 'aabb' )   , 'override: la original intacta' );

  %replay perezoso: N ediciones -> UN solo sync en el acceso
  Mz = msh( V , F );  B0 = Mz.BVH_;
  Mz.DEBUG = true;
  out = evalc( 'Mz.V = Mz.V + 0.01; Mz.V = Mz.V + 0.01; Bz = Mz.BVH_;' );
  assert( numel( strfind( out , 'EVENT' ) ) == 2 , 'replay: dos eventos anotados' );
  assert( numel( strfind( out , 'RPLAY' ) ) == 1 && contains( out , 'sync absoluto' ) , ...
          'replay: un solo sync para N ediciones' );
  assert( isequal( Bz.child4 , B0.child4 ) && ~isequal( Bz.bounds4 , B0.bounds4 ) , ...
          'replay: jerarquia intacta, bounds refrescados' );

  fprintf( 'ALL msh class tests passed.\n' );
end
