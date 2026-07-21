function bench_fanClosestElement
%BENCH_FANCLOSESTELEMENT  fan (aproximado) vs bvhClosestElement (exacto).
%
%   Por malla y regimen: velocidad (min-of-3), % de acierto del elemento y
%   error del resto (relativo al tamano LOCAL). Asevera la cota d_fan >= d_ex
%   en cada celda. La etapa 1 reutiliza un Bnodes precalculado (el caso de uso
%   honesto: quien itera, cachea); el fan CSR se recalcula en cada llamada
%   (parte del contrato de fanClosestElement).
%
%   NOTA historica: el antiguo approximateClosestElement (blob fusionado +
%   kernels pt4/fan4 AVX) era mas rapido en mallas puras de triangulos; se
%   retiro a cambio de la composabilidad (semillas externas, point-BVH
%   estandar). Este bench documenta lo que cuesta esa decision.

  addpath( fullfile( fileparts( mfilename('fullpath') ) , '..' , 'MESH' ) );
  rng(11);

  MESHES = {};
  V = randn( 26000 ,3);  V = V ./ sqrt( sum( V.^2 ,2) );
  MESHES{end+1} = { 'esfera 52k tri' , struct( 'xyz',V , 'tri',convhulln(V) ) };
  s = linspace( 0 , 40*pi , 30001 ).';
  MESHES{end+1} = { 'helice 30k seg' , struct( 'xyz',[cos(s),sin(s),0.05*s] , 'tri',[(1:3e4).',(2:3e4+1).'] ) };
  Xt = randn( 20000 ,3);
  MESHES{end+1} = { 'delaunay 135k tet' , struct( 'xyz',Xt , 'tri',delaunayn(Xt) ) };

  nQ = 1e5;
  for m = 1:numel( MESHES )
    name = MESHES{m}{1};  M = MESHES{m}{2};
    ext  = max( max(M.xyz) - min(M.xyz) );
    B    = BVH( M );
    Mn   = struct( 'xyz',M.xyz , 'tri',(1:size(M.xyz,1)).' );
    Bn   = BVH( Mn );
    %tamano local: media de las aristas incidentes al vertice mas cercano
    hloc = ext / sqrt( size( M.tri ,1) );

    REG = { 'surf' , 0 ; 'near' , 0.01 ; 'mid' , 0.3 ; 'far' , 3 };
    fprintf( '%-18s' , name );
    for r = 1:size( REG ,1)
      w = randi( size(M.xyz,1) , nQ ,1);
      P = M.xyz(w,:) + REG{r,2}*ext*randn( nQ ,3);

      tx = Inf;  for k = 1:3, tic; [ ~ , ~ , d0 ] = bvhClosestElement( {M,B} , P );   tx = min(tx,toc); end
      tf = Inf;  for k = 1:3, tic; [ ef , ~ , df ] = fanClosestElement( {M,Bn} , P ); tf = min(tf,toc); end
      assert( all( df >= d0 - 1e-12 ) , '%s/%s: COTA VIOLADA' , name , REG{r,1} );
      hit = abs( df - d0 ) <= 1e-12 + 1e-9*max(d0,1);
      err = max( ( df - d0 ) / hloc );
      fprintf( '  %s x%.2f (hit %.0f%%, err/h %.2g)' , REG{r,1} , tx/tf , 100*mean(hit) , err );
    end
    fprintf( '\n' );
  end
end
