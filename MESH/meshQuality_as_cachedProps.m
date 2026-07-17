function M = meshQuality_as_cachedProps( M , varargin )
%MESHQUALITY_AS_CACHEDPROPS  Las metricas de meshQuality como cachedProps de un msh.
%
%   M = meshQuality_as_cachedProps( M )                    todas las de su celltype
%   M = meshQuality_as_cachedProps( M , 'area' , ... )     solo las pedidas
%
%   EJEMPLO de definicion DINAMICA de cachedProps (no forma parte de la clase):
%   registra UNA ENTRADA POR METRICA, cada una con su computeFcn delegando en
%   meshQuality y, donde la matematica lo permite, un handler del evento
%   'transform' que ACTUALIZA el valor cacheado sin recomputar nada:
%
%     adimensionales (ratios, angulos)    INVARIANTES ante semejanzas
%     longitudes / areas                  escalan s / s^2 (semejanzas)
%     volume / signedvolume / orientation EXACTOS ante CUALQUIER afin: |det(A)|,
%                                         det(A), sign(det(A))
%     normal                              rota: N * (det(R)*R).'  (el factor det
%                                         cubre las reflexiones)
%
%   Si la T no cumple lo que el handler exige (p.ej. escala anisotropa sobre un
%   handler de semejanza), el handler LANZA y la cascada del replay recomputa —
%   nunca un valor silenciosamente incorrecto.
%
%   Acceso y gestion, como cualquier cachedProp:
%       q = M.aspectratio_;             % == M.cached.aspectratio
%       M.cached                        % tabla (estado + eventos)
%       M = M.cached.area.removeProp;   % desregistrar una
%
%   Invalidacion: cualquier edicion de coordenadas fuera de transform() y
%   cualquier edicion de conectividad tiran el valor (se recomputa al pedirlo).
%
%   NOTAS: (1) registra con el NOMBRE CANONICO de la metrica ('area',
%   'aspectratio', ...): si ya tenias una cachedProp con ese nombre, la
%   REDEFINE. (2) meshQuality comparte intermedios entre metricas pedidas en
%   UNA llamada; aqui cada metrica se computa por separado — la cache paga en
%   el acceso REPETIDO, no en el primer barrido (para un barrido unico llama a
%   meshQuality directamente). (3) el conjunto valido depende del celltype en
%   el momento del registro; si luego cambias el tipo de celdas, las
%   definiciones obsoletas erroran al recomputar (honesto).
%
% See also msh/defineCachedProp, meshQuality, msh_CLASS_TUTORIAL.md (secc. 6).

  %tabla por celltype: { nombre canonico , tipo de update ante transform }
  %  'inv'  invariante (semejanzas)      'len'  escala s (semejanzas)
  %  'area' escala s^2 (semejanzas)      'adet' factor |det(A)| (afin exacto)
  %  'det'  factor det(A) (afin exacto)  'sgn'  factor sign(det(A)) (afin exacto)
  %  'rotN' transporta normales          ''     sin handler (solo invalidar)
  switch M.ct
    case 3
      TBL = { 'length','len' };
    case 5
      TBL = { 'lengths','len' ; 'minlength','len' ; 'maxlength','len' ; ...
              'angles','inv' ; 'minangle','inv' ; 'maxangle','inv' ; ...
              'area','area' ; 'normal','rotN' ; ...
              'heights','len' ; 'minheight','len' ; 'maxheight','len' ; ...
              'inradius','len' ; 'circumradius','len' ; ...
              'aspectratio','inv' ; 'aspectfrobenius','inv' ; ...
              'edgeratio','inv' ; 'radiusratio','inv' ; 'relativesize','inv' };
    case 10
      TBL = { 'lengths','len' ; 'minlength','len' ; 'maxlength','len' ; ...
              'angles','inv' ; 'minangle','inv' ; 'maxangle','inv' ; ...
              'areas','area' ; ...
              'volume','adet' ; 'signedvolume','det' ; 'orientation','sgn' ; ...
              'inradius','len' ; 'circumradius','len' ; ...
              'edgeratio','inv' ; 'radiusratio','inv' ; 'aspectratio','inv' ; ...
              'relativesize','inv' ; ...
              'heights','len' ; 'minheight','len' ; 'maxheight','len' ; ...
              'dihedral','inv' ; 'mindihedral','inv' ; 'maxdihedral','inv' };
    otherwise
      error('meshQuality_as_cachedProps:celltype', ...
            'celltype %g no soportado (segmentos=3, triangulos=5, tets=10).', M.ct );
  end

  %subconjunto pedido (default: todas)
  if ~isempty( varargin )
    want = lower( cellfun( @char , varargin ,'uni',0) );
    [ ok , w ] = ismember( want , TBL(:,1) );
    if ~all( ok )
      error('meshQuality_as_cachedProps:metric', ...
            'metrica(s) desconocida(s) para celltype %g: %s (validas: %s).', ...
            M.ct , strjoin( want(~ok) , ', ' ) , strjoin( TBL(:,1).' , ', ' ) );
    end
    TBL = TBL( w ,:);
  end

  for i = 1:size( TBL ,1)
    name = TBL{i,1};
    cf   = @(m) meshQuality( toStruct( m ) , name );      %captura name por valor
    args = { 'changeCoords' , [] , 'changeConnectivity' , [] };
    switch TBL{i,2}
      case 'inv',  h = @(v,m,T) sameIfSimilarity( v , T , m.nsd );
      case 'len',  h = @(v,m,T) v * simScale( T , m.nsd );
      case 'area', h = @(v,m,T) v * simScale( T , m.nsd )^2;
      case 'adet', h = @(v,m,T) v * abs( det( linPart( T ) ) );
      case 'det',  h = @(v,m,T) v * det( linPart( T ) );
      case 'sgn',  h = @(v,m,T) v * sign( det( linPart( T ) ) );
      case 'rotN', h = @(v,m,T) v * normalRot( T , m.nsd ).';
      otherwise,   h = [];
    end
    if ~isempty( h ), args = [ args , { 'transform' , h } ]; end       %#ok<AGROW>
    M = M.defineCachedProp( name , cf , args{:} );
  end

end

%------------------------------------------------------------------ helpers
function A = linPart( T )
  %parte lineal 3x3 (tets viven en 3D)
  T = double( T );
  if size( T ,1) < 3 || size( T ,2) < 3
    error('meshQuality_as_cachedProps:T','T sin parte lineal 3x3.');
  end
  A = T( 1:3 , 1:3 );
end

function [ R , s ] = simParts( T , d )
  %descompone una SEMEJANZA en rotacion R (ortonormal, det=+-1) y escala s;
  %ERROR si T no es semejanza -> la cascada del replay recomputa
  T = double( T );
  if d == 2
    if size( T ,1) < 2 || size( T ,2) < 2, error('meshQuality_as_cachedProps:notSimilarity','T invalida.'); end
    A = T( 1:2 , 1:2 );
  else
    A = linPart( T );
    d = 3;
  end
  s = abs( det( A ) ) ^ ( 1/d );
  if s == 0, error('meshQuality_as_cachedProps:notSimilarity','T es singular.'); end
  R = A / s;
  if norm( R.' * R - eye( d ) , 'fro' ) > 1e-9
    error('meshQuality_as_cachedProps:notSimilarity','T no es una semejanza.');
  end
end

function s = simScale( T , d )
  [ ~ , s ] = simParts( T , d );
end

function v = sameIfSimilarity( v , T , d )
  simParts( T , d );      %solo valida (error si no es semejanza): v no cambia
end

function R = normalRot( T , d )
  %matriz que transporta NORMALES de cara: cross(Ra,Rb) = det(R) * R * cross(a,b)
  %-> el factor det(R) hace las reflexiones exactas
  R = simParts( T , d );
  R = R * sign( det( R ) );
end
