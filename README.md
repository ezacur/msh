# msh

Clase MATLAB para gestionar mallas (polilíneas, superficies de triángulos,
volúmenes de tetraedros y mallas mixtas) con **derivados cacheados
perezosamente** que se auto-invalidan — o se auto-actualizan barato — al
editar la malla, y un **motor propio de queries geométricas** (elemento más
cercano, intersección rayo-malla e intersección malla-malla) acelerado por BVH
con MEX optimizados.

```matlab
M = msh( V , F );                       % V: nV x 3 double, F: nF x k int32
[e,cp,d,bc] = M.ClosestElement( P );    % construye el BVH una vez, lo cachea
M2 = M.Transform( T );                  % semejanza: el BVH se PLIEGA en O(1)
M.V = V + delta;                        % deformación: refit perezoso, no rebuild
M.CP                                    % tabla de CPs: valor/estado y eventos
```

## La convención del API: tres espacios de nombres por CASO

| Caso | Qué es | Ejemplos |
|---|---|---|
| `MAYÚSCULAS` | propiedades | `M.V`, `M.F`, `M.VIZ`, `M.INFO`, `M.DEBUG`, `M.CP` |
| `Capitalized` | métodos | `M.Transform`, `M.ClosestElement`, `M.DefineCP`, `M.Plot` |
| `minúscula` | CPs (cached props) | `M.bvh`, `M.boundary`, `M.esup` |

`DefineCP` **exige** la minúscula inicial. Excepciones documentadas: los
contadores `M.nsd`/`M.nV`/`M.nF`/`M.ct`, los métodos del protocolo MATLAB
(`subsref`, `disp`, `loadobj`) y los alias transicionales de solo lectura
`xyz`/`tri`/`celltype`.

Ojo: las formas funcionales **no** despachan a la clase. `plot(M)` o
`transform(M,T)` caen a las funciones del path — usa `M.Plot()` y
`M.Transform(T)`.

## Estructura del repo

| Directorio | Contenido |
|---|---|
| `@msh\` | la clase (`msh.m`) + `private\` (funciones legadas en transición a msh-nativas) |
| `BVH\` | motor de queries, **independiente de la clase**: `BVH` (builder SAH BVH4, 6 tipos de volumen), `bvhClosestElement`, `bvhIntersectRay`, `bvhIntersectMesh`, `fanClosestElement` (aproximado), `plotBVH` (visor), MEXes `*_mx.cpp`, tests y benchmarks |
| `MESH\` | toolbox legado de structs (`.xyz`/`.tri`), en extinción progresiva |
| `tools\` | utilidades (`transform`, `IntersectSurfaceRay`, `miniball`, `parseargs`, ...) |
| `uiTools\` | interacción sobre figuras (`silhouette`, `silhouetteByShading`, `headlight`, `backFaceCulling`) |
| `extras\` | `@I3D` (imágenes 3D); lo necesitan `MeshOffset`, `MeshSmoothZip`, `MeshZipSmooth` |
| raíz | `cacheHandle.m`/`cacheProxy.m`/`cachedOwner.m` (patrón de cache genérico, independiente del dominio: hereda de `cachedOwner` y tu clase es puro dominio — ver los cuatro `*_TUTORIAL.html` listados abajo), tutoriales y documentos de diseño |
| `dev\` | ejemplos del patrón de cache (`cart`, `route`, `serie`, `facto` —todo-incremental—, `ivp` —EDO con solución cacheada incrementalmente—, `collatz` —cadenas 3n+1: el atajo gana 100× y la primera CP compuesta—, `linsys` —`A\b` con la factorización cacheada: el análogo del `bvh`— + demos), el benchmark `bench_incremental.m` y su batería de verificación (`dev\tests\test_cache.m` — correr antes y después de tocar cualquier `cache*.m`) |

## Conceptos clave

- **CPs (cached props)**: registro *nombre → { computeFcn , eventos }* por malla.

  ```matlab
  M.bvh                      % LEE (perezoso: HIT / replay / computa y guarda)
  M.bvh_                     % RECALCULA a la fuerza
  M.CP                       % tabla de definiciones y estados
  M.CP.bvh                   % el valor (igual que M.bvh)
  M.CP.bvh.delete            % borra el VALOR (la definición queda)
  M = M.CP.bvh.removeCP      % borra definición y valor (== M.RemoveCP('bvh'))
  M = M.CP.bvh.set( x )      % siembra un valor a mano (aislado, COW)
  M = M.DefineCP( 'foo' , @(m) ... , evento , handler|[] , ... )
  M = M.RemoveCP( 'foo' )
  ```

  Y dos herramientas de diagnóstico en `cachedOwner` (forma de función = informe
  impreso; forma de punto = los datos):

  ```matlab
  VerifyCPs( o )     % ¿la caché sirve lo mismo que el compute? (handler que miente,
                     % evento no declarado). No toca la caché.
  ProfileCPs( o )    % ¿compensa cachear? mide trabajo / MISS / HIT por CP y separa
                     % "no debería ser CP" de "no conviene guardarla".
  o.ADVISE = true    % el consejero pasivo: avisa (una vez) de la CP que no se
                     % sostiene con los datos de ahora. Apagado por defecto.
  InfoCP( o , 'x' )  % la ficha de UNA CP: estado, eventos y su política, de quién
  o.CP.x.info        % sale, quién sale de ella, historial (recálculos/HIT/replay/
                     % caídas) y coste. `S = o.InfoCP('x')` da el struct.
  ```

  Las ediciones disparan **eventos** de específico a general — `M.V = ...` lanza
  `[changeNodeCount] [changeDim] changeCoords`, `M.F = ...` lanza
  `[changeFaceCount] changeConnectivity`, y `Transform()` lanza
  `transform(T) + changeCoords`. Para cada CP: evento no declarado =
  **sobrevive**; declarado con `[]` = **invalida**; declarado con handler =
  queda **pendiente** y el handler lo actualiza en el próximo acceso
  (`transform` es incremental `@(v,m,T)`; el resto son absolutos `@(v,m)`). Si
  un handler falla se degrada en cascada: transform → sync absoluto →
  recompute, nunca un valor silenciosamente incorrecto.

  De fábrica nacen 9 CPs, todas sobreescribibles: `bvh` (transform → plegado
  O(1), changeCoords → refit O(n)), `boundary`, `edges`, `esup`, `psup`,
  `esue`, `bbox`, `surfCent` (`[area, centroide]`) y `triNormals` (transform →
  rotación O(n) de las filas).

- **Semántica de valor + cache compartida (copy-on-write)**: las copias
  comparten los valores cacheados hasta que una edita; las hermanas nunca
  pierden su cache. El handle (`cacheHandle`) es privado y `Transient`, así que
  MATLAB lo destruye por conteo de referencias en cuanto la malla que lo posee
  muere: no hay ni hace falta recolector propio.

- **Tipos obligatorios**: `M.V` siempre `double`, `M.F` siempre `int32`
  (0-padded; 4 nodos no-cero = tetraedro). Los atributos por nodo/cara
  (`M.AddField('xyzFOO',v)`, `M.GetField`) conservan el tipo que les des.

- **Motor BVH**: blob autocontenido y serializable, marco de semejanza global,
  refit persistente, batería de volúmenes (`aabb` por defecto, `sphere`, `obb`,
  `kdop`, `rss`, `lss`). Single-thread a la par o por delante de los MEX
  especializados históricos; multihilo (OpenMP, sigue `maxNumCompThreads`)
  varias veces más rápido.

- **`M.DEBUG = true`** narra por consola la cadena de procesos: `HIT`/`MISS`/
  `RPLAY`/`RECMP` de cada CP con tiempos, qué cae / queda pendiente /
  sobrevive en cada evento, `Transform` y queries. También
  `msh(V,F,'DEBUG',true)`.

CPs definidas por el usuario, ejemplo real listo para usar:

```matlab
M = meshQuality_as_cachedProps( M );   % registra una CP por métrica de calidad
q = M.aspectratio;                     % lee (perezoso)
q = M.aspectratio_;                    % recalcula a la fuerza
```

## Instalación

Requiere MATLAB **R2021b+** (desarrollado y probado en R2022a) y un compilador
C++ (MSVC) para los MEX.

```matlab
addpath C:\repos\msh          % @msh + cacheHandle + cacheProxy + cachedOwner
addpath C:\repos\msh\BVH      % motor de queries
addpath C:\repos\msh\MESH     % toolbox legado
addpath C:\repos\msh\tools
addpath C:\repos\msh\uiTools
addpath C:\repos\msh\extras   % @I3D (lo usan MeshOffset / MeshSmoothZip / MeshZipSmooth)
```

Compilar los MEX, una sola vez. `/openmp` habilita el multihilo y `-lut` el
soporte de Ctrl-C (`utIsInterruptPending`); los tres MEX sin esos flags son
mono-hilo por diseño:

```matlab
cd C:\repos\msh\BVH
mex BVH_mx.cpp                                                     % builder del blob
mex COMPFLAGS="$COMPFLAGS /openmp" -lut bvhClosestElement_mx.cpp   % elemento más cercano
mex COMPFLAGS="$COMPFLAGS /openmp" -lut bvhIntersectRay_mx.cpp     % rayos
mex COMPFLAGS="$COMPFLAGS /openmp" -lut fanClosestElement_mx.cpp   % cercano aproximado
mex bvhTriPairs_mx.cpp                                             % pares candidatos malla-malla
mex triTriPairs_mx.cpp                                             % intersección tri-tri exacta

cd ..\tools
mex IntersectSurfaceRay_mx.cpp                                     % rayos, camino legado
mex miniball_mx.cpp
```

### Git LFS

Los binarios (`*.mexw64`, `*.dll`) se versionan con **Git LFS**, así que hace
falta `git-lfs` instalado y en el PATH (`git lfs install`). Sin él los hooks
del repo salen con código 2 — **`git push` queda bloqueado** — y cualquier
binario que recommitees se guardaría en claro en vez de como puntero LFS.

### Dependencia externa: `mTools`

`MESH\` **no es autocontenido**: 58 de sus 129 funciones llaman a helpers del
toolbox `mTools` del autor, que vive fuera de este repo — `getPlane`,
`distance2Plane`, `ndmat`, `maxnorm`, `minv`, `argmin`, `hplot3d`,
`InterpolatingSplines`, `MatchPoints`, `FarthestPointSampling`, los
envoltorios VTK y la E/S `read_*`/`write_*`, entre otros. Sin `mTools` en el
path:

- `@msh\`, `@msh\private\`, `BVH\` y `uiTools\` funcionan **enteros**: el
  núcleo no toca nada de eso;
- `tools\transform.m` funciona con una matriz `T`, pero **falla** en cuanto la
  malla lleva un campo `NORMALS` (necesita `maketransform` y `cbrt`) — y con
  ello falla `M.Transform` sobre mallas con normales;
- de `MESH\` funciona lo que no dependa de esos helpers.

## Tests y documentación

```matlab
cd C:\repos\msh\MESH
test_msh                        % la clase: construcción, tipos, eventos, replay, COW, proxy
cd ..\BVH
test_BVH                        % closest-element, 6 volúmenes, refit, fold, 2-D, degenerados
test_bvhIntersectRay            % rayos vs oráculo
test_bvhClosestElement_tets     % semántica tsearchn en tets
test_bvhIntersectMesh           % intersección malla-malla vs fuerza bruta
test_fanClosestElement          % localizador aproximado: cota, formas de llamada, guards
```

Esas 6 suites pasan en R2022a. Los `test_meshNormals*` y `test_MeshSubdivide*`
de `MESH\` necesitan `mTools` (ver arriba); `test_meshNormals_fix` necesita
además `MESH\old_meshNormals.mat`, que no está en el repo.

Benchmarks en `BVH\`: `bench_BVH`, `bench_fanClosestElement`,
`bench_seedCeiling`, `bench_vertexSeed`.

- `msh_CLASS_TUTORIAL.md` — tutorial completo de la clase (secciones 0-18).
- `BVH_TUTORIAL.md` — tutorial del motor de queries.
- `msh_DESIGN.md`, `msh_QUERY_ENGINE_DESIGN.md` — documentos de diseño
  (históricos: conservan la nomenclatura antigua `msh*` del motor y la API
  anterior a la convención por caso).

El patrón de cache tiene los suyos, uno por pregunta (todos en la raíz):

| Tutorial | Responde a |
|---|---|
| `cache_TUTORIAL.html` | **cómo funciona** por dentro: estados, replay, COW, el porqué de cada pieza (ejemplo `cart`) |
| `route_TUTORIAL.html` | **cómo se usa**, línea a línea, en una clase nueva (`route`) — y qué cubre la batería |
| `serie_TUTORIAL.html` | **cómo se decide** la política: conservador primero, atajos después (`serie`) |
| `incremental_TUTORIAL.html` | **cuándo el atajo paga** — con los tres casos medidos: pierde (`facto`), empata (`ivp`), gana 100× (`collatz`) |

```matlab
cd C:\repos\msh\dev
cartDemo, routeDemo, serieDemo                       % las 7 sesiones narradas
factoDemo, ivpDemo, collatzDemo, linsysDemo
bench_incremental                                    % ¿ahorra el atajo? (mide)
cd tests
test_cache                                           % la batería del patrón
```

## Estado

En desarrollo activo. La clase y el motor están operativos y testeados; el
toolbox legado `MESH\` se irá absorbiendo en `@msh\private` hasta desaparecer.
