# msh

Clase MATLAB para gestionar mallas (polilíneas, superficies de triángulos,
volúmenes de tetraedros y mallas mixtas) con **derivados cacheados
perezosamente** que se auto-invalidan — o se auto-actualizan barato — al
editar la malla, y un **motor propio de queries geométricas** (elemento más
cercano e intersección rayo-malla) acelerado por BVH con MEX optimizados.

```matlab
M = msh( V , F );                 % V: nV x 3 double, F: nF x k int32
[e,cp,d,bc] = M.closestElement( P );   % construye el BVH una vez, lo cachea
M2 = M.transform( T );                 % semejanza: el BVH se PLIEGA en O(1)
M.V = V + delta;                       % deformación: refit perezoso, no rebuild
q = M.aspectratio_;                    % cachedProps definibles por el usuario
```

## Estructura del repo

| Directorio | Contenido |
|---|---|
| `@msh\` | la clase (`msh.m`) + `private\` (funciones legadas en transición a msh-nativas) |
| `BVH\` | motor de queries, **independiente de la clase**: `BVH` (builder SAH BVH4, 6 tipos de volumen), `bvhClosestElement`, `bvhIntersectRay`, `plotBVH` (visor), MEXes `*_mx.cpp` y tests |
| `MESH\` | toolbox legado de structs (`.xyz`/`.tri`), en extinción progresiva |
| `tools\` | utilidades (`transform`, `IntersectSurfaceRay`, ...) |
| raíz | `cacheHandle.m`/`cacheView.m` (maquinaria de cache genérica), tutoriales y documentos de diseño |

## Conceptos clave

- **cachedProps**: registro *nombre → { computeFcn, eventos }* por malla.
  Acceso `M.BVH_` o `M.cached.BVH`; las ediciones disparan eventos
  (`changeCoords`, `changeConnectivity`, `transform(T)`, ...) y cada entrada
  sobrevive, cae, o queda pendiente de un **replay perezoso** (p. ej. el BVH
  se refita en O(n) o se pliega en O(1) en vez de reconstruirse). Definibles
  por el usuario: `M = M.defineCachedProp(nombre, @(m)..., evento, handler)`.
- **Semántica de valor + cache compartida (copy-on-write)**: las copias
  comparten los valores cacheados hasta que una edita; las hermanas nunca
  pierden su cache.
- **Tipos obligatorios**: `M.V` siempre `double`, `M.F` siempre `int32`
  (0-padded; 4 nodos no-cero = tetraedro).
- **Motor BVH**: blob autocontenido y serializable, marco de semejanza global,
  refit persistente, batería de volúmenes (aabb/sphere/obb/kdop/rss/lss).
  Single-thread a la par o por delante de los MEX especializados históricos;
  multihilo (OpenMP, sigue `maxNumCompThreads`) varias veces más rápido.

## Instalación

Requiere MATLAB R2021b+ (probado en R2022a) y un compilador C++ (MSVC) para
los MEX:

```matlab
addpath C:\repos\msh          % @msh + cacheHandle + cacheView
addpath C:\repos\msh\BVH      % motor de queries
addpath C:\repos\msh\MESH     % toolbox legado
addpath C:\repos\msh\tools

cd C:\repos\msh\BVH           % compilar los MEX (una sola vez)
mex BVH_mx.cpp
mex COMPFLAGS="$COMPFLAGS /openmp" bvhClosestElement_mx.cpp
mex COMPFLAGS="$COMPFLAGS /openmp" bvhIntersectRay_mx.cpp
```

Los binarios (`*.mexw64`, `*.dll`) se versionan con **Git LFS**.

## Tests y documentación

```matlab
cd MESH; test_msh                      % la clase (construcción, eventos, replay, proxy...)
cd ..\BVH
test_BVH                               % closest-element, 6 volúmenes, refit, fold
test_bvhIntersectRay                   % rayos vs oráculo
test_bvhClosestElement_tets            % semántica tsearchn en tets
```

- `msh_CLASS_TUTORIAL.md` — tutorial completo de la clase (18 secciones).
- `BVH_TUTORIAL.md` — tutorial del motor de queries.
- `msh_DESIGN.md`, `msh_QUERY_ENGINE_DESIGN.md` — documentos de diseño
  (históricos; conservan la nomenclatura antigua `msh*` del motor).

## Estado

En desarrollo activo. La clase y el motor están operativos y testeados; el
toolbox legado `MESH\` se irá absorbiendo en `@msh\private` hasta desaparecer.
