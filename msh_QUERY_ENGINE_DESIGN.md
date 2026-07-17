# Motor unificado de queries geométricas — Diseño

> Estado: **propuesta aprobada en conversación (2026-07-16), pendiente de
> implementación por fases**. Sustituye/absorbe: el BVH de esferas actual
> (`mshBVH`/`mshClosestElement`, fase 0 ya operativa), el BVH4 de
> `tools\IntersectSurfaceRay_mx.cpp`, y el motor legado
> `vtkClosestElement`/`distanceFrom` (locator de una sola malla, no
> serializable, recalculado ante cualquier cambio — el problema que originó
> todo esto).

## 0. Objetivos y no-objetivos

**Objetivos** (performance como objetivo explícito del proyecto):
- Un solo motor para **closest-element** (todas las celltypes: puntos,
  segmentos, triángulos, quads, tets, mixtas) y **ray-mesh** (triángulos),
  compartiendo build, layout de nodos, hojas, marco global y caché.
- Estructuras de búsqueda **como valores**: serializables, copiables,
  cacheables en `msh` (dos claves, COW), `save`/`load`-ables, válidas en
  `parfor`. Nada de estado global ni handles opacos.
- Transformaciones: semejanza → **O(1)** (plegado en el marco); deformación →
  refit O(n) manteniendo jerarquía; afín anisótropa → update conservador O(n)
  sin coordenadas, o refit exacto con ellas.
- Dos volúmenes envolventes **a elección del usuario** (`'sphere'` | `'aabb'`,
  con `'auto'` heurístico), mismo árbol, mismos kernels de hoja.

**No-objetivos** (por ahora): motion blur, instancing multi-nivel, GPU,
curvas/superficies paramétricas.

## 1. Decisiones tomadas (con su porqué)

### 1.1 Disciplina de distancias cuadradas (la pregunta del sqrt)

El recorrido **ya es** todo-cuadrático en el camino caliente; el `sqrt` de las
esferas no es evitable *exactamente*, pero sí amortizable:

- **Esfera**: la cota inferior es `|p−c| − r`. Al elevarla al cuadrado aparece
  el término cruzado `2·r·|p−c|`: la resta del radio vive en el espacio de
  distancias *lineales* (offset de Minkowski), no hay identidad algebraica que
  lo evite. El truco usado: comparar `|p−c|² > (r + best)²`, manteniendo
  `best = √best2` **solo cuando el mínimo mejora** (5–15 veces por query, no
  por nodo). Coste real: despreciable, pero existe.
- **AABB**: la cota `Σ max(0, lo−x, x−hi)²` es **cuadrática nativa** — se
  compara con `best2` directamente, cero `sqrt` en toda la query. Punto
  (pequeño) extra para AABB.

### 1.2 Esferas vs AABB según el tipo de malla (la intuición del plano)

Confirmada y adoptada como heurística de `'auto'`: en mallas **superficiales**,
a partir de cierto nivel del árbol los clusters son localmente planos y la
esfera desperdicia una dimensión entera (esfera de un disco de radio R ≈ bola
R³; el box lo abraza: R×R×h). En mallas **volumétricas** (tets) los clusters
son 3D-isótropos y la esfera es competitiva con test más barato y branchless.
Heurística `'auto'` por nodo raíz o global: razón de anisotropía/planitud del
PCA (análoga al umbral ≥2 del ROOT-PCA existente) + celltype (10 → sphere,
5/3 → aabb). Se calibra con benchmark en la fase 2.

Nota templada por datos propios: con warm-start + Dmax, la estanqueidad del
volumen mueve el end-to-end menos de lo que la teoría sugiere (penalización
×1.02–1.10 con árbol obsoleto tras warp violento). Se esperan 10–25% entre
volúmenes, no ×2 — otra razón para ofrecer ambos y medir, no dogmatizar.

### 1.3 Marco global de semejanza (`T_frame`)

Todo árbol lleva un `T_frame` (4×4, semejanza). Reglas:

| Evento | Acción | Coste |
|---|---|---|
| rígida / escala uniforme | `T_frame ← S∘T_frame`; nodos intactos | **O(1)** |
| deformación (coords) | refit de volúmenes en el marco actual | O(n) |
| afín anisótropa, sin coords | esferas: radios×σmax; AABB: 8 esquinas transformadas | O(n), conservador |
| afín anisótropa, con coords | refit exacto | O(n) |
| degradación acumulada | re-marco (PCA fresco) + rebuild, por heurística de calidad (Σ áreas/volúmenes vs build) | O(n log n) |

Query: `p' = T_frame⁻¹(p)` (9 flops/punto), resultados des-transformados a la
salida (`cp`, `d·s`, `Dmax/s` a la entrada). El miedo a "AABB rotado 45°"
desaparece: los boxes viven en el marco material y las rotaciones nunca los
tocan. **No se puede** plegar una afín no-similar en el marco para queries de
distancia (no preserva el orden de las distancias) — por eso la fila 3.

Nota: el refit bottom-up de AABBs (merge de hijos) es **exacto** (min/max
distribuye sobre uniones) y O(n); el de esferas por merge acumula holgura, por
eso el refit de esferas se hace por rango de elementos (como hoy).

### 1.4 Estructura como blob serializable (la caché)

El árbol empaquetado se devuelve a MATLAB como **un array plano tipado con
cabecera versionada** (magic, versión, volumen, contadores, offsets, T_frame).
El mex de query lo valida (magic + tamaños) y lo usa sin reconstruir nada.

- Vive en la caché de `msh` con el esquema de dos claves ya operativo
  (`'BVH'` geometría / `'BVH:hierarchy'` topología) + `T_frame`.
- `save`/`load` y `parfor` funcionan gratis (es un valor).
- Mata el LRU-por-fingerprint interno de `IntersectSurfaceRay_mx`: la caché de
  `msh` sabe *exactamente* cuándo el árbol es válido, sin hashear O(n) bytes
  por llamada.
- Lección miniball: cero estado persistente en los mex.
- Alineación: cargas SIMD no alineadas (`loadu`) — coste ~nulo en CPUs
  modernas, sin dependencia del alineamiento de arrays MATLAB.

### 1.5 Embree: no (por ahora), y el porqué

Se parece en espíritu (SAH, wide-BVH, SIMD, paquetes), pero delegarlo todo
sería un mal trade para ESTE proyecto:

1. **Embree es de rayos-sobre-triángulos.** El closest-point va por
   `rtcPointQuery` con callback por primitiva (pierde el batching SIMD) y no
   cubre segmentos, nubes, tets ni mixtas — que son la mitad del valor del
   motor. Nuestro 7.4 µs/pt con warm-start+Dmax no quedaría muy atrás, si
   queda atrás.
2. **Sus escenas no son serializables**: handles opacos por sesión →
   reaparece el problema vtk (rebuild por sesión y por worker, nada de
   save/load, nada de "estructura como valor"). Anti-tesis de §1.4.
3. Dependencia de despliegue real (DLLs, TBB, versionado) para un toolbox de
   laboratorio en Windows.
4. Donde Embree ganaría claro (throughput masivo de rayos en mallas enormes),
   el workload actual (paridad de meshIsInterior, siluetas, oclusión) no lo
   exige; el BVH4 propio ya mide bien ahí.

**Puerta abierta**: el API de `msh` queda backend-agnóstico (`backend='embree'`
posible algún día para rayos-triángulo puros, como acelerador opcional). Gate
de decisión: si aparece un workload de rayos donde el motor propio quede >×3
por detrás de forma sostenida.

## 2. Arquitectura

```
                    ┌──────────────────────────────────────────┐
                    │  msh cache: 'BVH' (blob) + 'BVH:hierarchy'│
                    └───────────────┬──────────────────────────┘
                                    │ blob (valor serializable, T_frame dentro)
        ┌───────────────────────────┼─────────────────────────────┐
        │ mshBVH (build/refit/     │ mshClosestElement           │ IntersectSurfaceRay
        │  transform, .m + _mx)     │  (_mx kernel CP)             │  (_mx kernel RAY)
        └───────────┬───────────────┴──────────────┬──────────────┴───────┐
                    ▼                              ▼                      ▼
             BUILD compartido               TRAVERSAL compartido     KERNELS de hoja
             - SAH binned (16 bins)         - nodos wide (4-ary)     - PreTri4 (v0,e1,e2):
             - hojas adaptativas              SoA float conservador     * Möller-Trumbore 4-wide (ray)
               (min 4 / max 16)             - test 4 hijos/pase SIMD    * Ericson (CP: ab=e1, ac=e2 ¡gratis!)
             - Morton para orden            - orden por distancia     - Seg/Tet/Point/Quad packs (solo CP)
             - PCA/frame en build           - CP: best2, warm-start,  - bc / uv por demanda
             - volumen: sphere|aabb|auto      Dmax, stale-pop
                                            - RAY: slab/esfera-ray,
                                              modos first|last|all|any
```

**Nodo wide (BVH4), SoA, bounds en float redondeado conservador** (robado de
IntersectSurfaceRay_mx; float = mitad de memoria, poda válida, tests exactos en
double en las hojas):
- `aabb4`: lox[4],hix[4],loy[4],hiy[4],loz[4],hiz[4] float + 4 hijos int32 = **128 B**.
- `sphere4`: cx[4],cy[4],cz[4],r[4] float + 4 hijos int32 = **80 B**.
Un pase SIMD testea los 4 hijos (distancia² al box / a la esfera), máscara +
orden por distancia → push. El coste por box se divide ~×4: mejora con mejor
ROI que cualquier apretón de volumen (por eso **no** volúmenes duales por nodo:
esfera∩AABB = 2 tests + 80→128+ B para un recorte que el AABB solo ya casi da).

**Hojas compartidas**: bloques contiguos por rango, 4-wide con carriles
degenerados de relleno. Triángulos en PreTri4 sirven a ambos kernels sin
duplicar memoria. Celltypes no-triángulo: packs propios, solo kernel CP (los
rayos filtran a hojas de triángulos en mallas mixtas).

## 3. Fases y criterios de aceptación

| Fase | Contenido | Aceptación |
|---|---|---|
| **P0** ✅ | motor esferas actual (mex, Dmax, refit, tests 11 bloques) | hecho — es la referencia viva |
| **P1** ✅ | formato blob v1 + `T_frame` (plegado de semejanzas; query pipeline des/transforma; 2D genuino; spot-check anti-obsolescencia con auto-rebuild) | ✅ resultados idénticos a P0 (14 bloques); fold medido 0.017 ms @52k tris (vs refit 5 ms) = O(1); save/load idéntico. Blob v1 = struct autocontenido: `X/Tri` en espacio de construcción + `frame` 4×4 + `nsd` (de columnas, NUNCA del contenido: malla 3D con z≡0 sigue 3D y puede deformarse fuera de plano) + `tet4` + `version`. 2D: transforms 2×2/2×3/3×3-homogéneo elevados con z-scale = escala 2D (mantiene la semejanza 3D del padding). Afín anisótropa → bake (frame→geometría, radios×σmax, frame=eye) |
| **P2** ✅ | build SAH binned (mex `mshBVH_mx`, agnóstico a la geometría: recibe esferas/cajas/centroides de elementos) + colapso BVH2→BVH4 + nodos wide (`sphere4`/`aabb4`, floats conservadores, aritmética double) + opción `volume` + recorrido ancho en ambos backends + refit/bake v2 | ✅ == v1 en toda la suite (sphere+aabb × matlab+mex). **Tabla medida (monohilo): 52k tris → aabb 7.1 ≈ v1 6.9 µs/pt, sphere 10.4; 26k tets → aabb 2.4 (×3.5 sobre v1 8.5), sphere 8.7.** Hallazgo: **AABB gana también en tets** (la esfera de un tet es holgada; el test de caja es sqrt-free) → `'auto'` = aabb siempre; la hipótesis "volumétrica→esferas" quedó refutada por datos. Partición única SAH para ambos volúmenes (confirmada la decisión). Notas honestas: (a) el objetivo ≤5 µs/pt se cumple en tets (2.4), NO aún en superficies — el cuello ya son los tests exactos de hoja, la palanca es el kernel 4-wide de P3 (PreTri4); (b) hojas SAH por defecto [2 16], knob `mshBVH(M,[minL maxL],vol)` ([4 32] da ~6% en superficies); (c) el refit v2 en MATLAB (15 ms @2.6k) es hoy más LENTO que el rebuild v2 vía mex (3 ms) — su valor actual es la persistencia de jerarquía para la caché; refit en C va a P4 |
| **P3** ✅ | hojas empaquetadas EN EL BUILD (`pkV/pkS/pkT/pkE` en el blob, en orden perm: muere el gather por llamada; refit/bake las mantienen) + kernel Ericson **4-wide AVX branchless** (mismo orden de regiones que el escalar → paridad exacta; lanes degeneradas → redo escalar; detección CPUID + fallback) + **`mshIntersectRay`(+`_mx`): rayos sobre el MISMO blob** (Möller-Trumbore double con la tolerancia inclusiva 1e-9 del oráculo, modos `first|last|all|any` con su semántica exacta de línea no acotada y ventana de oclusión; sólo celdas triángulo, mixtas saltan el resto; rayos plegados por el marco — t invariante) | ✅ 25 bloques verdes; rayos == `IntersectSurfaceRay` (oráculo) en todos los modos. **Hallazgo clave de perfilado: el coste de closest-point lo domina la DISTANCIA del punto a la malla** — casquete tangente ~O(√n) — no los kernels: 52k tris → lejos 7.4 µs/pt, **cerca de superficie 0.80 µs/pt** (el workload real: registro/stick/offset); lejos ⇒ `Dmax` (0.22). Por eso el kernel 4-wide no movió el benchmark lejano (correcto pero no era el cuello). Rayos: 1.5 vs 0.4 µs/ray (`first`), 0.63 vs 0.33 (`any`) — el mex especializado gana ×2.5-3.7 (SIMD BVH4+MT4+poda ordenada): margen para P4; el valor entregado es la ARQUITECTURA (un solo blob cacheado/transformable por malla para ambas queries, sin LRU-fingerprint) |
| **P4** ✅ | recorrido de rayos ORDENADO por t con poda por mejor-t y culling de stack (correcto con t negativos: la clave es la cota del intervalo), `invd` hoisteado por rayo, slab-test de los 4 slots en un pase AVX (double sobre bounds float — intervalos idénticos al escalar), MT4 para `first/last/all` (escalar early-exit para `any`), **refit en C** (modo refit de `mshBVH_mx`: 2.0 ms vs rebuild 3.4 — antes 23 ms en MATLAB) | ✅ 25 bloques verdes. **Rayos monohilo: paridad con `IntersectSurfaceRay_mx`** (first 0.44-0.46 vs 0.39-0.40 µs/ray; any 0.38 vs 0.33). **Multihilo (el mex del oráculo es monohilo por diseño; el nuestro paraleliza por rayos): first 0.162 µs/ray con 16 hilos físicos = ×2.5 MEJOR que el especializado**; closest-point lejano 8.2 → 0.77 µs/pt (×10.6 con 32). Ojo operativo: primer parallel-region por proceso paga la creación del pool; hyperthreads (32) rinden peor que cores físicos (16) — el wrapper pasa `maxNumCompThreads` (= físicos) que es el óptimo. Pendiente menor trasladado a backlog: walk de point-location d==0 para tets, heurística de re-marco, gate Embree (con paridad+MT lograda, el gate se aleja) |

Regla transversal: la implementación MATLAB pura se conserva siempre como
oráculo y fallback (`backend='matlab'`), y cada fase corre la suite completa
(`test_mshBVH`, `test_mshClosestElement_tets` + los que traiga cada fase).

## 4. Riesgos y mitigaciones

- **Máquina de desarrollo inestable bajo carga** (access violations en MATLAB
  puro, ~1/4 runs pesados): protocolo diary/marcadores + reintento; los
  benchmarks se repiten ×3 y se toma mediana.
- **Bounds float conservadores**: redondeo SIEMPRE hacia fuera
  (nextafter/±ulp); test de estrés con geometrías en los límites de precisión.
- **Celltypes mixtas en SIMD**: carriles degenerados de relleno por tipo
  (como PreTri4); nunca mezclar tipos en un mismo bloque de 4.
- **Ambición vs regresión**: cada fase es shippeable por sí sola; P0 no se
  toca hasta que P2 lo iguale en tests.
