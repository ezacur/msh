# Clase `msh` — Diseño

> Estado: **diseño + esqueleto para revisión**. Nada implementado todavía
> (los cuerpos de los métodos en `@msh/msh.m` son *stubs* con `TODO`).
> La única pieza escrita "de verdad" es `mshCache.m`, porque es el mecanismo
> central y conviene verlo concreto antes de decidir.

## 1. Objetivo

Encapsular en una clase MATLAB el modelo de malla que hoy vive como `struct`
suelto (con las ~126 funciones `mesh*`/`Mesh*` de `MESH\`), añadiendo:

- **Renombrado** de los atributos principales: `.xyz → .nodes`, `.tri → .faces`.
- **Caché perezosa** de estructuras derivadas (EsuP, PsuP, EsuE, edges,
  normales, bbox…): se calculan la primera vez que se piden, se guardan, y
  quedan disponibles para llamadas posteriores.
- **Auto-invalidación consciente de dependencias**: al editar la malla, solo se
  tira lo que dejó de ser válido; lo que no depende del cambio se conserva.
- **Semántica por valor** (como cualquier variable MATLAB), pero con la caché
  compartida por referencia y con *copy-on-write*, de modo que instancias
  hermanas nunca se pisan la caché.

## 2. Modelo de datos

| Atributo (nuevo) | Antiguo | Contenido | Forma |
|---|---|---|---|
| `.nodes` | `.xyz` | coordenadas de los nodos | `nV × 2` o `nV × 3` (en 2D, Z ≡ 0) |
| `.faces` | `.tri` | conectividad | `nF × k`, `k∈{1,2,3,4}`; mixtas → rellenar con `0` |
| campos por nodo | `.xyzFIELD*` | escalares/vectores por vértice | `nV × …` |
| campos por cara | `.triATT*` | escalares/vectores por celda | `nF × …` |
| `.celltype` | `.celltype` | código VTK (der.) | escalar o `nF × 1` |

**Tipo de celda** (VTK, igual que hoy vía `meshCelltype`): `k` nodos/cara →
1 = punto, 2 = segmento(3), 3 = triángulo(5), 4 = quad(9, 2D) / tetra(10, 3D).
Por eso **se conserva el nº de columnas nativo de `.nodes`** (2 vs 3): es lo que
desambigua quad vs tetra. El padding a 3D (Z=0) se hace *on demand* donde haga
falta, no en el almacenamiento.

**Atributos especiales** (subconjunto conocido, con tratamiento propio):

| Atributo | Tipo | Notas |
|---|---|---|
| `.triNORMALS`, `.xyzNORMALS` | campo geométrico | norma euclídea = 1. Derivable (auto) **o** fijable por el usuario. |
| `.xyzUV`, `.xyzUV2`, `.xyzUV3` | campo por nodo | coords de textura. |
| `.texture` | dato | bitmap `H×W×3` asociado a la malla. |

**Presentación / metadatos** (datos puros, sin efecto en la caché de derivados;
se serializan; su edición no dispara invalidación):

| Atributo | Tipo | Notas |
|---|---|---|
| `.vizProp` | struct | preferencias de dibujo: `FaceColor`, `EdgeColor`, `EdgeAlpha`, `FaceAlpha`, `LineWidth`… Normalmente vacío. |
| `.info` | struct | metadatos libres: `info.description`, `info.filename`, … |

`plot(M)` fusiona opciones con precedencia **defaults de `plotMESH` < `M.vizProp`
< args explícitos**. Se apoya en que `plotMESH` ya lee opciones desde los campos
del struct y en que el `patch(... defOPTS, mOPTS, varargin)` deja **ganar al
último**; así `plot(M,'FaceColor','b')` pisa a `M.vizProp.FaceColor`. `vizProp` e
`info` NO se vuelcan en `toStruct()` (las 126 funciones no las entienden); `plot`
inyecta `vizProp` por su cuenta.

Las normales son un **híbrido**: si el usuario las fija, se respetan y se
serializan como dato; si no, se calculan perezosamente y se cachean. En ambos
casos participan en la invalidación (ver §5): editar geometría/conectividad las
vacía y se regeneran cuando se vuelven a pedir. Esto arregla el *footgun* actual
(`meshNormals` hoy confía a ciegas en `M.triNORMALS` y avisa "rebuild after
editing").

## 3. Semántica: value class + caché por referencia con copy-on-write

`msh` es una **value class** (copia al asignar, como los `struct` de hoy → no
rompe el estilo `M2 = op(M1)` de las 126 funciones). La caché vive en un objeto
**handle** privado (`mshCache`), es decir una *referencia* — el "puntero" del
enunciado. La regla de oro:

> **Leer** la caché (rellenar un derivado que faltaba) escribe en el handle
> compartido → beneficia a todas las copias que comparten esa referencia (mismos
> datos ⇒ mismos derivados). **Editar** la malla dispara *copy-on-write*: la
> instancia editada recibe un handle de caché **nuevo y propio** con solo los
> derivados que sobreviven al cambio; las instancias hermanas siguen con el suyo.

### 3.1 Traza del ejemplo de aceptación

```matlab
M1 = msh( nodes , faces );                 % cache H0 (vacía)
M2 = metodo_que_usa_EsuP( M1 );            % M2 comparte H0; al pedir EsuP se
                                           % rellena H0 -> M1.EsuP y M2.EsuP existen
M3 = rotate( M1 );                         % edita coords -> COW: M3 recibe H1
                                           % (copia de H0 sin lo geométrico); EsuP
                                           % depende solo de la topología -> SOBREVIVE
M1.faces(1:10,:) = [];                     % edita conectividad -> COW: M1 recibe H2
                                           % SIN EsuP (invalidada). M3 sigue con H1,
                                           % así que M3.EsuP intacta.
```

| Paso | `M1.cache` | `M2.cache` | `M3.cache` | `EsuP` |
|---|---|---|---|---|
| `msh(...)` | H0∅ | — | — | — |
| `metodo(M1)` | H0{EsuP} | H0{EsuP} | — | compartida en H0 |
| `rotate(M1)` | H0{EsuP} | H0{EsuP} | **H1{EsuP}** | H1 = copia de survivors |
| `M1.faces(…)=[]` | **H2{∅}** | H0{EsuP} | H1{EsuP} | M1 invalida; M3 conserva |

Esto es exactamente el comportamiento pedido. La clave es que **el `set` que
invalida crea un handle nuevo** (no muta el compartido).

### 3.2 Por qué no basta con un `struct` normal

En una value class, un *getter* recibe una **copia** del objeto; si memoizara en
una propiedad normal, la escritura se perdería. Para que `M.EsuP` calcule **y
guarde** de forma transparente hace falta que el almacén sea un **handle**
(la referencia persiste a través de la copia del getter). De ahí `mshCache`.

## 4. Derivados cacheados y sus dependencias

Cada derivado declara de qué depende. Al editar se calculan los *tags tocados* y
se invalida todo derivado cuyo conjunto de dependencias los interseque.

**Tags:** `Connectivity` (topología de `.faces`), `NodeCoords` (valores de
`.nodes`), `NodeCount`, `FaceCount`.

| Derivado | Delegación | Dependencias | ¿`rotate` lo invalida? |
|---|---|---|---|
| `EsuP` (elems por punto) | `meshEsuP` | Connectivity | **No** |
| `PsuP` (puntos por punto) | `meshPsuP` | Connectivity | **No** |
| `EsuE` (elems por elem) | `meshEsuE` | Connectivity | **No** |
| `edges` (aristas únicas) | `meshEdges` | Connectivity | **No** |
| `celltype` | `meshCelltype` | Connectivity | **No** |
| `boundary` | `meshBoundary*` | Connectivity | **No** |
| `triNORMALS`/`xyzNORMALS` | `meshNormals` | NodeCoords + Connectivity | **Sí** |
| `bbox` | `meshBB` | NodeCoords | **Sí** |
| `surface`/`centroid` | `meshSurface` | NodeCoords + Connectivity | **Sí** |
| longitudes de arista | `meshEdges` | NodeCoords + Connectivity | **Sí** |
| `quality` | `meshQuality` | NodeCoords + Connectivity | **Sí** |
| `BVH` (closest-element) | `mshBVH` | NodeCoords + NodeCount + Connectivity + FaceCount | especial, ver §4.1 |

### 4.1 BVH: actualizar en vez de invalidar

El BVH de **esferas** (`MESH\mshBVH.m` + `MESH\mshClosestElement.m`, reemplazo
in-house del motor `vtkClosestElement`, mismas salidas `[e,cp,d,bc]`) es el
derivado más caro de reconstruir y el único con una vía barata de *seguir* a la
malla: ante una transformación afín, `mshBVH(B,T)` actualiza centros y radios
en O(n) sin reconstruir (exacto para rígidas + escala uniforme; conservador pero
correcto para afines generales — radios × mayor valor singular).

Por eso `transform()` no lo invalida: rescata el BVH cacheado **antes** del
`set.nodes` (que dispara el COW normal), y re-siembra la versión actualizada en
la caché nueva.

**Jerarquía persistente (deformaciones).** Para ediciones no-afines de las
coordenadas (deformación de la malla) tampoco hace falta reconstruir: la
*jerarquía* del árbol (perm/child/range) solo depende de la conectividad,
mientras que las *esferas* dependen además de las coordenadas. Se cachean como
dos claves con dependencias distintas:

| clave | contenido | deps |
|---|---|---|
| `BVH` | árbol + esferas al día | NodeCoords + NodeCount + Connectivity + FaceCount |
| `BVH:hierarchy` | árbol reutilizable | Connectivity + FaceCount |

Al deformar cae `BVH` pero sobrevive `BVH:hierarchy`; el siguiente acceso hace
**refit** — `mshBVH(B, M)`: conserva el árbol y recalcula solo las esferas en
O(n), cada nodo desde su propio rango de elementos (cotas tan ajustadas como
las de una construcción fresca). Bajo deformación grande el árbol pierde
*calidad de partición* (los splits se eligieron en las coordenadas viejas),
es decir velocidad — nunca corrección: como las hojas hacen tests exactos, los
resultados son idénticos a los de un BVH reconstruido (asertado en el test).
Solo cambiar la conectividad tira las dos claves. TODO opcional: heurística de
re-build cuando la degradación acumulada supere un factor (p.ej. razón entre la
suma de radios actual y la de construcción).

A diferencia del locator persistente de VTK (una malla a la vez, con el patrón
`vtkClosestElement(M)` … `vtkClosestElement([],[])` + `onCleanup`), cada `msh`
lleva su BVH en su caché: tantas mallas simultáneas como se quiera, sin estado
global.

**Backend compilado.** `MESH\mshClosestElement_mx.cpp` (MSVC, OpenMP opcional)
implementa el recorrido con: nodos empaquetados de 48 B, bloques de elementos
de hoja contiguos en orden `perm` (vértices pre-gathered), poda sin `sqrt` vía
`(r+best)²`, *culling* de entradas obsoletas del stack al desapilar, hijo
cercano primero, warm-start con el ganador del punto anterior, lote de queries
en orden Morton (un chunk contiguo por hilo), stack fijo iterativo y entradas
100% validadas (un `B` corrupto da error, nunca access violation).
`mshClosestElement` lo usa automáticamente (`backend` `'auto'|'matlab'|'mex'`);
los hilos siguen `maxNumCompThreads`. El `leafSize` por defecto de `mshBVH` es
adaptativo: 16 con el MEX, 256 en MATLAB puro.

**Radio de búsqueda `Dmax`** (`mshClosestElement(M,P,B,Dmax)`, ambos backends):
siembra la cota best-so-far con `Dmax` en vez de `∞`, de modo que todo lo que
quede más lejos se poda desde la raíz — un punto fuera de radio cuesta ~una
visita de nodo (medido: 0.20 µs/pt vs 5.0 µs/pt de la búsqueda completa en nube
lejana, ×25). Sin elemento dentro de `Dmax` (estricto: `d < Dmax`) devuelve
`e = d = cp = bc = NaN`. Convención unificada: cualquier "sin respuesta"
(fuera de radio, punto no finito) es NaN, nunca 0.

Rendimiento medido (monohilo, 52k tris, 20k puntos): MATLAB puro 244 µs/pt
(×33 vs fuerza bruta), **MEX 7.4 µs/pt (×33 sobre el MATLAB, ×1000 sobre
fuerza bruta)**, con concordancia **bit-exacta** entre backends
(`max|Δd| = 0` en la suite). Construcción despreciable (2–20 ms).

Clasificación de una edición:

- `set.nodes` con mismas filas, distintos valores → `{NodeCoords}`.
  Con distinto nº de filas → `{NodeCoords, NodeCount}`.
- `set.faces` → `{Connectivity, FaceCount}`.
- `M.nodes(i,:) = …` / `M.faces(i,:) = []` → MATLAB llama al `set` completo, así
  que caen en los casos anteriores.

**Modos de un derivado.** `meshEsuP`/`EsuE` admiten `sparse|cell|index` (y `EsuE`
además `bynode|byedge|byface`). La **propiedad** `M.EsuP` devuelve la forma
canónica cacheada (`sparse`); el **método** `M.esup(mode)` / `M.esue(mode,by)`
ofrece las variantes (derivables baratas desde la sparse o cacheadas por clave
`EsuE:byedge`, etc.).

## 5. Datos vs. caché al editar

- **Caché (derivados)**: se *invalida* (se vacía) por dependencia. Se regenera al
  pedirse.
- **Datos (campos de usuario)**: se *mantienen*. Si cambia el nº de nodos/caras,
  las filas de los campos deben acompañar el cambio.

Problema conocido: `M.faces(1:10,:) = []` por la vía cruda **no sabe** qué filas
se quitaron, así que no puede remapear los campos por cara. Estrategia (igual de
conservadora que el `Mesh.m` actual, que recorta/rellena con aviso):

- `set.faces`/`set.nodes` reconcilian tamaños de campos (recortar/rellenar con
  `NaN`) **con `warning`** si el recuento cambió.
- Para edición estructural que preserva campos, usar métodos dedicados que
  delegan en las funciones existentes: `M.removeFaces(idx)`, `M.removeNodes(idx)`,
  `M.tidy(...)` (→ `MeshRemoveFaces`, `MeshRemoveNodes`, `MeshTidy`), que recortan
  **todos** los `xyz*`/`tri*` a la vez.

## 6. Integración con el toolbox (puente + delegar)

- **Alias de compatibilidad**: `M.xyz` ≡ `M.nodes`, `M.tri` ≡ `M.faces`
  (propiedades `Dependent` de lectura/escritura).
- **Puente struct**: `S = M.toStruct()` materializa el `struct` legado
  (`.xyz`, `.tri`, `.xyzFIELD*`, `.triATT*`, `.celltype`, `.triNORMALS`, …).
  `msh.fromStruct(S)` hace el camino inverso. El **constructor absorbe** la lógica
  de conversión del actual `MESH\Mesh.m` (handles de patch, `triangulation`,
  `delaunay`, contornos, etc.).
- **Métodos que delegan** (patrón de 3 líneas):

  ```matlab
  function M = tidy(M, varargin)
    M = msh.fromStruct( MeshTidy( M.toStruct(), varargin{:} ) );
  end
  ```

  Así se reutilizan las 126 funciones sin reescribir nada. Migración posterior
  (opcional) de las *core* a métodos nativos, sin urgencia.

Requisito: `MESH\` en el path (`addpath` del toolbox).

## 7. Persistencia (save/load)

- `nodes`, `faces`, campos y especiales fijados (normales/UV/textura del usuario)
  son propiedades normales → se serializan.
- La caché (`cache_`) es **`Transient`** → no se guarda. Tras `load` queda vacía;
  `loadobj` la reinicializa a un `mshCache` fresco y se repuebla perezosamente.

## 8. Layout de ficheros

```
C:\repos\msh\
├─ @msh\msh.m         ← classdef msh   (esqueleto: props + firmas + stubs)
├─ mshCache.m         ← classdef mshCache < handle  (mecanismo de caché, escrito)
├─ msh_DESIGN.md      ← este documento
└─ MESH\ …            ← las ~126 funciones existentes (sin tocar)
```

(Opcional a futuro) los métodos que delegan pueden moverse a ficheros sueltos
dentro de `@msh\` para no inflar `msh.m`.

## 9. Decisiones abiertas (a cerrar en la revisión)

1. **Ergonomía de campos `M.xyzFOO` / `M.triBAR`.**
   - *(A, recomendada)* sobrecargar `subsref`/`subsasgn` para que
     `M.xyzTemp = …` funcione como en el `struct` de hoy. Máxima fidelidad; es el
     punto de más riesgo de implementación (hay que enrutar `.`, `()`, `{}` y no
     romper el *dispatch* de métodos/propiedades).
   - *(B, simple)* API por métodos: `M = M.addField('xyzTemp', v)`,
     `v = M.getField('xyzTemp')`, + `nodeFields`/`faceFields` como `struct`.
     Robusto, pero se pierde el acceso directo tipo campo.
2. **Padding 2D→3D**: almacenar nativo (2 o 3 col) y padear on-demand
   *(recomendado, preserva quad vs tetra)*, o normalizar siempre a 3 col con un
   flag `is2D`.
3. **`rotate`/`transform` y normales**: por defecto **invalidar** y recomputar;
   optimización posible: rotar las normales cacheadas en vez de tirarlas.
4. **Nombre del constructor corto**: la clase es `msh`; ¿alias `Msh`/`MSH`?
5. **¿`celltype` como propiedad de solo lectura** (siempre derivada) o dato
   fijable para forzar quad vs tetra en el caso ambiguo de 4 nodos?
