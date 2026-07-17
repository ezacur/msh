# Tutorial: la clase `msh`

Contenedor de mallas con **semántica de valor** y un sistema de **cachedProps**:
derivados definidos por un registro *nombre → { cómo se computa, cómo reacciona
a los eventos de edición }*, calculados perezosamente, y que ante cada edición
**sobreviven**, **caen**, o quedan **pendientes de una actualización barata**
(replay perezoso). Envuelve el toolbox legado de structs (`MESH\`) sin
sustituirlo: cualquier función `mesh*`/`Mesh*` sigue funcionando a través del
puente `toStruct()`.

La idea en una frase: *pides `M.cached.BVH` y se calcula; lo vuelves a pedir y
es gratis; editas la malla y cada derivado hace lo más barato que su definición
permita — morir, sobrevivir, o actualizarse solo.*

> Documentos hermanos: `msh_DESIGN.md` (decisiones de diseño),
> `BVH_TUTORIAL.md` (el motor de queries `BVH` / `bvhClosestElement` /
> `bvhIntersectRay` en profundidad).

---

## 0. Instalación y verificación

```matlab
addpath C:\repos\msh          % @msh + cacheHandle + cacheView
addpath C:\repos\msh\BVH      % motor de queries: BVH / bvhClosestElement / ...
addpath C:\repos\msh\MESH     % toolbox legado
addpath C:\repos\msh\tools    % transform (usado por M.transform)
```

Para las queries hacen falta los MEX del motor (una sola vez, ver
`BVH_TUTORIAL.md` §0). Verificación completa de la clase:

```matlab
cd C:\repos\msh\MESH
test_msh          % 11 bloques: construcción, eventos, replay, proxy, queries...
```

---

## 1. Construcción y datos

```matlab
M = msh( V , F );             % coordenadas + conectividad
M = msh( S );                 % struct legado (.xyz/.tri + campos)
M = msh( M0 );                % copia (valor; comparte cache mientras no edites)
M = msh();                    % malla vacía
M = msh( V , F , 'xyzTemp' , t , 'DEBUG' , true );   % pares nombre/valor extra
```

Desde struct se aceptan varios dialectos (la **entrada** habla el idioma
legado; los nombres nuevos son de la clase):

| Campo de entrada | Interpretación |
|---|---|
| `.xyz` / `.tri` | formato nativo del toolbox |
| `.vertices` / `.faces` | estilo `patch`/FV; si `min(faces(:))==0` se asume 0-based y se corrige a +1 |
| `.xyzNOMBRE` / `.triNOMBRE` | atributos por nodo / por cara (crop/NaN-pad con warning si las filas no cuadran) |
| `.texture` | va a `M.INFO.texture` |
| cualquier otro campo | va a parar a `M.INFO` (no se pierde nada) |

**Ojo con el alias `.faces`**: la corrección 0-based se dispara con cualquier
cero, así que las mallas **mixtas 0-padded deben entrar por `.tri`** (el cero
ahí es relleno, no un índice).

### Los datos y sus reglas

| Propiedad | Qué es |
|---|---|
| `M.V` | coordenadas `nV × 2\|3`, **siempre `double`** |
| `M.F` | conectividad `nF × k`, k ∈ {1,2,3,4}, 0-padded, **siempre `int32`** |
| `M.nsd` | 2 ó 3 — el número de columnas de `V` |
| `M.nV`, `M.nF` | número de vértices / caras |
| `M.ct` | código VTK del tipo de celda (alias legado: `M.celltype`) |
| `M.viz` | preferencias de plot (datos, sin efecto en cache) |
| `M.INFO` | metadatos libres — la textura vive en `M.INFO.texture` |
| `M.DEBUG` | narrador de la cadena de procesos (§16) |

* **La dimensión `nsd` es el número de columnas de `V`, nunca se infiere de
  los valores**: una malla 3D con todas las z = 0 sigue siendo 3D (puede
  deformarse fuera del plano sin cambiar de tipo).
* **planar / flat** (descripción, no tipo): una malla con `nsd == 3` cuyos
  vértices caen todos en **un** plano es **planar** (`M.isPlanar()`); si ese
  plano es exactamente `z == 0` es **flat** (`M.isFlat()`, el caso "casi 2D").
  Toda flat es planar. El display las anuncia (§17).
* `M.F`: 4 nodos no-cero = **tetraedro**, siempre (los quads no existen aquí).
  Mixtas: ceros de relleno **al final**.
* **Tipos obligatorios**: los setters validan (enteros, no negativos, en rango
  int32, ≤ nV) y convierten cualquier entrada numérica. Los **atributos**
  por nodo/cara conservan el tipo que les des (`single`, `uint8`, ...). En la
  frontera de los MEX el contrato ya se cumplía por construcción: blobs con
  coordenadas `double` e índices `int32`; la conectividad nunca cruza al C++.

```matlab
V = randn(700,3);  V = V./sqrt(sum(V.^2,2));      % esfera unidad
M = msh( V , convhulln(V) );

M.nV          % 700
M.nF          % ~1396
M.nsd         % 3
M.ct          % 5 (VTK: triángulos)   == M.celltype (alias)
```

---

## 2. cachedProps: el modelo

Cada malla lleva un **registro** (parte del valor: viaja con las copias y se
guarda en el `.mat`) de derivados cacheables:

```
nombre  ->  computeFcn  @(m) ...        cómo producir el valor
            eventos     evento -> handler | []    cómo reaccionar al editar
```

y una **cache de valores** (un handle compartido entre copias, con
copy-on-write al editar: las copias hermanas nunca se pisan).

Acceso — dos puertas equivalentes:

```matlab
M.cached.BVH          % la forma "real"
M.BVH_                % el alias: sufijo '_' = "esto es una cachedProp"
```

El primer acceso computa y guarda (**MISS**); los siguientes son **HIT**
gratis. `M.cached` a secas imprime la tabla completa de definiciones y
estados:

```
  cachedProps (9) -- acceso M.cached.<nombre> o M.<nombre>_ :
    BVH         [blob BVH: 1,396 elems, aabb]
                  eventos: transform->handler, changeCoords->handler, changeConnectivity->invalida
    boundary    [0x2 double]
                  eventos: changeConnectivity->invalida
    bbox        (sin calcular)
                  eventos: changeCoords->invalida
    ...
```

### Convención de nombres

| Forma | Qué es |
|---|---|
| `M.V`, `M.F`, `M.ct`, `M.plot()`, ... | propiedades y métodos **reales** de la clase |
| `M.nombre_` — sufijo `_` | **reservado en exclusiva** para los alias de cachedProps: `M.nombre_` ≡ `M.cached.nombre` |
| `VERTICES`, `FACES`, `VATTS`, `FATTS`, `cachePROPS`, `CACHE` — MAYÚSCULAS | almacenamiento privado interno (invisible desde fuera; solo relevante si extiendes `msh.m`) |

**Alias legados** (interceptados en `subsref`, solo lectura): `M.celltype` →
`M.ct`, `M.xyz` → `M.V`, `M.tri` → `M.F`. Los dos últimos existen para la
transición (§14): permiten que código legado *de solo lectura* trate un `msh`
como si fuera el struct viejo; el saneamiento futuro los retirará. Como el
sufijo es namespace exclusivo de la cache, una cachedProp puede llamarse
**como quieras — incluso `nodes` o `V`** — sin ambigüedad:

```matlab
M = M.defineCachedProp( 'nodes' , @(m) m.nV , 'nodecount' , [] );
M.V              % la propiedad real (coordenadas)
M.nodes_         % tu cachedProp (== M.cached.nodes)
```

---

## 3. Eventos: qué dispara cada edición

Las ediciones no "invalidan por tags": **disparan eventos**, del más
específico al más general:

| Edición | Eventos disparados |
|---|---|
| `M.V = X` (deformación) | `changeCoords` |
| `M.V = X` (cambia nº de filas) | `changeNodeCount` + `changeCoords` |
| `M.V = X` (cambia 2D↔3D) | `changeDim` + `changeCoords` |
| `M.F = T` | `changeConnectivity` (+ `changeFaceCount` si cambian las filas) |
| `M = M.transform( T )` | `transform(T)` + `changeCoords` |

Y cada cachedProp declara, para cada evento, una de tres reacciones:

* **no declarado** → el valor **sobrevive** intacto (el evento no le afecta);
* **declarado con `[]`** → el valor **cae** (invalidar, recomputar al pedirse);
* **declarado con un handler** → el valor queda **pendiente**: el handler lo
  actualizará *barato* en el próximo acceso (§4).

El ejemplo canónico es el BVH de fábrica, que es pura declaración:

```matlab
% (así está registrado internamente)
M = M.defineCachedProp( 'BVH' , @(m) BVH( toStruct(m) ) , ...
      'transform'          , @(B,m,T) BVH( B , T )            , ... % plegado O(1)
      'changeCoords'       , @(B,m)   BVH( B , toStruct(m) )  , ... % refit O(n)
      'changeConnectivity' , [] );                                     % rebuild
```

Deformas → el blob queda pendiente y se *refita* al pedirlo. Transformas → se
*pliega* en O(1). Cambias la topología → cae y se reconstruye. Y `boundary`,
que solo declara `changeConnectivity`, **sobrevive a todas las deformaciones y
transformaciones sin inmutarse**.

---

## 4. Replay perezoso: los handlers corren en el acceso, no en el evento

Cuando un evento encuentra un handler, **no lo ejecuta**: anota el evento en
el log de la entrada y sigue. El trabajo se hace una sola vez, en el próximo
acceso:

```
M.V = X1;        % [msh] EVENT {changeCoords}: pendientes {BVH} | caen {surfCent} | ...
M.V = X2;        % [msh] EVENT {changeCoords}: pendientes {BVH} | ...
M.V = X100;      % ... 100 ediciones: coste ~0
B = M.BVH_;      % [msh] RPLAY 'BVH' -> sync absoluto via changeCoords en 4.4 ms
                 % UN solo refit, contra la malla FINAL
```

Reglas del replay:

* `transform` es **incremental**: su handler recibe cada `T` del log
  (`@(v,m,T)`) y se aplican en secuencia — dos transforms = dos pliegues O(1).
* Los demás eventos son **absolutos**: su handler recibe la malla *actual*
  (`@(v,m)`) que ya contiene todo lo ocurrido — así N deformaciones colapsan
  en un solo sync, y un log mixto (transform + edición ciega) lo subsume un
  único sync absoluto.
* **Cascada ante fallos**: si el handler de `transform` lanza (p.ej. la T no
  era una semejanza y el plegado del BVH da `BVH:notSimilarity`), se cae al
  sync absoluto (refit); si ese también falla, se recomputa desde cero. Nunca
  verificas qué tipo de T tienes — la cascada lo descubre sola.

Nota de diseño: esto generaliza (y sustituye) el viejo truco de las "dos
claves" del BVH — el valor pendiente con su log *es* la jerarquía
superviviente esperando el refit.

---

## 5. Las cachedProps de fábrica

Toda malla nace con estas definiciones (todas sobreescribibles, §6):

| Nombre | Qué es | Eventos declarados |
|---|---|---|
| `BVH` | blob de aceleración (`BVH`, aabb) | `transform`→plegado, `changeCoords`→refit, `changeConnectivity`→invalida |
| `boundary` | facetas del borde abierto (`MeshBoundary`; matriz cruda, vacía si es cerrada) | `changeConnectivity`→invalida |
| `edges` | aristas únicas (`meshEdges`) | `changeConnectivity`→invalida |
| `EsuP` / `PsuP` | adyacencias por punto (sparse) | `changeConnectivity`, `changeNodeCount`→invalida |
| `EsuE` | elems vecinos por elem (sparse) | `changeConnectivity`→invalida |
| `bbox` | `[min ; max]` (`meshBB`) | `changeCoords`→invalida |
| `surfCent` | `[área, centroide]` (una llamada a `meshSurface` sirve ambos) | `changeCoords`, `changeConnectivity`→invalida |
| `triNORMALS` | normales por cara (`meshNormals`) | `changeCoords`, `changeConnectivity`→invalida, **`transform`→rotación O(n)** |

**El acceso a estas entradas es SOLO por el sufijo o el proxy** (`M.BVH_`,
`M.boundary_`, `M.surfCent_`, `M.cached.EsuP`, ...) — el nombre desnudo queda
reservado para datos y métodos, **sin excepciones**. El único derivado con
nombre desnudo es `M.ct` (trivial, se computa directo, sin cache).

* `M.surfCent_` devuelve `[área, centroide]` en una fila (una llamada a
  `meshSurface` sirve los dos valores): `sc = M.surfCent_; a = sc(1);
  c = sc(2:4);`
* `M.triNORMALS_` es la normal **canónica** por cara (cross products crudos).
  El atributo `NORMALS` del usuario (p.ej. normales suavizadas con
  `meshNormals(M,k)`, o medidas) es un dato independiente: viaja en
  `toStruct()`, las funciones legadas lo honran, y `transform()` lo rota —
  pero no se mezcla con la cachedProp.
* Normales por nodo: **parametrizadas** (`meshNormals(M,'a'|'u'|'b'|...)` — no
  hay "el valor" que cachear), así que no hay entrada de fábrica: llama a
  `meshNormals(toStruct(M),'angle')` directamente, o define la tuya con el
  método fijado (`defineCachedProp('nAngle', @(m) meshNormals(toStruct(m),'angle'), ...)`).

El criterio general (cache vs atributo): **la cache es una consecuencia de la
malla — siempre correcta o inexistente; un atributo es una afirmación del
usuario — sobrevive a las ediciones y puede quedarse rancio en silencio.**

---

## 6. Define las tuyas: `defineCachedProp`

```matlab
M = M.defineCachedProp( nombre , @(m) ... [, evento , handler|[] , ...] );
M = M.removeCachedProp( nombre );
```

La `computeFcn` recibe **el msh** (usa `toStruct(m)` para las funciones
legadas, o cualquier otra cachedProp: componen). Los nombres de evento admiten
alias (`'coords'`, `'connectivity'`, `'nodecount'`, ...). Registrar **no
ejecuta nada** — solo enseña a la malla qué hacer cuando se lo pidas.

```matlab
% aristas del borde, solo sensible a la topología
M = M.defineCachedProp( 'freeEdges' , @(m) MeshBoundary( m.F ) , ...
                        'connectivity' , [] );

% una métrica de calidad por cara (meshQuality parametrizada -> una entrada por métrica)
M = M.defineCachedProp( 'aspectRatio' , ...
        @(m) meshQuality( toStruct(m) , 'aspectratio' ) , ...
        'coords' , [] , 'connectivity' , [] );

% con handler de transform: el centroide de vértices se transforma EXACTO bajo
% cualquier afín (la media conmuta con T) -> nunca hace falta recomputarlo
M = M.defineCachedProp( 'centro' , @(m) mean( m.V , 1 ) , ...
        'coords' , [] , 'nodecount' , [] , ...
        'transform' , @(c,m,T) [ c , 1 ] * T(1:3,:).' );

v = M.cached.aspectRatio;      % primer acceso: computa
v = M.aspectRatio_;            % alias; ahora es un hit
M2 = M.transform( T );  M2.centro_;   % replay incremental: sin recomputar
```

**Ejemplo completo y verificado**: `MESH\meshQuality_as_cachedProps.m` registra
*todas* las métricas de `meshQuality` (una cachedProp por métrica, según el
celltype) con la matemática de actualización de cada una — las adimensionales
son invariantes ante semejanzas, longitudes/áreas escalan `s`/`s²`, y
`volume`/`signedvolume`/`orientation` se actualizan **exactas ante cualquier
afín** vía `det(A)`:

```matlab
M = meshQuality_as_cachedProps( M );              % todas las de su celltype
M = meshQuality_as_cachedProps( M , 'area' , 'aspectratio' );   % o un subconjunto
q = M.aspectratio_;                % computa y cachea
M2 = M.transform( T );             % semejanza: TODAS quedan pendientes...
q2 = M2.aspectratio_;              % ...y se actualizan sin recomputar
```

**Sobreescribir una de fábrica** reemplaza la definición completa (y descarta
el valor viejo, que ya no corresponde) — redeclara los handlers que quieras
conservar:

```matlab
M = M.defineCachedProp( 'BVH' , @(m) BVH( toStruct(m) , [2 8] , 'sphere' ) , ...
      'transform'    , @(B,m,T) BVH( B , T ) , ...
      'changeCoords' , @(B,m)   BVH( B , toStruct(m) ) , ...
      'connectivity' , [] );
% M.BVH_, M.closestElement, M.intersectRay usan ahora el blob de esferas
```

Como `msh` es una value class, `defineCachedProp` **devuelve** la malla
modificada — sin el `M =` no pasa nada (así funciona el lenguaje; las
variables auxiliares son responsabilidad del programador).

Detalles: sin detección de ciclos (una computeFcn que se pida a sí misma
recursa); los function handles anónimos serializan en el `.mat` con la letra
pequeña de siempre (capturan su workspace).

---

## 7. El proxy `M.cached`: la API completa

| Expresión | Hace | Devuelve |
|---|---|---|
| `M.cached` | tabla de definiciones y estados | vista |
| `M.cached.BVH` | el valor (computa/replay si hace falta) | valor |
| `M.cached.BVH.frame` | indexa **dentro** del valor | `blob.frame` |
| `M.cached.BVH.delete` | borra el **valor** (la definición queda) | — |
| `M = M.cached.BVH.removeProp` | borra definición **y** valor | msh nuevo |
| `M = M.cached.BVH.set( x )` | siembra un valor a mano ("dangerous, but...") | msh nuevo |
| `M.cached.BVH.changeCoords` | el handler del evento (`[]` = invalida) | handle |
| `M.cached.BVH.changeCoords( B , M )` | lo ejecuta con tus argumentos | valor |

Todo funciona igual por el alias (`M.BVH_.delete`, `M.BVH_.nE`, ...).

**Semántica de compartición** — la distinción importante:

* `.delete` es un *statement* que muta el handle **compartido**: las copias
  hermanas también pierden ese valor. Benigno — como mucho, alguien recomputa.
* `.set(x)` es **conservador (COW)**: devuelve un msh nuevo con su propia
  cache sembrada; las hermanas ni se enteran. Sigue siendo peligroso en el
  sentido que importa — nadie verifica que `x` sea correcto.
* `.removeProp` y `defineCachedProp` tocan la *definición* (que vive en el
  valor): devuelven un msh nuevo, hermanas intactas.

Colisiones: tras el nombre, `delete`/`removeProp`/`set` y los nombres de
evento son operaciones reservadas; cualquier otro nombre se reenvía como
indexación del valor. Si un valor tuviera un campo llamado como una operación,
sácalo en dos pasos (`B = M.cached.BVH; B.delete` ya indexa normal).

---

## 8. El mecanismo valor + handle + copy-on-write

`msh` es una value class (copias independientes), pero los **valores**
cacheados viven en un handle compartido (COW = *copy-on-write*: se comparte
hasta que alguien escribe). La combinación da la semántica del ejemplo
fundacional:

```matlab
M1 = msh( V , F );
e  = M1.EsuP_;           % cachedProp (sufijo!): se calcula y queda en la
                         % cache COMPARTIDA
M2 = M1;                 % copia: mismos datos, misma cache
M3 = M1.transform( RT ); % rotar no toca la topología...
e3 = M3.EsuP_;           % ...así que M3 YA TIENE EsuP (hit, gratis)

M1.F(1:10,:) = [];       % edición de conectividad en M1
% -> M1 recibe una cache NUEVA donde EsuP/boundary/BVH... han caído
% -> M2 y M3 CONSERVAN la suya intacta
M2.EsuP_;                % sigue siendo un hit
```

Cada evento hace **copy-on-write del handle**: la instancia editada evoluciona
a un handle nuevo llevándose los supervivientes y los pendientes; las hermanas
conservan el suyo. Resolver un pendiente o depositar un MISS sí muta el handle
compartido — es benigno (mismo valor correcto para todas).

---

## 9. `transform()`

```matlab
M2 = M.transform( T );     % T homogénea 4×4 / 3×4 / 3×3 (2D: 3×3 homogénea 2D)
```

Dispara el evento semántico `transform(T)` — la única edición que *sabe qué
está pasando* (el setter de `V` es ciego: no puede distinguir una traslación
de una deformación). Con las definiciones de fábrica:

* `BVH` → pendiente de **plegado O(1)** (~14 µs) si T es semejanza; si no, la
  cascada lo baja a refit O(n).
* `triNORMALS` cacheadas → pendiente de **rotación O(n)** de las filas.
* Los atributos `NORMALS` del usuario se rotan aquí mismo (vía `tools\transform`).
* `bbox`, `surfCent` → caen (declaran `changeCoords` sin handler).
* `boundary`, `edges`, `EsuP`... → **sobreviven** (ni se enteran).

Si trasladas "a mano" (`M.V = M.V + t`) pagas la versión ciega:
`changeCoords`, refit en vez de plegado. Usa `transform` para transformar.

```matlab
B1 = M.BVH_;
M2 = M.transform( [0.9*RT , [1;2;3] ; 0 0 0 1] );
B2 = M2.BVH_;                      % replay: plegado O(1)
isequal( B2.child4 , B1.child4 )   % true  (misma jerarquía)
isequal( B2.X      , B1.X )        % true  (misma geometría empacada)
isequal( B2.frame  , B1.frame )    % false (el marco absorbió T)
```

---

## 10. Edición y guardias

```matlab
M.V = X2;                  % nV×2|3 numérico; error si deja caras colgando
M.F = T2;                  % nF×k, k∈1..4, enteros ≥ 0; error si referencia
                           % vértices inexistentes o índices no enteros
M.F(1:10,:) = [];          % la edición indexada también pasa por el setter
M.V(end+1,:) = [0 0 0];    % crecer vértices es legal
```

* Encoger `M.V` bajo el máximo índice usado → `error('msh:nodes',…)`;
  índices fuera de rango / no enteros en `M.F` → `error('msh:faces',…)`.
* Si cambias el número de vértices/caras con atributos asociados, se
  **reconcilian** (crop / NaN-pad) con `warning('msh:field',…)` — para
  remapear de verdad usa `removeNodes`/`removeFaces`/`tidy`.

---

## 11. Queries

```matlab
[e,cp,d,bc,F] = M.closestElement( P );          % elemento más cercano
[e,cp,d]      = M.closestElement( P , Dmax );   % con radio de corte
[xyz,c,t,rid] = M.intersectRay( rays );         % 'first' por defecto
[xyz,c,t,rid] = M.intersectRay( rays , 'all' ); % 'first'|'last'|'all'|'any'
```

Usan (y rellenan de paso) la cachedProp `BVH` — la que esté definida: si la
sobreescribiste con esferas u otro leaf size, es la que trabaja. Tras un
`transform` el blob llega plegado; tras una deformación, refitado. **Nunca hay
blob rancio**: la clase intercepta todas las ediciones. Detalles del motor
(salidas `bc`/`F`, semántica de `Dmax`, modos de rayos, rendimiento) en
`BVH_TUTORIAL.md`.

Visualizador del blob (↑/↓ profundidad, `a` todo, `l` hojas, `f` marco):

```matlab
M.plotBVH();
```

---

## 12. Atributos por nodo y por cara

```matlab
M = M.addField( 'xyzTemp'  , temp  );   % por nodo  (prefijo explícito)
M = M.addField( 'triLabel' , label );   % por cara
M = M.addField( 'Color'    , c     );   % sin prefijo: se infiere por filas
                                        % (error si nV == nF: sé explícito)
v = M.getField( 'xyzTemp' );            % con o sin prefijo
M.hasField( 'triLabel' )                % true
L = M.fieldNames();                     % L.node = {'xyzTemp',…}, L.face = …
M = M.rmField( 'xyzTemp' );
```

* Filas validadas contra `nV`/`nF` al añadir; reconciliación crop/NaN-pad al
  cambiar tamaños (§10). Cada atributo conserva su propio tipo.
* `NORMALS` es especial solo en una cosa: `transform()` **rota** los atributos
  `xyzNORMALS`/`triNORMALS` (vía `tools\transform`) en vez de dejarlos rancios.
  Son datos tuyos, independientes de la cachedProp `triNORMALS_` (§5).
* En `toStruct()` los atributos salen con su prefijo legado (`xyz*`/`tri*`).

---

## 13. `viz`, `INFO`, `plot`, textura

```matlab
M.viz.FaceColor = [0.8 0.2 0.2];       % preferencias persistentes de plot
M.viz.EdgeColor = 'none';
M.INFO.paciente = 'ID-042';            % metadatos libres
M.INFO.texture  = imread('piel.png');  % la textura es un metadato más

M.plot();                              % plotMESH con las preferencias
M.plot('FaceAlpha',0.5);               % precedencia: defaults < viz < args
```

`viz` e `INFO` son datos normales: se copian con la malla, se guardan en el
`.mat`, y no participan en la cache. Los campos desconocidos de un struct de
entrada acaban en `INFO`; `INFO.texture` sale como `.texture` en `toStruct()`
(y un `.texture` de entrada entra ahí).

---

## 14. Puente con el toolbox legado

```matlab
S = M.toStruct();            % struct .xyz/.tri/.xyzF*/.triA*/.texture
R = MeshFlipFaces( S );      % cualquier función legada funciona
M2 = msh( R );               % y de vuelta
```

El puente habla el idioma viejo: exporta `.xyz` (double) y `.tri` (int32).

**Transición a msh-nativas**: `@msh\private\` contiene copias de las legadas
que la clase usa internamente (`MeshBoundary`, `meshEdges`, `meshEsuP/PsuP/EsuE`,
`meshBB`, `meshSurface`, `meshCelltype`) — desde los métodos de `msh` esas
copias tienen precedencia sobre `MESH\`, así que se pueden ir saneando
(reescribir sobre `.V`/`.F`) sin tocar a los consumidores legados, que siguen
usando la copia de `MESH\`. Cuando una función legada se quede sin consumidores
fuera de la clase, se borra de `MESH\` y su vida sigue solo en `private`. Los
alias de lectura `M.xyz`/`M.tri` (§2) ayudan en el mientras tanto: una legada
de solo lectura acepta un `msh` directamente.

Delegaciones ya envueltas (devuelven un `msh` nuevo, con cache fresca):

```matlab
M = M.tidy( … );                 % MeshTidy
M = M.removeFaces( idx );        % MeshRemoveFaces
M = M.removeNodes( idx );        % MeshRemoveNodes (remapea caras y atributos)
M = M.append( M2 , S3 , … );     % MeshAppend (acepta msh y structs mezclados)
```

---

## 15. `save` / `load`

```matlab
save malla.mat M
L = load('malla.mat');
L.M.closestElement(P);     % funciona: el BVH se recomputa perezoso
```

Las **definiciones** (el registro, con sus handlers) viajan en el `.mat`; los
**valores** cacheados no (el handle es Transient): un `.mat` de `msh` pesa lo
que pesan sus datos y todo se recalcula al primer uso. Para transportar un BVH
precalculado, guárdalo aparte (`B = M.BVH_` es un struct-valor serializable).

---

## 16. `DEBUG`: ver la cadena de procesos

```matlab
M.DEBUG = true;                    % o msh(V,F,'DEBUG',true)
```

Narra por consola cada paso — la tabla de eventos de §3, en vivo:

```
[msh] MISS  'boundary' -> calculado en 1.84 ms
[msh] HIT   'boundary'
[msh] SET   V 700x3 -> 700x3
[msh] EVENT {changeCoords}: pendientes {BVH} | caen {surfCent} | sobreviven {boundary}
[msh] RPLAY 'BVH' -> sync absoluto via changeCoords en 4.42 ms
[msh] TRANS transform() sobre 700 nodos
[msh] EVENT {transform+changeCoords}: pendientes {BVH, triNORMALS} | caen {...} | ...
[msh] RPLAY 'BVH' -> 1 transform(s) incremental(es) en 0.05 ms
[msh] QUERY closestElement: 500 puntos
[msh] QUERY closestElement resuelta en 0.62 ms
[msh] DEF   cachedProp 'aspectRatio' definida (eventos: changeCoords, changeConnectivity)
```

`DEBUG` es un dato normal (viaja con las copias, no afecta a la cache).

---

## 17. Display

Mostrar la variable (o `disp(M)`) **no** evalúa ningún getter — solo enseña lo
que ya existe, pero con lo que existe es exhaustivo:

```
M =
  msh: 700 vertices (3D), 1,401 faces (mixed: 5 segments + 1,396 triangles)
    node fields: xyzTemp [700x2 double]
    viz:  FaceColor = [0.8 0.2 0.2]
    INFO: paciente = 'ID-042'
    cached:
      BVH       [blob BVH: 1,401 elems, aabb]
      boundary  [42x2 double]
      surfCent  (pendiente: changeCoords)
    definidas sin calcular: EsuE, EsuP, PsuP, bbox, edges, triNORMALS
```

* Cabecera con los conteos con **separadores de miles** (`2,000,000 vertices`
  se distingue de `200,000` de un vistazo), desglose de **tipos de celda**
  (mixtas incluidas) y la descripción **planar/flat**: `(3D, flat (z = 0))` si
  todas las z son 0 exactas, `(3D, planar: z = 2.5)` si el plano es
  axis-aligned, `(3D, planar)` si es un plano cualquiera. La dimensión sigue
  siendo el número de columnas.
* Atributos con dimensiones y clase, `viz`/`INFO` como pares `nombre = valor`
  (la textura aparece dentro de `INFO`).
* Entradas cacheadas con su **valor si es escalar o vector de ≤ 3 números**,
  su estado pendiente (con los eventos del log), o `[tamaño clase]`.
* `M.cached` da la vista completa con los eventos declarados por entrada.

---

## 18. Limitaciones y pendientes

* Sin detección de ciclos entre cachedProps (una computeFcn que se pida a sí
  misma recursa).
* Eventos semánticos futuros: `removeFaces(idx)` / `removeNodes(map)` con los
  mapas de índices que ya devuelven `MeshTidy`/`MeshRemove*` — permitirían
  updates quirúrgicos de las adyacencias (borrar filas de `EsuP` en vez de
  recomputarla).
* Las delegaciones (`tidy`, `append`, …) devuelven un `msh` con cache fresca
  y **registro de fábrica** (no heredan tus `defineCachedProp`).
* El acceso por el proxy paga un pequeño peaje de `subsref` (µs); los métodos
  y propiedades reales no se interceptan en la práctica.
* Los `.mat` guardados con esquemas de propiedades anteriores al renombrado
  (nodes/faces/vizProp/info) no cargan limpio — proyecto en desarrollo.
* No hay aritmética de mallas sobrecargada (`M1 + M2`…): usa `append`.
