# Tutorial: la clase `msh`

Contenedor de mallas con **semántica de valor** y un sistema de **CPs**
(*cached props*): derivados definidos por un registro *nombre → { cómo se
computa, cómo reacciona a los eventos de edición }*, calculados perezosamente,
y que ante cada edición **sobreviven**, **caen**, o quedan **pendientes de una
actualización barata** (replay perezoso). Envuelve el toolbox legado de structs
(`MESH\`) sin sustituirlo: cualquier función `mesh*`/`Mesh*` sigue funcionando
a través del puente `ToStruct()`.

La idea en una frase: *pides `M.bvh` y se calcula; lo vuelves a pedir y es
gratis; editas la malla y cada derivado hace lo más barato que su definición
permita — morir, sobrevivir, o actualizarse solo. Y si quieres forzar el
recálculo: `M.bvh_`.*

> Documentos hermanos: `msh_DESIGN.md` (decisiones de diseño),
> `BVH_TUTORIAL.md` (el motor de queries `BVH` / `bvhClosestElement` /
> `bvhIntersectRay` en profundidad).

---

## 0. Instalación y verificación

```matlab
addpath C:\repos\msh          % @msh + cacheHandle + cacheView
addpath C:\repos\msh\BVH      % motor de queries: BVH / bvhClosestElement / ...
addpath C:\repos\msh\MESH     % toolbox legado
addpath C:\repos\msh\tools    % transform (usado por M.Transform)
```

Para las queries hacen falta los MEX del motor (una sola vez, ver
`BVH_TUTORIAL.md` §0). Verificación completa de la clase:

```matlab
cd C:\repos\msh\MESH
test_msh          % 12 bloques: construcción, eventos, replay, proxy, queries...
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
| `M.VIZ` | preferencias de plot (datos, sin efecto en cache) |
| `M.INFO` | metadatos libres — la textura vive en `M.INFO.texture` |
| `M.DEBUG` | narrador de la cadena de procesos (§16) |

* **La dimensión `nsd` es el número de columnas de `V`, nunca se infiere de
  los valores**: una malla 3D con todas las z = 0 sigue siendo 3D (puede
  deformarse fuera del plano sin cambiar de tipo).
* **planar / flat** (descripción, no tipo): una malla con `nsd == 3` cuyos
  vértices caen todos en **un** plano es **planar** (`M.IsPlanar()`); si ese
  plano es exactamente `z == 0` es **flat** (`M.IsFlat()`, el caso "casi 2D").
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

## 2. La convención del API: tres espacios de nombres por CASO

El miembro que tocas se reconoce por cómo está escrito:

| Caso | Qué es | Ejemplos |
|---|---|---|
| `MAYÚSCULAS` | **propiedades** (datos) | `V`, `F`, `VIZ`, `INFO`, `DEBUG`, y el proxy `CP` |
| `Capitalized` | **métodos** públicos | `Plot`, `Transform`, `Tidy`, `ClosestElement`, `DefineCP`, ... |
| `minúscula` | **CPs** (derivados cacheados) | `bvh`, `esup`, `boundary`, `triNormals`, ... |

Excepciones documentadas: los **contadores** `nsd`/`nV`/`nF`/`ct` (minúsculos,
se leen mil veces al día), los métodos del **protocolo MATLAB** (`subsref`,
`disp`, `loadobj`, ... — el lenguaje exige esos nombres) y los **alias
transicionales** `xyz`/`tri`/`celltype` (§14).

### El modelo de las CPs

Cada malla lleva un **registro** (parte del valor: viaja con las copias y se
guarda en el `.mat`) de derivados cacheables:

```
nombre  ->  computeFcn  @(m) ...        cómo producir el valor
            eventos     evento -> handler | []    cómo reaccionar al editar
```

y una **cache de valores** (un handle compartido entre copias, con
copy-on-write al editar: las copias hermanas nunca se pisan — §8).

### Acceso a las CPs: desnudo LEE, sufijo `_` RECALCULA

```matlab
M.bvh          % LEE: HIT si está fresca, replay si está pendiente,
               %      computa y guarda si no existe (MISS)
M.bvh_         % RECALCULA a la fuerza: descarta el valor (y su log),
               %      computa fresco, guarda y devuelve
M.bvh.frame    % indexar dentro del valor (también M.bvh_.frame)
M.CP           % el plano de control: tabla, .delete, .set, eventos (§7)
```

El primer acceso desnudo computa y guarda (**MISS**); los siguientes son
**HIT** gratis. La cache es invisible en el uso corriente — escribes `M.esup`
como si fuera un dato más. `M.CP` a secas imprime la tabla completa:

```
  CPs (9) -- leer M.<nombre> | recalcular M.<nombre>_ | control M.CP.<nombre> :
    bvh         [blob BVH: 1,396 elems, aabb]
                  eventos: transform->handler, changeCoords->handler, changeConnectivity->invalida
    boundary    [0x2 double]
                  eventos: changeConnectivity->invalida
    bbox        (sin calcular)
                  eventos: changeCoords->invalida
    ...
```

La convención por caso hace los tres espacios **disjuntos por construcción**:
una CP no puede pisar una propiedad ni un método (empiezan distinto), y
`DefineCP` lo exige (nombre en minúscula, ni contadores ni alias). Así una CP
puede llamarse casi **como quieras**:

```matlab
M = M.DefineCP( 'nodes' , @(m) m.nV , 'nodecount' , [] );
M.V              % la propiedad real (coordenadas)
M.nodes          % tu CP (lectura perezosa)
M.nodes_         % tu CP, recalculada a la fuerza
```

**Ojo — formas funcionales**: los métodos son Capitalized, así que `plot(M)` o
`transform(M,T)` **ya no despachan a la clase** — caen a las funciones del
path (`transform(M,T)` "medio funciona" vía los alias de lectura pero
**devuelve un struct**, perdiendo la clase y la cache). Usa siempre la forma
punto: `M.Plot()`, `M.Transform(T)`.

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
| `M = M.Transform( T )` | `transform(T)` + `changeCoords` |

Y cada CP declara, para cada evento, una de tres reacciones:

* **no declarado** → el valor **sobrevive** intacto (el evento no le afecta);
* **declarado con `[]`** → el valor **cae** (invalidar, recomputar al pedirse);
* **declarado con un handler** → el valor queda **pendiente**: el handler lo
  actualizará *barato* en el próximo acceso (§4).

El ejemplo canónico es el `bvh` de fábrica, que es pura declaración:

```matlab
% (así está registrado internamente)
M = M.DefineCP( 'bvh' , @(m) BVH( ToStruct(m) ) , ...
      'transform'          , @(B,m,T) BVH( B , T )            , ... % plegado O(1)
      'changeCoords'       , @(B,m)   BVH( B , ToStruct(m) )  , ... % refit O(n)
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
M.V = X1;        % [msh] EVENT {changeCoords}: pendientes {bvh} | caen {surfCent} | ...
M.V = X2;        % [msh] EVENT {changeCoords}: pendientes {bvh} | ...
M.V = X100;      % ... 100 ediciones: coste ~0
B = M.bvh;       % [msh] RPLAY 'bvh' -> sync absoluto via changeCoords en 4.4 ms
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
  era una semejanza y el plegado del bvh da `BVH:notSimilarity`), se cae al
  sync absoluto (refit); si ese también falla, se recomputa desde cero. Nunca
  verificas qué tipo de T tienes — la cascada lo descubre sola.
* El **sufijo `_` se salta todo esto**: `M.bvh_` no consulta ni el valor ni el
  log — recomputa desde cero sí o sí. Es la herramienta de "no me fío / quiero
  uno fresco", no la de uso diario.

Nota de diseño: esto generaliza (y sustituye) el viejo truco de las "dos
claves" del BVH — el valor pendiente con su log *es* la jerarquía
superviviente esperando el refit.

---

## 5. Las CPs de fábrica

Toda malla nace con estas definiciones (todas sobreescribibles, §6):

| Nombre | Qué es | Eventos declarados |
|---|---|---|
| `bvh` | blob de aceleración (`BVH`, aabb) | `transform`→plegado, `changeCoords`→refit, `changeConnectivity`→invalida |
| `boundary` | facetas del borde abierto (`MeshBoundary`; matriz cruda, vacía si es cerrada) | `changeConnectivity`→invalida |
| `edges` | aristas únicas (`meshEdges`) | `changeConnectivity`→invalida |
| `esup` / `psup` | adyacencias por punto (sparse) | `changeConnectivity`, `changeNodeCount`→invalida |
| `esue` | elems vecinos por elem (sparse) | `changeConnectivity`→invalida |
| `bbox` | `[min ; max]` (`meshBB`) | `changeCoords`→invalida |
| `surfCent` | `[área, centroide]` (una llamada a `meshSurface` sirve ambos) | `changeCoords`, `changeConnectivity`→invalida |
| `triNormals` | normales por cara (`meshNormals`) | `changeCoords`, `changeConnectivity`→invalida, **`transform`→rotación O(n)** |

* `M.surfCent` devuelve `[área, centroide]` en una fila (una llamada a
  `meshSurface` sirve los dos valores): `sc = M.surfCent; a = sc(1);
  c = sc(2:4);`
* `M.triNormals` es la normal **canónica** por cara (cross products crudos).
  El atributo `NORMALS` del usuario (p.ej. normales suavizadas con
  `meshNormals(M,k)`, o medidas) es un dato independiente: viaja en
  `ToStruct()`, las funciones legadas lo honran, y `Transform()` lo rota —
  pero no se mezcla con la CP.
* Normales por nodo: **parametrizadas** (`meshNormals(M,'a'|'u'|'b'|...)` — no
  hay "el valor" que cachear), así que no hay entrada de fábrica: llama a
  `meshNormals(ToStruct(M),'angle')` directamente, o define la tuya con el
  método fijado (`DefineCP('nAngle', @(m) meshNormals(ToStruct(m),'angle'), ...)`).

El criterio general (cache vs atributo): **la cache es una consecuencia de la
malla — siempre correcta o inexistente; un atributo es una afirmación del
usuario — sobrevive a las ediciones y puede quedarse rancio en silencio.**

---

## 6. Define las tuyas: `DefineCP`

```matlab
M = M.DefineCP( nombre , @(m) ... [, evento , handler|[] , ...] );
M = M.RemoveCP( nombre );
```

La `computeFcn` recibe **el msh** (usa `ToStruct(m)` para las funciones
legadas, o cualquier otra CP: componen). Los nombres de evento admiten
alias (`'coords'`, `'connectivity'`, `'nodecount'`, ...). Registrar **no
ejecuta nada** — solo enseña a la malla qué hacer cuando se lo pidas. El
nombre debe empezar en **minúscula** (§2); `DefineCP` valida y rechaza
contadores y alias reservados.

```matlab
% aristas del borde, solo sensible a la topología
M = M.DefineCP( 'freeEdges' , @(m) MeshBoundary( m.F ) , ...
                'connectivity' , [] );

% una métrica de calidad por cara (meshQuality parametrizada -> una entrada por métrica)
M = M.DefineCP( 'aspectRatio' , ...
        @(m) meshQuality( ToStruct(m) , 'aspectratio' ) , ...
        'coords' , [] , 'connectivity' , [] );

% con handler de transform: el centroide de vértices se transforma EXACTO bajo
% cualquier afín (la media conmuta con T) -> nunca hace falta recomputarlo
M = M.DefineCP( 'centro' , @(m) mean( m.V , 1 ) , ...
        'coords' , [] , 'nodecount' , [] , ...
        'transform' , @(c,m,T) [ c , 1 ] * T(1:3,:).' );

v = M.aspectRatio;             % primer acceso: computa
v = M.aspectRatio;             % ahora es un hit
M2 = M.Transform( T );  M2.centro;    % replay incremental: sin recomputar
```

**Ejemplo completo y verificado**: `MESH\meshQuality_as_cachedProps.m` registra
*todas* las métricas de `meshQuality` (una CP por métrica, según el
celltype) con la matemática de actualización de cada una — las adimensionales
son invariantes ante semejanzas, longitudes/áreas escalan `s`/`s²`, y
`volume`/`signedvolume`/`orientation` se actualizan **exactas ante cualquier
afín** vía `det(A)`:

```matlab
M = meshQuality_as_cachedProps( M );              % todas las de su celltype
M = meshQuality_as_cachedProps( M , 'area' , 'aspectratio' );   % o un subconjunto
q = M.aspectratio;                 % computa y cachea
M2 = M.Transform( T );             % semejanza: TODAS quedan pendientes...
q2 = M2.aspectratio;               % ...y se actualizan sin recomputar
```

**Sobreescribir una de fábrica** reemplaza la definición completa (y descarta
el valor viejo, que ya no corresponde) — redeclara los handlers que quieras
conservar:

```matlab
M = M.DefineCP( 'bvh' , @(m) BVH( ToStruct(m) , [2 8] , 'sphere' ) , ...
      'transform'    , @(B,m,T) BVH( B , T ) , ...
      'changeCoords' , @(B,m)   BVH( B , ToStruct(m) ) , ...
      'connectivity' , [] );
% M.bvh, M.ClosestElement, M.IntersectRay usan ahora el blob de esferas
```

Como `msh` es una value class, `DefineCP` **devuelve** la malla modificada —
sin el `M =` no pasa nada (así funciona el lenguaje; las variables auxiliares
son responsabilidad del programador).

Detalles: sin detección de ciclos (una computeFcn que se pida a sí misma
recursa); los function handles anónimos serializan en el `.mat` con la letra
pequeña de siempre (capturan su workspace).

---

## 7. El proxy `M.CP`: el plano de control

| Expresión | Hace | Devuelve |
|---|---|---|
| `M.CP` | tabla de definiciones y estados | vista |
| `M.CP.bvh` | el valor (computa/replay si hace falta, como `M.bvh`) | valor |
| `M.CP.bvh.frame` | indexa **dentro** del valor | `blob.frame` |
| `M.CP.bvh.delete` | borra el **valor** (la definición queda) | — |
| `M = M.CP.bvh.removeProp` | borra definición **y** valor | msh nuevo |
| `M = M.CP.bvh.set( x )` | siembra un valor a mano ("dangerous, but...") | msh nuevo |
| `M.CP.bvh.changeCoords` | el handler del evento (`[]` = invalida) | handle |
| `M.CP.bvh.changeCoords( B , M )` | lo ejecuta con tus argumentos | valor |

Las **operaciones** viven SOLO en el proxy: `M.bvh.delete` NO borra nada
(indexa `delete` dentro del blob → error). La lectura desnuda y el sufijo son
las puertas de uso diario; `M.CP` es la sala de máquinas.

**Semántica de compartición** — la distinción importante:

* `.delete` es un *statement* que muta el handle **compartido**: las copias
  hermanas también pierden ese valor. Benigno — como mucho, alguien recomputa.
* `M.<nombre>_` (recalcular) también escribe al handle **compartido** — misma
  regla benigna: las que comparten cache no han divergido (el COW separa al
  editar), así que el valor fresco les vale a todas.
* `.set(x)` es **conservador (COW)**: devuelve un msh nuevo con su propia
  cache sembrada; las hermanas ni se enteran. Sigue siendo peligroso en el
  sentido que importa — nadie verifica que `x` sea correcto. (Y ojo: un
  `M.<nombre>_` posterior pisa el valor sembrado con el recomputado.)
* `.removeProp` y `DefineCP` tocan la *definición* (que vive en el
  valor): devuelven un msh nuevo, hermanas intactas.

Colisiones: tras el nombre, `delete`/`removeProp`/`set` y los nombres de
evento son operaciones reservadas del proxy; cualquier otro nombre se reenvía
como indexación del valor. Si un valor tuviera un campo llamado como una
operación, sácalo en dos pasos (`B = M.bvh; B.delete` ya indexa normal).

---

## 8. El mecanismo valor + handle + copy-on-write

`msh` es una value class (copias independientes), pero los **valores**
cacheados viven en un handle compartido (COW = *copy-on-write*: se comparte
hasta que alguien escribe). La combinación da la semántica del ejemplo
fundacional:

```matlab
M1 = msh( V , F );
e  = M1.esup;            % CP: se calcula y queda en la cache COMPARTIDA
M2 = M1;                 % copia: mismos datos, misma cache
M3 = M1.Transform( RT ); % rotar no toca la topología...
e3 = M3.esup;            % ...así que M3 YA TIENE esup (hit, gratis)

M1.F(1:10,:) = [];       % edición de conectividad en M1
% -> M1 recibe una cache NUEVA donde esup/boundary/bvh... han caído
% -> M2 y M3 CONSERVAN la suya intacta
M2.esup;                 % sigue siendo un hit
```

Cada evento hace **copy-on-write del handle**: la instancia editada evoluciona
a un handle nuevo llevándose los supervivientes y los pendientes; las hermanas
conservan el suyo. Resolver un pendiente, depositar un MISS o un recalculo
`_` sí mutan el handle compartido — es benigno (mismo valor correcto para
todas).

---

## 9. `Transform()`

```matlab
M2 = M.Transform( T );     % T homogénea 4×4 / 3×4 / 3×3 (2D: 3×3 homogénea 2D)
```

Dispara el evento semántico `transform(T)` — la única edición que *sabe qué
está pasando* (el setter de `V` es ciego: no puede distinguir una traslación
de una deformación). Con las definiciones de fábrica:

* `bvh` → pendiente de **plegado O(1)** (~14 µs) si T es semejanza; si no, la
  cascada lo baja a refit O(n).
* `triNormals` cacheadas → pendiente de **rotación O(n)** de las filas.
* Los atributos `NORMALS` del usuario se rotan aquí mismo (vía `tools\transform`).
* `bbox`, `surfCent` → caen (declaran `changeCoords` sin handler).
* `boundary`, `edges`, `esup`... → **sobreviven** (ni se enteran).

Si trasladas "a mano" (`M.V = M.V + t`) pagas la versión ciega:
`changeCoords`, refit en vez de plegado. Usa `Transform` para transformar.

```matlab
B1 = M.bvh;
M2 = M.Transform( [0.9*RT , [1;2;3] ; 0 0 0 1] );
B2 = M2.bvh;                       % replay: plegado O(1)
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
  remapear de verdad usa `RemoveNodes`/`RemoveFaces`/`Tidy`.

---

## 11. Queries

```matlab
[e,cp,d,bc,F] = M.ClosestElement( P );          % elemento más cercano
[e,cp,d]      = M.ClosestElement( P , Dmax );   % con radio de corte
[xyz,c,t,rid] = M.IntersectRay( rays );         % 'first' por defecto
[xyz,c,t,rid] = M.IntersectRay( rays , 'all' ); % 'first'|'last'|'all'|'any'
```

Usan (y rellenan de paso) la CP `bvh` — la que esté definida: si la
sobreescribiste con esferas u otro leaf size, es la que trabaja. Tras un
`Transform` el blob llega plegado; tras una deformación, refitado. **Nunca hay
blob rancio**: la clase intercepta todas las ediciones. Detalles del motor
(salidas `bc`/`F`, semántica de `Dmax`, modos de rayos, rendimiento) en
`BVH_TUTORIAL.md`.

Visualizador del blob (↑/↓ profundidad, `a` todo, `l` hojas, `f` marco):

```matlab
M.PlotBVH();
```

---

## 12. Atributos por nodo y por cara

```matlab
M = M.AddField( 'xyzTemp'  , temp  );   % por nodo  (prefijo explícito)
M = M.AddField( 'triLabel' , label );   % por cara
M = M.AddField( 'Color'    , c     );   % sin prefijo: se infiere por filas
                                        % (error si nV == nF: sé explícito)
v = M.GetField( 'xyzTemp' );            % con o sin prefijo
M.HasField( 'triLabel' )                % true
L = M.FieldNames();                     % L.node = {'xyzTemp',…}, L.face = …
M = M.RmField( 'xyzTemp' );
```

* Filas validadas contra `nV`/`nF` al añadir; reconciliación crop/NaN-pad al
  cambiar tamaños (§10). Cada atributo conserva su propio tipo.
* `NORMALS` es especial solo en una cosa: `Transform()` **rota** los atributos
  `xyzNORMALS`/`triNORMALS` (vía `tools\transform`) en vez de dejarlos rancios.
  Son datos tuyos, independientes de la CP `triNormals` (§5).
* En `ToStruct()` los atributos salen con su prefijo legado (`xyz*`/`tri*`).

---

## 13. `VIZ`, `INFO`, `Plot`, textura

```matlab
M.VIZ.FaceColor = [0.8 0.2 0.2];       % preferencias persistentes de plot
M.VIZ.EdgeColor = 'none';
M.INFO.paciente = 'ID-042';            % metadatos libres
M.INFO.texture  = imread('piel.png');  % la textura es un metadato más

M.Plot();                              % plotMESH con las preferencias
M.Plot('FaceAlpha',0.5);               % precedencia: defaults < VIZ < args
```

`VIZ` e `INFO` son datos normales: se copian con la malla, se guardan en el
`.mat`, y no participan en la cache. Los campos desconocidos de un struct de
entrada acaban en `INFO`; `INFO.texture` sale como `.texture` en `ToStruct()`
(y un `.texture` de entrada entra ahí).

---

## 14. Puente con el toolbox legado

```matlab
S = M.ToStruct();            % struct .xyz/.tri/.xyzF*/.triA*/.texture
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
M = M.Tidy( … );                 % MeshTidy
M = M.RemoveFaces( idx );        % MeshRemoveFaces
M = M.RemoveNodes( idx );        % MeshRemoveNodes (remapea caras y atributos)
M = M.Append( M2 , S3 , … );     % MeshAppend (acepta msh y structs mezclados)
```

---

## 15. `save` / `load`

```matlab
save malla.mat M
L = load('malla.mat');
L.M.ClosestElement(P);     % funciona: el bvh se recomputa perezoso
```

Las **definiciones** (el registro, con sus handlers) viajan en el `.mat`; los
**valores** cacheados no (el handle es Transient): un `.mat` de `msh` pesa lo
que pesan sus datos y todo se recalcula al primer uso. Para transportar un BVH
precalculado, guárdalo aparte (`B = M.bvh` es un struct-valor serializable).

---

## 16. `DEBUG`: ver la cadena de procesos

```matlab
M.DEBUG = true;                    % o msh(V,F,'DEBUG',true)
```

Narra por consola cada paso — la tabla de eventos de §3, en vivo:

```
[msh] MISS  'boundary' -> calculado en 1.84 ms
[msh] HIT   'boundary'
[msh] RECMP 'boundary' -> recalculado a la fuerza en 1.79 ms
[msh] SET   V 700x3 -> 700x3
[msh] EVENT {changeCoords}: pendientes {bvh} | caen {surfCent} | sobreviven {boundary}
[msh] RPLAY 'bvh' -> sync absoluto via changeCoords en 4.42 ms
[msh] TRANS Transform() sobre 700 nodos
[msh] EVENT {transform+changeCoords}: pendientes {bvh, triNormals} | caen {...} | ...
[msh] RPLAY 'bvh' -> 1 transform(s) incremental(es) en 0.05 ms
[msh] QUERY ClosestElement: 500 puntos
[msh] QUERY ClosestElement resuelta en 0.62 ms
[msh] DEF   CP 'aspectRatio' definida (eventos: changeCoords, changeConnectivity)
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
    VIZ:  FaceColor = [0.8 0.2 0.2]
    INFO: paciente = 'ID-042'
    CPs:
      bvh       [blob BVH: 1,401 elems, aabb]
      boundary  [42x2 double]
      surfCent  (pendiente: changeCoords)
    definidas sin calcular: bbox, edges, esue, esup, psup, triNormals
```

* Cabecera con los conteos con **separadores de miles** (`2,000,000 vertices`
  se distingue de `200,000` de un vistazo), desglose de **tipos de celda**
  (mixtas incluidas) y la descripción **planar/flat**: `(3D, flat (z = 0))` si
  todas las z son 0 exactas, `(3D, planar: z = 2.5)` si el plano es
  axis-aligned, `(3D, planar)` si es un plano cualquiera. La dimensión sigue
  siendo el número de columnas.
* Atributos con dimensiones y clase, `VIZ`/`INFO` como pares `nombre = valor`
  (la textura aparece dentro de `INFO`).
* Entradas cacheadas con su **valor si es escalar o vector de ≤ 3 números**,
  su estado pendiente (con los eventos del log), o `[tamaño clase]`.
* `M.CP` da la vista completa con los eventos declarados por entrada.

---

## 18. Limitaciones y pendientes

* Sin detección de ciclos entre CPs (una computeFcn que se pida a sí
  misma recursa).
* Eventos semánticos futuros: `RemoveFaces(idx)` / `RemoveNodes(map)` con los
  mapas de índices que ya devuelven `MeshTidy`/`MeshRemove*` — permitirían
  updates quirúrgicos de las adyacencias (borrar filas de `esup` en vez de
  recomputarla).
* Las delegaciones (`Tidy`, `Append`, …) devuelven un `msh` con cache fresca
  y **registro de fábrica** (no heredan tus `DefineCP`).
* El acceso a CPs paga un pequeño peaje de `subsref` (µs); los métodos
  y propiedades reales no se interceptan en la práctica.
* Las formas funcionales no despachan a la clase (§2): `plot(M)` /
  `transform(M,T)` caen al path — usa `M.Plot()` / `M.Transform(T)`.
* Los `.mat` guardados con el API anterior cargan sus datos, pero el registro
  serializado trae los nombres/handlers viejos (`BVH`, `toStruct`) — para
  mallas antiguas: reconstruye con `M = msh( struct viejo )`.
* No hay aritmética de mallas sobrecargada (`M1 + M2`…): usa `Append`.
```
