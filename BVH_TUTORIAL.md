# Tutorial: motor de queries geométricas (`BVH` · `bvhClosestElement` · `bvhIntersectRay` · `plotBVH`)

Motor in-house de **elemento más cercano** y **intersección rayo-malla** sobre
una única estructura de aceleración cacheable, transformable y deformable.
Sustituye a `vtkClosestElement` / `distanceFrom` (y cubre el rol de
`tsearchn` en tets) con una diferencia clave: la estructura de búsqueda es un
**valor** — se guarda, se copia, se transporta con la malla, y puedes tener
tantas como mallas, a la vez.

> Documentos hermanos: `msh_QUERY_ENGINE_DESIGN.md` (diseño e internals),
> `msh_DESIGN.md` (la futura clase `msh`, que cacheará todo esto sola).

---

## 0. Instalación (una sola vez)

Los tres MEX se compilan una vez (MSVC configurado vía `mex -setup`) y son
**obligatorios** — el motor no tiene versiones MATLAB-puro:

```matlab
cd C:\repos\msh\BVH
mex BVH_mx.cpp                                        % build SAH + refit
mex COMPFLAGS="$COMPFLAGS /openmp" bvhClosestElement_mx.cpp
mex COMPFLAGS="$COMPFLAGS /openmp" bvhIntersectRay_mx.cpp
```

**Hilos**: los MEX paralelizan siguiendo `maxNumCompThreads` de tu sesión.
En `matlab -singleCompThread` (o tras `maxNumCompThreads(1)`) corren en serie;
`maxNumCompThreads(k)` los limita a `k` hilos. Perfecto para workflows de
high-throughput: la sesión manda.

Verificación:

```matlab
test_BVH                    % closest-element completo (15 bloques)
test_bvhClosestElement_tets    % semántica tsearchn en tets
test_bvhIntersectRay           % rayos vs IntersectSurfaceRay (oráculo)
```

---

## 1. El formato de malla

Una malla es un `struct` con dos campos:

| Campo | Contenido |
|---|---|
| `.xyz` | vértices, `nV × 3` (3D) o `nV × 2` (2D genuino) |
| `.tri` | conectividad, `nF × k` con `k` = 1 (puntos), 2 (segmentos), 3 (triángulos), 4 (**tetraedros**) |

**4 nodos = tetraedro, siempre** (los quads no existen en este motor).

**Mallas mixtas**: filas rellenadas con `0` **al final** (`[a b 0]` = segmento
en una matriz de 3 columnas). El tipo de cada fila lo da su número de nodos
no-cero.

**Regla de dimensionalidad**: la dimensión es `size(xyz,2)` y **nunca** se
infiere del contenido. Una malla 3D con todas las z = 0 sigue siendo 3D —
puede deformarse fuera del plano sin invalidar nada. Una malla 2D de verdad
tiene 2 columnas.

Mallas de ejemplo del tutorial:

```matlab
% Superficie de triángulos (esfera unidad, ~1.4k triángulos)
V = randn(700,3);  V = V./sqrt(sum(V.^2,2));
Mtri = struct('xyz', V, 'tri', convhulln(V));

% Polilínea de segmentos (arco 3D)
s = linspace(0, 3*pi/2, 200).';
Mseg = struct('xyz', [cos(s), sin(s), 0.2*s], 'tri', [(1:199).', (2:200).']);

% Volumen de tetraedros (nube en el cubo unidad)
W = rand(300,3);
Mtet = struct('xyz', W, 'tri', delaunayn(W));

% Nube de puntos (1 nodo por fila)
Mpts = struct('xyz', rand(500,3), 'tri', (1:500).');

% Mixta: triángulos + segmentos (0-padded)
Mmix = struct('xyz', V, 'tri', [Mtri.tri ; [(1:20).', (21:40).', zeros(20,1)]]);

% 2D genuina (disco triangulado, vértices de DOS columnas)
X2  = rand(200,2);
M2d = struct('xyz', X2, 'tri', delaunayn(X2));
```

---

## 2. Arranque en 30 segundos

```matlab
P = randn(1000,3);                            % puntos de consulta

[e, cp, d, bc, F] = bvhClosestElement(Mtri, P);
%   e  : fila de Mtri.tri del elemento más cercano
%   cp : punto más cercano SOBRE la malla (1000 x 3)
%   d  : distancia euclídea
%   bc : coordenadas baricéntricas de cp en su elemento (robustas, ver §4)
%   F  : dónde cae cp (vértice/arista/cara/interior) y si es borde de la malla

[xyz, cell, t] = bvhIntersectRay(Mtri, [0 0 -5 ; 0 0 5]);   % un rayo
```

---

## 3. El blob `B`: construir una vez, consultar mil veces

`BVH` construye la estructura de aceleración — el **blob** — y las queries
la reutilizan. La forma recomendada empaqueta malla y blob juntos: *"sobre el
primer argumento, busca lo más cercano al segundo"*:

```matlab
B = BVH(Mtri);                                   % una vez (~ms)
S = {Mtri, B};                                      % la "malla acelerada"

[e1, cp1, d1] = bvhClosestElement(S, P1);           % queries baratas
[e2, cp2, d2] = bvhClosestElement(S, P2);
[xyz, cell, t] = bvhIntersectRay(S, rays);          % el MISMO blob sirve rayos
% (la forma clásica bvhClosestElement(M, P, B) sigue funcionando)
```

Propiedades que lo diferencian del locator persistente de VTK:

**a) Tantas mallas como quieras, a la vez** — cada malla lleva su blob, sin
estado global:

```matlab
Stri = {Mtri, BVH(Mtri)};   Sseg = {Mseg, BVH(Mseg)};   Stet = {Mtet, BVH(Mtet)};
[~,~,dA] = bvhClosestElement(Stri, P);
[~,~,dB] = bvhClosestElement(Sseg, P);
[~,~,dC] = bvhClosestElement(Stet, P);
[~, quien] = min([dA, dB, dC], [], 2);              % ¿a qué malla está más cerca?
```

**b) Es un valor** — `save`/`load` con resultados idénticos, viaja a workers
de `parfor`, se copia con una asignación.

**c) Se auto-protege** — cada query hace un *chequeo de frescura*: compara 4
vértices testigo de `M` contra la reconstrucción `frame(B.X)` (coste:
~microsegundos, independiente del tamaño). Si editaste la malla sin actualizar
el blob, la query **da error** (`bvhClosestElement:staleBVH`) — a propósito:
como `B` es un valor, un rebuild interno se descartaría y se repetiría
silenciosamente en cada llamada posterior. Recuperación explícita:

```matlab
try
    [e, cp, d] = bvhClosestElement({M, B}, P);
catch
    B = BVH(B, M);                 % refit (misma conectividad) o BVH(M)
    [e, cp, d] = bvhClosestElement({M, B}, P);
end
```

No es un hash (una edición quirúrgica que esquive los 4 testigos no se
detecta); es una red de seguridad para el uso suelto — cuando el blob viva en
la caché de la clase `msh`, toda edición queda interceptada por los setters y
el chequeo será redundante.

### Opciones de construcción

```matlab
B = BVH(M);                       % por defecto: SAH, hojas [2 16], volumen AABB
B = BVH(M, [4 32]);               % hojas SAH adaptativas: minLeaf=4, maxLeaf=32
                                     %   (n <= minLeaf siempre es hoja; el SAH puede
                                     %    agrupar hasta maxLeaf si dividir no compensa)
B = BVH(M, 8);                    % un escalar s equivale a [s s] (hoja fija)
B = BVH(M, Inf);                  % una única hoja = fuerza bruta (referencias)
```

**Batería de volúmenes** — todos bajo la MISMA partición SAH (la jerarquía es
idéntica; solo cambia el "cage" de cada slot):

```matlab
B = BVH(M, [], 'aabb');           % (defecto) cajas en el marco de construcción
B = BVH(M, [], 'sphere');         % esferas envolventes
B = BVH(M, [], 'obb');            % cajas ORIENTADAS: cada slot PCA-alineado
                                     %   a SUS nodos contenidos
B = BVH(M, [], 'kdop');           % 14-DOP: AABB + los 4 pares de planos
                                     %   diagonales (1,±1,±1)
B = BVH(M, [], 'rss');            % rectángulo-barrido-por-esfera (estilo PQP):
                                     %   rect PCA por slot + radio
B = BVH(M, [], 'lss');            % cápsulas (alias 'capsule'): segmento PCA
                                     %   por slot + radio
B = BVH(M, [], 'aabb','noframe'); % AABB alineado al MUNDO (sin marco de
                                     %   centroide/PCA) -- pieza de comparación
```

Guía medida (52k tris / 26k tets / hélice de 30k segmentos; monohilo, µs/pt):

| volumen | superficie (lejos) | tets | wireframe (hélice) |
|---|---|---|---|
| `aabb` (defecto) | 6.5 | **2.4** | 1.07 |
| `sphere` | 9.8 | 9.0 | — |
| `obb` | 2.3 | 3.9 | 0.91 |
| `kdop` | 8.4 | 3.0 | — |
| **`rss`** | **2.0** (×3.3 sobre aabb) | 4.4 | **0.75** |
| `lss` | 10.5 | 6.7 | 1.06 |

Lectura: los volúmenes **orientados/barridos (rss, obb) conquistan el caso
lejano** en superficies y curvas (abrazan los clusters localmente planos o
tubulares en su propia orientación, encogiendo el "casquete tangente" —
lo que dábamos por física resultó ser holgura del AABB); el **AABB sigue
mandando en tets** y en queries cerca de la superficie; el k-DOP no paga sus
slabs extra; la cápsula decepciona incluso en su terreno (una hélice curva
rompe la premisa de segmentos rectos por slot — el rect del RSS la abraza
mejor). Todos dan resultados idénticos (asertado): elegir volumen solo cambia
velocidad. En rayos, `rss`/`lss` usan tests conservadores de caja expandida
(idénticos resultados, poda algo menor) — para rayos usa `aabb`.

### El marco de construcción: centrado + PCA

El blob guarda su geometría en un **marco propio**: siempre centrado en el
centroide, y **alineado a los ejes principales (PCA)** cuando la nube es
anisótropa (razón de ejes ≥ 2). Dos consecuencias importantes:

- **Un blob elíptico en diagonal SÍ alinea sus cajas**: los AABBs viven en el
  marco material de la malla, no en los ejes del mundo. Y como las rotaciones
  posteriores se pliegan en el marco (§6), la calidad nunca se degrada por
  orientación.
- **Precisión float garantizada lejos del origen**: los bounds de nodo son
  `single` redondeados hacia fuera; el redondeo es *relativo* (1 ulp a la
  magnitud del valor). Sin centrado, una malla de tamaño 1 situada a 1e8 del
  origen tendría cajas infladas ±8 unidades (ulp de single a 1e8) — correctas
  pero inútiles. Centrada, las coordenadas internas son del tamaño de la malla
  y el inflado es ~1e-7 relativo. Para esferas, el radio se redondea hacia
  ARRIBA e incluye el error de redondeo del centro (`|c - single(c)|`), así
  que la magnitud está contemplada también ahí. Testeado con la malla a 1e7.

### Visualizar el blob

```matlab
plotBVH(B)          % la geometría del blob + los volúmenes de nodo
plotBVH(B, Mtri)    % ídem, dibujando tu malla
```

Interactivo: **↑/↓** cambia la profundidad del árbol; **a** superpone todos
los niveles (color por profundidad); **l** solo hojas; **f** muestra los ejes
del marco (verás el centrado y la rotación PCA); **r** resetea la vista. Hojas
en rojo, nodos internos en azul; el título muestra profundidad, nº de nodos y
slots, tipo de volumen y tamaños de hoja. Todo se dibuja en coordenadas de
mundo (a través del marco): ves exactamente lo que ven las queries — pruébalo
tras un `BVH(B,T)` o un refit.

---

## 4. `bvhClosestElement` a fondo

```matlab
[e, cp, d, bc, F] = bvhClosestElement({M,B}, P, Dmax)     % forma recomendada
[e, cp, d, bc, F] = bvhClosestElement(M, P, B, Dmax)      % forma clásica
```

- `e`: fila de `M.tri` del ganador · `cp`: punto más cercano · `d`: distancia.
- `bc`: baricéntricas de `cp` en su elemento, **robustas**: calculadas con la
  forma de productos cruzados (como `vtkClosestElement/calcular_barycentric` —
  la forma de Gram cancela catastróficamente en *slivers*; ésta es exacta
  hasta aspecto 1e10), recortadas a `[0,1]` y renormalizadas a suma 1.
  Testeado con un triángulo de aspecto 1e8.
- `F` (5º output): **clasificación del punto más cercano**:

  | `F.type` | significado |
  |---|---|
  | 0 | sin respuesta (fuera de `Dmax` / punto no finito) |
  | 1 | sobre un VÉRTICE del elemento |
  | 2 | sobre una ARISTA (interior de ella) |
  | 3 | interior de una CARA (triángulo, o cara de tet) |
  | 4 | INTERIOR del elemento (dentro de un tet) |

  `F.onBoundary`: `true` si `cp` cae sobre el **borde abierto** de la malla —
  para superficies de triángulos, una arista de `MeshBoundary` o uno de sus
  vértices; para wireframes, un extremo libre (nodo de grado 1). Siempre
  `false` en mallas cerradas, tets, nubes y mixtas. El patrón de `bc > tol`
  te dice *cuál* vértice/arista es. (Es la misma información con la que
  `distanceFrom` negaba `d` en el borde; aquí sale explícita.)

  **Coste**: `bc` y `F` se calculan **solo si los pides** (4º/5º output).
  `F.onBoundary` ejecuta `MeshBoundary(M.tri)` en cada llamada con 5 outputs
  (O(nE·log nE), ~ms) — para llamadas repetidas sobre la misma malla es el
  candidato obvio a cachear (lo hará la clase `msh`; `distanceFrom` lo
  memoiza con `persistent`). La clasificación `F.type` en sí es vectorizada
  y despreciable.

### Sin respuesta (convención)

Punto sin elemento dentro de `Dmax`: **`e = 0`, `d = Inf`, `cp`/`bc` = NaN,
`F.type = 0`**. Punto de query no finito: igual pero `d = NaN`.

```matlab
[e, cp, d] = bvhClosestElement(S, P, 0.1);    % radio de búsqueda
cerca = e > 0;                                % o isfinite(d)
```

`Dmax` hace *early-exit* real: un punto más lejos que `Dmax` cuesta ~una
visita de nodo (medido: 0.11 µs/punto, ×65 sobre la búsqueda completa).

`Dmax` también acepta un **vector nP×1**: cota/radio **por punto**, sembrada
como best-so-far inicial de cada query. Sirve para radios de búsqueda
heterogéneos o para sembrar cotas superiores de heurísticas (si la cota es
*alcanzable*, ínflala `(1+1e-9)` para que el elemento que la produce se
re-encuentre; un miss bajo cota alcanzable significa "el candidato era el
ganador"). Nota empírica (`bench_vertexSeed`): sembrar la cota del
vértice-más-cercano / 1-ring **no acelera** este motor — el warm-start y la
travesía ordenada ya auto-siembran tras el primer descenso, y el coste
far-field está dominado por la cáscara tangente, que ninguna cota superior
puede podar.

### Por tipo de celda

```matlab
[e,cp,d,bc]   = bvhClosestElement({Mtri,Btri}, P);   % triángulos: bc = [u v w]
[e,cp,d,bc]   = bvhClosestElement({Mseg,Bseg}, P);   % segmentos: bc = [1-t, t]
[e,cp,d]      = bvhClosestElement(Mpts, P);          % nube: vecino más cercano
[e,cp,d,bc,F] = bvhClosestElement({Mtet,Btet}, P);   % tets: interior -> d=0,
                                                     %   F.type=4, bc del punto
[e,cp,d]      = bvhClosestElement(Mmix, P);          % mixta: cada fila su métrica
[e,cp,d]      = bvhClosestElement(M2d, rand(400,2)); % 2D: cp sale con z=0 exacto
```

### Reemplazo de `tsearchn` (point-location en tets)

Mismo tet contenedor y mismas baricéntricas (verificado al 100%), ×7–×60 más
rápido, y funciona en mallas de tets **deformadas** (no-Delaunay):

```matlab
[e, cp, d, bc] = bvhClosestElement({Mtet,Btet}, P);
tid = e;  tid(d > 0) = NaN;          % == tsearchn(W, Mtet.tri, P)
```

### Localizador APROXIMADO: `approximateClosestElement`

Misma firma y convenciones, pero heurístico: vértice más cercano (BVH de
puntos) + distancia **exacta** a su abanico de elementos (EsuP, fusionado en
el MEX). Garantías: `d_apx >= d_exacta` **siempre** (cota superior, asertada
en `bench_approximate`); el elemento es el correcto el ~95–99 % de las veces
en mallas razonables, con error acotado por la escala local. Todos los
celltypes; en **nubes de puntos y polilíneas es exacto** en la práctica
(100 % medido).

```matlab
Ba = approximateClosestElement( M );                 % blob propio (~= coste del exacto)
[e,cp,d,bc,F] = approximateClosestElement( {M,Ba} , P [, Dmax] );
% OJO: Dmax corta por la distancia AL VERTICE (etapa 1), no al elemento
```

El MEX lleva su propia artillería: hojas de vértices como bloques SoA con
kernel AVX de 4 puntos (`pt4`), abanicos de triángulos pre-empacados como
bloques PreTri4 barridos con el kernel 4-wide del motor exacto (`fan4`, solo
mallas puras de triángulos), hojas grandes `[32 128]` (barrido medido), Morton
+ warm-start + pool fusionado. Cuándo pagar la imprecisión (medido, 1 hilo,
`bench_approximate`): **far-field ×2.7–6.3** (22→3.6 µs/pt en 200k tris),
**caja media ×1.1–5.1**, near-surface ×1.3–2.1. Evítalo en mallas
**anisótropas/sliver** (hit cae al 69–82 %) y para point-location interior en
tets (hit 83 % en el régimen interior, y además ahí es ×0.85 — más lento: los
abanicos de Delaunay tienen ~20 tets por vértice).

---

## 5. `bvhIntersectRay` a fondo

```matlab
[xyz, cell_id, t, ray_id] = bvhIntersectRay({M,B}, rays, MODE)
```

- `rays`: `2×3` (`[p0; p1]`), `N×6` o `2×3×N`. Impacto = `p0 + t·(p1−p0)`,
  `t` **no acotado** (negativos = detrás de `p0`). Solo celdas triángulo
  (en mixtas el resto se ignora).

```matlab
[xyz, cell, t] = bvhIntersectRay(S, rays, 'first');  % menor t (picking)
[xyz, cell, t] = bvhIntersectRay(S, rays, 'last');   % mayor t (cara de salida)
[~, cell]      = bvhIntersectRay(S, rays, 'any');    % ¿segmento p0->p1 ocluido?
ocluido = cell > 0;                                  %   (early-exit, el más barato)
[xyz, cell, t, rid] = bvhIntersectRay(S, rays, 'all');  % todos, orden (rayo, t)
```

`'any'` acepta impactos con `1e-9 < t < 1-1e-5`: las bandas de guarda dejan
que `p1` esté exactamente SOBRE la superficie sin auto-reportarse.

**Tuning para trabajo intensivo de rayos**: el default de hojas del blob
(`[2 16]`) está optimizado para closest-point; el kernel de rayos 4-wide rinde
mejor con hojas grandes. Medido a 52k tris (best-of-3, monohilo):
`[16 64]` → `first` **0.315 µs/ray vs 0.347 de `IntersectSurfaceRay_mx`**
(×1.10 más rápido) y `any` 0.30 vs 0.29 (paridad). Si tu carga es mayormente
rayos: `B = BVH(M, [16 64])`. Con hilos (`maxNumCompThreads`), ×2.5
adicional.

### Receta: interior/exterior por paridad de cruces

```matlab
lejos = [10 10 10];                                  % punto fuera seguro
rays  = [P, repmat(lejos, size(P,1), 1)];
[~, ~, t, rid] = bvhIntersectRay(S, rays, 'all');
w = t > 0 & t < 1;
n = accumarray(rid(w), 1, [size(P,1), 1]);
dentro = mod(n, 2) == 1;                             % impar = interior
```

*Caveat*: la paridad pura puede fallar en casos rozados — un rayo que pasa
exactamente por una arista compartida puede contarse en las dos caras
incidentes (la tolerancia inclusiva evita "colarse" entre triángulos, al
precio de posibles dobles conteos), y uno tangente toca sin cruzar. Para uso
robusto: varios rayos en direcciones aleatorias + voto por mayoría. Nota:
`meshIsInterior` NO usa paridad de rayos — sus métodos son etapas de poda
(bbox, miniball, elipse, icosaedro, seeds) + closest-element con signo de
pseudonormal ('mesh'/'pseudonormal') + winding numbers + tetgen; la receta de
paridad de arriba es una capacidad nueva de este motor.

---

## 6. Transformaciones: mover la malla sin tocar la estructura

Las semejanzas (rotación + traslación + escala uniforme, reflexiones
incluidas) se **pliegan** en el marco global en **O(1)** — no se toca ni un
nodo (medido: 0.014 ms en 52k triángulos):

```matlab
B = BVH(Mtri);
ang = 0.7;  sc = 1.8;
R = [cos(ang) -sin(ang) 0; sin(ang) cos(ang) 0; 0 0 1] * sc;
T = [R, [10; -2; 5]; 0 0 0 1];

B2 = BVH(B, T);                                   % O(1): pliega el marco
M2 = Mtri;  M2.xyz = M2.xyz * T(1:3,1:3).' + T(1:3,4).';   % la malla sí se mueve

[e, cp, d] = bvhClosestElement({M2, B2}, P);         % exacto
[xyz, cell, t] = bvhIntersectRay({M2, B2}, rays2);   % t invariante ante semejanzas
```

- **Pasa siempre la malla transformada** junto al blob plegado (el chequeo de
  frescura compara ambos). Los plegados **componen**.
- **¿Y si `T(1:3,1:3)` no es una semejanza?** (escala anisótropa, cizalla): no
  se puede plegar para queries de distancia, y el motor **da error**
  (`BVH:notSimilarity`) — tú decides, y lo natural es el **refit** contra
  la malla ya transformada (conserva la jerarquía, O(n)):

  ```matlab
  Mt = transform(M, T);              % la malla siempre se transforma
  try,   Bt = BVH(B, T);          % semejanza: plegado O(1)
  catch, Bt = BVH(B, Mt);         % afín general: refit O(n)
  end
  ```

- **2D**: transforms `2×2`, `2×3` y `3×3` homogénea 2D (elevadas con
  z-escala = escala 2D).

---

## 7. Deformaciones: refit con jerarquía persistente

Cuando cambian las coordenadas (no la conectividad), `BVH(B, M2)` conserva
el árbol Y el marco, y recalcula la geometría en C — 1.8 ms vs 3.3 ms del
rebuild a 52k triángulos:

```matlab
M2 = Mtri;
M2.xyz = Mtri.xyz + 0.03*sin(4*Mtri.xyz(:,[2 3 1]));   % deformación pequeña
B = BVH(B, M2);                     % REFIT: mismo árbol, geometría al día
[e, cp, d] = bvhClosestElement({M2, B}, P);
```

- **Corrección absoluta**: resultados idénticos a un blob reconstruido
  (asertado); con deformaciones grandes solo se degrada la *velocidad*
  (medido ×1.02–1.10 tras warps violentos).
- **¿Y la orientación global en un refit AABB?** El marco se **conserva** —
  es *material*, no espacial: una deformación suave no rota la malla, y las
  cajas se recalculan (ajustadas) en ese marco. Si acumulas muchísima
  deformación y quieres re-estimar el PCA: rebuild (`BVH(M2)`), que
  re-encuadra. Míralo con `plotBVH(B, M2)` antes y después.
- Requiere la **misma conectividad**; si cambió, error explícito → rebuild.

### Patrón de bucle de deformación

```matlab
B = BVH(M);
for it = 1:100
    M.xyz = M.xyz + paso_de_deformacion(M);
    B = BVH(B, M);                          % refit O(n) por iteración
    [e, cp, d] = bvhClosestElement({M, B}, Pobjetivo);
end
```

### Escena: una malla rígida + una deformándose

```matlab
BA = BVH(MA);   BB = BVH(MB);
for frame = 1:nFrames
    T  = pose_del_frame(frame);
    BA = BVH(BA, T);                        % O(1)
    MA.xyz = MA.xyz * T(1:3,1:3).' + T(1:3,4).';
    MB.xyz = MB.xyz + latido(frame, MB);
    BB = BVH(BB, MB);                       % O(n) en C
    [~, cpB, dAB] = bvhClosestElement({MB, BB}, MA.xyz, 5.0);   % con Dmax
end
```

---

## 8. Interpolación con `bc`

`bc` convierte el closest-point en un interpolador de campos por vértice
(la misma mecánica que usa `MeshQuery` con el motor viejo, ahora a nivel de
motor y con `bc` robustas):

```matlab
[e, cp, d, bc] = bvhClosestElement({Mtri, Btri}, P);
idx = Mtri.tri(e, :);
Fp  = bc(:,1).*Fld(idx(:,1),:) + bc(:,2).*Fld(idx(:,2),:) + bc(:,3).*Fld(idx(:,3),:);

% Tets: 4 pesos (interpolación volumétrica exacta en interiores)
% "Pegar" puntos a la superficie:  P_pegados = cp;
```

---

## 9. Rendimiento: chuleta

Monohilo en la máquina de desarrollo (52k tris / 26k tets, 20k queries):

| Query | µs por query |
|---|---|
| closest-point, punto **cerca** de la superficie (workload típico) | **0.8** (×5-6 vs `vtkClosestElement`) |
| closest-point, punto lejano (aabb / rss) | 6.5 / **2.0** (×11-38 vs vtk) |
| closest-point lejano con `Dmax` | **0.11** |
| point-location / closest en **tets** | **2.3–2.6** (×7–60 vs `tsearchn`) |
| rayo `first` / `any` (blob `[16 64]`) | **0.315 / 0.29** (×1.10 / paridad vs `IntersectSurfaceRay_mx`) |
| build / refit / plegado de semejanza | 40 ms / 2 ms / **0.014 ms** |

Reglas:

1. **Reutiliza el blob**; actualízalo barato: semejanza → plegado O(1);
   deformación → refit O(n); conectividad → rebuild.
2. **Hilos**: `maxNumCompThreads` manda (medido con 16 hilos físicos: rayos
   0.16 µs/rayo — ×2.5 sobre el mex especializado, que es monohilo; closest
   lejano ×10). La primera región paralela del proceso paga la creación del
   pool (~decenas de ms, una vez); los hyperthreads rinden peor que los cores
   físicos.
3. Lotes coherentes: el motor ya reordena internamente (Morton + warm-start).

---

## 10. Limitaciones y avisos honestos

- Los MEX son obligatorios (no hay rutas MATLAB-puro).
- `bvhIntersectRay` solo intersecta **triángulos**.
- El **refit** exige la misma conectividad.
- 0-padding con ceros **al final** de la fila.
- El chequeo de frescura son 4 vértices testigo — detecta ediciones normales,
  no sabotajes quirúrgicos.
- `d` es distancia **sin signo**; para el signo: paridad de rayos (§5),
  `meshIsInterior`, o el flag `F.onBoundary` según el caso de uso.
- `F.onBoundary` solo aplica a superficies de triángulos puros y wireframes.

---

## 11. Chuleta de firmas

```matlab
B  = BVH(M)                            % build (SAH [2 16], AABB, marco PCA+centroide)
B  = BVH(M, [minL maxL], 'sphere')     % hojas y volumen a medida
B  = BVH(M, s)                         % escalar == [s s]; Inf = fuerza bruta
B  = BVH(B, T)                         % TRANSFORMAR: semejanza O(1) / afín bake
B  = BVH(B, M2)                        % REFIT a malla deformada (C, O(n))

[e,cp,d,bc,F] = bvhClosestElement({M,B}, P)         % forma recomendada
[e,cp,d,bc,F] = bvhClosestElement({M,B}, P, Dmax)   % radio de búsqueda
[e,cp,d,bc,F] = bvhClosestElement(M, P)             % auto-build
[e,cp,d,bc,F] = bvhClosestElement(M, P, B, Dmax)    % forma clásica

[xyz,cell,t,rid] = bvhIntersectRay({M,B}, rays, MODE)   % first|last|all|any
[xyz,cell,t,rid] = bvhIntersectRay(M, rays)             % auto-build, 'first'

h = plotBVH(B, M)                      % visualizador interactivo del blob
```
