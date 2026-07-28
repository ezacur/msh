/* bvhKernels.h  --  shared closest-point kernels for the BVH query MEXes.
 *
 *   Single source of truth for the geometric primitives (point/segment/triangle
 *   /tetrahedron distance), the REGION-EXACT barycentric emitters, the 4-wide
 *   AVX Ericson kernel, and the Morton bit-spread -- included by
 *   bvhClosestElement_mx.cpp and fanClosestElement_mx.cpp so a fix or a
 *   precision change lives in ONE place (they used to be copied verbatim).
 *
 *   All symbols are `static`/`static inline` -> internal linkage, one private
 *   copy per translation unit (each MEX is its own TU), no ODR concerns.
 *
 * See also bvhClosestElement_mx.cpp, fanClosestElement_mx.cpp.
 */
#ifndef BVH_KERNELS_H
#define BVH_KERNELS_H

#include <cmath>
#include <cstdint>
#include <limits>
#include <immintrin.h>
#if defined(_MSC_VER)
#include <intrin.h>
#endif

static const double INF = std::numeric_limits<double>::infinity();

struct Elem { double cx, cy, cz, r; };     /* element sphere, pkS column */

static inline double sq( double v ) { return v*v; }

/* ------------------------------------------------- predicados de orientacion
 * orient3d: SIGNO de det[ a-d ; b-d ; c-d ] (>0 si a,b,c se ven antihorarios
 * desde d, 0 si los cuatro son coplanares). Es el predicado del que cuelga toda
 * la interseccion malla-malla: si el signo de un vertice respecto al plano del
 * triangulo contrario es EXACTO, entonces "que aristas cruzan" es una decision
 * combinatoria exacta que depende solo de (arista, triangulo) -- nunca de por
 * que triangulo se llego -- y por tanto los dos triangulos que comparten una
 * arista coinciden SIEMPRE. Eso es lo que permite soldar los segmentos por clave
 * canonica (arista padre + parametro) en vez de por tolerancia geometrica.
 *
 * Dos etapas:
 *   1) doble con FILTRO de error dinamico (Shewchuk): si |det| supera la cota,
 *      el signo es cierto y se devuelve. Es el camino rapido y el 99.9% de los
 *      casos, sin coste extra apreciable sobre un det3 normal.
 *   2) si no, se recalcula en DOUBLE-DOUBLE (~32 digitos): los productos son
 *      EXACTOS via twoProd, y las sumas acumulan error ~2^-106 relativo.
 *
 * ALCANCE HONESTO: la etapa 2 resuelve el signo salvo que el determinante
 * verdadero sea menor que ~2^-106 relativo, en cuyo caso se informa 0
 * (coplanar). No es la aritmetica de expansiones completa de Shewchuk, que es
 * exacta tambien para el cero; es ~2^80 veces mas fina que un umbral 1e-8 y
 * suficiente para decidir topologia.
 *
 * COMO SUBIR EL LISTON (si algun dia hace falta el cero BIT-EXACTO): sustituir
 * la etapa 2 por EXPANSIONES DE SHEWCHUK, no por aritmetica racional. El motivo
 * es que un predicado solo necesita el SIGNO de un determinante, y eso se calcula
 * con +, - y x -- sin una sola division. El resultado exacto de esas operaciones
 * sobre doubles es un ENTERO (por una potencia de dos), asi que no hay
 * denominadores que representar: una expansion (suma no evaluada de doubles no
 * solapados, con twoSum/twoProd que ya estan aqui arriba) lo representa exacto,
 * es autocontenida (sin GMP ni bignum) y en el camino filtrado cuesta lo mismo
 * que ahora. Los racionales cobrarian la contabilidad de denominadores a cambio
 * de nada, y ademas CASCADEAN: en cuanto un punto construido alimenta otra
 * decision, los denominadores se multiplican y la longitud crece
 * exponencialmente con la profundidad -- que es justo por lo que el kernel de
 * construcciones exactas de CGAL es lento.
 * Los denominadores solo aparecen de verdad en las CONSTRUCCIONES (el punto de
 * corte, t = num/den), y este diseno las evita por otra via: la identidad de cada
 * punto es COMBINATORIA (ver triTriSegment mas abajo) y las coordenadas se
 * calculan una sola vez al final, solo para salida, sin realimentar ninguna
 * decision. Asi el redondeo afecta a la geometria del resultado pero NUNCA a su
 * topologia. Si algun dia hay que construir sobre puntos construidos (fase de
 * booleanas), la respuesta moderna son los PREDICADOS INDIRECTOS (guardar el
 * punto simbolicamente como "corte de la arista e con el triangulo t" y evaluar
 * los predicados exactamente desde las entradas originales), tampoco racionales. */

static inline void twoSum( double a, double b, double* s, double* e )
{ const double t = a + b; const double z = t - a; *e = ( a - ( t - z ) ) + ( b - z ); *s = t; }

static inline void twoProd( double a, double b, double* p, double* e )
{ const double t = a*b; *e = std::fma( a, b, -t ); *p = t; }   /* exacto: fma */

/* d = d + (p,pe), acumulador double-double */
static inline void ddAdd( double* hi, double* lo, double p, double pe )
{
  double s, e;
  twoSum( *hi, p, &s, &e );
  e += *lo + pe;
  twoSum( s, e, hi, lo );
}

/* etapa 2 en funcion APARTE y NO inline: el camino frio no debe inlinearse, y
 * ademas es necesario -- con la etapa double-double inline MSVC revienta con un
 * internal compiler error al optimizar las funciones que llaman al predicado
 * decenas de veces (triTriSegment: 24, triTriCoplanar: hasta 54). */
static int orient3dSlow( double ax, double ay, double az,
                         double bx, double by, double bz,
                         double cx, double cy, double cz )
{
  /* double-double. Cada producto de dos dobles es exacto (twoProd); los seis
   * terminos triples se acumulan con twoSum. */
  double hi = 0.0, lo = 0.0;
  const double X[3] = { ax, ay, az };
  const double U[6] = { by, bz, bz, bx, bx, by };
  const double W[6] = { cz, cy, cx, cz, cy, cx };
  const int    sg[6] = { 1, -1, 1, -1, 1, -1 };
  for( int k = 0; k < 6; ++k ) {
    double p, pe;
    twoProd( U[k], W[k], &p, &pe );          /* exacto */
    /* multiplicar por X[k/2] mantiene la estructura: (p+pe)*x */
    double q, qe, r, re;
    twoProd( p , X[k/2], &q, &qe );
    twoProd( pe, X[k/2], &r, &re );
    const double s = (double)sg[k];
    ddAdd( &hi, &lo, s*q, s*qe );
    ddAdd( &hi, &lo, s*r, s*re );
  }
  const double v = hi + lo;
  if( v > 0.0 ) return  1;
  if( v < 0.0 ) return -1;
  return 0;
}

static inline int orient3d( const double* a, const double* b, const double* c,
                            const double* d )
{
  const double ax=a[0]-d[0], ay=a[1]-d[1], az=a[2]-d[2];
  const double bx=b[0]-d[0], by=b[1]-d[1], bz=b[2]-d[2];
  const double cx=c[0]-d[0], cy=c[1]-d[1], cz=c[2]-d[2];

  const double bzcy=bz*cy, bycz=by*cz;
  const double bxcz=bx*cz, bzcx=bz*cx;
  const double bycx=by*cx, bxcy=bx*cy;
  const double det = ax*(bycz-bzcy) + ay*(bzcx-bxcz) + az*(bxcy-bycx);

  /* filtro dinamico: cota del error de redondeo del calculo de arriba */
  const double perm = fabs(ax)*( fabs(bycz)+fabs(bzcy) )
                    + fabs(ay)*( fabs(bzcx)+fabs(bxcz) )
                    + fabs(az)*( fabs(bxcy)+fabs(bycx) );
  const double bound = 7.7715611723761058e-16 * perm;      /* (7+56eps)*eps */
  if( det >  bound ) return  1;
  if( det < -bound ) return -1;
  if( perm == 0.0 )  return  0;              /* degenerado exacto: todo cero */
  return orient3dSlow( ax,ay,az, bx,by,bz, cx,cy,cz );
}

/* orient2d: SIGNO de det[ b-a ; c-a ] en 2D (>0 si a,b,c giran antihorario).
 * Misma estructura que orient3d -- filtro dinamico y escalada a double-double --
 * y hace falta por lo mismo: en el caso COPLANAR todos los orient3d son 0, asi
 * que las decisiones (dentro/fuera, cruce de aristas) pasan a ser 2D en el plano
 * comun, y tienen que seguir siendo exactas o la soldadura se rompe igual.
 * La escalada vive tambien en una funcion APARTE y NO inline, por lo mismo. */
static int orient2dSlow( double bax, double bay, double cax, double cay )
{
  double hi = 0.0, lo = 0.0, p, pe;
  twoProd(  bax, cay, &p, &pe );  ddAdd( &hi, &lo,  p,  pe );
  twoProd(  bay, cax, &p, &pe );  ddAdd( &hi, &lo, -p, -pe );
  const double v = hi + lo;
  if( v > 0.0 ) return  1;
  if( v < 0.0 ) return -1;
  return 0;
}
static inline int orient2d( double ax, double ay, double bx, double by,
                            double cx, double cy )
{
  const double bax = bx-ax, bay = by-ay;
  const double cax = cx-ax, cay = cy-ay;
  const double l = bax*cay, r = bay*cax;
  const double det = l - r;
  const double bound = 6.6613381477509392e-16 * ( fabs(l) + fabs(r) );  /* (3+8eps)eps */
  if( det >  bound ) return  1;
  if( det < -bound ) return -1;
  if( l == 0.0 && r == 0.0 ) return 0;
  return orient2dSlow( bax, bay, cax, cay );
}

/* Para trabajar en el plano de un triangulo coplanar se PROYECTA sobre el par de
 * ejes con mayor area, o sea se descarta la coordenada de mayor |normal|. El
 * orden ciclico (drop x -> (y,z), drop y -> (z,x), drop z -> (x,y)) preserva la
 * orientacion, asi que los signos de orient2d siguen significando lo mismo.
 * Solo se compara MAGNITUD de la normal: no decide topologia, solo elige la
 * proyeccion menos degenerada, asi que no necesita ser exacto. */
static inline void planeAxes( const double* p, const double* q, const double* r,
                              int* i0, int* i1 )
{
  const double ux=q[0]-p[0], uy=q[1]-p[1], uz=q[2]-p[2];
  const double vx=r[0]-p[0], vy=r[1]-p[1], vz=r[2]-p[2];
  const double nx = fabs( uy*vz - uz*vy );
  const double ny = fabs( uz*vx - ux*vz );
  const double nz = fabs( ux*vy - uy*vx );
  if( nx >= ny && nx >= nz )      { *i0 = 1; *i1 = 2; }   /* fuera x */
  else if( ny >= nz )             { *i0 = 2; *i1 = 0; }   /* fuera y */
  else                            { *i0 = 0; *i1 = 1; }   /* fuera z */
}

/* area con signo x2 del triangulo (a,b,c) proyectado sobre los ejes (i0,i1):
 * LINEAL en cada punto, que es lo que permite sacar los parametros de un cruce
 * como cociente de dos de ellas. */
static inline double area2d( const double* a, const double* b, const double* c,
                             int i0, int i1 )
{
  return ( b[i0]-a[i0] )*( c[i1]-a[i1] ) - ( b[i1]-a[i1] )*( c[i0]-a[i0] );
}

/* determinante CON SIGNO det[p-x ; q-x ; r-x] -- la magnitud de la que orient3d
 * toma el signo. Proporcional a la distancia con signo de x al plano (p,q,r), y
 * LINEAL en x, que es lo que permite sacar el parametro del corte como un
 * cociente de dos de ellos. */
static inline double det3sub( const double* p, const double* q, const double* r,
                              const double* x )
{
  const double ax=p[0]-x[0], ay=p[1]-x[1], az=p[2]-x[2];
  const double bx=q[0]-x[0], by=q[1]-x[1], bz=q[2]-x[2];
  const double cx=r[0]-x[0], cy=r[1]-x[1], cz=r[2]-x[2];
  return ax*(by*cz-bz*cy) + ay*(bz*cx-bx*cz) + az*(bx*cy-by*cx);
}

/* ------------------------------------------- interseccion triangulo-triangulo
 * Devuelve los extremos del segmento de corte por CLAVE CANONICA (arista padre
 * + parametro), NUNCA por coordenadas -- para que dos triangulos que comparten
 * una arista produzcan claves BIT-IDENTICAS y la soldadura posterior sea un
 * `unique` por filas sin tolerancia (la tecnica de MESH\MeshZeroContour.m).
 *
 * FORMULACION. En vez de la cascada de permutaciones de Devillers-Guigue, se usa
 * la caracterizacion combinatoria: los extremos del segmento de corte son
 * SIEMPRE puntos donde una ARISTA de un triangulo atraviesa el OTRO triangulo.
 * Enumerar las 6 aristas es un pelo mas caro (24 orient3d en el caso completo
 * frente a ~12) pero cada resultado ES una clave canonica por construccion, y
 * todo el predicado se decide con orient3d: sin coordenadas y sin divisiones.
 * El rechazo, que es el caso comun, cuesta 3 orient3d.
 *
 * EL PARAMETRO. t = Da/(Da-Db) con Da = det(p,q,r,a): depende SOLO del triangulo
 * contrario y del extremo de la arista, jamas de por que triangulo se llego. Se
 * mide desde el vertice de id GLOBAL menor, asi que sale bit-identico desde
 * cualquiera de los dos triangulos que comparten la arista. Un signo 0 significa
 * que el vertice esta EN el plano: la clave colapsa a la forma de vertice, que es
 * el mismo snapping que hace MeshZeroContour con t==0 y t==1.
 *
 * ALCANCE: el caso COPLANAR (los seis signos 0) se DETECTA exactamente y se
 * informa con su propio codigo; recortar el poligono comun es una rama aparte,
 * todavia sin implementar.
 *
 * HUECO CONOCIDO -- IDENTIDAD ENTRE MALLAS. Las claves son canonicas DENTRO de
 * cada malla, pero los ids de A y de B viven en espacios de nombres distintos.
 * Un punto que es vertice de A Y vertice de B recibe por tanto DOS claves
 * distintas: (onB=0, ia=idA) y (onB=1, ia=idB). El dedupe combinatorio de aqui no
 * puede verlo, asi que un contacto puntual vertice-sobre-vertice se informa como
 * IX_SEGMENT con dos claves en vez de IX_TOUCH con una. Verificado con
 * T1=[-1 0 0;1 0 0;0 2 0] / T2=[1 0 0;2 1 1;2 -1 1].
 * RESOLUCION (en la fase de soldadura, no aqui): las claves de VERTICE se sueldan
 * entre mallas por igualdad EXACTA de coordenadas -- y eso es licito porque la
 * coordenada de una clave de vertice es una coordenada de ENTRADA, sin ninguna
 * aritmetica encima, asi que comparar bits es exacto y no es una tolerancia. Las
 * claves de ARISTA se sueldan por clave canonica, como ya hacen. Las dos vias son
 * exactas; ninguna introduce un epsilon. */

#define IX_NONE      0    /* no se cortan */
#define IX_SEGMENT   1    /* se cortan en un segmento (2 claves) */
#define IX_TOUCH     2    /* se tocan en un solo punto (1 clave) */
#define IX_COPLANAR  3    /* coplanares: pendiente de la rama de recorte 2D */

/* una clave canonica de punto de interseccion. TRES tipos: el caso no coplanar
 * solo produce los dos primeros; la rama COPLANAR introduce el tercero (cruce de
 * una arista de T1 con una de T2), que es la razon de fijar este contrato antes
 * de construir el recorrido. */
#define IXK_VERTEX  1   /* un vertice: ia. ib = -1                              */
#define IXK_EDGE    2   /* punto sobre la arista canonica (ia<ib), parametro t  */
#define IXK_EDGEX   3   /* cruce de la arista (ia<ib) de T1 con la (ja<jb) de T2 */
struct IxKey {
  int    kind;  /* IXK_*                                                        */
  int    onB;   /* kinds 1 y 2: 0 = la arista/vertice es de T1, 1 = de T2.
                 * kind 3: no aplica (la clave nombra las DOS aristas)          */
  int    ia;    /* id GLOBAL menor (o el unico, si es un vertice)               */
  int    ib;    /* id GLOBAL mayor; -1 si es un VERTICE                         */
  int    ja;    /* kind 3: arista de T2, id menor; si no, -1                    */
  int    jb;    /* kind 3: arista de T2, id mayor; si no, -1                    */
  double t;     /* kinds 2 y 3: parametro desde ia hacia ib. kind 1: 0          */
  double s;     /* kind 3: parametro desde ja hacia jb. si no, 0                */
};
static inline IxKey ixVertex( int onB , int g )
{ IxKey k; k.kind=IXK_VERTEX; k.onB=onB; k.ia=g; k.ib=-1; k.ja=-1; k.jb=-1; k.t=0.0; k.s=0.0; return k; }

/* igualdad de claves: comparacion EXACTA de enteros y de los parametros. Es la
 * operacion sobre la que descansa toda la soldadura, asi que no lleva tolerancia
 * a proposito -- dos triangulos vecinos producen el mismo parametro BIT a BIT
 * porque se calcula desde la arista canonica. */
static inline bool sameKey( const IxKey& a , const IxKey& b )
{
  if( a.kind != b.kind ) return false;
  if( a.kind == IXK_VERTEX ) return a.onB == b.onB && a.ia == b.ia;
  if( a.kind == IXK_EDGE   ) return a.onB == b.onB && a.ia == b.ia &&
                                    a.ib  == b.ib  && a.t  == b.t;
  return a.ia == b.ia && a.ib == b.ib && a.ja == b.ja && a.jb == b.jb &&
         a.t  == b.t  && a.s  == b.s;
}

/* ¿el segmento (a,b) atraviesa el triangulo (p,q,r)? sa/sb son los signos ya
 * calculados de a y b respecto al plano de (p,q,r). Exacto: solo orient3d.
 * Devuelve en `zed` CUALES de los tres tetraedros salieron 0, que es lo que dice
 * si el corte cae sobre una ARISTA o un VERTICE de (p,q,r) -- imprescindible para
 * darle la clave canonica correcta (ver la nota en edgeVsTriKey). */
static inline bool segCrossesTri( const double* a, const double* b,
                                  const double* p, const double* q, const double* r,
                                  int sa, int sb, int zed[3] )
{
  zed[0] = zed[1] = zed[2] = 0;
  /* sa != sb cubre los dos casos vivos: lados opuestos, o un extremo EN el plano.
   * sa == sb incluye "ambos 0" = arista contenida en el plano, que es coplanar y
   * lo trata la rama de coplanares, no esta. */
  if( sa == sb ) return false;
  /* el punto de corte cae dentro de (p,q,r) si los tres tetraedros que forma la
   * recta ab con las tres aristas del triangulo giran igual */
  const int o1 = orient3d( a, b, p, q );
  const int o2 = orient3d( a, b, q, r );
  const int o3 = orient3d( a, b, r, p );
  if( o1 == 0 && o2 == 0 && o3 == 0 ) return false;   /* ab en el plano: coplanar */
  const int s = ( o1 != 0 ) ? o1 : ( ( o2 != 0 ) ? o2 : o3 );
  if( o1 != 0 && o1 != s ) return false;
  if( o2 != 0 && o2 != s ) return false;
  if( o3 != 0 && o3 != s ) return false;
  zed[0] = ( o1 == 0 );  zed[1] = ( o2 == 0 );  zed[2] = ( o3 == 0 );
  return true;
}

/* parametro del cruce de la recta (la,lb) con la recta (ma,mb), sabiendo que son
 * COPLANARES (es lo que garantiza orient3d(la,lb,ma,mb)==0). Se resuelve en el
 * plano que ambas generan, con las mismas areas con signo que la rama coplanar,
 * de modo que el valor depende SOLO de los cuatro extremos -- no del triangulo
 * por el que se llego. Eso es lo que lo hace canonico. */
static inline double crossParam( const double* la, const double* lb,
                                 const double* ma, const double* mb )
{
  int i0, i1;  planeAxes( la, lb, ma, &i0, &i1 );
  const double A1 = area2d( ma, mb, la, i0, i1 );
  const double A2 = area2d( ma, mb, lb, i0, i1 );
  if( A1 == A2 ) {                      /* proyeccion degenerada: prueba otro par */
    planeAxes( la, lb, mb, &i0, &i1 );
    const double B1 = area2d( ma, mb, la, i0, i1 );
    const double B2 = area2d( ma, mb, lb, i0, i1 );
    if( B1 == B2 ) return 0.0;
    double t = B1/(B1-B2);
    return ( t < 0.0 ) ? 0.0 : ( ( t > 1.0 ) ? 1.0 : t );
  }
  double t = A1/(A1-A2);
  return ( t < 0.0 ) ? 0.0 : ( ( t > 1.0 ) ? 1.0 : t );
}

/* clave canonica del corte de la arista (ga,gb) con el plano de (p,q,r).
 * Da/Db son los determinantes con signo; sa/sb sus signos exactos. */
/* Clave canonica del corte de la arista (ga,gb) con el triangulo (p,q,r), cuyos
 * vertices tienen ids gp,gq,gr. `zed` dice cuales de los tres tetraedros salieron
 * cero, o sea si el corte cae sobre una arista o un vertice del triangulo.
 *
 * POR QUE HACE FALTA. Si el corte cae en el INTERIOR del triangulo, la clave
 * "punto en la arista (ga,gb) con parametro t" es canonica: t se calcula contra
 * ese triangulo y los dos triangulos que comparten la arista dan el mismo valor.
 * Pero si el corte cae sobre una ARISTA del triangulo, esa arista la comparten
 * DOS triangulos, y t medido contra cada uno de sus planos difiere en los ultimos
 * bits -- dos representaciones del mismo lugar, dos redondeos. La clave deja de
 * ser canonica y la soldadura abre la curva. Medido: la misma arista de A recibia
 * t = 0.61263911317434949 y t = 0.6126391131743496 segun por que cara de B se
 * llegaba, y el lazo quedaba con 4 extremos sueltos.
 * La clave correcta en ese caso es el CRUCE ARISTA x ARISTA (IXK_EDGEX), cuyos
 * parametros dependen solo de los cuatro extremos; y si caen dos ceros, el corte
 * esta en un VERTICE del triangulo y la clave es ese vertice. */
static inline IxKey ixKeyEdgeTri( int onB, int ga, int gb,
                                  const double* A, const double* B,
                                  const double* p, const double* q, const double* r,
                                  int gp, int gq, int gr,
                                  int sa, int sb, const int zed[3] )
{
  IxKey k;  k.kind = IXK_EDGE;  k.onB = onB;  k.t = 0.0;  k.s = 0.0;
  k.ib = -1;  k.ja = -1;  k.jb = -1;
  if( sa == 0 ) { k.kind = IXK_VERTEX; k.ia = ga; return k; }  /* el corte ES a */
  if( sb == 0 ) { k.kind = IXK_VERTEX; k.ia = gb; return k; }  /* ... o b       */

  const int nz = zed[0] + zed[1] + zed[2];
  if( nz >= 2 ) {                     /* dos ceros: el corte es un VERTICE de pqr */
    int gv;
    if(      zed[0] && zed[1] ) gv = gq;      /* aristas pq y qr -> q */
    else if( zed[1] && zed[2] ) gv = gr;      /* qr y rp -> r */
    else                        gv = gp;      /* rp y pq -> p */
    k.kind = IXK_VERTEX;  k.onB = 1 - onB;    /* el vertice es del OTRO triangulo */
    k.ia = gv;  return k;
  }
  if( nz == 1 ) {                     /* un cero: el corte esta en una ARISTA de pqr */
    const double *ma, *mb;  int fa, fb;
    if(      zed[0] ) { ma=p; mb=q; fa=gp; fb=gq; }
    else if( zed[1] ) { ma=q; mb=r; fa=gq; fb=gr; }
    else              { ma=r; mb=p; fa=gr; fb=gp; }
    const double *la = A, *lb = B;  int ea = ga, eb = gb;
    if( eb < ea ) { const int z=ea; ea=eb; eb=z;  la=B; lb=A; }   /* canonicas */
    if( fb < fa ) { const int z=fa; fa=fb; fb=z;  const double* w=ma; ma=mb; mb=w; }
    /* la clave EDGEX nombra SIEMPRE (arista de A , arista de B) en ese orden, para
     * que salga identica se llegue recorriendo las aristas de A o las de B */
    k.kind = IXK_EDGEX;  k.onB = 0;
    if( onB == 0 ) { k.ia=ea; k.ib=eb; k.ja=fa; k.jb=fb;
                     k.t = crossParam( la, lb, ma, mb );
                     k.s = crossParam( ma, mb, la, lb ); }
    else           { k.ia=fa; k.ib=fb; k.ja=ea; k.jb=eb;
                     k.t = crossParam( ma, mb, la, lb );
                     k.s = crossParam( la, lb, ma, mb ); }
    return k;
  }
  /* arista canonica: siempre del id menor al mayor, para que t sea el mismo
   * venga de donde venga */
  const double* lo = A;  const double* hi = B;
  int gl = ga, gh = gb;
  if( gb < ga ) { lo = B; hi = A; gl = gb; gh = ga; }
  const double Dl = det3sub( p, q, r, lo );
  const double Dh = det3sub( p, q, r, hi );
  const double den = Dl - Dh;
  k.ia = gl;  k.ib = gh;
  k.t  = ( den != 0.0 ) ? ( Dl / den ) : 0.0;
  if( k.t < 0.0 ) k.t = 0.0; else if( k.t > 1.0 ) k.t = 1.0;
  return k;
}

/* Kernel: T1 y T2 son 9 dobles (3 vertices xyz); g1 y g2 los ids GLOBALES de sus
 * vertices. Devuelve IX_* y hasta 2 claves canonicas en out.
 * Coste: 3 orient3d para rechazar (el caso comun), 24 en el caso completo. */
static int triTriSegment( const double* T1, const int* g1,
                          const double* T2, const int* g2,
                          IxKey out[2] , int* nOut )
{
  const double *p1=T1, *q1=T1+3, *r1=T1+6;
  const double *p2=T2, *q2=T2+3, *r2=T2+6;
  *nOut = 0;

  const int s1[3] = { orient3d(p2,q2,r2,p1), orient3d(p2,q2,r2,q1), orient3d(p2,q2,r2,r1) };
  if( ( s1[0]>0 && s1[1]>0 && s1[2]>0 ) || ( s1[0]<0 && s1[1]<0 && s1[2]<0 ) )
    return IX_NONE;                                    /* T1 entero a un lado */
  const int s2[3] = { orient3d(p1,q1,r1,p2), orient3d(p1,q1,r1,q2), orient3d(p1,q1,r1,r2) };
  if( ( s2[0]>0 && s2[1]>0 && s2[2]>0 ) || ( s2[0]<0 && s2[1]<0 && s2[2]<0 ) )
    return IX_NONE;
  if( s1[0]==0 && s1[1]==0 && s1[2]==0 )
    return IX_COPLANAR;                                /* deteccion EXACTA */

  const double* V1[3] = { p1,q1,r1 };
  const double* V2[3] = { p2,q2,r2 };
  static const int E[3][2] = { {0,1},{1,2},{2,0} };
  IxKey acc[6];  int n = 0;
  int zed[3];
  for( int e = 0; e < 3; ++e ) {
    const int a=E[e][0], b=E[e][1];
    if( segCrossesTri( V1[a],V1[b], p2,q2,r2, s1[a],s1[b], zed ) )
      acc[n++] = ixKeyEdgeTri( 0, g1[a],g1[b], V1[a],V1[b], p2,q2,r2,
                               g2[0],g2[1],g2[2], s1[a],s1[b], zed );
  }
  for( int e = 0; e < 3; ++e ) {
    const int a=E[e][0], b=E[e][1];
    if( segCrossesTri( V2[a],V2[b], p1,q1,r1, s2[a],s2[b], zed ) )
      acc[n++] = ixKeyEdgeTri( 1, g2[a],g2[b], V2[a],V2[b], p1,q1,r1,
                               g1[0],g1[1],g1[2], s2[a],s2[b], zed );
  }
  /* dedupe EXACTO por clave: un vertice que esta en el plano contrario lo emiten
   * las DOS aristas incidentes, con claves bit-identicas. Es el mismo argumento
   * que hace funcionar la soldadura global. */
  int m = 0;
  for( int i = 0; i < n; ++i ) {
    bool dup = false;
    for( int j = 0; j < m; ++j )
      if( sameKey( acc[i] , out[j] ) ) { dup = true; break; }
    if( !dup ) { if( m < 2 ) out[m] = acc[i]; ++m; }
  }
  *nOut = ( m > 2 ) ? 2 : m;
  if( m == 0 ) return IX_NONE;
  if( m == 1 ) return IX_TOUCH;
  return IX_SEGMENT;
}

/* ------ rama COPLANAR: solape de dos triangulos en el mismo plano ------------
 * La interseccion es un poligono convexo de hasta 6 vertices, y esos vertices
 * son de exactamente tres clases: vertice de T1 dentro de T2, vertice de T2
 * dentro de T1, y cruce de una arista de T1 con una de T2 (la clave IXK_EDGEX,
 * que el caso no coplanar nunca produce). Se enumeran EXPLICITAMENTE en vez de
 * recortar con Sutherland-Hodgman, porque asi cada punto sale ya con su clave
 * canonica en lugar de haber que rastrear la procedencia a traves del recorte.
 * Cuesta hasta 54 orient2d, pero los pares coplanares son raros y aqui manda la
 * correccion.
 *
 * El ORDEN del poligono se resuelve con orient2d (envoltura convexa por gift
 * wrapping), NO ordenando por angulo: un atan2 sobre coordenadas redondeadas
 * seria una decision tomada sobre una construccion, que es justo lo que este
 * diseno evita. Las coordenadas que devuelve son solo para salida.
 *
 * Devuelve cuantos vertices tiene el solape, EN ORDEN: 0 = sin solape, 1 = toque
 * en un punto, 2 = contacto en un segmento, >=3 = solape con area (el llamante
 * emite los segmentos consecutivos, cerrando el ultimo con el primero). */

/* ¿(x,y) dentro del triangulo 2D? 1 dentro, 0 en el borde, -1 fuera. Exacto.
 * `ccw` normaliza el sentido de giro del triangulo. */
static inline int inTri2d( const double* A, const double* B, const double* C,
                           double x, double y, int i0, int i1, int ccw )
{
  const int s1 = ccw * orient2d( A[i0],A[i1], B[i0],B[i1], x,y );
  const int s2 = ccw * orient2d( B[i0],B[i1], C[i0],C[i1], x,y );
  const int s3 = ccw * orient2d( C[i0],C[i1], A[i0],A[i1], x,y );
  if( s1 < 0 || s2 < 0 || s3 < 0 ) return -1;
  if( s1 == 0 || s2 == 0 || s3 == 0 ) return  0;
  return 1;
}

static int triTriCoplanar( const double* T1, const int* g1,
                           const double* T2, const int* g2,
                           IxKey out[6], double xyz[18] )
{
  const double* V1[3] = { T1, T1+3, T1+6 };
  const double* V2[3] = { T2, T2+3, T2+6 };
  int i0, i1;  planeAxes( T1, T1+3, T1+6, &i0, &i1 );
  const int c1 = orient2d( V1[0][i0],V1[0][i1], V1[1][i0],V1[1][i1], V1[2][i0],V1[2][i1] );
  const int c2 = orient2d( V2[0][i0],V2[0][i1], V2[1][i0],V2[1][i1], V2[2][i0],V2[2][i1] );
  if( c1 == 0 || c2 == 0 ) return 0;            /* degenerado en el plano */

  IxKey  K[15];  double P[15][3];  int n = 0;
  static const int E[3][2] = { {0,1},{1,2},{2,0} };

  for( int v = 0; v < 3; ++v )                  /* vertices de T1 dentro de T2 */
    if( inTri2d( V2[0],V2[1],V2[2], V1[v][i0],V1[v][i1], i0,i1, c2 ) >= 0 ) {
      K[n] = ixVertex( 0, g1[v] );
      P[n][0]=V1[v][0]; P[n][1]=V1[v][1]; P[n][2]=V1[v][2];  ++n;
    }
  for( int v = 0; v < 3; ++v )                  /* vertices de T2 dentro de T1 */
    if( inTri2d( V1[0],V1[1],V1[2], V2[v][i0],V2[v][i1], i0,i1, c1 ) >= 0 ) {
      K[n] = ixVertex( 1, g2[v] );
      P[n][0]=V2[v][0]; P[n][1]=V2[v][1]; P[n][2]=V2[v][2];  ++n;
    }
  for( int e = 0; e < 3 && n < 15; ++e ) {      /* cruces arista x arista */
    const double *a = V1[E[e][0]], *b = V1[E[e][1]];
    for( int f = 0; f < 3 && n < 15; ++f ) {
      const double *c = V2[E[f][0]], *d = V2[E[f][1]];
      const int o1 = orient2d( a[i0],a[i1], b[i0],b[i1], c[i0],c[i1] );
      const int o2 = orient2d( a[i0],a[i1], b[i0],b[i1], d[i0],d[i1] );
      const int o3 = orient2d( c[i0],c[i1], d[i0],d[i1], a[i0],a[i1] );
      const int o4 = orient2d( c[i0],c[i1], d[i0],d[i1], b[i0],b[i1] );
      if( o1 == 0 || o2 == 0 || o3 == 0 || o4 == 0 ) continue;  /* toque: ya lo
                                * aporta la clase de vertice; no duplicar aqui */
      if( o1 == o2 || o3 == o4 ) continue;                      /* no se cruzan */
      /* CANONICALIZAR PRIMERO y medir el parametro desde los extremos canonicos.
       * Calcularlo en el orden local y luego hacer t = 1-t al ordenar los ids
       * deja el parametro referido a una direccion y el punto a la otra: el punto
       * sale mal (lo estuvo, en los 2 de 6 cruces cuyos ids venian decrecientes,
       * y casi lo achaco a la referencia en vez a mi codigo). */
      int ea = g1[E[e][0]], eb = g1[E[e][1]];
      const double *la = a, *lb = b;
      if( eb < ea ) { const int z=ea; ea=eb; eb=z;  la=b; lb=a; }
      int fa = g2[E[f][0]], fb = g2[E[f][1]];
      const double *ma = c, *mb = d;
      if( fb < fa ) { const int z=fa; fa=fb; fb=z;  ma=d; mb=c; }
      const double A1 = area2d( ma,mb,la, i0,i1 ), A2 = area2d( ma,mb,lb, i0,i1 );
      const double B1 = area2d( la,lb,ma, i0,i1 ), B2 = area2d( la,lb,mb, i0,i1 );
      const double tt = ( A1 != A2 ) ? A1/(A1-A2) : 0.0;   /* de ea hacia eb */
      const double ss = ( B1 != B2 ) ? B1/(B1-B2) : 0.0;   /* de fa hacia fb */
      IxKey k;  k.kind=IXK_EDGEX;  k.onB=0;
      k.ia=ea; k.ib=eb; k.ja=fa; k.jb=fb;  k.t=tt;  k.s=ss;
      K[n] = k;
      for( int q = 0; q < 3; ++q ) P[n][q] = la[q] + tt*( lb[q]-la[q] );
      ++n;
    }
  }
  if( n == 0 ) return 0;

  int m = 0;                                    /* dedupe EXACTO por clave */
  for( int i = 0; i < n; ++i ) {
    bool dup = false;
    for( int j = 0; j < m; ++j ) if( sameKey( K[i] , out[j] ) ) { dup = true; break; }
    if( dup ) continue;
    if( m < 6 ) {
      out[m] = K[i];
      xyz[3*m]=P[i][0]; xyz[3*m+1]=P[i][1]; xyz[3*m+2]=P[i][2];
      ++m;
    }
  }
  /* COLAPSAR puntos con coordenadas BIT-IDENTICAS. Hace falta de verdad: un
   * vertice compartido por las dos mallas entra dos veces (una clave por malla,
   * porque los ids viven en espacios de nombres distintos -- el hueco documentado
   * arriba), y dos puntos coincidentes hacen degenerar la envoltura convexa. La
   * comparacion exacta es licita porque estas son coordenadas de ENTRADA (claves
   * de vertice) o construidas de forma identica, no un epsilon geometrico. */
  {
    int k2 = 0;
    for( int i = 0; i < m; ++i ) {
      bool same = false;
      for( int j = 0; j < k2; ++j )
        if( xyz[3*i] == xyz[3*j] && xyz[3*i+1] == xyz[3*j+1] && xyz[3*i+2] == xyz[3*j+2] )
          { same = true; break; }
      if( same ) continue;
      if( k2 != i ) {
        out[k2] = out[i];
        xyz[3*k2]=xyz[3*i]; xyz[3*k2+1]=xyz[3*i+1]; xyz[3*k2+2]=xyz[3*i+2];
      }
      ++k2;
    }
    m = k2;
  }
  if( m < 3 ) return m;

  /* ORDEN convexo por gift wrapping, decidido solo con orient2d */
  int  ord[6];  bool used[6] = {false,false,false,false,false,false};
  int  start = 0;                               /* el de menor (x,y) proyectado */
  for( int i = 1; i < m; ++i )
    if( xyz[3*i+i0] <  xyz[3*start+i0] ||
      ( xyz[3*i+i0] == xyz[3*start+i0] && xyz[3*i+i1] < xyz[3*start+i1] ) ) start = i;
  int cur = start, h = 0;
  while( h < m ) {
    ord[h++] = cur;  used[cur] = true;
    int nxt = -1;
    for( int cand = 0; cand < m; ++cand ) {
      if( cand == cur ) continue;
      if( used[cand] && cand != start ) continue;
      if( nxt < 0 ) { nxt = cand; continue; }
      const int o = orient2d( xyz[3*cur+i0],xyz[3*cur+i1],
                              xyz[3*nxt+i0],xyz[3*nxt+i1],
                              xyz[3*cand+i0],xyz[3*cand+i1] );
      if( o < 0 ) nxt = cand;                   /* mas a la derecha: mejor apoyo */
    }
    if( nxt < 0 || nxt == start ) break;        /* cerrado (o sin candidatos) */
    cur = nxt;
  }
  /* h es cuantos entraron DE VERDAD en la envoltura: leer ord[] mas alla seria
   * usar basura como indice (lo era, y reventaba con dos puntos coincidentes) */
  IxKey  ko[6];  double po[18];
  for( int i = 0; i < h; ++i ) {
    ko[i] = out[ ord[i] ];
    po[3*i]=xyz[3*ord[i]]; po[3*i+1]=xyz[3*ord[i]+1]; po[3*i+2]=xyz[3*ord[i]+2];
  }
  for( int i = 0; i < h; ++i ) {
    out[i] = ko[i];
    xyz[3*i]=po[3*i]; xyz[3*i+1]=po[3*i+1]; xyz[3*i+2]=po[3*i+2];
  }
  return h;
}

/* ---------------------------------------------------------------- primitives */

static inline double d2Point( const double* a, const double* p, double* cp )
{
  cp[0]=a[0]; cp[1]=a[1]; cp[2]=a[2];
  return sq(p[0]-a[0]) + sq(p[1]-a[1]) + sq(p[2]-a[2]);
}

/* closest point on a segment. Los EXTREMOS se copian, no se interpolan: con
 * cp = a + t*(b-a) y t==1 el resultado es fl(a + fl(b-a)), que NO es b (falla
 * hasta 1 ulp de la magnitud de la coordenada). Y la baricentrica de ese caso es
 * [0,1], que reconstruye b EXACTO -- asi que cp y bc*V discrepaban un ulp
 * justamente en el caso mas facil. Copiar el extremo lo hace exacto por los dos
 * lados y cuesta lo mismo: la rama del clamp ya estaba ahi. (a se copiaba bien
 * por accidente: a + 0*v == a.) Importa en mallas de SEGMENTOS/polilineas, donde
 * el extremo es un caso corriente y no lo captura ninguna region de vertice. */
static inline double d2Seg( const double* a, const double* b, const double* p, double* cp )
{
  const double v0=b[0]-a[0], v1=b[1]-a[1], v2=b[2]-a[2];
  const double vv=v0*v0+v1*v1+v2*v2;
  double t = (p[0]-a[0])*v0 + (p[1]-a[1])*v1 + (p[2]-a[2])*v2;
  t = ( vv > 0.0 ) ? t/vv : 0.0;
  if( t <= 0.0 )      { cp[0]=a[0]; cp[1]=a[1]; cp[2]=a[2]; }
  else if( t >= 1.0 ) { cp[0]=b[0]; cp[1]=b[1]; cp[2]=b[2]; }
  else                { cp[0]=a[0]+t*v0; cp[1]=a[1]+t*v1; cp[2]=a[2]+t*v2; }
  return sq(p[0]-cp[0]) + sq(p[1]-cp[1]) + sq(p[2]-cp[2]);
}

/* closest point on a triangle: Ericson, "Real-Time Collision Detection" 5.1.5 */
static double d2Tri( const double* A, const double* B, const double* C,
                     const double* p, double* cp )
{
  const double ab0=B[0]-A[0], ab1=B[1]-A[1], ab2=B[2]-A[2];
  const double ac0=C[0]-A[0], ac1=C[1]-A[1], ac2=C[2]-A[2];
  const double ap0=p[0]-A[0], ap1=p[1]-A[1], ap2=p[2]-A[2];

  const double d1 = ab0*ap0 + ab1*ap1 + ab2*ap2;
  const double d2 = ac0*ap0 + ac1*ap1 + ac2*ap2;
  if( d1 <= 0.0 && d2 <= 0.0 ) return d2Point( A, p, cp );          /* vertex A */

  const double bp0=p[0]-B[0], bp1=p[1]-B[1], bp2=p[2]-B[2];
  const double d3 = ab0*bp0 + ab1*bp1 + ab2*bp2;
  const double d4 = ac0*bp0 + ac1*bp1 + ac2*bp2;
  if( d3 >= 0.0 && d4 <= d3 ) return d2Point( B, p, cp );           /* vertex B */

  const double vc = d1*d4 - d3*d2;
  if( vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 ) {                       /* edge AB */
    const double den = d1 - d3;
    const double t = ( den != 0.0 ) ? d1/den : 0.0;
    cp[0]=A[0]+t*ab0; cp[1]=A[1]+t*ab1; cp[2]=A[2]+t*ab2;
    return sq(p[0]-cp[0]) + sq(p[1]-cp[1]) + sq(p[2]-cp[2]);
  }

  const double cq0=p[0]-C[0], cq1=p[1]-C[1], cq2=p[2]-C[2];
  const double d5 = ab0*cq0 + ab1*cq1 + ab2*cq2;
  const double d6 = ac0*cq0 + ac1*cq1 + ac2*cq2;
  if( d6 >= 0.0 && d5 <= d6 ) return d2Point( C, p, cp );           /* vertex C */

  const double vb = d5*d2 - d1*d6;
  if( vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 ) {                       /* edge AC */
    const double den = d2 - d6;
    const double t = ( den != 0.0 ) ? d2/den : 0.0;
    cp[0]=A[0]+t*ac0; cp[1]=A[1]+t*ac1; cp[2]=A[2]+t*ac2;
    return sq(p[0]-cp[0]) + sq(p[1]-cp[1]) + sq(p[2]-cp[2]);
  }

  const double va = d3*d6 - d5*d4;
  if( va <= 0.0 && (d4-d3) >= 0.0 && (d5-d6) >= 0.0 ) {             /* edge BC */
    const double den = (d4-d3) + (d5-d6);
    const double t = ( den != 0.0 ) ? (d4-d3)/den : 0.0;
    cp[0]=B[0]+t*(C[0]-B[0]); cp[1]=B[1]+t*(C[1]-B[1]); cp[2]=B[2]+t*(C[2]-B[2]);
    return sq(p[0]-cp[0]) + sq(p[1]-cp[1]) + sq(p[2]-cp[2]);
  }

  const double den = va + vb + vc;                                  /* interior */
  if( den != 0.0 && std::isfinite( den ) ) {
    const double v = vb/den, w = vc/den;
    cp[0]=A[0]+v*ab0+w*ac0; cp[1]=A[1]+v*ab1+w*ac1; cp[2]=A[2]+v*ab2+w*ac2;
    return sq(p[0]-cp[0]) + sq(p[1]-cp[1]) + sq(p[2]-cp[2]);
  }

  /* fully degenerate triangle: closest of the three edges */
  double c2[3], c3[3];
  double q1 = d2Seg( A, B, p, cp );
  double q2 = d2Seg( A, C, p, c2 );
  double q3 = d2Seg( B, C, p, c3 );
  if( q2 < q1 ) { q1=q2; cp[0]=c2[0]; cp[1]=c2[1]; cp[2]=c2[2]; }
  if( q3 < q1 ) { q1=q3; cp[0]=c3[0]; cp[1]=c3[1]; cp[2]=c3[2]; }
  return q1;
}

static inline double det3( const double* u, const double* v, const double* w )
{
  return u[0]*(v[1]*w[2]-v[2]*w[1])
       - u[1]*(v[0]*w[2]-v[2]*w[0])
       + u[2]*(v[0]*w[1]-v[1]*w[0]);
}

static double d2Tet( const double* A, const double* B, const double* C,
                     const double* D, const double* p, double* cp )
{
  const double BA[3]={B[0]-A[0],B[1]-A[1],B[2]-A[2]};
  const double CA[3]={C[0]-A[0],C[1]-A[1],C[2]-A[2]};
  const double DA[3]={D[0]-A[0],D[1]-A[1],D[2]-A[2]};
  const double d0 = det3( BA, CA, DA );
  if( d0 != 0.0 && std::isfinite( d0 ) ) {
    const double Bp[3]={B[0]-p[0],B[1]-p[1],B[2]-p[2]};
    const double Cp[3]={C[0]-p[0],C[1]-p[1],C[2]-p[2]};
    const double Dp[3]={D[0]-p[0],D[1]-p[1],D[2]-p[2]};
    const double pA[3]={p[0]-A[0],p[1]-A[1],p[2]-A[2]};
    const double l1 = det3( Bp, Cp, Dp ) / d0;
    const double l2 = det3( pA, CA, DA ) / d0;
    const double l3 = det3( BA, pA, DA ) / d0;
    const double l4 = 1.0 - l1 - l2 - l3;
    const double tol = -1e-12;
    if( l1 >= tol && l2 >= tol && l3 >= tol && l4 >= tol ) {        /* inside */
      cp[0]=p[0]; cp[1]=p[1]; cp[2]=p[2];
      return 0.0;
    }
  }
  double c2[3];
  double q  = d2Tri( A, B, C, p, cp );
  double q2 = d2Tri( A, B, D, p, c2 );
  if( q2 < q ) { q=q2; cp[0]=c2[0]; cp[1]=c2[1]; cp[2]=c2[2]; }
  q2 = d2Tri( A, C, D, p, c2 );
  if( q2 < q ) { q=q2; cp[0]=c2[0]; cp[1]=c2[1]; cp[2]=c2[2]; }
  q2 = d2Tri( B, C, D, p, c2 );
  if( q2 < q ) { q=q2; cp[0]=c2[0]; cp[1]=c2[1]; cp[2]=c2[2]; }
  return q;
}

/* exact test against packed element j (verts at v = pkV + 12*j).
 * 4 nonzero nodes ALWAYS mean a tetrahedron. */
static inline double d2Elem( const double* v, int k, const double* p, double* cp )
{
  switch( k ) {
    case 1: return d2Point( v, p, cp );
    case 2: return d2Seg( v, v+3, p, cp );
    case 3: return d2Tri( v, v+3, v+6, p, cp );
    case 4: return d2Tet( v, v+3, v+6, v+9, p, cp );
  }
  return INF;   /* k == 0: empty 0-padded row */
}

/* --------------------------------------------- region-exact barycentrics
 * The closest-point search classifies each hit into a REGION (vertex / edge /
 * face / interior). Emitting the barycentric coordinates FROM that region --
 * instead of reverse-engineering them from the (rounded) closest point -- makes
 * edge and vertex hits give EXACT zeros for the off-components and keeps the
 * whole thing well-conditioned on slivers (needle triangles), where the
 * from-cp cross-product forms divide by a vanishing area and lose 5-15 digits.
 * bcTri3 mirrors d2Tri's cascade bit-for-bit so the region matches the cp the
 * traversal already found. */
static void bcTri3( const double* A, const double* B, const double* C,
                    const double* p, double bc[3] )
{
  const double ab0=B[0]-A[0], ab1=B[1]-A[1], ab2=B[2]-A[2];
  const double ac0=C[0]-A[0], ac1=C[1]-A[1], ac2=C[2]-A[2];
  const double ap0=p[0]-A[0], ap1=p[1]-A[1], ap2=p[2]-A[2];
  const double d1 = ab0*ap0 + ab1*ap1 + ab2*ap2;
  const double d2 = ac0*ap0 + ac1*ap1 + ac2*ap2;
  if( d1 <= 0.0 && d2 <= 0.0 ) { bc[0]=1.0; bc[1]=0.0; bc[2]=0.0; return; }   /* A */

  const double bp0=p[0]-B[0], bp1=p[1]-B[1], bp2=p[2]-B[2];
  const double d3 = ab0*bp0 + ab1*bp1 + ab2*bp2;
  const double d4 = ac0*bp0 + ac1*bp1 + ac2*bp2;
  if( d3 >= 0.0 && d4 <= d3 ) { bc[0]=0.0; bc[1]=1.0; bc[2]=0.0; return; }    /* B */

  const double vc = d1*d4 - d3*d2;
  if( vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 ) {                                 /* AB */
    const double den = d1 - d3;
    const double t = ( den != 0.0 ) ? d1/den : 0.0;
    bc[0]=1.0-t; bc[1]=t; bc[2]=0.0; return;
  }

  const double cq0=p[0]-C[0], cq1=p[1]-C[1], cq2=p[2]-C[2];
  const double d5 = ab0*cq0 + ab1*cq1 + ab2*cq2;
  const double d6 = ac0*cq0 + ac1*cq1 + ac2*cq2;
  if( d6 >= 0.0 && d5 <= d6 ) { bc[0]=0.0; bc[1]=0.0; bc[2]=1.0; return; }    /* C */

  const double vb = d5*d2 - d1*d6;
  if( vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 ) {                                 /* AC */
    const double den = d2 - d6;
    const double t = ( den != 0.0 ) ? d2/den : 0.0;
    bc[0]=1.0-t; bc[1]=0.0; bc[2]=t; return;
  }

  const double va = d3*d6 - d5*d4;
  if( va <= 0.0 && (d4-d3) >= 0.0 && (d5-d6) >= 0.0 ) {                       /* BC */
    const double den = (d4-d3) + (d5-d6);
    const double t = ( den != 0.0 ) ? (d4-d3)/den : 0.0;
    bc[0]=0.0; bc[1]=1.0-t; bc[2]=t; return;
  }

  const double den = va + vb + vc;                                           /* interior */
  if( den != 0.0 && std::isfinite( den ) ) {
    const double v = vb/den, w = vc/den;
    bc[0]=1.0-v-w; bc[1]=v; bc[2]=w; return;
  }

  /* fully degenerate: nearest of the three edges, emit its 1-D parametrization */
  double c1[3], c2[3], c3[3];
  const double q1 = d2Seg( A, B, p, c1 );
  const double q2 = d2Seg( A, C, p, c2 );
  const double q3 = d2Seg( B, C, p, c3 );
  if( q1 <= q2 && q1 <= q3 ) {
    const double L=(ab0*ab0+ab1*ab1+ab2*ab2); double t=(L>0)?(d1/L):0.0; if(t<0)t=0; else if(t>1)t=1;
    bc[0]=1.0-t; bc[1]=t; bc[2]=0.0;
  } else if( q2 <= q3 ) {
    const double L=(ac0*ac0+ac1*ac1+ac2*ac2); double t=(L>0)?(d2/L):0.0; if(t<0)t=0; else if(t>1)t=1;
    bc[0]=1.0-t; bc[1]=0.0; bc[2]=t;
  } else {
    const double bcx=C[0]-B[0], bcy=C[1]-B[1], bcz=C[2]-B[2];
    const double L=bcx*bcx+bcy*bcy+bcz*bcz;
    double t=(L>0)?(((p[0]-B[0])*bcx+(p[1]-B[1])*bcy+(p[2]-B[2])*bcz)/L):0.0; if(t<0)t=0; else if(t>1)t=1;
    bc[0]=0.0; bc[1]=1.0-t; bc[2]=t;
  }
}

/* enforce the barycentric invariants as hard as floating point allows, WITHOUT
 * disturbing the region-exact structure. Three steps, each chosen so it cannot
 * degrade what the region cascade already got right:
 *
 *   1) clamp negatives to 0 (rounding at a region boundary, or the tet interior
 *      tolerance, can dip a hair below 0). Exact zeros survive untouched.
 *
 *   2) normalize ONLY if the weights do not already sum to 1, and normalize by
 *      DIVISION. Dividing is ONE rounding per component; multiplying by a
 *      precomputed 1/s is TWO (the reciprocal, then the product). This runs once
 *      per query on the winner only, so buying speed with accuracy here is the
 *      wrong trade. And the `s != 1.0` guard means every vertex hit and most
 *      edge hits are left BIT-EXACT instead of being rescaled by rounding noise.
 *
 *   3) make the sum EXACTLY 1 by pushing the residue into the LARGEST weight.
 *      Step 2 alone leaves the sum ~1 ulp off, because each quotient rounds
 *      independently. Two details matter and both were measured:
 *        - the residue goes to the LARGEST weight, never to a small one and
 *          never to a zero: that is where a 1-ulp fix is relatively smallest,
 *          and it is the only choice that cannot destroy a structural zero.
 *        - the residue is measured against the sum in the SAME left-to-right
 *          grouping the caller will use ( ((b0+b1)+b2)+b3 ), and corrected in a
 *          tiny fixed-point loop. Computing "1 - sum of the others" in a
 *          different grouping looks equivalent but is not: it fixed triangles
 *          and made 4-component tet weights WORSE, because the two groupings
 *          round differently. Targeting the actual sum fixes both.
 *      It is a no-op whenever the structure was already consistent.
 *
 * A pathological all-nonpositive vector (never happens for a real closest
 * point) snaps to the first vertex. */
static inline void finalizeBC( double bc[4] )
{
  if( bc[0] < 0.0 ) bc[0] = 0.0;   if( bc[1] < 0.0 ) bc[1] = 0.0;
  if( bc[2] < 0.0 ) bc[2] = 0.0;   if( bc[3] < 0.0 ) bc[3] = 0.0;
  const double s = bc[0] + bc[1] + bc[2] + bc[3];
  if( !( s > 0.0 ) ) { bc[0]=1.0; bc[1]=0.0; bc[2]=0.0; bc[3]=0.0; return; }
  if( s != 1.0 ) { bc[0]/=s; bc[1]/=s; bc[2]/=s; bc[3]/=s; }
  int k = 0;
  if( bc[1] > bc[k] ) k = 1;
  if( bc[2] > bc[k] ) k = 2;
  if( bc[3] > bc[k] ) k = 3;
  for( int it = 0; it < 3; ++it ) {
    const double t = ( ( bc[0] + bc[1] ) + bc[2] ) + bc[3];   /* el orden del que suma */
    if( t == 1.0 ) break;
    const double nk = bc[k] + ( 1.0 - t );
    if( !( nk > 0.0 ) ) break;           /* no empujar la mayor a <=0: inalcanzable */
    bc[k] = nk;
  }
}

/* region-exact barycentrics of p's closest point on packed element (type k),
 * padded to 4 components (unused = 0). Caller applies finalizeBC. */
static void bcElem( const double* v, int k, const double* p, double bc[4] )
{
  bc[0]=bc[1]=bc[2]=bc[3]=0.0;
  switch( k ) {
    case 1: bc[0]=1.0; return;
    case 2: {
      const double* A=v; const double* B=v+3;
      const double abx=B[0]-A[0], aby=B[1]-A[1], abz=B[2]-A[2];
      const double L2=abx*abx+aby*aby+abz*abz;
      double t = (p[0]-A[0])*abx + (p[1]-A[1])*aby + (p[2]-A[2])*abz;
      t = ( L2 > 0.0 ) ? t/L2 : 0.0;
      if( t<0.0 ) t=0.0; else if( t>1.0 ) t=1.0;
      bc[0]=1.0-t; bc[1]=t; return;
    }
    case 3: { double b3[3]; bcTri3( v, v+3, v+6, p, b3 );
              bc[0]=b3[0]; bc[1]=b3[1]; bc[2]=b3[2]; return; }
    case 4: {
      const double* A=v; const double* B=v+3; const double* C=v+6; const double* D=v+9;
      const double BA[3]={B[0]-A[0],B[1]-A[1],B[2]-A[2]};
      const double CA[3]={C[0]-A[0],C[1]-A[1],C[2]-A[2]};
      const double DA[3]={D[0]-A[0],D[1]-A[1],D[2]-A[2]};
      const double d0 = det3( BA, CA, DA );
      if( d0 != 0.0 && std::isfinite( d0 ) ) {
        const double Bp[3]={B[0]-p[0],B[1]-p[1],B[2]-p[2]};
        const double Cp[3]={C[0]-p[0],C[1]-p[1],C[2]-p[2]};
        const double Dp[3]={D[0]-p[0],D[1]-p[1],D[2]-p[2]};
        const double pA[3]={p[0]-A[0],p[1]-A[1],p[2]-A[2]};
        const double l1 = det3( Bp, Cp, Dp ) / d0;
        const double l2 = det3( pA, CA, DA ) / d0;
        const double l3 = det3( BA, pA, DA ) / d0;
        const double l4 = 1.0 - l1 - l2 - l3;
        const double tol = -1e-12;
        if( l1 >= tol && l2 >= tol && l3 >= tol && l4 >= tol ) {
          bc[0]=l1; bc[1]=l2; bc[2]=l3; bc[3]=l4; return;
        }
      }
      /* outside/degenerate: nearest face, map its triangle bc onto the 4 verts */
      double c[3], cb[3], b3[3];  int face = 0;
      double q  = d2Tri( A, B, C, p, c  );
      double q2 = d2Tri( A, B, D, p, cb ); if( q2 < q ) { q=q2; face=1; }
      q2        = d2Tri( A, C, D, p, cb ); if( q2 < q ) { q=q2; face=2; }
      q2        = d2Tri( B, C, D, p, cb ); if( q2 < q ) { q=q2; face=3; }
      switch( face ) {
        case 0: bcTri3(A,B,C,p,b3); bc[0]=b3[0]; bc[1]=b3[1]; bc[2]=b3[2]; break;
        case 1: bcTri3(A,B,D,p,b3); bc[0]=b3[0]; bc[1]=b3[1]; bc[3]=b3[2]; break;
        case 2: bcTri3(A,C,D,p,b3); bc[0]=b3[0]; bc[2]=b3[1]; bc[3]=b3[2]; break;
        case 3: bcTri3(B,C,D,p,b3); bc[1]=b3[0]; bc[2]=b3[1]; bc[3]=b3[2]; break;
      }
      return;
    }
  }
}

/* ------------------------------------------------- 4-wide triangle kernel */

static int g_avx = -1;
static bool useAVX( void )
{
  if( g_avx < 0 ) {
#if defined(_MSC_VER)
    int ci[4];  __cpuid( ci, 1 );
    const bool osxsave = ( ci[2] >> 27 ) & 1;
    const bool avx     = ( ci[2] >> 28 ) & 1;
    g_avx = ( osxsave && avx && ( ( _xgetbv(0) & 6 ) == 6 ) ) ? 1 : 0;
#else
    g_avx = __builtin_cpu_supports( "avx" ) ? 1 : 0;
#endif
  }
  return g_avx == 1;
}

static inline __m256d mm_dot3( __m256d ax, __m256d ay, __m256d az,
                               __m256d bx, __m256d by, __m256d bz )
{
  return _mm256_add_pd( _mm256_add_pd( _mm256_mul_pd(ax,bx), _mm256_mul_pd(ay,by) ),
                        _mm256_mul_pd(az,bz) );
}

/* branchless Ericson over one PreTri4 block: (v0,e1,e2) in 4-wide SoA --
 * ALIGNED sequential loads, no marshalling, edges precomputed at build.
 * Same region order and same per-region arithmetic as the scalar d2Tri, so
 * selected lanes match it to the last bit; genuinely degenerate lanes are
 * redone with the scalar fallback. Null-padding lanes (id 0) produce a finite
 * garbage result that the CALLER filters by id. */
static void tri4blk( const double* blk, const double* p,
                     double d2o[4], double cpo[12] )
{
  const __m256d Ax  = _mm256_loadu_pd( blk      ), Ay  = _mm256_loadu_pd( blk +  4 ), Az  = _mm256_loadu_pd( blk +  8 );
  const __m256d abx = _mm256_loadu_pd( blk + 12 ), aby = _mm256_loadu_pd( blk + 16 ), abz = _mm256_loadu_pd( blk + 20 );
  const __m256d acx = _mm256_loadu_pd( blk + 24 ), acy = _mm256_loadu_pd( blk + 28 ), acz = _mm256_loadu_pd( blk + 32 );
  const __m256d Bx = _mm256_add_pd( Ax, abx ), By = _mm256_add_pd( Ay, aby ), Bz = _mm256_add_pd( Az, abz );
  const __m256d Cx = _mm256_add_pd( Ax, acx ), Cy = _mm256_add_pd( Ay, acy ), Cz = _mm256_add_pd( Az, acz );
  const __m256d px = _mm256_set1_pd( p[0] ), py = _mm256_set1_pd( p[1] ), pz = _mm256_set1_pd( p[2] );
  const __m256d zero = _mm256_setzero_pd(), one = _mm256_set1_pd( 1.0 );

  const __m256d apx = _mm256_sub_pd(px,Ax), apy = _mm256_sub_pd(py,Ay), apz = _mm256_sub_pd(pz,Az);
  const __m256d d1 = mm_dot3( abx,aby,abz, apx,apy,apz );
  const __m256d d2 = mm_dot3( acx,acy,acz, apx,apy,apz );
  const __m256d bpx = _mm256_sub_pd(px,Bx), bpy = _mm256_sub_pd(py,By), bpz = _mm256_sub_pd(pz,Bz);
  const __m256d d3 = mm_dot3( abx,aby,abz, bpx,bpy,bpz );
  const __m256d d4 = mm_dot3( acx,acy,acz, bpx,bpy,bpz );
  const __m256d cqx = _mm256_sub_pd(px,Cx), cqy = _mm256_sub_pd(py,Cy), cqz = _mm256_sub_pd(pz,Cz);
  const __m256d d5 = mm_dot3( abx,aby,abz, cqx,cqy,cqz );
  const __m256d d6 = mm_dot3( acx,acy,acz, cqx,cqy,cqz );
  const __m256d vc = _mm256_sub_pd( _mm256_mul_pd(d1,d4), _mm256_mul_pd(d3,d2) );
  const __m256d vb = _mm256_sub_pd( _mm256_mul_pd(d5,d2), _mm256_mul_pd(d1,d6) );
  const __m256d va = _mm256_sub_pd( _mm256_mul_pd(d3,d6), _mm256_mul_pd(d5,d4) );

#define LE(a,b) _mm256_cmp_pd( a, b, _CMP_LE_OQ )
#define GE(a,b) _mm256_cmp_pd( a, b, _CMP_GE_OQ )
  const __m256d mA  = _mm256_and_pd( LE(d1,zero), LE(d2,zero) );
  const __m256d mB  = _mm256_and_pd( GE(d3,zero), LE(d4,d3) );
  const __m256d mAB = _mm256_and_pd( _mm256_and_pd( LE(vc,zero), GE(d1,zero) ), LE(d3,zero) );
  const __m256d mC  = _mm256_and_pd( GE(d6,zero), LE(d5,d6) );
  const __m256d mAC = _mm256_and_pd( _mm256_and_pd( LE(vb,zero), GE(d2,zero) ), LE(d6,zero) );
  const __m256d d43 = _mm256_sub_pd( d4, d3 ), d56 = _mm256_sub_pd( d5, d6 );
  const __m256d mBC = _mm256_and_pd( _mm256_and_pd( LE(va,zero), GE(d43,zero) ), GE(d56,zero) );
#undef LE
#undef GE

  /* safe divisions (unselected lanes discard the garbage) */
  const __m256d eqz1 = _mm256_cmp_pd( _mm256_sub_pd(d1,d3), zero, _CMP_EQ_OQ );
  const __m256d denAB = _mm256_blendv_pd( _mm256_sub_pd(d1,d3), one, eqz1 );
  const __m256d tAB = _mm256_div_pd( d1, denAB );
  const __m256d eqz2 = _mm256_cmp_pd( _mm256_sub_pd(d2,d6), zero, _CMP_EQ_OQ );
  const __m256d denAC = _mm256_blendv_pd( _mm256_sub_pd(d2,d6), one, eqz2 );
  const __m256d tAC = _mm256_div_pd( d2, denAC );
  const __m256d sBC = _mm256_add_pd( d43, d56 );
  const __m256d eqz3 = _mm256_cmp_pd( sBC, zero, _CMP_EQ_OQ );
  const __m256d tBC = _mm256_div_pd( d43, _mm256_blendv_pd( sBC, one, eqz3 ) );
  const __m256d denI = _mm256_add_pd( _mm256_add_pd( va, vb ), vc );
  const __m256d eqzI = _mm256_cmp_pd( denI, zero, _CMP_EQ_OQ );
  const __m256d rden = _mm256_div_pd( one, _mm256_blendv_pd( denI, one, eqzI ) );
  const __m256d vI = _mm256_mul_pd( vb, rden ), wI = _mm256_mul_pd( vc, rden );

  /* priority cascade: A, B, AB, C, AC, BC, interior (Ericson order) */
  __m256d cx = Ax, cy = Ay, cz = Az;
  __m256d done = mA;
  __m256d sel;
#define PICK(m,qx,qy,qz) \
  sel = _mm256_andnot_pd( done, m ); \
  cx = _mm256_blendv_pd( cx, qx, sel ); \
  cy = _mm256_blendv_pd( cy, qy, sel ); \
  cz = _mm256_blendv_pd( cz, qz, sel ); \
  done = _mm256_or_pd( done, m );
  PICK( mB, Bx, By, Bz );
  PICK( mAB, _mm256_add_pd(Ax,_mm256_mul_pd(tAB,abx)),
             _mm256_add_pd(Ay,_mm256_mul_pd(tAB,aby)),
             _mm256_add_pd(Az,_mm256_mul_pd(tAB,abz)) );
  PICK( mC, Cx, Cy, Cz );
  PICK( mAC, _mm256_add_pd(Ax,_mm256_mul_pd(tAC,acx)),
             _mm256_add_pd(Ay,_mm256_mul_pd(tAC,acy)),
             _mm256_add_pd(Az,_mm256_mul_pd(tAC,acz)) );
  PICK( mBC, _mm256_add_pd(Bx,_mm256_mul_pd(tBC,_mm256_sub_pd(Cx,Bx))),
             _mm256_add_pd(By,_mm256_mul_pd(tBC,_mm256_sub_pd(Cy,By))),
             _mm256_add_pd(Bz,_mm256_mul_pd(tBC,_mm256_sub_pd(Cz,Bz))) );
  /* interior: everything not yet done */
  sel = _mm256_andnot_pd( done, _mm256_castsi256_pd( _mm256_set1_epi64x( -1 ) ) );
  cx = _mm256_blendv_pd( cx, _mm256_add_pd(Ax,_mm256_add_pd(_mm256_mul_pd(vI,abx),_mm256_mul_pd(wI,acx))), sel );
  cy = _mm256_blendv_pd( cy, _mm256_add_pd(Ay,_mm256_add_pd(_mm256_mul_pd(vI,aby),_mm256_mul_pd(wI,acy))), sel );
  cz = _mm256_blendv_pd( cz, _mm256_add_pd(Az,_mm256_add_pd(_mm256_mul_pd(vI,abz),_mm256_mul_pd(wI,acz))), sel );

  const __m256d dx = _mm256_sub_pd( px, cx ), dy = _mm256_sub_pd( py, cy ), dz = _mm256_sub_pd( pz, cz );
  const __m256d dd = mm_dot3( dx,dy,dz, dx,dy,dz );

  _mm256_storeu_pd( d2o, dd );
  double cxx[4], cyy[4], czz[4];
  _mm256_storeu_pd( cxx, cx );  _mm256_storeu_pd( cyy, cy );  _mm256_storeu_pd( czz, cz );
  for( int l = 0; l < 4; ++l ) { cpo[3*l]=cxx[l]; cpo[3*l+1]=cyy[l]; cpo[3*l+2]=czz[l]; }

  /* degenerate lanes (interior selected with a zeroed denominator): redo scalar */
  const int degm = _mm256_movemask_pd( _mm256_and_pd( sel, eqzI ) );
  if( degm ) {
    for( int l = 0; l < 4; ++l )
      if( degm & (1<<l) ) {
        const double A3[3] = { blk[   l], blk[ 4+l], blk[ 8+l] };
        const double B3[3] = { A3[0]+blk[12+l], A3[1]+blk[16+l], A3[2]+blk[20+l] };
        const double C3[3] = { A3[0]+blk[24+l], A3[1]+blk[28+l], A3[2]+blk[32+l] };
        d2o[l] = d2Tri( A3, B3, C3, p, cpo+3*l );
      }
  }
#undef PICK
}

/* ---------------------------------------------------------------- Morton */

static inline uint32_t expandBits( uint32_t v )
{
  v = ( v * 0x00010001u ) & 0xFF0000FFu;
  v = ( v * 0x00000101u ) & 0x0F00F00Fu;
  v = ( v * 0x00000011u ) & 0xC30C30C3u;
  v = ( v * 0x00000005u ) & 0x49249249u;
  return v;
}

#endif /* BVH_KERNELS_H */
