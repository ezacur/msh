/* BVH_mx  --  binned-SAH BVH4 builder (geometry-agnostic core).
 *
 *   [bounds4, child4, srange, perm, pk4, pk4id, s4] = ...
 *        BVH_mx( eC , eR , eLo , eHi , minLeaf , maxLeaf , vol , eV , eN )
 *
 *   pk4/pk4id/s4: the PreTri4 leaf pool -- every ALL-TRIANGLE leaf slot gets
 *   its triangles preprocessed as (v0, e1, e2) in 4-wide SoA blocks of
 *   36 doubles (v0x*4, v0y*4, v0z*4, e1x*4 ... e2z*4), PADDED to multiples of
 *   4 with null lanes (id 0, zero triangle: a natural miss for Moller-Trumbore
 *   and id-filtered in closest-point). pk4id holds the 1-based PACKED element
 *   position per lane (0 = padding); s4 (int32 [8 x nN4]) holds per slot
 *   [blockStart; blockCount] (1-based, count 0 = no blocks: non-triangle or
 *   internal slot).
 *
 *     eC   nE x 3  double  element bounding-sphere centers (also the SAH
 *                          centroids used for binning)
 *     eR   nE x 1  double  element bounding-sphere radii
 *     eLo  nE x 3  double  element AABB lower corners
 *     eHi  nE x 3  double  element AABB upper corners
 *     minLeaf, maxLeaf     adaptive leaves: n <= minLeaf is always a leaf; the
 *                          SAH may keep up to maxLeaf together when splitting
 *                          does not pay
 *     vol  1..4            node volume type: 1 = spheres, 2 = AABBs,
 *                          3 = OBBs (per-slot PCA over the contained element
 *                          vertices), 4 = 14-DOPs (AABB + the 4 diagonal
 *                          direction pairs (1,±1,±1))
 *     eV   12 x nE double  element vertices, zero-padded (REQUIRED for 3|4)
 *     eN   nE x 1  double  nonzero-vertex count per element (REQUIRED for 3|4)
 *
 *   The PARTITION is the same for every volume type (single hierarchy, see
 *   msh_QUERY_ENGINE_DESIGN.md): binned SAH (16 bins on the largest
 *   centroid-extent axis, AABB surface-area metric), then the binary tree is
 *   COLLAPSED into a 4-ary BVH (each node adopts up to 4 grandchildren,
 *   expanding the largest child first). Only the requested node VOLUMES are
 *   emitted, computed exactly from the elements in each slot's range.
 *
 *   Outputs (nN4 nodes, column-major, one column per node; every slot field is
 *   a GROUP OF 4 consecutive floats, one per slot):
 *     bounds4  single [S x nN4], S by volume:
 *       spheres (16): cx,cy,cz,r
 *       AABB    (24): lox,hix,loy,hiy,loz,hiz
 *       OBB     (60): a1(3),a2(3),a3(3) axes + lo(3),hi(3) axis-coords bounds.
 *                     The axes are rounded to float FIRST and the bounds are
 *                     then computed against the ROUNDED axes, so containment
 *                     in the test space is exact.
 *       14-DOP  (56): the 24 AABB rows (same layout as vol 2) + lo,hi pairs
 *                     for the diagonals (1,1,1),(1,-1,1),(1,1,-1),(1,-1,-1)
 *                     (unnormalized projections x±y±z).
 *              CONSERVATIVELY rounded outward (float bounds, exact double
 *              geometry inside them). Empty slots: never-hit markers.
 *     child4   int32 [4 x nN4]   > 0 internal child (1-based), -1 leaf, 0 empty
 *     srange   int32 [8 x nN4]   per-slot element range [lo;hi] (1-based) into
 *                                perm -- present for ALL slots (refits use it)
 *     perm     int32 [nE x 1]    element permutation (leaf ranges contiguous)
 *
 *   Nodes are emitted parent-before-children (DFS preorder): reverse id order
 *   is a valid bottom-up sweep.
 *
 *   Compile:  mex BVH_mx.cpp
 *
 * See also BVH, bvhClosestElement_mx.
 */

#include "mex.h"
#include <cmath>
#include <cstdint>
#include <cstring>
#include <vector>
#include <algorithm>
#include <limits>

static const double INF = std::numeric_limits<double>::infinity();

#define SAH_BINS  16
#define SAH_CTRAV 1.0
#define SAH_CPRIM 1.0

struct BN { int lo, hi, l, r; };   /* binary node: range into perm, children (-1 = leaf) */

/* outward float rounding: the float interval must CONTAIN the double value */
static inline float fdown( double v )
{
  float f = (float)v;
  if( (double)f > v ) f = std::nextafterf( f, -std::numeric_limits<float>::infinity() );
  return f;
}
static inline float fup( double v )
{
  float f = (float)v;
  if( (double)f < v ) f = std::nextafterf( f,  std::numeric_limits<float>::infinity() );
  return f;
}

static inline double boxArea( const double lo[3], const double hi[3] )
{
  const double dx = hi[0]-lo[0], dy = hi[1]-lo[1], dz = hi[2]-lo[2];
  if( dx < 0.0 || dy < 0.0 || dz < 0.0 ) return 0.0;   /* empty */
  return 2.0*( dx*dy + dx*dz + dy*dz );
}

/* eigenvectors of a symmetric 3x3 (cyclic Jacobi): columns of V, sorted by
 * eigenvalue DESCENDING, right-handed. Robust for rank-deficient inputs. */
static void eig3( const double Cin[9], double V[9] )
{
  double A[9];  for( int i = 0; i < 9; ++i ) A[i] = Cin[i];
  for( int i = 0; i < 9; ++i ) V[i] = ( i%4 == 0 ) ? 1.0 : 0.0;
  for( int sweep = 0; sweep < 50; ++sweep ) {
    /* largest off-diagonal */
    int p = 0, q = 1;  double mx = std::fabs( A[1] );
    if( std::fabs( A[2] ) > mx ) { mx = std::fabs( A[2] ); p = 0; q = 2; }
    if( std::fabs( A[5] ) > mx ) { mx = std::fabs( A[5] ); p = 1; q = 2; }
    const double scale = std::fabs(A[0]) + std::fabs(A[4]) + std::fabs(A[8]);
    if( mx <= 1e-15 * ( scale > 0.0 ? scale : 1.0 ) ) break;
    const double app = A[ p + 3*p ], aqq = A[ q + 3*q ], apq = A[ p + 3*q ];
    const double theta = 0.5*( aqq - app ) / apq;
    const double t = ( theta >= 0.0 ? 1.0 : -1.0 ) /
                     ( std::fabs(theta) + std::sqrt( theta*theta + 1.0 ) );
    const double c = 1.0/std::sqrt( t*t + 1.0 ), s = t*c;
    for( int k = 0; k < 3; ++k ) {                     /* rotate A (sym) */
      const double akp = A[ k + 3*p ], akq = A[ k + 3*q ];
      A[ k + 3*p ] = c*akp - s*akq;
      A[ k + 3*q ] = s*akp + c*akq;
    }
    for( int k = 0; k < 3; ++k ) {
      const double apk = A[ p + 3*k ], aqk = A[ q + 3*k ];
      A[ p + 3*k ] = c*apk - s*aqk;
      A[ q + 3*k ] = s*apk + c*aqk;
    }
    for( int k = 0; k < 3; ++k ) {                     /* accumulate V */
      const double vkp = V[ k + 3*p ], vkq = V[ k + 3*q ];
      V[ k + 3*p ] = c*vkp - s*vkq;
      V[ k + 3*q ] = s*vkp + c*vkq;
    }
  }
  double ev[3] = { A[0], A[4], A[8] };
  int ord[3] = { 0, 1, 2 };                            /* sort descending */
  if( ev[ord[0]] < ev[ord[1]] ) std::swap( ord[0], ord[1] );
  if( ev[ord[1]] < ev[ord[2]] ) std::swap( ord[1], ord[2] );
  if( ev[ord[0]] < ev[ord[1]] ) std::swap( ord[0], ord[1] );
  double W[9];
  for( int j = 0; j < 3; ++j )
    for( int k = 0; k < 3; ++k ) W[ k + 3*j ] = V[ k + 3*ord[j] ];
  /* right-handed: a3 = a1 x a2 if needed */
  const double det =
      W[0]*( W[4]*W[8] - W[7]*W[5] )
    - W[3]*( W[1]*W[8] - W[7]*W[2] )
    + W[6]*( W[1]*W[5] - W[4]*W[2] );
  if( det < 0.0 ) { W[6] = -W[6]; W[7] = -W[7]; W[8] = -W[8]; }
  for( int i = 0; i < 9; ++i ) V[i] = W[i];
}

/* ------------------------------------------------------------------ builder */
struct Builder {
  const double *eC, *eR, *eLo, *eHi;
  const double *eV;      /* 12 x nE element vertices (vol 3|4), zero-padded */
  const double *eN;      /* nE     nonzero-vertex counts (vol 3|4)          */
  mwSize nE;
  int minLeaf, maxLeaf, vol;

  std::vector<int> perm;
  std::vector<BN>  bn;

  /* BVH4 output pools */
  std::vector<float>   b4;    /* S per node (16 or 24)  */
  std::vector<int32_t> c4;    /* 4 per node             */
  std::vector<int32_t> r4;    /* 8 per node             */
  int S;                      /* floats per node column */

  int buildBin( int lo, int hi )
  {
    const int nid = (int)bn.size();
    bn.push_back( { lo, hi, -1, -1 } );
    const int n = hi - lo + 1;
    if( n <= minLeaf ) return nid;

    /* centroid bounds over the range */
    double cmn[3] = {  INF,  INF,  INF };
    double cmx[3] = { -INF, -INF, -INF };
    for( int i = lo; i <= hi; ++i ) {
      const int e = perm[i];
      for( int a = 0; a < 3; ++a ) {
        const double c = eC[ e + a*nE ];
        if( c < cmn[a] ) cmn[a] = c;
        if( c > cmx[a] ) cmx[a] = c;
      }
    }
    int axis = 0;
    double ext = cmx[0]-cmn[0];
    if( cmx[1]-cmn[1] > ext ) { axis = 1; ext = cmx[1]-cmn[1]; }
    if( cmx[2]-cmn[2] > ext ) { axis = 2; ext = cmx[2]-cmn[2]; }

    int mid = -1;
    if( !( ext > 0.0 ) ) {                    /* coincident centroids */
      if( n <= maxLeaf ) return nid;
      mid = lo + n/2;                         /* arbitrary median split */
    } else {
      /* --- binned SAH on `axis` --- */
      int    bcnt[SAH_BINS];
      double blo[SAH_BINS][3], bhi[SAH_BINS][3];
      for( int b = 0; b < SAH_BINS; ++b ) {
        bcnt[b] = 0;
        for( int a = 0; a < 3; ++a ) { blo[b][a] = INF; bhi[b][a] = -INF; }
      }
      const double scale = (double)SAH_BINS * ( 1.0 - 1e-9 ) / ext;
      for( int i = lo; i <= hi; ++i ) {
        const int e = perm[i];
        int b = (int)( ( eC[ e + axis*nE ] - cmn[axis] ) * scale );
        if( b < 0 ) b = 0;  if( b >= SAH_BINS ) b = SAH_BINS-1;
        bcnt[b]++;
        for( int a = 0; a < 3; ++a ) {
          const double l = eLo[ e + a*nE ], h = eHi[ e + a*nE ];
          if( l < blo[b][a] ) blo[b][a] = l;
          if( h > bhi[b][a] ) bhi[b][a] = h;
        }
      }
      /* suffix (right) accumulations */
      double rlo[SAH_BINS][3], rhi[SAH_BINS][3];  int rcnt[SAH_BINS];
      {
        double alo[3] = {  INF,  INF,  INF }, ahi[3] = { -INF, -INF, -INF };
        int    ac = 0;
        for( int b = SAH_BINS-1; b >= 0; --b ) {
          ac += bcnt[b];
          for( int a = 0; a < 3; ++a ) {
            if( blo[b][a] < alo[a] ) alo[a] = blo[b][a];
            if( bhi[b][a] > ahi[a] ) ahi[a] = bhi[b][a];
            rlo[b][a] = alo[a];  rhi[b][a] = ahi[a];
          }
          rcnt[b] = ac;
        }
      }
      /* prefix sweep, evaluate split after bin b (left = bins 0..b) */
      double bestCost = INF;  int bestB = -1;
      {
        const double saP = boxArea( rlo[0], rhi[0] );
        double alo[3] = {  INF,  INF,  INF }, ahi[3] = { -INF, -INF, -INF };
        int    ac = 0;
        if( saP > 0.0 ) {
          for( int b = 0; b < SAH_BINS-1; ++b ) {
            ac += bcnt[b];
            for( int a = 0; a < 3; ++a ) {
              if( blo[b][a] < alo[a] ) alo[a] = blo[b][a];
              if( bhi[b][a] > ahi[a] ) ahi[a] = bhi[b][a];
            }
            if( ac == 0 || rcnt[b+1] == 0 ) continue;
            const double cost = SAH_CTRAV + SAH_CPRIM *
              ( boxArea(alo,ahi)*ac + boxArea(rlo[b+1],rhi[b+1])*rcnt[b+1] ) / saP;
            if( cost < bestCost ) { bestCost = cost; bestB = b; }
          }
        }
      }
      if( n <= maxLeaf && (double)n * SAH_CPRIM <= bestCost ) return nid;  /* leaf pays */

      if( bestB >= 0 ) {                      /* partition by bin id */
        int i = lo, j = hi;
        while( i <= j ) {
          const int e = perm[i];
          int b = (int)( ( eC[ e + axis*nE ] - cmn[axis] ) * scale );
          if( b < 0 ) b = 0;  if( b >= SAH_BINS ) b = SAH_BINS-1;
          if( b <= bestB ) ++i;
          else { std::swap( perm[i], perm[j] ); --j; }
        }
        mid = i;
      }
      if( mid <= lo || mid > hi ) {           /* degenerate: median fallback */
        mid = lo + n/2;
        std::nth_element( perm.begin()+lo, perm.begin()+mid, perm.begin()+hi+1,
          [&]( int a, int b ){ return eC[ a + axis*nE ] < eC[ b + axis*nE ]; } );
      }
    }

    const int L = buildBin( lo,  mid-1 );
    const int R = buildBin( mid, hi    );
    bn[nid].l = L;  bn[nid].r = R;
    return nid;
  }

  /* exact slot volumes over a perm range, conservatively rounded to float */
  void volumize( int col, int k, int lo, int hi )
  {
    if( vol == 3 ) {                         /* OBB: PCA over the slot vertices */
      double mu[3] = { 0.0, 0.0, 0.0 };
      double nv = 0.0;
      for( int i = lo; i <= hi; ++i ) {
        const int e = perm[i];
        const int n = (int)eN[e];
        const double* q = eV + (size_t)e*12;
        for( int j = 0; j < n; ++j )
          for( int a = 0; a < 3; ++a ) mu[a] += q[3*j+a];
        nv += n;
      }
      if( nv > 0.0 ) for( int a = 0; a < 3; ++a ) mu[a] /= nv;
      double Cv[9] = {0,0,0,0,0,0,0,0,0};
      for( int i = lo; i <= hi; ++i ) {
        const int e = perm[i];
        const int n = (int)eN[e];
        const double* q = eV + (size_t)e*12;
        for( int j = 0; j < n; ++j ) {
          const double dx = q[3*j]-mu[0], dy = q[3*j+1]-mu[1], dz = q[3*j+2]-mu[2];
          Cv[0]+=dx*dx; Cv[4]+=dy*dy; Cv[8]+=dz*dz;
          Cv[1]+=dx*dy; Cv[2]+=dx*dz; Cv[5]+=dy*dz;
        }
      }
      Cv[3]=Cv[1]; Cv[6]=Cv[2]; Cv[7]=Cv[5];
      double Vv[9];
      eig3( Cv, Vv );
      /* round the AXES to float FIRST, then bound projections on THOSE axes */
      float af[9];
      double ad[9];
      for( int i = 0; i < 9; ++i ) { af[i] = (float)Vv[i]; ad[i] = (double)af[i]; }
      double mn[3] = { INF, INF, INF }, mx[3] = { -INF, -INF, -INF };
      for( int i = lo; i <= hi; ++i ) {
        const int e = perm[i];
        const int n = (int)eN[e];
        const double* q = eV + (size_t)e*12;
        for( int j = 0; j < n; ++j )
          for( int a = 0; a < 3; ++a ) {
            const double pr = q[3*j]*ad[3*a] + q[3*j+1]*ad[3*a+1] + q[3*j+2]*ad[3*a+2];
            if( pr < mn[a] ) mn[a] = pr;
            if( pr > mx[a] ) mx[a] = pr;
          }
      }
      float* B = &b4[ (size_t)col*S ];
      for( int a = 0; a < 3; ++a ) {
        B[ 4*(3*a  ) + k ] = af[3*a  ];      /* axis a: 3 components */
        B[ 4*(3*a+1) + k ] = af[3*a+1];
        B[ 4*(3*a+2) + k ] = af[3*a+2];
        B[ 4*( 9+a ) + k ] = fdown( mn[a] );
        B[ 4*(12+a ) + k ] = fup(   mx[a] );
      }
      return;
    }
    if( vol == 4 ) {                         /* 14-DOP: AABB + 4 diagonal pairs */
      double mn[7], mx[7];
      for( int a = 0; a < 7; ++a ) { mn[a] = INF; mx[a] = -INF; }
      for( int i = lo; i <= hi; ++i ) {
        const int e = perm[i];
        const int n = (int)eN[e];
        const double* q = eV + (size_t)e*12;
        for( int j = 0; j < n; ++j ) {
          const double x = q[3*j], y = q[3*j+1], z = q[3*j+2];
          const double pr[7] = { x, y, z, x+y+z, x-y+z, x+y-z, x-y-z };
          for( int a = 0; a < 7; ++a ) {
            if( pr[a] < mn[a] ) mn[a] = pr[a];
            if( pr[a] > mx[a] ) mx[a] = pr[a];
          }
        }
      }
      float* B = &b4[ (size_t)col*S ];
      for( int a = 0; a < 3; ++a ) {         /* first 24 rows: AABB layout */
        B[ 4*(2*a  ) + k ] = fdown( mn[a] );
        B[ 4*(2*a+1) + k ] = fup(   mx[a] );
      }
      for( int a = 3; a < 7; ++a ) {         /* diagonals */
        B[ 4*( 2*a ) + k ] = fdown( mn[a] );
        B[ 4*(2*a+1) + k ] = fup(   mx[a] );
      }
      return;
    }
    if( vol == 5 ) {                         /* RSS: PCA rect + swept radius   */
      double mu[3] = { 0.0, 0.0, 0.0 };
      double nv = 0.0;
      for( int i = lo; i <= hi; ++i ) {
        const int e = perm[i];  const int n = (int)eN[e];
        const double* q = eV + (size_t)e*12;
        for( int j = 0; j < n; ++j )
          for( int a = 0; a < 3; ++a ) mu[a] += q[3*j+a];
        nv += n;
      }
      if( nv > 0.0 ) for( int a = 0; a < 3; ++a ) mu[a] /= nv;
      double Cv[9] = {0,0,0,0,0,0,0,0,0};
      for( int i = lo; i <= hi; ++i ) {
        const int e = perm[i];  const int n = (int)eN[e];
        const double* q = eV + (size_t)e*12;
        for( int j = 0; j < n; ++j ) {
          const double dx = q[3*j]-mu[0], dy = q[3*j+1]-mu[1], dz = q[3*j+2]-mu[2];
          Cv[0]+=dx*dx; Cv[4]+=dy*dy; Cv[8]+=dz*dz;
          Cv[1]+=dx*dy; Cv[2]+=dx*dz; Cv[5]+=dy*dz;
        }
      }
      Cv[3]=Cv[1]; Cv[6]=Cv[2]; Cv[7]=Cv[5];
      double Vv[9];  eig3( Cv, Vv );
      float af[9];  double ad[9];
      for( int i = 0; i < 9; ++i ) { af[i] = (float)Vv[i]; ad[i] = (double)af[i]; }
      double mn[3] = { INF, INF, INF }, mx[3] = { -INF, -INF, -INF };
      for( int i = lo; i <= hi; ++i ) {
        const int e = perm[i];  const int n = (int)eN[e];
        const double* q = eV + (size_t)e*12;
        for( int j = 0; j < n; ++j )
          for( int a = 0; a < 3; ++a ) {
            const double pr = q[3*j]*ad[3*a] + q[3*j+1]*ad[3*a+1] + q[3*j+2]*ad[3*a+2];
            if( pr < mn[a] ) mn[a] = pr;
            if( pr > mx[a] ) mx[a] = pr;
          }
      }
      const float w0f = (float)( ( mn[2] + mx[2] )/2 );      /* rect plane      */
      const double rw = std::max( std::fabs( mn[2] - (double)w0f ) ,
                                  std::fabs( mx[2] - (double)w0f ) );
      float* B = &b4[ (size_t)col*S ];
      for( int a = 0; a < 3; ++a ) {
        B[ 4*(3*a  ) + k ] = af[3*a  ];
        B[ 4*(3*a+1) + k ] = af[3*a+1];
        B[ 4*(3*a+2) + k ] = af[3*a+2];
      }
      B[ 4* 9 + k ] = fdown( mn[0] );  B[ 4*10 + k ] = fup( mx[0] );
      B[ 4*11 + k ] = fdown( mn[1] );  B[ 4*12 + k ] = fup( mx[1] );
      B[ 4*13 + k ] = w0f;
      B[ 4*14 + k ] = fup( rw );
      return;
    }
    if( vol == 6 ) {                         /* LSS/capsule: PCA segment + r   */
      double mu[3] = { 0.0, 0.0, 0.0 };
      double nv = 0.0;
      for( int i = lo; i <= hi; ++i ) {
        const int e = perm[i];  const int n = (int)eN[e];
        const double* q = eV + (size_t)e*12;
        for( int j = 0; j < n; ++j )
          for( int a = 0; a < 3; ++a ) mu[a] += q[3*j+a];
        nv += n;
      }
      if( nv > 0.0 ) for( int a = 0; a < 3; ++a ) mu[a] /= nv;
      double Cv[9] = {0,0,0,0,0,0,0,0,0};
      for( int i = lo; i <= hi; ++i ) {
        const int e = perm[i];  const int n = (int)eN[e];
        const double* q = eV + (size_t)e*12;
        for( int j = 0; j < n; ++j ) {
          const double dx = q[3*j]-mu[0], dy = q[3*j+1]-mu[1], dz = q[3*j+2]-mu[2];
          Cv[0]+=dx*dx; Cv[4]+=dy*dy; Cv[8]+=dz*dz;
          Cv[1]+=dx*dy; Cv[2]+=dx*dz; Cv[5]+=dy*dz;
        }
      }
      Cv[3]=Cv[1]; Cv[6]=Cv[2]; Cv[7]=Cv[5];
      double Vv[9];  eig3( Cv, Vv );
      const double a1[3] = { Vv[0], Vv[1], Vv[2] };          /* dominant axis   */
      double tmn = INF, tmx = -INF;
      for( int i = lo; i <= hi; ++i ) {
        const int e = perm[i];  const int n = (int)eN[e];
        const double* q = eV + (size_t)e*12;
        for( int j = 0; j < n; ++j ) {
          const double t = ( q[3*j]-mu[0] )*a1[0] + ( q[3*j+1]-mu[1] )*a1[1] + ( q[3*j+2]-mu[2] )*a1[2];
          if( t < tmn ) tmn = t;
          if( t > tmx ) tmx = t;
        }
      }
      float P0f[3], P1f[3];
      for( int a = 0; a < 3; ++a ) {
        P0f[a] = (float)( mu[a] + tmn*a1[a] );
        P1f[a] = (float)( mu[a] + tmx*a1[a] );
      }
      /* radius against the ROUNDED float segment: exact containment */
      const double Q0[3] = { P0f[0], P0f[1], P0f[2] };
      const double Q1[3] = { P1f[0], P1f[1], P1f[2] };
      const double sv[3] = { Q1[0]-Q0[0], Q1[1]-Q0[1], Q1[2]-Q0[2] };
      const double svv = sv[0]*sv[0] + sv[1]*sv[1] + sv[2]*sv[2];
      double r2 = 0.0;
      for( int i = lo; i <= hi; ++i ) {
        const int e = perm[i];  const int n = (int)eN[e];
        const double* q = eV + (size_t)e*12;
        for( int j = 0; j < n; ++j ) {
          double t = ( q[3*j]-Q0[0] )*sv[0] + ( q[3*j+1]-Q0[1] )*sv[1] + ( q[3*j+2]-Q0[2] )*sv[2];
          t = ( svv > 0.0 ) ? t/svv : 0.0;
          if( t < 0.0 ) t = 0.0;  if( t > 1.0 ) t = 1.0;
          const double cx = Q0[0]+t*sv[0], cy = Q0[1]+t*sv[1], cz = Q0[2]+t*sv[2];
          const double dd = ( q[3*j]-cx )*( q[3*j]-cx ) + ( q[3*j+1]-cy )*( q[3*j+1]-cy ) + ( q[3*j+2]-cz )*( q[3*j+2]-cz );
          if( dd > r2 ) r2 = dd;
        }
      }
      float* B = &b4[ (size_t)col*S ];
      for( int a = 0; a < 3; ++a ) {
        B[ 4*a     + k ] = P0f[a];
        B[ 4*(3+a) + k ] = P1f[a];
      }
      B[ 4*6 + k ] = fup( std::sqrt( r2 ) );
      return;
    }
    if( vol == 1 ) {                         /* sphere: bbox center + max(d+r) */
      double mn[3] = { INF, INF, INF }, mx[3] = { -INF, -INF, -INF };
      for( int i = lo; i <= hi; ++i ) {
        const int e = perm[i];
        for( int a = 0; a < 3; ++a ) {
          const double c = eC[ e + a*nE ];
          if( c < mn[a] ) mn[a] = c;
          if( c > mx[a] ) mx[a] = c;
        }
      }
      const double c[3] = { (mn[0]+mx[0])/2, (mn[1]+mx[1])/2, (mn[2]+mx[2])/2 };
      double r = 0.0;
      for( int i = lo; i <= hi; ++i ) {
        const int e = perm[i];
        const double dx = eC[e]-c[0], dy = eC[e+nE]-c[1], dz = eC[e+2*nE]-c[2];
        const double d  = std::sqrt( dx*dx + dy*dy + dz*dz ) + eR[e];
        if( d > r ) r = d;
      }
      const float cf[3] = { (float)c[0], (float)c[1], (float)c[2] };
      const double ddx = c[0]-(double)cf[0], ddy = c[1]-(double)cf[1], ddz = c[2]-(double)cf[2];
      const float rf = fup( r + std::sqrt( ddx*ddx + ddy*ddy + ddz*ddz ) );
      float* B = &b4[ (size_t)col*S ];
      B[ 0 + k ] = cf[0];  B[ 4 + k ] = cf[1];  B[ 8 + k ] = cf[2];  B[ 12 + k ] = rf;
    } else {                                 /* aabb: min/max of element boxes */
      double mn[3] = { INF, INF, INF }, mx[3] = { -INF, -INF, -INF };
      for( int i = lo; i <= hi; ++i ) {
        const int e = perm[i];
        for( int a = 0; a < 3; ++a ) {
          const double l = eLo[ e + a*nE ], h = eHi[ e + a*nE ];
          if( l < mn[a] ) mn[a] = l;
          if( h > mx[a] ) mx[a] = h;
        }
      }
      float* B = &b4[ (size_t)col*S ];
      B[ 0 + k ] = fdown( mn[0] );  B[ 4 + k ] = fup( mx[0] );
      B[ 8 + k ] = fdown( mn[1] );  B[ 12 + k ] = fup( mx[1] );
      B[ 16 + k ] = fdown( mn[2] ); B[ 20 + k ] = fup( mx[2] );
    }
  }

  void emptySlot( int col, int k )
  {
    const float FM = std::numeric_limits<float>::max();
    float* B = &b4[ (size_t)col*S ];
    if( vol == 1 ) {
      B[ 0 + k ] = B[ 4 + k ] = B[ 8 + k ] = FM;
      B[ 12 + k ] = -FM;                                    /* r < 0: never hit */
    } else if( vol == 3 || vol == 5 ) {
      for( int a = 0; a < 3; ++a ) {                        /* identity axes    */
        for( int c = 0; c < 3; ++c ) B[ 4*(3*a+c) + k ] = ( a == c ) ? 1.f : 0.f;
        B[ 4*( 9+a) + k ] =  FM;
        B[ 4*(12+a) + k ] = -FM;
      }
      if( vol == 5 ) { B[ 4*13 + k ] = 0.f;  B[ 4*14 + k ] = -FM; }
    } else if( vol == 6 ) {
      for( int a = 0; a < 6; ++a ) B[ 4*a + k ] = FM;
      B[ 4*6 + k ] = -FM;                                   /* r < 0            */
    } else {                                                /* aabb / kdop      */
      const int npair = ( vol == 4 ) ? 7 : 3;
      for( int a = 0; a < npair; ++a ) {
        B[ 4*(2*a  ) + k ] =  FM;
        B[ 4*(2*a+1) + k ] = -FM;
      }
    }
    c4[ (size_t)col*4 + k ] = 0;
    r4[ (size_t)col*8 + 2*k ] = 0;  r4[ (size_t)col*8 + 2*k + 1 ] = 0;
  }

  /* PreTri4 pool: every all-triangle leaf slot preprocessed as (v0,e1,e2) in
   * 4-wide SoA blocks of 36 doubles, padded with null lanes (id 0). Walk order
   * is deterministic, so a refit reproduces the exact same layout. */
  void buildPreTri4( std::vector<double>& pk4, std::vector<int32_t>& pk4id,
                     std::vector<int32_t>& s4 ) const
  {
    const size_t nN4 = c4.size()/4;
    s4.assign( 8*nN4, 0 );
    for( size_t col = 0; col < nN4; ++col )
      for( int k = 0; k < 4; ++k ) {
        if( c4[ col*4 + k ] != -1 ) continue;            /* leaf slots only */
        const int lo = (int)r4[ col*8 + 2*k ] - 1, hi = (int)r4[ col*8 + 2*k + 1 ] - 1;
        bool allTri = true;
        for( int i = lo; i <= hi && allTri; ++i )
          if( (int)eN[ perm[i] ] != 3 ) allTri = false;
        if( !allTri ) continue;
        const int nb = ( hi - lo + 1 + 3 )/4;
        s4[ col*8 + 2*k     ] = (int32_t)( pk4.size()/36 ) + 1;
        s4[ col*8 + 2*k + 1 ] = (int32_t)nb;
        for( int b = 0; b < nb; ++b ) {
          const size_t base = pk4.size();
          pk4.resize( base + 36, 0.0 );
          for( int l = 0; l < 4; ++l ) {
            const int i = lo + b*4 + l;
            if( i > hi ) { pk4id.push_back( 0 ); continue; }
            const int e = perm[i];
            const double* q = eV + (size_t)e*12;
            double* B = &pk4[ base ];
            B[ 0*4+l ] = q[0];        B[ 1*4+l ] = q[1];        B[ 2*4+l ] = q[2];
            B[ 3*4+l ] = q[3]-q[0];   B[ 4*4+l ] = q[4]-q[1];   B[ 5*4+l ] = q[5]-q[2];
            B[ 6*4+l ] = q[6]-q[0];   B[ 7*4+l ] = q[7]-q[1];   B[ 8*4+l ] = q[8]-q[2];
            pk4id.push_back( (int32_t)( i + 1 ) );
          }
        }
      }
  }

  int emit4( int b )     /* b: INTERNAL binary node -> BVH4 node index */
  {
    int slots[4];  int ns = 2;
    slots[0] = bn[b].l;  slots[1] = bn[b].r;
    while( ns < 4 ) {
      int best = -1, bc = -1;
      for( int k = 0; k < ns; ++k ) {
        const int s = slots[k];
        if( bn[s].l < 0 ) continue;                        /* leaf: not expandable */
        const int cnt = bn[s].hi - bn[s].lo + 1;
        if( cnt > bc ) { bc = cnt; best = k; }
      }
      if( best < 0 ) break;
      const int s = slots[best];
      slots[best]  = bn[s].l;
      slots[ns++]  = bn[s].r;
    }

    const int col = (int)( c4.size()/4 );
    b4.resize( b4.size() + S,  0.f );
    c4.resize( c4.size() + 4,  0 );
    r4.resize( r4.size() + 8,  0 );

    for( int k = 0; k < 4; ++k ) {
      if( k >= ns ) { emptySlot( col, k ); continue; }
      const int s = slots[k];
      r4[ (size_t)col*8 + 2*k     ] = bn[s].lo + 1;
      r4[ (size_t)col*8 + 2*k + 1 ] = bn[s].hi + 1;
      volumize( col, k, bn[s].lo, bn[s].hi );
      c4[ (size_t)col*4 + k ] = ( bn[s].l < 0 ) ? -1 : 0;  /* internal: fixed below */
    }
    for( int k = 0; k < ns; ++k ) {                        /* recurse AFTER reserving col */
      const int s = slots[k];
      if( bn[s].l >= 0 ) {
        /* NOTE: emit4 grows c4 (reallocation) -- evaluate it BEFORE indexing
         * c4, never inside the assignment (pre-C++17 evaluation order UB). */
        const int cid = emit4( s );
        c4[ (size_t)col*4 + k ] = cid + 1;
      }
    }
    return col;
  }
};

/* ------------------------------------------------------------------ gateway */
void mexFunction( int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[] )
{
  /* REFIT mode (5th arg is the int32 srange):
   *   bounds4 = BVH_mx( eC , eR , eLo , eHi , srange , child4 , perm , vol )
   * keeps the hierarchy, recomputes every slot volume from its element range
   * (exact for both volume types), conservative floats out. */
  if( nrhs >= 8 && mxIsInt32( prhs[4] ) ) {
    for( int i = 0; i < 4; ++i )
      if( !mxIsDouble(prhs[i]) || mxIsComplex(prhs[i]) || mxIsSparse(prhs[i]) )
        mexErrMsgIdAndTxt( "BVH_mx:type", "element arrays must be full real doubles." );
    if( !mxIsInt32(prhs[5]) || !mxIsInt32(prhs[6]) )
      mexErrMsgIdAndTxt( "BVH_mx:refit", "srange/child4/perm must be int32." );

    Builder W;
    W.nE  = mxGetM( prhs[0] );
    W.eC  = mxGetPr( prhs[0] );  W.eR  = mxGetPr( prhs[1] );
    W.eLo = mxGetPr( prhs[2] );  W.eHi = mxGetPr( prhs[3] );
    W.vol = (int)mxGetScalar( prhs[7] );
    if( W.vol < 1 || W.vol > 6 )
      mexErrMsgIdAndTxt( "BVH_mx:vol", "vol must be 1..6." );
    W.S = ( W.vol == 1 ) ? 16 : ( W.vol == 2 ) ? 24 : ( W.vol == 3 ) ? 60 : ( W.vol == 4 ) ? 56 : ( W.vol == 5 ) ? 60 : 28;
    if( nrhs < 10 || !mxIsDouble(prhs[8]) || !mxIsDouble(prhs[9]) ||
        mxGetM(prhs[8]) != 12 || mxGetN(prhs[8]) != W.nE ||
        mxGetNumberOfElements(prhs[9]) != W.nE )
      mexErrMsgIdAndTxt( "BVH_mx:refit", "refit requires eV (12 x nE) and eN (nE)." );
    W.eV = mxGetPr( prhs[8] );
    W.eN = mxGetPr( prhs[9] );

    const mwSize nN = mxGetN( prhs[4] );
    if( mxGetM(prhs[4]) != 8 || mxGetM(prhs[5]) != 4 || mxGetN(prhs[5]) != nN ||
        mxGetNumberOfElements(prhs[6]) != W.nE ||
        mxGetN(prhs[0]) != 3 || mxGetNumberOfElements(prhs[1]) != W.nE ||
        mxGetM(prhs[2]) != W.nE || mxGetM(prhs[3]) != W.nE )
      mexErrMsgIdAndTxt( "BVH_mx:refit", "inconsistent refit array sizes." );

    const int32_t* sr = (const int32_t*)mxGetData( prhs[4] );
    const int32_t* ch = (const int32_t*)mxGetData( prhs[5] );
    const int32_t* pm = (const int32_t*)mxGetData( prhs[6] );
    W.perm.resize( W.nE );
    for( mwSize i = 0; i < W.nE; ++i ) {
      const int e = (int)pm[i] - 1;
      if( e < 0 || e >= (int)W.nE )
        mexErrMsgIdAndTxt( "BVH_mx:refit", "corrupt perm." );
      W.perm[i] = e;
    }
    W.b4.assign( (size_t)W.S*nN, 0.f );
    W.c4.assign( (size_t)4*nN, 0 );     /* scratch for emptySlot */
    W.r4.assign( (size_t)8*nN, 0 );
    for( mwSize c = 0; c < nN; ++c )
      for( int k = 0; k < 4; ++k ) {
        if( ch[ c*4 + k ] == 0 ) { W.emptySlot( (int)c, k ); continue; }
        const int lo = (int)sr[ c*8 + 2*k ] - 1, hi = (int)sr[ c*8 + 2*k + 1 ] - 1;
        if( lo < 0 || hi < lo || hi >= (int)W.nE )
          mexErrMsgIdAndTxt( "BVH_mx:refit", "corrupt slot range." );
        W.volumize( (int)c, k, lo, hi );
      }
    plhs[0] = mxCreateNumericMatrix( W.S, nN, mxSINGLE_CLASS, mxREAL );
    memcpy( mxGetData(plhs[0]), W.b4.data(), sizeof(float)*W.b4.size() );
    if( nlhs > 1 ) {           /* refreshed PreTri4 pool (same layout, new geometry) */
      W.c4.assign( ch, ch + 4*nN );          /* real tree (scratch was zeroed) */
      W.r4.assign( sr, sr + 8*nN );
      std::vector<double> PK4;  std::vector<int32_t> PKID, S4;
      W.buildPreTri4( PK4, PKID, S4 );
      plhs[1] = mxCreateDoubleMatrix( 36, (mwSize)( PK4.size()/36 ), mxREAL );
      memcpy( mxGetPr(plhs[1]), PK4.data(), sizeof(double)*PK4.size() );
    }
    return;
  }

  if( nrhs < 7 )
    mexErrMsgIdAndTxt( "BVH_mx:nrhs",
                       "expected eC, eR, eLo, eHi, minLeaf, maxLeaf, vol." );
  for( int i = 0; i < 4; ++i )
    if( !mxIsDouble(prhs[i]) || mxIsComplex(prhs[i]) || mxIsSparse(prhs[i]) )
      mexErrMsgIdAndTxt( "BVH_mx:type", "inputs must be full real doubles." );

  Builder W;
  W.nE = mxGetM( prhs[0] );
  if( W.nE == 0 ) mexErrMsgIdAndTxt( "BVH_mx:empty", "no elements." );
  if( mxGetN(prhs[0]) != 3 || mxGetNumberOfElements(prhs[1]) != W.nE ||
      mxGetM(prhs[2]) != W.nE || mxGetN(prhs[2]) != 3 ||
      mxGetM(prhs[3]) != W.nE || mxGetN(prhs[3]) != 3 )
    mexErrMsgIdAndTxt( "BVH_mx:size", "inconsistent element array sizes." );

  W.eC  = mxGetPr( prhs[0] );  W.eR  = mxGetPr( prhs[1] );
  W.eLo = mxGetPr( prhs[2] );  W.eHi = mxGetPr( prhs[3] );
  W.minLeaf = (int)mxGetScalar( prhs[4] );  if( W.minLeaf < 1 ) W.minLeaf = 1;
  W.maxLeaf = (int)mxGetScalar( prhs[5] );  if( W.maxLeaf < W.minLeaf ) W.maxLeaf = W.minLeaf;
  W.vol     = (int)mxGetScalar( prhs[6] );
  if( W.vol < 1 || W.vol > 6 )
    mexErrMsgIdAndTxt( "BVH_mx:vol", "vol must be 1..6 (sphere, aabb, obb, kdop, rss, lss)." );
  W.S = ( W.vol == 1 ) ? 16 : ( W.vol == 2 ) ? 24 : ( W.vol == 3 ) ? 60 : ( W.vol == 4 ) ? 56 : ( W.vol == 5 ) ? 60 : 28;
  if( nrhs < 9 || !mxIsDouble(prhs[7]) || !mxIsDouble(prhs[8]) ||
      mxGetM(prhs[7]) != 12 || mxGetN(prhs[7]) != W.nE ||
      mxGetNumberOfElements(prhs[8]) != W.nE )
    mexErrMsgIdAndTxt( "BVH_mx:eV", "eV (12 x nE) and eN (nE) are required." );
  W.eV = mxGetPr( prhs[7] );
  W.eN = mxGetPr( prhs[8] );

  W.perm.resize( W.nE );
  for( mwSize i = 0; i < W.nE; ++i ) W.perm[i] = (int)i;
  W.bn.reserve( 2*W.nE );

  const int root = W.buildBin( 0, (int)W.nE - 1 );

  if( W.bn[root].l < 0 ) {   /* whole mesh fits in one leaf: single node4 */
    W.b4.resize( W.S, 0.f );  W.c4.resize( 4, 0 );  W.r4.resize( 8, 0 );
    W.r4[0] = 1;  W.r4[1] = (int32_t)W.nE;
    W.volumize( 0, 0, 0, (int)W.nE - 1 );
    W.c4[0] = -1;
    for( int k = 1; k < 4; ++k ) W.emptySlot( 0, k );
  } else {
    W.emit4( root );
  }

  const mwSize nN4 = (mwSize)( W.c4.size()/4 );

  plhs[0] = mxCreateNumericMatrix( W.S, nN4, mxSINGLE_CLASS, mxREAL );
  memcpy( mxGetData(plhs[0]), W.b4.data(), sizeof(float)*W.b4.size() );
  if( nlhs > 1 ) {
    plhs[1] = mxCreateNumericMatrix( 4, nN4, mxINT32_CLASS, mxREAL );
    memcpy( mxGetData(plhs[1]), W.c4.data(), sizeof(int32_t)*W.c4.size() );
  }
  if( nlhs > 2 ) {
    plhs[2] = mxCreateNumericMatrix( 8, nN4, mxINT32_CLASS, mxREAL );
    memcpy( mxGetData(plhs[2]), W.r4.data(), sizeof(int32_t)*W.r4.size() );
  }
  if( nlhs > 3 ) {
    plhs[3] = mxCreateNumericMatrix( W.nE, 1, mxINT32_CLASS, mxREAL );
    int32_t* p = (int32_t*)mxGetData( plhs[3] );
    for( mwSize i = 0; i < W.nE; ++i ) p[i] = (int32_t)( W.perm[i] + 1 );
  }
  if( nlhs > 4 ) {
    std::vector<double> PK4;  std::vector<int32_t> PKID, S4;
    W.buildPreTri4( PK4, PKID, S4 );
    const mwSize nB = (mwSize)( PK4.size()/36 );
    plhs[4] = mxCreateDoubleMatrix( 36, nB, mxREAL );
    memcpy( mxGetPr(plhs[4]), PK4.data(), sizeof(double)*PK4.size() );
    if( nlhs > 5 ) {
      plhs[5] = mxCreateNumericMatrix( 4, nB, mxINT32_CLASS, mxREAL );
      memcpy( mxGetData(plhs[5]), PKID.data(), sizeof(int32_t)*PKID.size() );
    }
    if( nlhs > 6 ) {
      plhs[6] = mxCreateNumericMatrix( 8, nN4, mxINT32_CLASS, mxREAL );
      memcpy( mxGetData(plhs[6]), S4.data(), sizeof(int32_t)*S4.size() );
    }
  }
}
