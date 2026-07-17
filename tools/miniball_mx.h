#include <vector>
// #include <cassert>
// #include <iomanip>
// #include <iostream>
// #include <cmath>
#include <numeric>
// #include <ostream>
// #include <algorithm>
#if defined( __GNUC__ ) && __GNUC__==2 && \
  __GNUC_MINOR__==95 && __GNUC_PATCHLEVEL__ <= 2
  #include <gcc2-95-2_fix.h>
// #else
//   #include <sstream>
#endif

namespace Seb { //'Seb' stands for "smalles enclosing ball"
  template<typename Float> class Point {
    public:
      typedef typename std::vector<Float>::const_iterator Const_iterator;
      typedef typename std::vector<Float>::iterator Iterator;
      Point( int d ): c( d ){}
      template<typename InputIterator> Point( int d,InputIterator first ): c( first,first+d ){}
      const Float& operator[]( unsigned int i ) const {
        return( c[i] );
      }
      Float& operator[]( unsigned int i ){
        return( c[i] );
      }
      Const_iterator begin() const {
        return( c.begin() );
      }
      Const_iterator end() const {
        return( c.end() );
      }
    private:
      std::vector<Float> c;
  };
}
namespace Seb {
  template<typename Float> inline Float sqr( const Float x ){
    return( x * x );
  }
  template<typename Float, class Pt, class PointAccessor> class Subspan {
    public:
      Subspan( unsigned int dim, const PointAccessor& S, int i );
      ~Subspan();
      void add_point( int global_index );
      void remove_point( unsigned int local_index );
      unsigned int size() const {
        return( r+1 );
      }
      bool is_member( unsigned int i ) const {
        return( membership[i] );
      }
      unsigned int global_index( unsigned int i ) const {
        return( members[i] );
      }
      unsigned int any_member() const {
        return( members[r] );
      }
      template<typename RandomAccessIterator1, typename RandomAccessIterator2> Float shortest_vector_to_span( RandomAccessIterator1 p, RandomAccessIterator2 w );
      template<typename RandomAccessIterator1,typename RandomAccessIterator2> void find_affine_coefficients( RandomAccessIterator1 c, RandomAccessIterator2 coeffs );
      Float representation_error();
    private:
      void append_column();
      void hessenberg_clear( unsigned int start );
      void special_rank_1_update();
      const PointAccessor &S;
      std::vector<bool> membership;
      const unsigned int dim;
      std::vector<unsigned int> members;
      Float **Q, **R;
      Float *u,*w;
      unsigned int r;
  };
}
namespace Seb {
  template<typename Float> inline void givens( Float& c, Float& s, const Float a, const Float b ){
    using std::abs;
    using std::sqrt;
    if( b == 0 ){
      c = 1;
      s = 0;
    } else if( abs( b ) > abs( a ) ){
      const Float t = a / b;
      s = 1 / sqrt( 1 + sqr( t ) );
      c = s * t;
    } else {
      const Float t = b / a;
      c = 1 / sqrt( 1 + sqr( t ) );
      s = c * t;
    }
  }
  template<typename Float, class Pt, class PointAccessor> Subspan<Float, Pt, PointAccessor>::Subspan( unsigned int dim, const PointAccessor& S, int index ): S( S ), membership( S.size() ), dim( dim ), members( dim+1 ){
    Q = new Float *[dim];
    R = new Float *[dim];
    for( unsigned int i=0; i<dim; ++i ){
      Q[i] = new Float[dim];
      R[i] = new Float[dim];
    }
    u = new Float[dim];
    w = new Float[dim];
    for( unsigned int i=0; i<dim; ++i ){
      for( unsigned int j=0; j<dim; ++j ){
        Q[i][j] =( i==j )? 1: 0;
      }
    }
    members[r = 0] = index;
    membership[index] = true;
  }
  template<typename Float, class Pt, class PointAccessor> Subspan<Float, Pt, PointAccessor>::~Subspan(){
    for( unsigned int i=0; i<dim; ++i ){
      delete[] Q[i];
      delete[] R[i];
    }
    delete[] Q;
    delete[] R;
    delete[] u;
    delete[] w;
  }
  template<typename Float, class Pt, class PointAccessor> void Subspan<Float, Pt, PointAccessor>::add_point( int index ){
    for( unsigned int i=0; i<dim; ++i ){
      u[i] = S[index][i] - S[members[r]][i];
    }
    append_column();
    membership[index] = true;
    members[r+1] = members[r];
    members[r]   = index;
    ++r;
  }
  template<typename Float, class Pt, class PointAccessor> void Subspan<Float, Pt, PointAccessor>::remove_point( const unsigned int local_index ){
    membership[global_index( local_index )] = false;
    if( local_index == r ){
      for( unsigned int i=0; i<dim; ++i ){
        u[i] = S[members[r]][i] - S[global_index( r-1 )][i];
      }
      --r;
      special_rank_1_update();
    } else {
      Float *dummy = R[local_index];
      for( unsigned int j = local_index+1; j < r; ++j ){
        R[j-1] = R[j];
        members[j-1] = members[j];
      }
      members[r-1] = members[r];
      R[--r] = dummy;
      hessenberg_clear( local_index );
    }
  }
  template<typename Float, class Pt, class PointAccessor>
  template<typename RandomAccessIterator1,
  typename RandomAccessIterator2>
  Float Subspan<Float, Pt, PointAccessor>::
  shortest_vector_to_span( RandomAccessIterator1 p, RandomAccessIterator2 w ){
    using std::inner_product;
    for( unsigned int i=0; i<dim; ++i ){
      w[i] = S[members[r]][i] - p[i];
    }
    for( unsigned int j = 0; j < r; ++j ){
      const Float scale = inner_product( w,w+dim,Q[j],Float( 0 ) );
      for( unsigned int i = 0; i < dim; ++i ){
        w[i] -= scale * Q[j][i];
      }
    }
    return( inner_product( w,w+dim,w,Float( 0 ) ) );
  }
  template<typename Float, class Pt, class PointAccessor> Float Subspan<Float, Pt, PointAccessor>::representation_error(){
    using std::abs;
    std::vector<Float> lambdas( size() );
    Float max = 0;
    Float error;
    for( unsigned int j = 0; j < size(); ++j ){
      find_affine_coefficients( S[global_index( j )],lambdas.begin() );
      error = abs( lambdas[j] - 1.0 );
      if( error > max ){ max = error; }
      for( unsigned int i = 0; i < j; ++i ){
        error = abs( lambdas[i] - 0.0 );
        if( error > max ){ max = error; }
      }
      for( unsigned int i = j+1; i < size(); ++i ){
        error = abs( lambdas[i] - 0.0 );
        if( error > max ){ max = error; }
      }
    }
    return( max );
  }
  template<typename Float, class Pt, class PointAccessor>
  template<typename RandomAccessIterator1,
  typename RandomAccessIterator2>
  void Subspan<Float, Pt, PointAccessor>::
  find_affine_coefficients( RandomAccessIterator1 p,
                           RandomAccessIterator2 lambdas ){
    for( unsigned int i=0; i<dim; ++i ){
      u[i] = p[i] - S[members[r]][i];
    }
    for( unsigned int i = 0; i < dim; ++i ){
      w[i] = 0;
      for( unsigned int k = 0; k < dim; ++k ){
        w[i] += Q[i][k] * u[k];
      }
    }
    Float origin_lambda = 1;
    for( int j = r-1; j>=0; --j ){
      for( unsigned int k=j+1; k<r; ++k ){
        w[j] -= *( lambdas+k ) * R[k][j];
      }
      origin_lambda -= *( lambdas+j ) = w[j] / R[j][j];
    }
    *( lambdas+r ) = origin_lambda;
  }
  template<typename Float, class Pt, class PointAccessor>
  void Subspan<Float, Pt, PointAccessor>::append_column(){
    for( unsigned int i = 0; i < dim; ++i ){
      R[r][i] = 0;
      for( unsigned int k = 0; k < dim; ++k ){
        R[r][i] += Q[i][k] * u[k];
      }
    }
    for( unsigned int j = dim-1; j > r; --j ){
      Float c, s;
      givens( c,s,R[r][j-1],R[r][j] );
      R[r][j-1] = c * R[r][j-1] + s * R[r][j];
      for( unsigned int i = 0; i < dim; ++i ){
        const Float a = Q[j-1][i];
        const Float b = Q[j][i];
        Q[j-1][i] =  c * a + s * b;
        Q[j][i]   =  c * b - s * a;
      }
    }
  }
  template<typename Float, class Pt, class PointAccessor>
  void Subspan<Float, Pt, PointAccessor>::hessenberg_clear( unsigned int pos ){
    for( ; pos < r; ++pos ){
      Float c, s;
      givens( c,s,R[pos][pos],R[pos][pos+1] );
      R[pos][pos] = c * R[pos][pos] + s * R[pos][pos+1];
      for( unsigned int j = pos+1; j < r; ++j ){
        const Float a = R[j][pos];
        const Float b = R[j][pos+1];
        R[j][pos]   =  c * a + s * b;
        R[j][pos+1] =  c * b - s * a;
      }
      for( unsigned int i = 0; i < dim; ++i ){
        const Float a = Q[pos][i];
        const Float b = Q[pos+1][i];
        Q[pos][i]   =  c * a + s * b;
        Q[pos+1][i] =  c * b - s * a;
      }
    }
  }
  template<typename Float, class Pt, class PointAccessor>
  void Subspan<Float, Pt, PointAccessor>::special_rank_1_update(){
    for( unsigned int i = 0; i < dim; ++i ){
      w[i] = 0;
      for( unsigned int k = 0; k < dim; ++k ){
        w[i] += Q[i][k] * u[k];
      }
    }
    for( unsigned int k = dim-1; k > 0; --k ){
      Float c, s;
      givens( c,s,w[k-1],w[k] );
      w[k-1] = c * w[k-1] + s * w[k];
      R[k-1][k]    = -s * R[k-1][k-1];
      R[k-1][k-1] *=  c;
      for( unsigned int j = k; j < r; ++j ){
        const Float a = R[j][k-1];
        const Float b = R[j][k];
        R[j][k-1] =  c * a + s * b;
        R[j][k]   =  c * b - s * a;
      }
      for( unsigned int i = 0; i < dim; ++i ){
        const Float a = Q[k-1][i];
        const Float b = Q[k][i];
        Q[k-1][i] =  c * a + s * b;
        Q[k][i]   =  c * b - s * a;
      }
    }
    for( unsigned int j = 0; j < r; ++j ){
      R[j][0] += w[0];
    }
    hessenberg_clear( 0 );
  }
}
namespace Seb {
  template<typename Float, class Pt = Point<Float>, class PointAccessor = std::vector<Pt> >
  class Smallest_enclosing_ball {
  public:
    typedef Float *Coordinate_iterator;
  public:
    Smallest_enclosing_ball( unsigned int d, const PointAccessor &P ): dim( d ), S( P ), up_to_date( true ), support( NULL ){
      allocate_resources();
      update();
    }
    ~Smallest_enclosing_ball(){
      deallocate_resources();
    }
  public:
    void invalidate(){
      up_to_date = false;
    }
  public:
    bool is_empty(){
      return( S.size() == 0 );
    }
    Float squared_radius(){
      if( !up_to_date ){ update(); }
      return( radius_square );
    }
    Float radius(){
      if( !up_to_date ){ update(); }
      return( radius_ );
    }
    Coordinate_iterator center_begin(){
      if( !up_to_date ){ update(); }
      return( center );
    }
    Coordinate_iterator center_end(){
      if( !up_to_date ){ update(); }
      return center+dim;
    }
  private:
    void allocate_resources();
    void deallocate_resources();
  private:
    void init_ball();
    Float find_stop_fraction( int& hinderer );
    bool successful_drop();
    void update();
  private:
    Smallest_enclosing_ball( const Smallest_enclosing_ball& );
    Smallest_enclosing_ball& operator=( const Smallest_enclosing_ball& );
  private:
    unsigned int dim;
    const PointAccessor &S;
    bool up_to_date;
    Float *center;
    Float radius_, radius_square;
    Subspan<Float, Pt, PointAccessor> *support;
  private:
    Float *center_to_aff;
    Float *center_to_point;
    Float *lambdas;
    Float  dist_to_aff, dist_to_aff_square;
  private:
    static const Float Eps;
  };
}
namespace Seb {
  template<typename Float, class Pt, class PointAccessor> void Smallest_enclosing_ball<Float, Pt, PointAccessor>::allocate_resources(){
    center            = new Float[dim];
    center_to_aff     = new Float[dim];
    center_to_point   = new Float[dim];
    lambdas           = new Float[dim+1];
  }
  template<typename Float, class Pt, class PointAccessor> void Smallest_enclosing_ball<Float, Pt, PointAccessor>::deallocate_resources(){
    delete[] center;
    delete[] center_to_aff;
    delete[] center_to_point;
    delete[] lambdas;
    if( support != NULL ){
      delete support;
    }
  }
  template<typename Float, class Pt, class PointAccessor>
  void Smallest_enclosing_ball<Float, Pt, PointAccessor>::init_ball(){
    for( unsigned int i = 0; i < dim; ++i ){
      center[i] = S[0][i];
    }
    radius_square = 0;
    unsigned int farthest = 0;
    for( unsigned int j = 1; j < S.size(); ++j ){
      Float dist = 0;
      for( unsigned int i = 0; i < dim; ++i ){
        dist += sqr( S[j][i] - center[i] );
      }
      if( dist >= radius_square ){
        radius_square = dist;
        farthest = j;
      }
      radius_ = sqrt( radius_square );
    }
    if( support != NULL ){
      delete support;
    }
    support = new Subspan<Float, Pt, PointAccessor>( dim,S,farthest );
  }
  template<typename Float, class Pt, class PointAccessor>
  bool Smallest_enclosing_ball<Float, Pt, PointAccessor>::successful_drop(){
    support->find_affine_coefficients( center,lambdas );
    unsigned int smallest = 0;
    Float minimum( 1 );
    for( unsigned int i=0; i<support->size(); ++i ){
      if( lambdas[i] < minimum ){
        minimum = lambdas[i];
        smallest = i;
      }
    }
    if( minimum <= 0 ){
      support->remove_point( smallest );
      return true;
    }
    return false;
  }
  template<typename Float, class Pt, class PointAccessor>
  Float Smallest_enclosing_ball<Float, Pt, PointAccessor>::find_stop_fraction( int& stopper ){
    using std::inner_product;
    Float scale =  1;
    stopper     = -1;
    for( unsigned int j = 0; j < S.size(); ++j ){
      if( !support->is_member( j ) ){
        for( unsigned int i = 0; i < dim; ++i ){
          center_to_point[i] = S[j][i] - center[i];
        }
        const Float dir_point_prod
        = inner_product( center_to_aff,center_to_aff+dim,
                        center_to_point,Float( 0 ) );
        if( dist_to_aff_square - dir_point_prod < Eps * radius_ * dist_to_aff ){
          continue;
        }
        Float bound = radius_square;
        bound -= inner_product( center_to_point,center_to_point+dim,
                               center_to_point,Float( 0 ) );
        bound /= 2 *( dist_to_aff_square - dir_point_prod );
        if( bound > 0 && bound < scale ){
          scale   = bound;
          stopper = j;
        }
      }
    }
    return scale;
  }
  template<typename Float, class Pt, class PointAccessor>
  void Smallest_enclosing_ball<Float, Pt, PointAccessor>::update(){
    up_to_date = true;
    init_ball();
    while( true ){
      while( ( dist_to_aff
              = sqrt( dist_to_aff_square
                     = support->shortest_vector_to_span( center,
                                                        center_to_aff ) ) )
             <= Eps * radius_ )
        if( !successful_drop() ){
          return;
        }
      int stopper;
      Float scale = find_stop_fraction( stopper );
      if( stopper >= 0 && support->size() <= dim ){
        for( unsigned int i = 0; i < dim; ++i ){
          center[i] += scale * center_to_aff[i];
        }
        const Pt& stop_point = S[support->any_member()];
        radius_square = 0;
        for( unsigned int i = 0; i < dim; ++i ){
          radius_square += sqr( stop_point[i] - center[i] );
        }
        radius_ = sqrt( radius_square );
        support->add_point( stopper );
      } else {
        for( unsigned int i=0; i<dim; ++i ){
          center[i] += center_to_aff[i];
        }
        const Pt& stop_point = S[support->any_member()];
        radius_square = 0;
        for( unsigned int i = 0; i < dim; ++i ){
          radius_square += sqr( stop_point[i] - center[i] );
        }
        radius_ = sqrt( radius_square );
        if( !successful_drop() ){
          return;
        }
      }
    }
  }
  template<typename Float, class Pt, class PointAccessor>
  const Float Smallest_enclosing_ball<Float, Pt, PointAccessor>::Eps = Float( 1e-14 );
}
