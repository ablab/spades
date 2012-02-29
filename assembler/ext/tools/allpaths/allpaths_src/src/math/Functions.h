///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef FUNCTIONS
#define FUNCTIONS

#include "Vec.h"
#include "math/Arith.h"
#include <cmath>
#include <numeric>

// ===========================================================================
//
// Min functions
//
// ==========================================================================

template <class T> inline T Min(const T& t1,const T& t2)                          
{ return (t1<t2)?t1:t2;             }

template <class T> inline T Min(const T& t1,const T& t2,const T& t3)              
{ return Min(t1,Min(t2,t3));        }

template <class T> inline T Min(const T& t1,const T& t2,const T& t3, const T& t4) 
{ return Min(Min(t1,t2),Min(t3,t4));}

template<class T> inline T Min( const vec<T>& v )
{    if ( v.size( ) == 0 ) FatalErr( "Minimum of zero length vector requested." );
     T answer = v[0];
     for ( unsigned int i = 1; i < v.size( ); i++ )
          answer = Min( answer, v[i] );
     return answer;    }

/// Compute a vector of length last-first which represents the minimum,
/// as computed by op, in a window of size nbhd on each side, of the
/// entries of the sequence [first, last).  Window entries that extend
/// outside the bounds of the input sequence are ignored.  For
/// instance, MinWindow of the 8 integers 01234567 with nbhd 3 is
/// 00001234.
///
/// The op should be similar to std::less in that op(a,b)==op(b,a)
/// implies a==b.
template<class FwdIt, class OutIt, class Op>
void MinWindow(FwdIt first, FwdIt last, OutIt out, int nbhd, Op op)
{
  // We consider a window [first, winLast) "centered" on center.  
  FwdIt winLast = first, center=first;
  // w is the distance from first to center.
  typename FwdIt::difference_type w=0;
  // Initialize the window: initially of size nbhd+1 but of size up to
  // 2*nbhd+1 when possible.  *smallest is the smallest value in the
  // window, according to operator<, and the leftmost such if there is
  // more than one.
  FwdIt smallest = first, nextSmallest = first;
  for (int i=0; i<=nbhd && winLast<last; ++winLast, ++i) {
    if (op(*winLast,*smallest)) {
      smallest = winLast;
    }
  }
  while (center!=last) {
    // Output an entry and move on to next
    *out = *smallest;
    ++out;
    ++center;
    ++w;
    // Now move the window along, moving first and winLast as necessary
    // to keep the window size as close to 2*nbhd+1 as possible while
    // keeping both in bounds.
    if (winLast<last) {
      if (op(*winLast,*smallest)) {
	smallest = winLast;
      }
      ++winLast;
    }
    if (w>nbhd) {
      --w;
      if (smallest==first++) {
	smallest = min_element(first, winLast, op);
      }
    }
  }
}

/// Compute a vector of length last-first which represents the minimum
/// in a window of size nbhd (on each side) of the entries of the
/// sequence [first, last).  Window entries that extend outside the
/// bounds of the input sequence are ignored.  For instance, MinWindow
/// of the 8 integers 01234567 with nbhd 3 is 00001234.
template<class FwdIt, class OutIt>
inline void MinWindow(FwdIt first, FwdIt last, OutIt out, int nbhd)
{
  MinWindow(first, last, out, nbhd,
	    std::less<typename std::iterator_traits<FwdIt>::value_type>());
}

// ==========================================================================
//
// Max functions
//
// ==========================================================================

template <class T> inline T Max(const T& t1,const T& t2)                          
{ return (t1>t2)?t1:t2;             }

template <class T> inline T Max(const T& t1,const T& t2,const T& t3)              
{ return (t1>t2) ? ((t1>t3) ? t1 : t3) : ((t2>t3) ? t2 : t3); }

template <class T> inline T Max(const T& t1,const T& t2,const T& t3,const T& t4)  
{ return Max(Max(t1,t2),Max(t3,t4));}

template <class T> inline T Max( const T& t1, const T& t2, const T& t3, const T& t4,
     const T& t5 )
{ return Max( Max(t1,t2), Max(t3,t4), t5 ); }

template <class T> inline T Max( const T& t1, const T& t2, const T& t3, const T& t4,
     const T& t5, const T& t6 )
{ return Max( Max(t1,t2), Max(t3,t4), t5, t6 ); }

template <class T> inline T Max( const T& t1, const T& t2, const T& t3, const T& t4,
     const T& t5, const T& t6, const T& t7 )
{ return Max( Max(t1,t2), Max(t3,t4), Max(t5,t6), t7 ); }

template <class T> inline T Max( const T& t1, const T& t2, const T& t3, const T& t4,
     const T& t5, const T& t6, const T& t7, const T& t8 )
{ return Max( Max(t1,t2), Max(t3,t4), Max(t5,t6), t7, t8 ); }

template<class T> inline T Max( const vector<T>& v )
{    if ( v.size( ) == 0 ) FatalErr( "Maximum of zero length vector requested." );
     T answer = v[0];
     for ( unsigned int i = 1; i < v.size( ); i++ )
          answer = Max( answer, v[i] );
     return answer;    }

// MaxMax: return maximum entry in a vec<vec>.

template<class T> inline T MaxMax( const vec< vec<T> >& v )
{    Bool first = True;
     T answer = 0;
     for ( int i = 0; i < v.isize( ); i++ )
     {    for ( int j = 0; j < v[i].isize( ); j++ )
          {    if (first)
               {    answer = v[i][j];
                    first = False;    }
               else answer = Max( answer, v[i][j] );    }    }
     if (first) FatalErr( "MaxMax called on empty vec<vec>." );
     return answer;    }            

/// Like MinWindow but for maximum.  For instance, MaxWindow
/// of the 8 integers 01234567 with nbhd 3 is 34567777.
template<class FwdIt, class OutIt>
inline void MaxWindow(FwdIt first, FwdIt last, OutIt out, int nbhd)
{
  MinWindow(first, last, out, nbhd,
	    std::greater<typename std::iterator_traits<FwdIt>::value_type>());
}

// ===========================================================================
//
// Sum functions
//
// ===========================================================================

///Templatized sum over all the elements of a container of numeric values.
///Prerequisites: V has begin(), end(), and a value_type that is
///convertible to longlong.
template<class V>
longlong Sum (const V & v) {
  longlong ret=0;
  return accumulate(v.begin(), v.end(), ret);
}

///Templatized cumulative sum over all the elements of a container of numeric values.
///Prerequisites: V has begin() and end() methods that return valid STL iterators
template<class V>
void CumulativeSum (const V & in, V & out ) {
  out.resize( in.size() );
  partial_sum( in.begin(), in.end(), out.begin() );
}

inline int Sum( const vec<int>& v )
{    int answer = 0;
     for ( unsigned int i = 0; i < v.size( ); i++ )
          answer += v[i];
     return answer;    }

inline longlong BigSum( const vec<int>& v )
{    longlong answer = 0;
     for ( unsigned int i = 0; i < v.size( ); i++ )
          answer += v[i];
     return answer;    }

inline int Sum( const vec<unsigned int>& v )
{    int answer = 0;
     for ( unsigned int i = 0; i < v.size( ); i++ )
          answer += v[i];
     return answer;    }

inline longlong BigSum( const vec<unsigned int>& v )
{    longlong answer = 0;
     for ( unsigned int i = 0; i < v.size( ); i++ )
          answer += v[i];
     return answer;    }

inline longlong Sum( const vec<longlong>& v )
{    longlong answer = 0;
     for ( unsigned int i = 0; i < v.size( ); i++ )
          answer += v[i];
     return answer;    }

inline longlong BigSum( const vec<longlong>& v )
{    longlong answer = 0;
     for ( unsigned int i = 0; i < v.size( ); i++ )
          answer += v[i];
     return answer;    }

inline uint64_t BigSum( const vec<uint64_t>& v )
{    uint64_t answer = 0;
     for ( size_t i = 0; i < v.size( ); i++ )
          answer += v[i];
     return answer;    }

inline float Sum( const vec<float>& v )
{    float answer = 0;
     for ( unsigned int i = 0; i < v.size( ); i++ )
          answer += v[i];
     return answer;    }

inline double Sum( const vec<double>& v )
{    double answer = 0;
     for ( unsigned int i = 0; i < v.size( ); i++ )
          answer += v[i];
     return answer;    }

inline long double Sum( const vec<long double>& v )
{    long double answer = 0;
     for ( unsigned int i = 0; i < v.size( ); i++ )
          answer += v[i];
     return answer;    }

inline void CumulativeSum( const vec<double>& in, vec<double> & out )
{    double value = 0.0;
  out.resize(in.size());

  for ( unsigned int i = 0; i < in.size( ); i++ ) {
          value += in[i];
	  out[i] = value;   
  }
}

inline int Sum( const vec<Bool>& v )
{    int answer = 0;
     for ( unsigned int i = 0; i < v.size( ); i++ )
          answer += v[i];
     return answer;    }

/// Compute a vector of length last-first which represents the sums in
/// a window of size nbhd (on each side) of the entries of the sequence
/// [first, last).  Window entries that extend outside the bounds of
/// the input sequence are ignored.  For instance, SumWindow of the 8
/// integers 01231011 with nbhd 3 is 67789863.
template<class FwdIt, class OutIt>
void SumWindow(FwdIt first, FwdIt last, OutIt out, int nbhd)
{
  // We consider a window [first, winLast) "centered" on center.  
  FwdIt winLast = first, center=first;
  // w is the distance from first to center.
  typename FwdIt::difference_type w=0;
  // Initialize the window: initially of size nbhd+1 but of size up to
  // 2*nbhd+1 when possible.  s is the current sum.
  typename FwdIt::value_type s=0;
  for (int i=0; i<=nbhd && winLast<last; ++winLast, ++i) {
    s += *winLast;
  }
  while (center!=last) {
    // Output an entry and move on to next
    *out = s;
    ++out;
    ++center;
    ++w;
    // Now move the window along, moving first and winLast as
    // necessary to keep w==nbhd if possible while keeping both in
    // bounds.  Update s whenever either pointer changes.
    if (winLast<last) {
      s += *winLast;
      ++winLast;
    }
    if (w>nbhd) {
      --w;
      s -= *first;
      ++first;
    }
  }
}

// ===========================================================================
//
// Mean and standard deviation
//
// ===========================================================================

class normal_distribution 
{
public:
  float mu_;
  float sigma_;
  
  normal_distribution(const float mu = 0.0, const float sigma = 0.0)
    : mu_(mu), sigma_(sigma) {}

  inline float Zscore(const float x) const 
  { return (0 == sigma_) ? 0 : (x - mu_) / sigma_;  }


  inline float ProbabilityDensity(const float x) const 
  { if (0 == sigma_) return 0;
    const float z = Zscore(x);
    return exp(-0.5 * z * z) * sqrt(0.5 / M_PI) / sigma_;
  }

  inline float ProbabilityLT(const float x) const 
  { return 0.5 * erfc(-Zscore(x) * sqrt(0.5)); }

  inline float ProbabilityGT(const float x) const 
  { return 1.0 - ProbabilityLT(x); } 

  //    X += Y
  inline normal_distribution & operator+=(const normal_distribution & Y)
  { 
    mu_ += Y.mu_;
    sigma_ = sqrt(sigma_ * sigma_ + Y.sigma_ * Y.sigma_);
    return *this;
  }
  
  //    X -= Y
  inline normal_distribution & operator-=(const normal_distribution & Y)
  { 
    mu_ -= Y.mu_;
    sigma_ = sqrt(sigma_ * sigma_ + Y.sigma_ * Y.sigma_);
    return *this;
  }

};


//      X == Y
inline bool operator==(const normal_distribution & X, 
                       const normal_distribution & Y) 
{ return X.mu_ == Y.mu_ && X.sigma_ == Y.sigma_; }

//      X < Y
inline bool operator<(const normal_distribution & X, 
                      const normal_distribution & Y) 
{ return X.mu_ < Y.mu_; }

//      X <= Y
inline bool operator<=(const normal_distribution & X, 
                       const normal_distribution & Y) 
{ return X.mu_ <= Y.mu_; }

//      X > Y
inline bool operator>(const normal_distribution & X, 
                      const normal_distribution & Y) 
{ return X.mu_ > Y.mu_; }

//      X >= Y
inline bool operator>=(const normal_distribution & X, 
                       const normal_distribution & Y) 
{ return X.mu_ >= Y.mu_; }


//      -X
inline normal_distribution operator-(const normal_distribution & X)
{ return normal_distribution(-X.mu_, X.sigma_); }


//      X + dx
inline normal_distribution operator+(const normal_distribution & X, 
                                     const float dx) 
{ return normal_distribution(X.mu_ + dx, X.sigma_); }


//      dx + X 
inline normal_distribution operator+(const float dx, 
                                     const normal_distribution & X)
{ return normal_distribution(dx + X.mu_, X.sigma_); }


//      X - dx
inline normal_distribution operator-(const normal_distribution & X, 
                                     const float dx) 
{ return normal_distribution(X.mu_ - dx, X.sigma_); }


//      dx - X 
inline normal_distribution operator-(const float dx, 
                                     const normal_distribution & X)
{ return normal_distribution(dx - X.mu_, X.sigma_); }



inline ostream & operator<<(ostream & os, const 
                            normal_distribution & X) 
{ return os << X.mu_ << "\t" << X.sigma_ << endl; }

inline istream &operator>>( istream &in, normal_distribution &nd ) {
  in >> nd.mu_ >> nd.sigma_;
  return in;
}


/*  CombineNormalDistribution                Filipe Ribeiro 2009-06-22
 *
 *   When you have different measurements (normaly distributed) of the same quantity, 
 *   you can compute the resultant normal distribution.   
 * 
 *   The resultant normal distribution (mu, sigma) satisfies:
 * 
 *      1 / sigma^2  =  sum(     1 / sigma[i]^2 )
 * 
 *     mu / sigma^2  =  sum( mu[i] / sigma[i]^2 )
 *
 */
normal_distribution CombineNormalDistributions(const vec<normal_distribution> & ds,
                                               float * score = 0);



/* ClusterNormalDistributions             Filipe Ribeiro 2010-04-16
 *  
 *  Extract clusters based on overlaps of n_sigma x sigma.
 */
vec< vec<normal_distribution> > 
ClusterNormalDistributions(const vec<normal_distribution> & ds,
                           const float n_sigma = 1);



///Calculate mean and stdev of a vector, return {-1,-1} if vector empty.
normal_distribution SafeMeanStdev(const vec<float> & points);

///Calculate mean and stdev based on count, sum and sumsq.
///return {-1,-1} if vector empty.
normal_distribution SafeMeanStdev(int count, double sum, double sumsq);

///Return optimal cutoff between distributions to minimize assignment error.
///Return nan if we cannot find one.
///ratio is the ratio of the prior probabilities: that is, we do not assume
///that both distributions have a total probability of 1. For example,
///if we know that there are 1000 items in d1 and only 100 in d2, we will
///put the cutoff closer to d1.mu_ than if there were equal numbers.
///Precondition: d1.mu_ < d2.mu_.
double OptimalCutoff(const normal_distribution & d1,
		     const normal_distribution & d2, 
		     const double ratio = 1);


template<class V>
inline double Mean( const V & v )
{
  ForceAssert ( v.size() > 0 );
  double sum = 0.0;
  for ( int i = 0; i < int(v.size()); i++ )
    sum += v[i];
  return sum / static_cast<double>(v.size());
}

template<class V>
inline double StdDev( const V & v , const double mean )
{
  ForceAssert ( v.size() > 0 );
  double dev, sumdevsq = 0.0;
  for (int i = 0; i < int(v.size()); i++) {
    dev = v[i] - mean;
    sumdevsq += dev*dev;
  }
  return sqrt( sumdevsq/static_cast<double>(v.size()) );
}



template <class T> inline T Mean(const T& t1,const T& t2)                         
{ return (t1+t2)/2;                 }

template <class T> inline T Mean(const T& t1,const T& t2, const T& t3)            
{ return (t1+t2+t3)/3;              }

// ===========================================================================
//
// Weighted mean: return sum(v_i^2)/sum(v_i).
//
// ===========================================================================

inline Float WeightedMean( const vector<int>& v )
{    Float sum1 = 0, sum2 = 0;
     for ( unsigned int i = 0; i < v.size(); i++ )
     {    sum1 += static_cast<Float>(v[i]);
          sum2 += Float(v[i]) * Float(v[i]);    }
     ForceAssert( sum1 != 0 );
     return sum2 / sum1;    }
 

/// Compute a vector of length last-first which represents the means in
/// a window of size nbhd (on each side) of the entries of the sequence
/// [first, last).  Window entries that extend outside the bounds of
/// the input sequence are not included in the mean.  For instance,
/// MeanWindow of the 8 integers 01231011 with nbhd 3 is the sequence
/// of doubles 6/4 7/5 7/6 8/7 9/7 8/6 6/5 3/4.
template<class FwdIt, class OutIt>
void MeanWindow(FwdIt first, FwdIt last, OutIt out, int nbhd)
{
  // We consider a window [first, winLast) "centered" on center.  
  FwdIt winLast = first, center=first;
  // w is the distance from first to center.
  typename FwdIt::difference_type w=0;
  // wid is the number of elements in the window.
  typename FwdIt::difference_type wid=0;
  // Initialize the window: initially of size nbhd+1 but of size up to
  // 2*nbhd+1 when possible.  s is the current sum.
  typename FwdIt::value_type s=0;
  for ( ; wid<=nbhd && winLast<last; ++winLast, ++wid) {
    s += *winLast;
  }
  while (center!=last) {
    // Output an entry and move on to next
    *out = double(s)/wid;
    ++out;
    ++center;
    ++w;
    // Now move the window along, moving first and winLast as
    // necessary to keep w==nbhd if possible while keeping both in
    // bounds.  Update s whenever either pointer changes.
    if (winLast<last) {
      s += *winLast;
      ++winLast;
      ++wid;
    }
    if (w>nbhd) {
      --w;
      s -= *first;
      ++first;
      --wid;
    }
  }
}






///Probability of n events for a Poisson dist. characterized by lambda.
double Poisson(double lambda, unsigned int n);

///Probability of at most n events for a Poisson dist. characterized by lambda.
double PoissonCdf(double lambda, unsigned int n);

///Minimum number of events with cumulative probability greater than p
///for a Poisson dist. characterized by lambda.
int InversePoissonCdf(double lambda, double p);
     
// ===========================================================================
//
/// N50: given a list v of positive integers, suppose that each entry n is replaced
/// by n copies of itself.  Then the median of the resulting enlarged list is the
/// "N50" of the original list.
///
/// It is assumed that v is sorted, prior to calling N50.
///
// ==========================================================================

template<class T> T N50( const vec<T>& v );

// N50_size: Returns the N50 of the vector sizes in the input vec<vec<T> >.
// Calls N50, above.
// N50_size sorts internally; it does not require that the input be sorted.
template<class T> int N50_size( const vec<vec<T> >& v );

// ===========================================================================
//
/// NStatistics: given a sorted list v of positive integers, suppose
/// that each entry n is replaced by n copies of itself, resulting in m
/// total entries.  Reverse the order of this list, so that the largest
/// element comes first, and call this expanded, reversed list w..
/// Return a vector containing N10, N20, ... N90 for v, where NX is
/// the m*X/100'th element of w.
///
/// For example, if 
///  v = 1, 2, 3, 4
/// then
///  w = 4, 4, 4, 4, 3, 3, 3, 2, 2, 1
/// and
///  N10 = 4
///  N20 = 4
///  N30 = 4
///  N40 = 4
///  N50 = 3
///  N60 = 3
///  N70 = 3
///  N80 = 2
///  N90 = 2
///  
// ===========================================================================

template<class T> vec<T> NStatistics( const vec<T>& v );

// ===========================================================================
//
// Median
//
// Calculate the median of vector invec from invec[beg] to invec.size( )
// (class T needs operator+, operator/, and 2 to be defined).
// 
// WARNING - WARNING - WARNING: Median expects (but does not check for)
// invec to be sorted.
//
// ===========================================================================

template <class T> inline T Median(const T& t1,const T& t2, const T& t3) {
  if ((t1 >= t2 || t3 >= t2) && (t1 <= t2 || t3 <= t2)) return t2;
  if ((t2 >= t3 || t1 >= t3) && (t2 <= t3 || t1 <= t3)) return t3;
  if ((t3 >= t1 || t2 >= t1) && (t3 <= t1 || t2 <= t1)) return t1;
  Assert(0==1);
}

template <class T> inline T Median( const vec<T> &v )
{    Assert( v.nonempty( ) );
     return v[ v.size( )/2 ];    }

template <class T> inline T Median( const vec<T> &invec, unsigned int beg ) {
  unsigned int end = invec.size( );
  ForceAssert( beg < end );
  if ( end - beg == 1 )
    return invec[beg];
  T result = 0;
  if ( ( end - beg ) % 2 == 0 ) {
    unsigned int posB = beg + ( ( end - beg ) / 2 );
    unsigned int posA = posB - 1;
    result = ( invec[posA] + invec[posB] ) / (T) 2;
  }
  else
    result = invec[ beg + ( ( end - beg - 1 ) / 2 ) ];
  return result;
}

// ===========================================================================
//
// Absolute value
//
// ===========================================================================

template <class T> inline T Abs(const T &a) { return a >= 0 ? a : -a; }
template <> inline Bool Abs<Bool>(const Bool &a) { return a; }

///round() is not declared for alpha
#if __alpha
#define round(X) nint(X)
#endif

// ===========================================================================
//
// Exclusive or
//
// ===========================================================================

template <class T> Bool Xor(T a, T b) { return a? (!(Bool)b) : ((Bool)b); }
template <class T> bool XOR(T a, T b) { return a? (!(bool)b) : ((bool)b); }

/// ===========================================================================
///
/// IntervalOverlap:  Find the overlap between intervals [a,b) and [c,d), where
/// a,b,c,d are integers.  Return 0 if they don't overlap.
///
/// ===========================================================================

inline int IntervalOverlap( int a, int b, int c, int d )
{    return Max( 0, Min(b, d) - Max(a, c) );    }

inline unsigned int IntervalOverlap( unsigned int a, unsigned int b, 
     unsigned int c, unsigned int d )
{    unsigned int m = Max(a, c), M = Min(b, d);
     if ( m <= M ) return M - m;
     return 0;    }

inline longlong IntervalOverlap( longlong a, longlong b, longlong c, longlong d )
{    longlong m = Max(a, c), M = Min(b, d);
     if ( m <= M ) return M - m;
     return 0;    }

inline double IntervalOverlap( double a, double b, double c, double d )
{    double m = Max(a, c), M = Min(b, d);
     if ( m <= M ) return M - m;
     return 0;    }


/// Class to store a quadratic function.
/// operator() returns the function's value at x.
/// solutions() returns solutions to the equation (*this)(x) = 0;

struct QuadraticFunction { 
  double a, b, c; 
  QuadraticFunction(double a=0.0, double b=0.0, double c=0.0): 
    a(a), b(b), c(c) {}
  double operator()(double x) const {
    return (((a*x)+b)*x)+c;  // two mults, not three
  }
  /// Return the solution(s) of the equation, or quiet_NaN() for both if there 
  /// aren't any.
  pair<double,double> solutions() const;
};

inline int ifloor( double x ) { return int( floor(x) ); }
inline int iceil( double x ) { return int( ceil(x) ); }


/// Compute the cosine of the angle between two vectors [begin1,end1),
/// and [begin2,end2),  
/// making sure to shrink the larger vector if the two are not the 
/// same size.
/// Note that we have to put the data into a vec<double> because otherwise
/// the multiplication inside inner_product can easily overrun the limits
/// of the data type pointed to by ForwardIterX.

template<class ForwardIter1, class ForwardIter2>
float Cosine( ForwardIter1 begin1, ForwardIter1 end1,
	      ForwardIter2 begin2, ForwardIter2 end2) {
  vec<double> d1(begin1, end1);
  vec<double> d2(begin2, end2);
  int size = static_cast<int>(min(d1.size(), d2.size()));
  double norm1 = 0, norm2 = 0, ret = 0;
  norm1 = sqrt(inner_product(d1.begin(), d1.begin() + size, d1.begin(), 0.0));
  norm2 = sqrt(inner_product(d2.begin(), d2.begin() + size, d2.begin(), 0.0));
  ret = inner_product(d1.begin(), d1.begin() + size, d2.begin(), 0.0)
    / (norm1*norm2);
  if (ret > 1.0001) {
    cerr << "cosine too high: " << ret << "\n";
    for (int i=0; i != size; ++i) {
      cerr << *(begin1+i) << " " << (*begin2+i) << "\n";
    }
  }
  return static_cast<float>(ret);
}

inline int gcd( int x, int y )
{    ForceAssert( x >= 1 && y >= 1 );
     if ( x > y ) swap( x, y );
     if ( y % x == 0 ) return x; 
     return gcd( y % x, x );    }

Bool CombineMeasurements( const double g1, const double g2, const double d1,
     const double d2, const double dmult, double& g, double& d );

// Inverse cdf of normal distribution

double InverseNormalCDF(double p, double mu=0.0, double sigma=1.0);

// Point-estimate of binomial parameter. That is, given the number of
// passes observed in a set of trials, determine a good single value
// for the binomial probability parameter p.

double EstProbability(int passes, int trials);

// Confidence interval for binomial parameter. That is, given the
// number of passes observed in a set of trials, determine a range of
// values for the binomial probability parameter p, such that an event
// with probability <= error must have occurred if the true
// probability is not in the returned range.

pair<double, double> BinomialConfidenceInterval(int passes, int trials, double error=0.05);

///Phred quality score, with small number correction.

float QualScoreFloat(int correct, int wrong);

///Phred quality score, with small number correction, rounded to int.

inline int QualScore(int correct, int wrong){
  return int(round(QualScoreFloat(correct, wrong)));
}

/// Inlined, templatized square function.

template <class T> T Square(const T & t) { return t*t; }

/// Pow, integer power n^k.

template<class T> T IPow( T n, int k )
{    T p = 1;
     for ( int i = 0; i < k; i++ )
          p *= n;
     return p;    }

#endif
