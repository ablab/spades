///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/// Wraps an STL vector, adds asserts for debugging and useful methods.
/// \class Vec
/// Vec.h defines class vec, which wraps the STL class vector, in such a way
/// that if compiled with NDEBUG on, it is (or should be) the same as a vector
/// (with some added functionality -- see below), but otherwise does run-time 
/// checking of each vector reference to make sure it is in range.
///
/// In debug mode, resize and reserve will fail if you (in effect) ask for more
/// than 100GB (or 1.5 GB on 32-bit systems).
///
/// PLEASE help keep Vec.h tidy by avoiding placing functions that operate on
/// vec in this file - consider using VecUtilities.h instead.
///
/// See Also: VecUtilities.h

#ifndef VEC_H
#define VEC_H

#include <unistd.h>

#include <algorithm>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <memory>
#include <strstream>
#include <vector>
#include <functional>

#include "String.h"
#include "system/Assert.h"
#include "system/StaticAssert.h"
#include "STLExtensions.h"
#include "system/Types.h"
#include "system/TraceVal.h"
#include "feudal/BinaryStream.h"
#include "system/SortInPlace.h"
#include "system/System.h"
#include "Compare.h"
#include "system/file/FileReader.h"

// Size checker: identity operation unless it FatalErr's and never returns.
template <class T>
inline size_type ValidatedSize(ulonglong i)
{
  //  PRINT(i);
#ifndef NDEBUG
  static const ulonglong TOO_BIG = (4 == sizeof(char *)) 
    ? (ulonglong)(INT_MAX) // 2 GB
    : ULLCONST(1000 * 1000) * ULLCONST(1000 * 1000); // ~1TB

  if ( ulonglong(i) * ulonglong(sizeof(T)) > TOO_BIG ) {
    FatalErr( "Attempt to resize vec<T> object to " << i 
	      << ", where sizeof(T) = " << sizeof(T) << "."
	      << "\nTOO_BIG=" << TOO_BIG 
	      << ", request=" << i << " T's" );
  }
#endif
  return i;
}

/////////////////////////////////////////////////////////////////////////////
//
//  vec Class Declaration and Template Definitions
//

template <class T> class vec : public vector<T> {
    typedef std::vector<T> BaseT;

 public:

// ===========================================================================
//
// CONSTRUCTORS
//
// ===========================================================================

  vec()            : vector< T >() {}
  
  explicit vec(size_type n) : vector< T >(ValidatedSize<T>(n)) {}
  
  vec(size_type n, const T& value) : vector< T >( ValidatedSize<T>(n), value ) {}

  enum ConstructorBehavior { IDENTITY };

  vec( size_type n, const ConstructorBehavior constructor_type  ) : 
    vector<T>(ValidatedSize<T>(n))
  {    ForceAssert( constructor_type == IDENTITY );
       for ( size_type i = 0; i < n; i++ )
            (*this)[i] = i;    }

  vec(const vector<T>& v) : vector<T>( v ) {}

  // vec<double> can be constructed e.g. from "{1.3,5.2,9}".

  explicit vec( const String& s );

  template<class ForwardIterator>
  vec(ForwardIterator first, ForwardIterator last) : 
    vector<T>( first, last ) {}

  ///Asserts index within bounds.
  typename vector<T>::reference operator[]( size_type i ) {
    AssertLt( i,  vector<T>::size() ); // Asserts index within bounds
    return vector<T>::operator[](i);   // ... and returns the element
  }
  
  ///Asserts index within bounds.
  typename vector<T>::const_reference operator[](size_type i) const {
    AssertLt( i, vector<T>::size() );  // Asserts index within bounds
    return vector<T>::operator[](i);   // ... and returns the element
  }

  void resize( size_type i, T c = T( ) ) { 
    vector<T>::resize( ValidatedSize<T>(i), c ); 
  }

  void reserve( size_type i ) { vector<T>::reserve(ValidatedSize<T>(i)); }

  typename vector<T>::reference front( ) {
    AssertGt( vector<T>::size( ), 0u ); // Asserts index within bounds
    return vector<T>::front( );   // ... and returns the element
  }

  typename vector<T>::const_reference front( ) const {
    AssertGt( vector<T>::size( ), 0u ); // Asserts index within bounds
    return vector<T>::front( );   // ... and returns the element
  }

  typename vector<T>::reference back( ) {
    AssertGt( vector<T>::size( ), 0u ); // Asserts index within bounds
    return vector<T>::back( );   // ... and returns the element
  }

  typename vector<T>::const_reference back( ) const {
    AssertGt( vector<T>::size( ), 0u ); // Asserts index within bounds
    return vector<T>::back( );   // ... and returns the element
  }

// ===========================================================================
//
// FUNCTIONS TO PUSH ELEMENTS ONTO VECTORS
//
// ===========================================================================

  /// push_front (insert item before first element of vector)

  void push_front( const T& t1 )
  {    vector<T>::insert( vector<T>::begin( ), t1 );    }

  /// Generalized push_back, allowing up to 12 items pushed back at a time.

  void push_back( const T& t1 )
  {    vector<T>::push_back(t1);    }
  void push_back( const T& t1, const T& t2 )
  {    vector<T>::push_back(t1); vector<T>::push_back(t2);    }
  void push_back( const T& t1, const T& t2, const T& t3 )
  {    vector<T>::push_back(t1); vector<T>::push_back(t2);
       vector<T>::push_back(t3);    }
  void push_back( const T& t1, const T& t2, const T& t3, const T& t4 )
  {    vector<T>::push_back(t1); vector<T>::push_back(t2);
       vector<T>::push_back(t3); vector<T>::push_back(t4);    }
  void push_back( const T& t1, const T& t2, const T& t3, const T& t4, const T& t5 )
  {    vector<T>::push_back(t1); vector<T>::push_back(t2);
       vector<T>::push_back(t3); vector<T>::push_back(t4);
       vector<T>::push_back(t5);    }
  void push_back( const T& t1, const T& t2, const T& t3, const T& t4, const T& t5,
       const T& t6 )
  {    push_back( t1, t2, t3, t4, t5 );
       push_back( t6 );     }
  void push_back( const T& t1, const T& t2, const T& t3, const T& t4, const T& t5,
       const T& t6, const T& t7 )
  {    push_back( t1, t2, t3, t4, t5 );
       push_back( t6, t7 );     }
  void push_back( const T& t1, const T& t2, const T& t3, const T& t4, const T& t5,
       const T& t6, const T& t7, const T& t8 )
  {    push_back( t1, t2, t3, t4, t5 );
       push_back( t6, t7, t8 );     }
  void push_back( const T& t1, const T& t2, const T& t3, const T& t4, const T& t5,
       const T& t6, const T& t7, const T& t8, const T& t9 )
  {    push_back( t1, t2, t3, t4, t5 );
       push_back( t6, t7, t8, t9 );     }
  void push_back( const T& t1, const T& t2, const T& t3, const T& t4, const T& t5,
       const T& t6, const T& t7, const T& t8, const T& t9, const T& t10 )
  {    push_back( t1, t2, t3, t4, t5 );
       push_back( t6, t7, t8, t9, t10 );     }
  void push_back( const T& t1, const T& t2, const T& t3, const T& t4, const T& t5,
       const T& t6, const T& t7, const T& t8, const T& t9, const T& t10,
       const T& t11 )
  {    push_back( t1, t2, t3, t4, t5 );
       push_back( t6, t7, t8, t9, t10, t11 );     }
  void push_back( const T& t1, const T& t2, const T& t3, const T& t4, const T& t5,
       const T& t6, const T& t7, const T& t8, const T& t9, const T& t10,
       const T& t11, const T& t12 )
  {    push_back( t1, t2, t3, t4, t5 );
       push_back( t6, t7, t8, t9, t10, t11, t12 );     }

  // push: construct object from up to ten arbitrary arguments, then push it back.

  template<class X1>
  void push( const X1& x1 )
  {    vector<T>::push_back( T(x1) );    }
  template<class X1, class X2>
  void push( const X1& x1, const X2& x2 )
  {    vector<T>::push_back( T(x1, x2) );    }
  template<class X1, class X2, class X3>
  void push( const X1& x1, const X2& x2, const X3& x3 )
  {    vector<T>::push_back( T(x1, x2, x3) );    }
  template<class X1, class X2, class X3, class X4>
  void push( const X1& x1, const X2& x2, const X3& x3, const X4& x4 )
  {    vector<T>::push_back( T(x1, x2, x3, x4) );    }
  template<class X1, class X2, class X3, class X4, class X5>
  void push( const X1& x1, const X2& x2, const X3& x3, const X4& x4, const X5& x5 )
  {    vector<T>::push_back( T(x1, x2, x3, x4, x5) );    }
  template<class X1, class X2, class X3, class X4, class X5, class X6>
  void push( const X1& x1, const X2& x2, const X3& x3, const X4& x4, const X5& x5,
       const X6& x6 )
  {    vector<T>::push_back( T(x1, x2, x3, x4, x5, x6) );    }
  template<class X1, class X2, class X3, class X4, class X5, class X6, class X7>
  void push( const X1& x1, const X2& x2, const X3& x3, const X4& x4, const X5& x5,
       const X6& x6, const X7& x7 )
  {    vector<T>::push_back( T(x1, x2, x3, x4, x5, x6, x7) );    }
  template<class X1, class X2, class X3, class X4, class X5, class X6, class X7,
       class X8>
  void push( const X1& x1, const X2& x2, const X3& x3, const X4& x4, const X5& x5,
       const X6& x6, const X7& x7, const X8& x8 )
  {    vector<T>::push_back( T(x1, x2, x3, x4, x5, x6, x7, x8) );    }
  template<class X1, class X2, class X3, class X4, class X5, class X6, class X7,
       class X8, class X9>
  void push( const X1& x1, const X2& x2, const X3& x3, const X4& x4, const X5& x5,
       const X6& x6, const X7& x7, const X8& x8, const X9& x9 )
  {    vector<T>::push_back( T(x1, x2, x3, x4, x5, x6, x7, x8, x9) );    }
  template<class X1, class X2, class X3, class X4, class X5, class X6, class X7,
       class X8, class X9, class X10>
  void push( const X1& x1, const X2& x2, const X3& x3, const X4& x4, const X5& x5,
       const X6& x6, const X7& x7, const X8& x8, const X9& x9, const X10& x10 )
  {    vector<T>::push_back( T(x1, x2, x3, x4, x5, x6, x7, x8, x9, x10) );    }
  template<class X1, class X2, class X3, class X4, class X5, class X6, class X7,
       class X8, class X9, class X10, class X11>
  void push( const X1& x1, const X2& x2, const X3& x3, const X4& x4, const X5& x5,
       const X6& x6, const X7& x7, const X8& x8, const X9& x9, const X10& x10,
       const X11& x11 )
  {    vector<T>::push_back( T(x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11) );    }
  template<class X1, class X2, class X3, class X4, class X5, class X6, class X7,
       class X8, class X9, class X10, class X11, class X12>
  void push( const X1& x1, const X2& x2, const X3& x3, const X4& x4, const X5& x5,
       const X6& x6, const X7& x7, const X8& x8, const X9& x9, const X10& x10,
       const X11& x11, const X12& x12 )
  {    vector<T>::push_back( 
            T(x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12) );    }

  void push_back_copies( const T& t1, const size_type n )
  {    for ( size_type i = 0; i < n; i++ )
            vector<T>::push_back(t1);    }

  template <class U>
  void append( const vec<U>& y ) 
  {    insert( this->end( ), y.begin( ), y.end( ) );    }
  
  void append( const vec<T>& y, size_type i, size_type j ) {
    if ( j == y.size( ) ) insert( this->end( ), y.begin( ) + i, y.end( ) );
    else insert( this->end( ), y.begin( ) + i, y.begin( ) + j );   
  }

  // appends values in y, but only those whose indices are in entries
  // IDX should be either (unsigned) int or longlong depending on the size of y
  template<typename IDX> 
  void append( const vec<T>& y, const vec<IDX>& entries ) {
    AssertLe(y.size(), numeric_limits<IDX>::max());
    const size_type n = this->size( );
    resize( n + entries.size( ) );
    for ( size_type j = 0; j < entries.size( ); j++ )
      (*this)[ n + j ] = y[ entries[j] ];
  }

// ===========================================================================
//
// FUNCTIONS TO TEST IF VECTOR HAS AN ATTRIBUTE
//
// ===========================================================================

  Bool nonempty( ) const { return ! this->empty( ); }
  Bool solo( ) const { return this->size( ) == 1; }

  Bool Ordered( ) const
  {    for ( size_type i = 1; i < this->size( ); i++ )
            if ( (*this)[i] < (*this)[i-1] ) return False;
       return True;    }

  Bool UniqueOrdered( ) const
  {    for ( size_type i = 1; i < this->size( ); i++ )
            if ( (*this)[i] <= (*this)[i-1] ) return False;
       return True;    }

// ===========================================================================
//
// FUNCTIONS TO EXTRACT PART OF A VECTOR
//
// ===========================================================================

  /// SetToSubOf: Set *this to the len entries of that, starting at pos.
  /// It is OK for that to be the same as this.
  void SetToSubOf( vec const& that, size_type pos, size_type len )
  { AssertLe(pos,that.size());
    AssertLe(len,that.size()-pos);
    if ( this != &that )
    { BaseT::assign(that.begin()+pos,that.begin()+pos+len); }
    else
    { BaseT::resize(pos+len);
      BaseT::erase(BaseT::begin(),BaseT::begin()+pos); } }
  // Implementation note:  the STL carefully accounts for self-assignment, so
  // a simpler implementation would be simply
  // BaseT::assign(that.begin()+pos,that.begin()+pos+len);
  // regardless of whether we are slicing ourselves or someone else.
  // However, the STL implementation appears to assume that the elements of the
  // vector also handle self-assignment correctly, which in the case of our
  // stuff hardly seems like a good assumption.  The resize/erase technique
  // implemented above is almost as efficient, and only assumes that element
  // copying to non-self works correctly.

  inline friend vec<T> SubOf( const vec<T>& x, size_type start, size_type n )  {
    vec<T> s(n);
    copy(x.begin() + start, x.begin() + start + n, s.begin());
    return s;
  }


  // SetToSubOf: Set *this to the entries in x pointed to by indices.
  // IDX should be either (unsigned) int or longlong depending on the size of x
  template<typename IDX>
  void SetToSubOf( const vec<T>& x, const vec<IDX>& indices ) {
    AssertLe(x.size(), numeric_limits<IDX>::max());
    const size_type n = indices.size( );
    resize(n);
    for ( size_type i = 0; i < n; i++ )
      (*this)[i] = x[ indices[i] ];
  }

  // IDX should be either (unsigned) int or longlong depending on the size of x
  template<typename IDX>
  inline friend vec<T> SubOf( const vec<T>& x, const vec<IDX>& indices ) {
    AssertLe(x.size(), numeric_limits<IDX>::max());
    const size_type n = indices.size( );
    vec<T> s(n);
    for ( size_type i = 0; i < n; i++ )
      s[i] = x[ indices[i] ];
    return s;
  }
    
  void SetToRangeOf( const vec<T>& v, size_type i, size_type j ) {
    AssertLe( i, j );
    assign(v.begin() + i, v.begin() + j);
  }

  inline friend vec<T> RangeOf( const vec<T>& v, size_type i, size_type j ) {
    AssertLe( i, j );
    vec<T> r(v.begin() + i, v.begin() + j);
    return r;
  }

  int isize( ) const { return this->size( ); }

  void SetCat( const vec<T>& v1, const vec<T>& v2 ) {
    *this = v1;
    append(v2);
  }

  void clear_and_resize( size_type n ) {
    this->clear( );
    resize(n);
  }

  void resize_and_set( size_type n, const T& x ) {
    this->clear();
    resize(n, x);
  }

  void SetToReverseOf( const vec<T>& v ) {
    resize(v.size());
    reverse_copy(v.begin(), v.end(), this->begin());
  }

  void ReverseMe( ) {
    reverse(this->begin(), this->end());
  }

  Bool Contains( const vec<T>& v ) const {
    return (search(this->begin(), this->end(), v.begin(), v.end()) != this->end());
  }

  Bool Contains( const vec<T>& v, size_type pos ) const {
    if ( v.size( ) + pos > this->size( ) )
      return False;
    size_type j;
    for ( j = 0; j < v.size( ); j++ )
      if ( (*this)[pos+j] != v[j] )
	break;
    return j == v.size( );
  }
  
// ===========================================================================
//
// FUNCTIONS TO ERASE ELEMENTS FROM VECTORS - PART 1 (MEMBER FUNCTIONS)
//
// ===========================================================================

  /// Erase: erase range of elements, where range is given by half-open interval.

  void Erase( size_type start, size_type stop ) {
    erase( this->begin( ) + start, this->begin( ) + stop );
  }

  /// EraseValue: erase all entries having the given value.
  void EraseValue( const T& x ) {
    erase(remove(this->begin(), this->end(), x), this->end());
  }

  /// print values to ostream, separated by sep.
  void Print(ostream & os, const char * sep = " ") const {
    copy(this->begin(), this->end(), ostream_iterator<T>(os, sep));
  }

  /// print values to ostream, separated by sep, with newline at end.
  void Println(ostream & os, const char * sep = " ") const {
    Print(os, sep); os << endl;
  }

  ///Set myself from text stream containing list of values of unknown length.
  void ReadFromTextStream(istream & is) {
    this->clear();
    T t;
    while (true) {
      is >> t;
      if (!is) break;
      push_back(t);
    }
  }

  /// CountValue: count all entries having the given value.
  size_type CountValue( const T& x ) const {
    return count(this->begin(), this->end(), x);
  }


// ===========================================================================
//
// MEMBER FUNCTIONS TO DO ARITHMETIC ON VECTORS.
//
// ===========================================================================

  /// Multiply a vector by a constant of type X.

  template<class X>
  vec & operator*=(const X & x) {
    const size_type S = this->size();
    for (size_type i = 0; i != S; ++i) {
      (*this)[i] = static_cast<T>( x * (*this)[i]);
    }
    return *this;
  }

  /// Divide a vector by a constant of type X.

  template<class X>
  vec & operator/=(const X & x) {
    return this->operator*=(1.0/x);
  }

  /// Add two vectors together.
  /// If vx is longer, add up to this vector's size.
  /// if vx is shorter, add up to vx's size only.

  template<class X> vec & operator+=(const vec<X> & vx) {
    const size_type S = min(this->size(), vx.size());
    for (size_type i = 0; i != S; ++i) {
      (*this)[i] += static_cast<T>(vx[i]);
    }
    return *this;
  }

  //stand-alone operators are implemented in terms of op=
  //See meyers, more effective C++, item 22 for reasons.

  /// Multiply a vector by a constant of type X.

  template<class X> friend vec operator*(const vec & v, const X & x) {
    return vec(v) *= x;
  }

  /// Multiply a vector by a constant of type X.

  template<class X> friend vec operator*(const X & x, const vec & v) {
    return vec(v) *= x;
  }
  /// Divide a vector by a constant of type X.

  template<class X> friend vec operator/(const vec & v, const X & x) {
    return vec(v) /= x;
  }

  /// Add two vectors together.
  /// If vx is longer, add up to this vector's size.
  /// if vx is shorter, add up to vx's size only.

  template<class X> friend vec operator+(const vec & v, const vec<X> & vx) {
    return vec(v) += vx;
  }

  /// NextDiff(i): return index of next element after i that is different from
  /// the ith element.

  inline size_type NextDiff( size_type i )
  {    size_type j;
       for ( j = i + 1; j < this->size( ); j++ )
            if ( (*this)[j] != (*this)[i] ) break;
       return j;    }

  friend int compare( vec const& v1, vec const& v2 )
  {   typedef typename std::vector<T>::const_iterator Itr;
      Itr itr1 = v1.begin();
      Itr itr2 = v2.begin();
      using std::min; Itr end = itr1 + min(v1.size(),v2.size());
      int result = 0;
      while ( !result && itr1 != end )
      { result = ::compare(*itr1,*itr2); ++itr1; ++itr2; }
      if ( !result ) result = ::compare(v1.size(),v2.size());
      return result;
  }
};

template <class T>
struct Serializability< vec<T> > : public ExternallySerializable {};

//
//  End of vec Class Declaration and Template Definitions
//
/////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////
//
// Mutating Functions
//

template<class T> vector<T> Reverse( const vector<T>& v ) {
  return vector<T>( v.rbegin( ), v.rend( ) );
}

template<class T> vec<T> Reverse( const vec<T>& v ) {
  return vec<T>( v.rbegin( ), v.rend( ) );
}

template<class T> void ReverseThis( vec<T>& v ) {
  v.ReverseMe();
}

template<class T> void RandomShuffle( vec<T>& v ) {
  random_shuffle( v.begin( ), v.end( ) );
}


/////////////////////////////////////////////////////////////////////////////
//
// Search Functions
//

template<class T> typename vec<T>::size_type LowerBound( const vec<T>& v, const T& x ) {
  return lower_bound( v.begin( ), v.end( ), x ) - v.begin( );    
}

template<class T> typename vec<T>::size_type UpperBound( const vec<T>& v, const T& x ) {
  return upper_bound( v.begin( ), v.end( ), x ) - v.begin( );   
}

template<class T> inline bool Member( const vector<T>& v, const T& x ) {
  return (find(v.begin(), v.end(), x) != v.end());
}

/// Return the position of an element in a vector, else -1.
template<class T> 
inline typename vec<T>::difference_type Position( const vector<T>& v, const T& x ) {
  typename vec<T>::const_iterator pos = find(v.begin(), v.end(), x);
  if (pos != v.end())
    return (pos - v.begin());
  else
    return -1;
}

/// BinPosition.  Return the position of an element in a sorted vector, else -1.
/// If the element appears more than once, the position of one of its instances
/// is returned.
template<class T, class U> 
inline typename vec<T>::difference_type BinPosition( const vector<T>& v, const U& x1 )
{    if ( v.size( ) == 0 ) return -1;
     T const& x(x1);
     typename vec<T>::size_type first = 0, last = v.size( ) - 1, next;
     while (1)
     {    if (first == last) return ( !(x < v[last]) && !(v[last] < x) ) ? last : -1;
          next = first + (last - first) / 2;
          if ( x < v[next] ) last = next;
          else if ( v[next] < x ) first = next + 1;
          else return next;    }    }

template<class T, class U>
inline Bool BinMember( const vector<T>& v, const U& x ) {
  return BinPosition( v, x ) >= 0;
}

/// BinSubset: determine if v is a subset of w; assumes w only is sorted and that there
/// is no repetition.
template<class T> inline Bool BinSubset( const vector<T>& v, const vector<T>& w )
{    for ( typename vec<T>::size_type i = 0; i < v.size( ); i++ )
          if ( !BinMember( w, v[i] ) ) return False;
     return True;    }

/// Determine if v is a subset of w; assumes that there is no repetition.
template<class T> inline Bool Subset( const vec<T>& v, const vec<T>& w )
{    for ( typename vec<T>::size_type i = 0; i < v.size( ); i++ )
          if ( !Member( w, v[i] ) ) return False;
     return True;    }



/////////////////////////////////////////////////////////////////////////////
//
// Numerical Functions
//

template<class T> int SizeSum( const vec< vec<T> >& v )
{    int sum = 0;
  for ( typename vec<T>::size_type i = 0; i < v.size( ); i++ )
          sum += v[i].size( );
     return sum;    }

template<class T> longlong SizeSumLong( const vec< vec<T> >& v )
{    longlong sum = 0;
  for ( typename vec<T>::size_type i = 0; i < v.size( ); i++ )
          sum += v[i].size( );
     return sum;    }

inline bool Nonnegative( const vec<int>& v )
{    for ( vec<int>::size_type i = 0; i < v.size( ); i++ )
          if ( v[i] < 0 ) return false;
     return true;    }


/////////////////////////////////////////////////////////////////////////////
//
// Destroy - returns the memory used by a vector, more or less
//

// Method given in Stroustrup's book, The C++ Programming Language, Special 
// Edition (2000), p. 457.

template<class T> inline void Destroy( vec<T>& v )
{    v.clear( );
     vec<T> tmp = v;
     v.swap(tmp);    }
template<class T1, class T2> inline void Destroy( vec<T1>& v1, vec<T2>& v2 )
{    Destroy(v1), Destroy(v2);    }
template<class T1, class T2, class T3> 
inline void Destroy( vec<T1>& v1, vec<T2>& v2, vec<T3>& v3 )
{    Destroy(v1), Destroy(v2), Destroy(v3);    }
template<class T1, class T2, class T3, class T4> 
inline void Destroy( vec<T1>& v1, vec<T2>& v2, vec<T3>& v3, vec<T4>& v4 )
{    Destroy(v1), Destroy(v2), Destroy(v3), Destroy(v4);    }
template<class T1, class T2, class T3, class T4, class T5> 
inline void Destroy( vec<T1>& v1, vec<T2>& v2, vec<T3>& v3, vec<T4>& v4, 
     vec<T5>& v5 )
{    Destroy(v1), Destroy(v2), Destroy(v3), Destroy(v4), Destroy(v5);    }



/////////////////////////////////////////////////////////////////////////////
//
// Sort Functions
//

template<class T>
inline void Sort( vec<T>& v )
{
  TRACEVAL_STOP_TRACING_COPIES;
  using std::sort; sort( v.begin( ), v.end( ) );
  TRACEVAL_START_TRACING_COPIES;
}

template<class T, class StrictWeakOrdering >
inline void Sort( vec<T>& v, StrictWeakOrdering comp )
{
  TRACEVAL_STOP_TRACING_COPIES;
  using std::sort; sort( v.begin( ), v.end( ), comp );
  TRACEVAL_START_TRACING_COPIES;
}

template<class T>
inline void ReverseSort( vec<T>& v )
{
  TRACEVAL_STOP_TRACING_COPIES;
  using std::sort; sort( v.rbegin( ), v.rend( ) );
  TRACEVAL_START_TRACING_COPIES;
}

template<class T, class StrictWeakOrdering>
inline void ReverseSort( vec<T>& v, StrictWeakOrdering comp )
{
  TRACEVAL_STOP_TRACING_COPIES;
  using std::sort; sort( v.rbegin( ), v.rend( ), comp );
  TRACEVAL_START_TRACING_COPIES;
}

template <class T, class Comp>
void Unique( vec<T>& v, Comp isEqual )
{
    if ( v.size() <= 1 )
        return;

    TRACEVAL_STOP_TRACING_COPIES;

    typedef typename vec<T>::iterator Itr;
    using std::iter_swap;

    Itr dest(v.begin());
    Itr end(v.end());
    for ( Itr itr(dest+1); itr != end; ++itr )
        if ( !isEqual(*itr,*dest) )
            iter_swap(itr, ++dest);
    v.erase(dest+1,end);

    TRACEVAL_START_TRACING_COPIES;
}

/// Leaves only unique elements in a vector \c v; these elements
/// will be also sorted.
template <class T, class StrictWeakOrdering, class EqualPredicate> 
inline void UniqueSort(vec<T> & v, StrictWeakOrdering comp, EqualPredicate equal)
{
    if ( v.size() <= 1 )
        return;
    Sort(v, comp);
    Unique(v, equal);
}

/// Leaves only unique elements in a vector \c v; these elements
/// will be also sorted.
template <class T> 
inline void UniqueSort(vec<T> & v)
{
  UniqueSort(v,less<T>(),equal_to<T>());
}



/////////////////////////////////////////////////////////////////////////////
//
// Erasing Functions
//

/// EraseIf: wrapper around erase-remove_if idiom.
template<class T> void EraseIf( vec<T>& v, bool (T::*f)( ) const )
{
  TRACEVAL_STOP_TRACING_COPIES;
  v.erase( remove_if( v.begin( ), v.end( ), mem_fun_ref(f) ), v.end( ) );
  TRACEVAL_START_TRACING_COPIES;
}

/// Another version of EraseIf: erase v[x] if erase[x] = True.
template<class T> void EraseIf( vec<T>& v, const vec<Bool>& erase )
{
  TRACEVAL_STOP_TRACING_COPIES;
  typename vec<T>::size_type count = 0;
  for ( typename vec<T>::size_type i = 0; i < v.size( ); i++ ) {
    if ( ! erase[i] ) {
      if ( count != i ) 
        v[count] = v[i];
      ++count;
    }
  }
  v.resize(count);
  TRACEVAL_START_TRACING_COPIES;
}

/// Another version of EraseIf: erase v[x] if erase[x] = True.
template<class T> void EraseUnless( vec<T>& v, const vec<Bool>& keep )
{
  TRACEVAL_STOP_TRACING_COPIES;
  typename vec<T>::size_type count = 0;
  for ( typename vec<T>::size_type i = 0; i < v.size( ); i++ ) {
    if ( keep[i] ) {
      if ( count != i ) 
        v[count] = v[i];
      ++count;
    }
  }
  v.resize(count);
  TRACEVAL_START_TRACING_COPIES;
}

/// EraseTheseIndices: Erase some elements of a vector, as determined by a sorted
/// list of indices.  Not efficiently implemented.
template<class T, typename IDX> void EraseTheseIndices( vec<T>& v, const vec<IDX>& these )
{
  AssertLe(v.size(), static_cast<size_t>(numeric_limits<IDX>::max()));
  TRACEVAL_STOP_TRACING_COPIES;
  typename vec<T>::size_type count = 0;
  for ( typename vec<T>::size_type i = 0; i < v.size( ); i++ )
    {    if ( !BinMember( these, i ) )
      {    if ( count != i ) v[count] = v[i];
      ++count;    }    }
  v.resize(count);
  TRACEVAL_START_TRACING_COPIES;
}


/////////////////////////////////////////////////////////////////////////////
//
// Input and Output Functions
//

/// IsAsciiVec() returns whether or not filename is a saved vec.  It
/// does this by checking whether the first line of the file is an
/// ASCII representation of a number.  If it encounters a non-digit,
/// non-whitespace character before it finds a newline, it returns
/// false.  Otherwise, it returns true.

bool IsAsciiVec( const String &filename );

/// Get the number of elements in saved vec in either ASCII or
/// Binary0 format.
longlong AsciiOrBinary0VecSize( const String& filename );

/// ***** The following is outmoded: please use BinaryRead2, *****
/// ***** etc. (below) except for backward compatibility.    *****
///
/// BinaryRead and BinaryWrite allow one to efficiently read and write a vector
/// of simple objects.  
/// WARNING: This will not work if members of the objects are discontiguous in 
/// memory, as can happen if members are forcibly aligned to word boundaries.
/// WARNING: this format is limited to vectors of size INT_MAX or less

template<class T> void BinaryWrite( int fd, const vec<T>& v ) {    
  ForceAssertLe(v.size(), static_cast<unsigned>(INT_MAX));
  int n = v.size( );
  WriteBytes( fd, &n, sizeof(int) );
  if ( n > 0 ) WriteBytes( fd, &v[0], (longlong) sizeof(T) * (longlong) n );    
}

/// BinaryWriteComplex is a different version of BinaryWrite.  The intention 
/// is that it not be called directly from a .cc file, but instead that the 
/// relevant .h file define BinaryWrite for a given class to be 
/// BinaryWriteComplex.  Ditto for BinaryReadComplex.
/// WARNING: this format is limited to vectors of size INT_MAX or less
template<class T> void BinaryWriteComplex( int fd, const vec<T>& v )
{    ForceAssertLe(v.size(), static_cast<unsigned>(INT_MAX));
     int n = v.size( );
     WriteBytes( fd, &n, sizeof(int) );
     for ( int i = 0; i < n; i++ )
          BinaryWrite( fd, v[i] );    }

template<class T> void BinaryRead( int fd, vec<T>& v ) {    
  int n;
  ReadBytes( fd, &n, sizeof(int) );
  v.resize(n);
  if ( n > 0 ) ReadBytes( fd, &v[0], (longlong) sizeof(T) * (longlong) n );    
}

template<class T> void BinaryReadComplex( int fd, vec<T>& v )
{    int n;
     ReadBytes( fd, &n, sizeof(int) );
     v.resize(n);
     for ( int i = 0; i < n; i++ )
          BinaryRead( fd, v[i] );    }

/// BinaryReadSubset: read selected entries from file written with BinaryWrite.
template<class T> void BinaryReadSubset( const String& filename, 
     const vec<int>& ids, vec<T>& v )
{    FileReader fr(filename.c_str());
     v.resize( ids.size( ) );
     for ( int i = 0; i < ids.isize( ); i++ )
     {    fr.seek( sizeof(int) + ids[i] * sizeof(T) );
          fr.read( &v[i], sizeof(T) );    }   }

/// ***** The following is outmoded: please use BinaryRead2, *****
/// ***** etc. (below) except for backward compatibility.    *****

template<class T> void BinaryWrite0( const String& filename, const vec<T>& v )
{    int fd = OpenForWrite(filename);
     longlong n = v.size( );
     String length = ToString(n) + "\n";
     int k = length.size( );
     WriteBytes( fd, length.c_str( ), k );
     if ( n > 0 ) WriteBytes( fd, &v[0], (longlong) sizeof(T) * (longlong) n );
     close(fd);    }

template<class T> void BinaryRead0( const String& filename, vec<T>& v )
{    String ns;
     {    Ifstream( in, filename );
          in >> ns;    }
     ForceAssert( ns.IsInt( ) );
     longlong n = ns.Int( );
     FileReader fr(filename.c_str());
     fr.seek( ns.size( ) + 1 );
     v.resize(n);
     if ( n > 0 ) fr.read( &v[0], sizeof(T)*n );    }

template<class T> void BinaryReadSubset0( const String& filename, 
     const vec<int>& ids, vec<T>& v )
{    String ns;
     {    Ifstream( in, filename );
          in >> ns;    }
     ForceAssert( ns.IsInt( ) );
     longlong n = ns.Int( );
     int k = ns.size( );
     FileReader fr(filename.c_str());
     v.resize( ids.size( ) );
     for ( int i = 0; i < ids.isize( ); i++ )
     {    ForceAssertGe( ids[i], 0 );
          ForceAssertLt( ids[i], n );
          fr.seek( k + 1 + ids[i] * sizeof(T) );
          fr.read( &v[i], sizeof(T) );    }    }

//=============================================================================
//=============================================================================

/// BinaryRead2, BinaryWrite2, BinaryReadSubset2: same as above, but smarter:
///
/// - first 34 bytes = "binary format 2, header = 3 lines\n";
/// - next  13 bytes = number of entries, as ascii, left justified and blank-padded;
/// - next  15 bytes = "\nlittle endian\n" or "\nbig endian   \n";
/// - the objects.
///
/// This format should NOT be modified.  If we ever have need to modify it, we 
/// should create BinaryRead3, etc.
///
/// BREAD2 and BREADX2 mirror READ and READX in System.h.

template<class T> void BinaryWrite2( const String& filename, const vec<T>& v );
template<class T> void BinaryRead2( const String& filename, vec<T>& v,
     bool strict = false, const Bool append = false );
template<class T> void BinaryReadSubset2( const String& filename,
     const vec<int>& ids, vec<T>& v, Bool append = False, bool strict = false );
template<class T> void BinaryReadRange2( const String& filename,
     longlong from, longlong to, vec<T>& v, bool strict = false );
template<class T> longlong BinarySize2( const String& filename, bool strict = false );


///Append source to target.
void BinaryCat2(const String & target, const String & source);

/** Get the size of the data in a Binary2 file, return -1 if fail.
  That is, it will return 1 for a vec<char> file and 4 for a vec<int> file.
 Will return -1 but will not fail if the file is not in Binary2 format,
 so it can be used to test the format of a file.
*/
int GetBinary2ElementSize(const String & filename, bool strict = false);

#define BREAD2( FILE, TYPE, DATA )   \
     TYPE DATA;                      \
     BinaryRead2( FILE, DATA );

#define BREADX2( FILE, DATA )    \
     BinaryRead2( FILE, DATA );

//=============================================================================
//=============================================================================

/// BinaryRead3, BinaryWrite3, BinaryReadSubset3: same as above, but even smarter:
///
/// - first 34 bytes = "binary format 3, header = 4 lines\n";
/// - next  13 bytes = number of entries, as ascii, left justified and blank-padded;
/// - next  15 bytes = "\nlittle endian\n" or "\nbig endian   \n";
/// - next  34 bytes = "padding to make long word aligned\n";
/// - the objects.
///
/// This format should NOT be modified.  If we ever have need to modify it, we 
/// should create BinaryRead4, etc.
///
/// BREAD3 and BREADX3 mirror READ and READX in System.h.

template<class T> void BinaryWrite3( const String& filename, const vec<T>& v );
template<class T> void BinaryRead3( const String& filename, vec<T>& v, 
     bool strict = false, const Bool append = False );
template<class T> void BinaryReadSubset3( const String& filename,
     const vec<int>& ids, vec<T>& v, Bool append = False, bool strict = false );
template<class T> void BinaryReadRange3( const String& filename,
     longlong from, longlong to, vec<T>& v, bool strict = false );
template<class T> longlong BinarySize3( const String& filename, bool strict = false );
int GetBinary3ElementSize(const String & filename, bool strict = false );

///Append source to target.
void BinaryCat3(const String & target, const String & source);

#define BREAD3( FILE, TYPE, DATA )   \
     TYPE DATA;                      \
     BinaryRead3( FILE, DATA );

#define BREADX3( FILE, DATA )    \
     BinaryRead3( FILE, DATA );


/// This class allows the progressive writing of a binary format 3
/// file.  Implemented in VecTemplate.h.
///
/// If you want a lot of these at the same time, set the optional arg
/// keep_open to false, so that you don't run out of file descriptors.
///
/// Objects to be written should be simple objects that can be written
/// directly with WriteBytes.  Otherwise bad things will happen.

template <typename T> 
class Binary3Writer {
 public:
  Binary3Writer( ) : m_fd (-1) {}
  Binary3Writer( const String& filename, bool keep_open = true )
  { this->Open( filename, keep_open ); }
  ~Binary3Writer();

  void Write( const T& object );
  void WriteMultiple( const vec<T>& objects );

  void Open( const String& filename, bool keep_open = true );
  void Close();

  const String& Filename() { return m_filename; }
  // This will return the filename, even if the file has been Close()d.

 private:
  // Disallow copy constructor and assignment operator.  Not implemented.
  Binary3Writer( const Binary3Writer<T>& );
  Binary3Writer<T>& operator=( const Binary3Writer<T>& );

  int m_fd;           // used if m_keep_open
  String m_filename;  // used if ! m_keep_open
  longlong m_objectCount;
  bool m_keep_open;
};


/// \class Binary3Iter
/// This class lets you iterate through a binary format 3 vector file
/// stored on disk.  Use this if you just want to handle each element,
/// in order, in one pass through a loop.  This is more efficient in
/// terms of both time and memory (!) than loading the whole vector.
///
/// Usage for, eg, a vec<int>:
///
///   int s;
///   for( Binary3Iter<int> iter(filename, &s); iter.More(); iter.Next(&s) )
///
/// The constructor and each call to Next fill v with the next item.
/// Also, iter.N() = size of vector, and iter.I() = current index.

template <typename T>
class Binary3Iter {
public:
  Binary3Iter( const String& filename, T* p_to_fill, 
	       longlong max_memory = 128 * 1024 /* empirically good */ );

  void Next( T* p_to_fill );
  bool More() const { return m_globalIndex < m_globalSize; }
  longlong N() const { return m_globalSize; }
  longlong I() const { return m_globalIndex; }

private:
  void FillBuffer();

  FileReader mFR;
  longlong m_globalSize, m_globalIndex;
  unsigned int m_localIndex, m_maxsize;
  vec<T> m_data;
};



// ============================================================================
// ========================================================================

/// Helper for functions that load Binary format 2 or 3, whichever the file is.
int WhichBinaryFormat( const String& filename );

/// Return the number of elements of a binary format 2 or 3 file. 
/// Assert if file is not binary 2 or 3 format.
/// Note that this is not templatized, unlike the BinarySize*** functions.
longlong BinaryNumElements( const String & filename);




// ============================================================================
// ========================================================================

void PrettyPrint( ostream& o, const vec<int>& v, int max_items = 0,
     String terminator = "\n" );

void PrettyPrint( ostream& o, const vec<longlong>& v, int max_items = 0,
     String terminator = "\n" );

void PrettyPrint( ostream& o, const vec<double>& v, int max_items = 0,
     String terminator = "\n" );

void PrettyPrint( ostream& o, const vec<TraceInt>& v, int max_items = 0,
     String terminator = "\n" );


// CompactPrint: print a vector of objects, separated by blanks (default) or
// a user-specified separator.

template<class T> void CompactPrint( ostream& out, const vec<T>& v,
    String separator = " " )
{    for ( typename vec<T>::size_type i = 0; i < v.size( ); i++ )
     {    if ( i > 0 ) out << separator;
          out << v[i];    }    }


/// WriteAppend has an implementation for T = alignment_plus in Alignment.{h,cc}.

template<class T> void WriteAppend( const String& f, const vec<T>& v )
{    ForceAssert( !IsRegularFile( f + ".gz" ) );
     static longlong max_size_bound = longlong(10000000) * longlong(100000000);
     if ( !IsRegularFile(f) )
     {    ForceAssertLt( (longlong) v.size( ), max_size_bound );
          Ofstream( out, f );
          out << setfill( '0' ) << setw(15) << v.size( ) << "\n";
          for ( typename vec<T>::size_type i = 0; i <  v.size( ); i++ )
               out << v[i];    }
     else
     {    longlong n;
          {    Ifstream( in, f );
               in >> n;    }
          ForceAssertLt( n + (longlong) v.size( ), max_size_bound );
          /* TODO:
           * All this crap with the deprecated ostrstream class could be simply:
           * char buf[17]; sprintf(buf,"%015ld\n",n+v.size());
           * Unfortunately, killing the strstream include in this file breaks a
           * dozen cc files that refer to strstreams without doing their
           * own include. Note that osize is not unfrozen, so there's a memory
           * leak here, as well. */
          ostrstream osize;
          osize << setfill( '0' ) << setw(15) << n + v.size( ) << "\n";
          int fd = Open( f, O_WRONLY );
          WriteBytes( fd, osize.str( ), 16 );
          close(fd);
          ofstream out( f.c_str( ), ios::app );
          for ( typename vec<T>::size_type i = 0; i <  v.size( ); i++ )
               out << v[i];    }    }

/// a specialized version of WriteAppend for String

template <>
void WriteAppend( const String& f, const vec<String>& v );

template<class T> ostream& operator<<(ostream& s, const vec<T>& v)
{    s << v.size( ) << "\n";
     for ( typename vec<T>::size_type i = 0; i < v.size( ); i++ )
          s << v[i];
     return s;    }

#ifdef __DECCXX_VER
#pragma do_not_instantiate ostream& operator<<(ostream&, const vec<int>&)
#pragma do_not_instantiate ostream& operator<<(ostream&, const vec<longlong>&)
#pragma do_not_instantiate ostream& operator<<(ostream&, const vec<float>&)
#pragma do_not_instantiate ostream& operator<<(ostream&, const vec<String>&)
#endif

ostream& operator<<(ostream& s, const vec<unsigned short>& v);
ostream& operator<<(ostream& s, const vec<int>& v);
ostream& operator<<(ostream& s, const vec<longlong>& v);
ostream& operator<<(ostream& s, const vec<float>& v);
ostream& operator<<(ostream& s, const vec<double>& v);
ostream& operator<<(ostream& s, const vec<String>& v);

template<class T> istream& operator>>(istream& s, vec<T>& v)
{    typename vec<T>::size_type n;
     s >> n;
     v.resize(n);
     char c;
     s.get(c);
     for ( typename vec<T>::size_type i = 0; i < v.size( ); i++ )
       s >> v[i];   // Breaks cxx
     return s;    }


#ifdef __DECCXX_VER
#pragma do_not_instantiate istream& operator>>(istream&, vec<String>&)
#endif

istream& operator>>(istream& s, vec<String>& v);

/// Print out a matrix, with left-justified entries, and given separation between
/// columns.  (Justification may be changed by supplying an optional argument
/// consisting of a string of l's and r's.)

void PrintTabular( ostream& out, const vec< vec<String> >& rows, int sep,
     String justify = String( ) );

void PrintCSV(ostream& out, const vec< vec<String> >& rows);


void BinaryWrite( int fd, const vec< String >& v );
void BinaryRead( int fd, vec< String >& v );

#define For_(T,x,v) for( vec< T >::const_iterator x = v.begin(); x != v.end(); ++x )
#define ForMut_(T,x,v) for( vec< T >::iterator x = v.begin(); x != v.end(); ++x )

#endif
