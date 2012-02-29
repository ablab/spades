///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef STLEXTENSIONS_H
#define STLEXTENSIONS_H

using namespace std;

#include <functional>
#include <vector>
#include <bitset>
#include <algorithm>
#include <iostream>

#include <math.h>
#include "feudal/BinaryStream.h"
#include "system/StaticAssert.h"

/// minimum<T> is a function object.  
///
/// If f is an object of class minimum<T> and x and y are objects of
/// class T, then f(x,y) returns a copy of whichever of x or y is lower
/// as defined by T::operator<().
///
/// Example:
///
/// vector<int> V(10);
/// iota( V.begin(), V.end(), 1 );
///
/// copy( V.begin(), V.end(), ostream_iterator<int>( cout, " " ) );
/// cout << endl;
///
/// transform( V.begin(), V.end(), V.begin(), bind2nd(minimum<int>(), 5) );
///  
/// copy( V.begin(), V.end(), ostream_iterator<int>( cout, " " ) );
/// cout << endl;
///
/// The preceding code produces the following output:
///
/// 1 2 3 4 5 6 7 8 9 10
/// 1 2 3 4 5 5 5 5 5 5
 
template <class T>
struct minimum : public binary_function<T, T, T> 
{
  T operator() (const T& x, const T& y ) const { return ( x<y ? x : y ); }
};

/// maximum<T> is a function object.  
///
/// If f is an object of class maximum<T> and x and y are objects of
/// class T, then f(x,y) returns a copy of whichever of x or y has the
/// highest value as defined by T::operator<().
///
/// Example:
///
/// vector<int> V(10);
/// iota( V.begin(), V.end(), 1 );
///
/// copy( V.begin(), V.end(), ostream_iterator<int>( cout, " " ) );
/// cout << endl;
///
/// transform( V.begin(), V.end(), V.begin(), bind2nd(maximum<int>(), 5) );
///  
/// copy( V.begin(), V.end(), ostream_iterator<int>( cout, " " ) );
/// cout << endl;
///
/// The preceding code produces the following output:
///
/// 1 2 3 4 5 6 7 8 9 10
/// 5 5 5 5 5 6 7 8 9 10
 
template <class T>
struct maximum : public binary_function<T, T, T> 
{
  T operator() (const T& x, const T& y ) const { return ( x<y ? y : x ); }
};

/// modulus<float> is a function object that implements the template modulus<T>
/// modulus<double> is a function object that implements the template modulus<T>
///
/// The default modulus<T> defines its operator() as
///
/// { return x % y; }
///
/// floats and doubles do not have a % operator defined.  They require specific
/// math functions to perform modulus operations.  These implementations are provided
/// below.

namespace std {
    template <>
    struct modulus<float> : public binary_function<float, float, float>
    {
        float operator() (const float& x, const float& y) const { return fmodf( x, y ); }
    };

    template <>
    struct modulus<double> : public binary_function<double, double, double>
    {
        double operator() (const double& x, const double& y) const { return fmod( x, y ); }
    };
}

/// address_of<T> and dereference<T> are function objects.
///
/// If f is an object of class address_of<T> and x is an object of
/// class T, then f(x) returns the address of x, i.e. a T*.
///
/// If f is an object of class dereference<T> and x is a pointer to
/// an object of class T, then f(x) returns a copy of *x.
///
/// Example:
///
/// vector<int> V1(10);
/// iota( V1.begin(), V1.end(), 1 );
///
/// copy( V1.begin(), V1.end(), ostream_iterator<int>( cout, " " ) );
/// cout << endl;
///
/// vector<int*> V_ptrs(10);
/// transform( V1.begin(), V1.end(), V_ptrs.begin(), address_of<int>() );
///
/// vector<int> V2(10);
/// transform( V_ptrs.begin(), V_ptrs.end(), dereference<int>() );
///
/// copy( V2.begin(), V2.end(), ostream_iterator<int>( cout, " " ) );
/// cout << endl;
///  
/// The preceding code will produce the following output:
///
/// 1 2 3 4 5 6 7 8 9 10
/// 1 2 3 4 5 6 7 8 9 10
///
/// Note that address_of<T> is nearly always inexpensive, while 
/// dereference<T> calls the copy constructor of T.
 
template <class T>
struct address_of : public unary_function<T, T*> 
{
  const T* operator() (const T& x ) const { return &x; }
  T* operator() ( T& x ) const { return &x; }
};

template <class T>
struct dereference : public unary_function<T*, T> 
{
  T operator() (const T* x ) const { return *x; }
};


// g++ 3.x defines select1st and select2nd in <ext/functional>

#if __GNUC__ > 2
#include <ext/functional>
using __gnu_cxx::select1st;
using __gnu_cxx::select2nd;
#endif

// cxx does not define select1st and select2nd, so we provide them here.

#if defined( __DECCXX_VER )

/*
 *
 * Copyright (c) 1994
 * Hewlett-Packard Company
 *
 * Permission to use, copy, modify, distribute and sell this software
 * and its documentation for any purpose is hereby granted without fee,
 * provided that the above copyright notice appear in all copies and
 * that both that copyright notice and this permission notice appear
 * in supporting documentation.  Hewlett-Packard Company makes no
 * representations about the suitability of this software for any
 * purpose.  It is provided "as is" without express or implied warranty.
 *
 *
 * Copyright (c) 1996
 * Silicon Graphics Computer Systems, Inc.
 *
 * Permission to use, copy, modify, distribute and sell this software
 * and its documentation for any purpose is hereby granted without fee,
 * provided that the above copyright notice appear in all copies and
 * that both that copyright notice and this permission notice appear
 * in supporting documentation.  Silicon Graphics makes no
 * representations about the suitability of this software for any
 * purpose.  It is provided "as is" without express or implied warranty.
 */

template <class T_Pair>
struct Select1st : public unary_function<T_Pair, typename T_Pair::first_type> {
  const typename T_Pair::first_type& operator()(const T_Pair& __x) const {
    return __x.first;
  }
};

template <class T_Pair>
struct Select2nd : public unary_function<T_Pair, typename T_Pair::second_type>

  const typename T_Pair::second_type& operator()(const T_Pair& __x) const {
    return __x.second;
  }
};

template <class T_Pair> 
struct select1st
    : public Select1st<T_Pair> {};
template <class T_Pair> 
struct select2nd 
    : public Select2nd<T_Pair> {};

/*
 * End HP/SGI copyrighted material.
 */

#endif

// g++ 3.x defines is_sorted in <ext/algorithm>
// g++ 3.x defines iota in <ext/numeric>

#if __GNUC__ > 2
#include <ext/algorithm>
using __gnu_cxx::is_sorted;
#include <ext/numeric>
using __gnu_cxx::iota;
#else
#if __GNUC__ <= 2
#include <algorithm>
#include <numeric>
#endif
#endif

// cxx does not have iota or is_sorted, so we provide it here.

#if defined( __DECCXX_VER )

template <class ForwardIterator, class T>
void iota(ForwardIterator first, ForwardIterator last, T value) {
  while (first != last) *first++ = value++;
}

/*
 *
 * Copyright (c) 1994
 * Hewlett-Packard Company
 *
 * Permission to use, copy, modify, distribute and sell this software
 * and its documentation for any purpose is hereby granted without fee,
 * provided that the above copyright notice appear in all copies and
 * that both that copyright notice and this permission notice appear
 * in supporting documentation.  Hewlett-Packard Company makes no
 * representations about the suitability of this software for any
 * purpose.  It is provided "as is" without express or implied warranty.
 *
 *
 * Copyright (c) 1996
 * Silicon Graphics Computer Systems, Inc.
 *
 * Permission to use, copy, modify, distribute and sell this software
 * and its documentation for any purpose is hereby granted without fee,
 * provided that the above copyright notice appear in all copies and
 * that both that copyright notice and this permission notice appear
 * in supporting documentation.  Silicon Graphics makes no
 * representations about the suitability of this software for any
 * purpose.  It is provided "as is" without express or implied warranty.
 */

/// is_sorted, a predicated testing whether a range is sorted in
/// nondescending order.  This is an extension, not part of the C++
/// standard.

template <class _ForwardIter>
bool is_sorted(_ForwardIter __first, _ForwardIter __last)
{
  if (__first == __last)
    return true;

  _ForwardIter __next = __first;
  for (++__next; __next != __last; __first = __next, ++__next) {
    if (*__next < *__first)
      return false;
  }

  return true;
}

/// is_sorted, a predicated testing whether a range is sorted in
/// nondescending order.  This is an extension, not part of the C++
/// standard.
template <class _ForwardIter, class _StrictWeakOrdering>
bool is_sorted(_ForwardIter __first, _ForwardIter __last,
               _StrictWeakOrdering __comp)
{
  if (__first == __last)
    return true;

  _ForwardIter __next = __first;
  for (++__next; __next != __last; __first = __next, ++__next) {
    if (__comp(*__next, *__first))
      return false;
  }

  return true;
}


/*
 * End HP/SGI copyrighted material.
 */

#endif

///returns true if the order is strictly ascending.
template <class _ForwardIter>
bool is_sorted_strict(_ForwardIter __first, _ForwardIter __last)
{
  if (__first == __last)
    return true;

  _ForwardIter __next = __first;
  for (++__next; __next != __last; __first = __next, ++__next) {
    if (*__next <= *__first)
      return false;
  }

  return true;
}

///returns true if the order is strictly ascending.
template <class _ForwardIter, class _StrictWeakOrdering>
bool is_sorted_strict(_ForwardIter __first, _ForwardIter __last,
               _StrictWeakOrdering __comp)
{
  if (__first == __last)
    return true;

  _ForwardIter __next = __first;
  for (++__next; __next != __last; __first = __next, ++__next) {
    if (!(__comp(*__first, *__next) ) )
      return false;
  }

  return true;
}

/// sort_unique combines the STL algorithms sort and unique.  The uniquification
/// part of the algorithm uses iter_swap and operator< to do its dirtywork --
/// the same operations as are used by the sorting part of the algorithm --
/// unlike the STL version of unique, which uses operator== and assignment.
template <class RAItr>
RAItr sort_unique( RAItr start, RAItr const& stop )
{ using std::distance; using std::sort; using std::iter_swap;
  if ( distance(start,stop) < 2 ) return stop;
  sort(start,stop);
  RAItr dest(start);
  for ( ++start; start != stop; ++start )
    if ( *dest < *start && ++dest != start ) iter_swap(dest,start);
  return ++dest; }

template <class RAItr, class Comp>
RAItr sort_unique( RAItr start, RAItr const& stop, Comp comp )
{ using std::distance; using std::sort; using std::iter_swap;
  if ( distance(start,stop) < 2 ) return stop;
  sort(start,stop,comp);
  RAItr dest(start);
  for ( ++start; start != stop; ++start )
    if ( comp(*dest,*start) && ++dest != start ) iter_swap(dest,start);
  return ++dest; }

/// convenience function to uniquely sort a vector-like container that has
/// begin, end, and resize methods.
template <class V>
V& sort_unique( V& container )
{ typename V::iterator start(container.begin());
  container.resize(sort_unique(start,container.end())-start);
  return container; }

/// copy_if is an algorithm.
///
/// copy_if is so tremendously useful, it's hard to understand why it's
/// not in the standard.  For notes on this implementation, see
/// Effective STL (Meyers, 2001), Item 36, pp.154-156.

template<typename InputIterator, typename OutputIterator, typename Predicate>
OutputIterator copy_if( InputIterator begin, InputIterator end,
			OutputIterator destBegin, Predicate p )
{
  while ( begin != end )
  {
    if ( p(*begin) )
      *destBegin++ = *begin;
    ++begin;
  }

  return destBegin;
}


/// dereference_compare is a functor which allows comparison of
/// iterators using the operator< of the objects they point to.
///
/// usage example: suppose Foo::operator< exists.  Then:
///   list<Foo> foo_list;
///   vec< list<Foo>::iterator > foo_iters;
///   sort( foo_iters.begin(), foo_iters.end(), 
///         dereferenced_compare<list<Foo>::iterator>() );

template<typename Iterator>
struct dereferenced_compare : 
  public binary_function<Iterator, Iterator, bool>
{
  bool operator()(Iterator lhs, Iterator rhs) const {
    return( *lhs < *rhs );
  }
};

/*
   Template: update_min

   Update a running minimum: compare a value to the current minimum and update the minimum if the value is smaller.
 */
template <typename NumericType>
void update_min( NumericType& currentMin, const NumericType& val ) {
  if (val < currentMin )
    currentMin = val;
}

/*
   Template: update_max

   Update a running maximum: compare a value to the current maximum and update the maximum if the value is smaller.
 */
template <typename NumericType>
void update_max( NumericType& currentMax, const NumericType& val ) {
  if (val > currentMax )
    currentMax = val;
}

#define ForEach_Mut(x,containerType,c) \
    for ( containerType::iterator x = c.begin(); x != c.end(); x ++ ) 

#define ForEach(x,containerType,c) \
    for ( containerType::const_iterator x = c.begin() ; x != c.end(); x ++ )

template <class CONTAINER, class T> inline
  bool STLContains( const CONTAINER& c, const T& x ) {
  return c.find(x) != c.end();
}

/**
   Class: ShallowViewOf

   A shallow view of a given class: comparison operators are taken from the class,
   but copy constructor and assignment are the default.

   You would not normally create instances of ShallowViewOf<T>, but you can cast
   a (T *) to a (ShallowViewOf<T> *) before calling std::sort(), to avoid unnecessary
   deep copying.

   The template parameter T must be a model of STL concept LessThanComparable.
*/
template <class T>
class ShallowViewOf {
  char data_[sizeof(T)];
};

template <class T>
inline int operator<( const ShallowViewOf<T>& v1, const ShallowViewOf<T>& v2) { 
  return ((T&)v1) < ((T&)v2); 
}

template <class T>
inline int operator>( const ShallowViewOf<T>& v1, const ShallowViewOf<T>& v2) { 
  return ((T&)v1) > ((T&)v2); 
}

template <class T>
inline int operator<=( const ShallowViewOf<T>& v1, const ShallowViewOf<T>& v2) { 
  return ((T&)v1) <= ((T&)v2); 
}
template <class T>
inline int operator>=( const ShallowViewOf<T>& v1, const ShallowViewOf<T>& v2) { 
  return ((T&)v1) >= ((T&)v2); 
}
template <class T>
inline int operator==( const ShallowViewOf<T>& v1, const ShallowViewOf<T>& v2) { 
  return ((T&)v1) == ((T&)v2); 
}
template <class T>
inline int operator!=( const ShallowViewOf<T>& v1, const ShallowViewOf<T>& v2) { 
  return ((T&)v1) != ((T&)v2); 
}

/**
   Functor: order_ShallowView

   An ordering of ShallowViewOf<T> that just calls the given ordering on T.
*/
template <typename T, typename StrictWeakOrdering>
  struct order_ShallowView: public binary_function< T, T, bool > {
    private:
     StrictWeakOrdering orderOnT_;
    public:
     order_ShallowView<T, StrictWeakOrdering>( StrictWeakOrdering orderOnT ): orderOnT_(orderOnT) { }
    
    bool operator() ( const ShallowViewOf<T>& v1, const ShallowViewOf<T>& v2 ) const {
      return orderOnT_( (const T&)v1, (const T&)v2 );
    }
  };


/**
   Macro: DEFINE_SWAP

   Define std::swap() and std::iter_swap() routines for the given class to do a bitwise
   swap, bypassing any copy constructors or assignment operators defined for the class.
   This can speed up the execution of STL routines on containers of this class's items.
*/
#define DEFINE_SWAP(T)                                                                  \
 namespace std {                                                                        \
    template <> inline void swap< T >( T& v1, T& v2 ) {                                 \
      std::swap( (ShallowViewOf<T>&)v1, (ShallowViewOf<T>&)v2 );                        \
      STATIC_ASSERT_M( sizeof( ShallowViewOf<T> ) == sizeof(T), ShallowViewAlignmentProblem ); \
    }                                                                                   \
    template <> inline                                                                  \
      void iter_swap< vector<T>::iterator,                                              \
                      vector<T>::iterator > ( vector<T>::iterator i1,                   \
 	 				      vector<T>::iterator i2 ) {                \
      std::swap( *i1, *i2 );                                                            \
    }                                                                                   \
 }                                                                                      \
                                                                                        \
 typedef int define_swap_and_iter_swap__ ## T ## _


/**
   Functor: cmp_functor

   Turn a C++ comparator function into a binary functor that always
   calls this specific comparator function.  Since the comparator function
   is a template parameter of this functor, the compiler may have an easier
   time inlining the function body than if you had passed around
   a pointer to the function itself.
*/
template<class T, bool myComparator(const T&, const T&) >
  struct cmp_functor: public binary_function< T, T, bool > {
    bool operator() ( const T& a, const T& b ) const {
      return myComparator(a,b);
    }
  };

#define COMPARE_BY2(T,F1,F2)  \
 friend bool operator< ( const T& v1, const T& v2 ) { \
   return v1.F1 < v2.F1 ? true : \
     ( v1.F1 > v2.F1 ? false : ( v1.F2 < v2.F2 ) );  \
 }

#define COMPARE_BY3(T,F1,F2,F3)  \
 friend bool operator< ( const T& v1, const T& v2 ) { \
   return  \
     (v1.F1) < (v2.F1) ? true : \
     (v1.F1) > (v2.F1) ? false: \
     (v1.F2) < (v2.F2) ? true : \
     (v1.F2) > (v2.F2) ? false : \
       (v1.F3) < (v2.F3);        \
 }

// Return the canonical version of the two pairs
// (a,b) and (b,a).
template <class _T1, class _T2>
inline pair< _T1, _T2 > CanonPair( const pair< _T1, _T2 >& p ) {
  return p.first <= p.second ? p : make_pair( p.second, p.first );
}

template <class _T1, class _T2>
inline pair< _T1, _T2 > CanonPair( const _T1& v1, const _T2& v2  ) {
  return CanonPair( make_pair( v1, v2 ) );
}



/// Class triple is just like class pair, but has three elements.

template <class _T1, class _T2, class _T3>
struct triple {
  typedef _T1 first_type;
  typedef _T2 second_type;
  typedef _T3 third_type;
  _T1 first;
  _T2 second;
  _T3 third;
  triple() : first(), second(), third() {}
  triple(const _T1& __a, const _T2& __b, const _T3& __c) 
       : first(__a), second(__b), third(__c) {}
  template <class _U1, class _U2, class _U3>
  triple(const triple<_U1, _U2, _U3>& __p) 
       : first(__p.first), second(__p.second), third(__p.third) {}

     size_t writeBinary( BinaryWriter& writer ) const
     {    size_t count = writer.write(first);
          count += writer.write(second);
          count += writer.write(third);
          return count;    }

     void readBinary( BinaryReader& reader )
     {    reader.read( &first );
          reader.read( &second );
          reader.read( &third );    }

     static size_t externalSizeof() { return 0; }
};

template<class T1, class T2, class T3>
struct Serializability< triple<T1,T2,T3> > : public SelfSerializable {};

template <class _T1, class _T2, class _T3>
inline bool operator==(const triple<_T1, _T2, _T3>& __x, 
     const triple<_T1, _T2, _T3>& __y)
{ 
  return __x.first == __y.first && __x.second == __y.second
       && __x.third == __y.third;
}

template <class _T1, class _T2, class _T3>
inline bool operator<(const triple<_T1, _T2, _T3>& __x, 
     const triple<_T1, _T2, _T3>& __y)
{ 
     if ( __x.first < __y.first ) return true;
     if ( __x.first > __y.first ) return false;
     if ( __x.second < __y.second ) return true;
     if ( __x.second > __y.second ) return false;
     if ( __x.third < __y.third ) return true;
     return false;

}

template <class _T1, class _T2, class _T3>
inline bool operator!=(const triple<_T1, _T2, _T3>& __x, 
     const triple<_T1, _T2, _T3>& __y) {
  return !(__x == __y);
}

template <class _T1, class _T2, class _T3>
inline bool operator>(const triple<_T1, _T2, _T3>& __x, 
     const triple<_T1, _T2, _T3>& __y) {
  return __y < __x;
}

template <class _T1, class _T2, class _T3>
inline bool operator<=(const triple<_T1, _T2, _T3>& __x, 
     const triple<_T1, _T2, _T3>& __y) {
  return !(__y < __x);
}

template <class _T1, class _T2, class _T3>
inline bool operator>=(const triple<_T1, _T2, _T3>& __x, 
     const triple<_T1, _T2, _T3>& __y) {
  return !(__x < __y);
}

template <class _T1, class _T2, class _T3>
inline triple<_T1, _T2, _T3> make_triple(const _T1& __x, const _T2& __y, const _T3& __z)
{
  return triple<_T1, _T2, _T3>(__x, __y, __z);
}

// BinPosition1.  Using .first, return the position of an element in a sorted 
// vector, else -1.  If the element appears more than once, the position of 
// one of its instances is returned.

template<class T1, class T2, class U>
int64_t BinPosition1( const vector< pair<T1,T2> >& v, const U& x1 )
{    if ( v.size( ) == 0 ) return -1;
     T1 const& x(x1);
     size_t first = 0, last = v.size( ) - 1, next;
     while(1)
     {    if (first == last) return ( !(x < v[last].first) && !(v[last].first < x) ) ? last : -1;
          next = first + (last - first) / 2;
          if ( x < v[next].first ) last = next;
          else if ( v[next].first < x ) first = next + 1;
          else return next;    }    }

template<class T1, class T2, class T3, class U>
int64_t BinPosition1( const vector< triple<T1,T2,T3> >& v, const U& x1 )
{    if ( v.size( ) == 0 ) return -1;
     T1 const& x(x1);
     size_t first = 0, last = v.size( ) - 1, next;
     while(1)
     {    if (first == last) return ( !(x < v[last].first) && !(v[last].first < x) ) ? last : -1;
          next = first + (last - first) / 2;
          if ( x < v[next].first ) last = next;
          else if ( v[next].first < x ) first = next + 1;
          else return next;    }    }

template <class _T1, class _T2, class _T3>
  inline ostream& operator<< ( ostream& out, const triple<_T1, _T2, _T3>& __x ) {
  out << "(" << __x.first << ", " << __x.second << ", " << __x.third << ")";
  return out;
}

template < size_t N >
bool operator< ( const bitset< N >& b1, const bitset< N >& b2 ) {
  for ( size_t i = 0; i < N; i++ )
    if ( !b1.test( i )  &&  b2.test( i ) )
      return true;
  return false;
}

template < size_t N >
bool operator<= ( const bitset< N >& b1, const bitset< N >& b2 ) {
  return b1 < b2  ||  b1 == b2;
}

template < size_t N >
bool operator> ( const bitset< N >& b1, const bitset< N >& b2 ) {
  return !( b1 <= b2 );
}

template < size_t N >
bool operator>= ( const bitset< N >& b1, const bitset< N >& b2 ) {
  return !( b1 < b2 );
}

template <class T>
class LtBySize : public std::binary_function<T,T,bool>
{    public:
     bool operator()( T const& t1, T const& t2 )
     { return t1.size() < t2.size(); }
};

#endif
