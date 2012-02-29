///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/* VecUtilities.h
 *
 * Collection of vec class utilities and helpers. 
 *
 * A grab bag of functions that deal with vecs, moved from Vec.h
 * PLEASE help keep Vec.h tidy by placing any new functions here.
 *
 * See also: Vec.h
 *
 * Contains:
 *
 * MkVec
 * JoinVec
 * UniqueSortAndCount
 * SortSync
 * ReverseSortSync
 * UniqueSortSync
 * SortIndex
 * Intersection
 * Meet
 * WhatPermutation
 * PermuteVec
 * indirect_compare
 */

#ifndef VEC_UTILITIES_H
#define VEC_UTILITIES_H

#include "Vec.h"
#include <cstddef>


/////////////////////////////////////////////////////////////////////////////
//
//  MkVec - creates a new vector from a list of entries
//
/////////////////////////////////////////////////////////////////////////////

template<class T> vec<T> MkVec( const T& v1 ) {
  vec<T> v;
  v.push_back( v1 );
  return v;
}

template<class T> vec<T> MkVec( const T& v1, const T& v2 ) {
  vec<T> v;
  v.push_back( v1, v2 );
  return v;
}

template<class T> vec<T> MkVec( const T& v1, const T& v2, const T& v3 ) {
  vec<T> v;
  v.push_back( v1, v2, v3 );
  return v;
}

template<class T> vec<T> MkVec( const T& v1, const T& v2, const T& v3,
				const T& v4 ) {
  vec<T> v;
  v.push_back( v1, v2, v3, v4 );
  return v;
}

template<class T> vec<T> MkVec( const T& v1, const T& v2, const T& v3,
				const T& v4, const T& v5 ) {
  vec<T> v;
  v.push_back( v1, v2, v3, v4, v5 );
  return v;
}

template<class T> vec<T> MkVec( const T& v1, const T& v2, const T& v3,
				const T& v4, const T& v5, const T& v6 ) {
  vec<T> v;
  v.push_back( v1, v2, v3, v4, v5, v6 );
  return v;
}

template<class T> vec<T> MkVec( const T& v1, const T& v2, const T& v3,
				const T& v4, const T& v5, const T& v6,
				const T& v7 ) {
  vec<T> v;
  v.push_back( v1, v2, v3, v4, v5, v6, v7 );
  return v;
}

template<class T> vec<T> MkVec( const T& v1, const T& v2, const T& v3,
				const T& v4, const T& v5, const T& v6,
				const T& v7, const T& v8 ) {
  vec<T> v;
  v.push_back( v1, v2, v3, v4, v5, v6, v7, v8 );
  return v;
}

template<class T> vec<T> MkVec( const T& v1, const T& v2, const T& v3,
				const T& v4, const T& v5, const T& v6,
				const T& v7, const T& v8, const T& v9 ) {
  vec<T> v;
  v.push_back( v1, v2, v3, v4, v5, v6, v7, v8, v9 );
  return v;
}


/////////////////////////////////////////////////////////////////////////////
//
//  JoinVec - creates a new vector from a list of vectors
//
/////////////////////////////////////////////////////////////////////////////

template<class T> vec<T> JoinVecs( const vec<T>& v1, const vec<T>& v2 ) {
  vec<T> v( v1 );
  v.append( v2 );
  return v;
}

template<class T> vec<T> JoinVecs( const vec<T>& v1, const vec<T>& v2, const vec<T>& v3 ) {
  vec<T> v( v1 );
  v.append( v2 );
  v.append( v3 );
  return v;
}

template<class T> vec<T> JoinVecs( const vec<T>& v1, const vec<T>& v2, const vec<T>& v3,
				   const vec<T>& v4 ) {
  vec<T> v( v1 );
  v.append( v2 );
  v.append( v3 );
  v.append( v4 );
  return v;
}

template<class T> vec<T> JoinVecs( const vec<T>& v1, const vec<T>& v2, const vec<T>& v3,
				   const vec<T>& v4, const vec<T>& v5 ) {
  vec<T> v( v1 );
  v.append( v2 );
  v.append( v3 );
  v.append( v4 );
  v.append( v5 );
  return v;
}

template<class T> vec<T> JoinVecs( const vec<T>& v1, const vec<T>& v2, const vec<T>& v3,
				   const vec<T>& v4, const vec<T>& v5, const vec<T>& v6 ) {
  vec<T> v( v1 );
  v.append( v2 );
  v.append( v3 );
  v.append( v4 );
  v.append( v5 );
  v.append( v6 );
  return v;
}


/////////////////////////////////////////////////////////////////////////////
//
//  UniqueSortAndCount
//
/////////////////////////////////////////////////////////////////////////////

/// After this method is applied to a vector \c v (first argument), this vector
/// will contain only unique elements (as defined by EqualPredicate)
/// found in its original content. These elements will be sorted in ascending order
/// (according to the StrictWeakOrdering). The number of times each 
/// returned unique element v[i] occured in the original content of the vector v will
/// be returned in counts[i] (old content of counts vector, if any, will be destroyed). 

template <class T, class StrictWeakOrdering, class EqualPredicate> 
void UniqueSortAndCount(vec<T> & v, StrictWeakOrdering comp, EqualPredicate equal, 
			vec<unsigned int> & counts) {
  counts.clear();
  if ( v.empty() ) return;

  TRACEVAL_STOP_TRACING_COPIES;
  
  Sort( v, comp );
  
  typename vec<T>::iterator iter = v.begin();
  typename vec<T>::iterator iter_end = v.end();
  typename vec<T>::iterator last_unique_elem_iter = iter; 
  typename vec<T>::iterator move_to = v.begin();

  // we know that v is non-empty, so it's ok to increase iter
  // right away in the next line:
  for ( ++iter; iter != iter_end ; ++iter ) {
    // we follow the standard assumptions of STL here.
    // the latters assume that operator== is defined
    // (but not necessarily !=):
    if ( equal(*iter,*last_unique_elem_iter ) ) {
      continue;
    }
    counts.push_back( distance(last_unique_elem_iter, iter) );
    *(++move_to) = *(last_unique_elem_iter = iter);
  }
  counts.push_back( distance(last_unique_elem_iter, iter) );
  v.erase(++move_to,iter_end);

  TRACEVAL_START_TRACING_COPIES;
}


/// When this method is applied to a vector \c v (first argument), it will be 
/// first shrunk to contain only unique elements (as defined by EqualPredicate)
/// found in its original content. These elements will be also sorted in 
/// ascending order (according to the StrictWeakOrdering) and the number of 
/// times each unique element occured in the original content of the vector v 
/// will be counted. The FilterPredicate will be then applied to each such 
/// count and upon return from this method <v> and <counts> will contain only 
/// unique elements (still sorted) and counts, respectively, such that for each 
/// counts[i] the value of passes(counts[i]) is <true>

template <class T, class StrictWeakOrdering, class EqualPredicate, class FilterPredicate> 
void UniqueSortAndCount(vec<T> & v, StrictWeakOrdering comp, EqualPredicate equal, 
			      vec<unsigned int> & counts, FilterPredicate pass) {
  counts.clear();
  if ( v.empty() ) return;

  TRACEVAL_STOP_TRACING_COPIES;
  
  Sort( v, comp );
  
  typename vec<T>::iterator iter = v.begin();
  typename vec<T>::iterator iter_end = v.end();
  typename vec<T>::iterator last_unique_elem_iter = iter; 
  typename vec<T>::iterator move_to = v.begin();
  --move_to; // point before the first element - in this
             // method we may not need to keep the first element at all!

  // we know that v is non-empty, so it's ok to increase iter
  // right away in the next line:
  for ( ++iter; iter != iter_end ; ++iter ) {
    // we follow the standard assumptions of STL here.
    // the latters assume that operator== is defined
    // (but not necessarily !=):
    if ( equal(*iter,*last_unique_elem_iter ) ) {
      continue; // we are inside a stretch of same element repetions, 
                // keep counting...
    }
    unsigned int cnt = distance(last_unique_elem_iter, iter); 
    if ( pass(cnt) ) {
      counts.push_back(cnt); // record count and store element only if passes()
      *(++move_to) = *(last_unique_elem_iter) ;
    }
    last_unique_elem_iter = iter; // we just encountered new unique element,
                                  // record its position (but not store yet!)
  }
  unsigned int cnt = distance(last_unique_elem_iter, iter_end); 
  if ( pass(cnt) ) { // take care of the last unique elem - not stored yet!
      counts.push_back( cnt );
      *(++move_to) = *(last_unique_elem_iter) ;
  }      
  v.erase(++move_to,iter_end);

  TRACEVAL_START_TRACING_COPIES;
}

/// After this method is applied to a vector \c v (first argument), this vector
/// will contain only unique elements
/// found in its original content. These elements will be sorted in ascending order.
/// The number of times each returned unique element v[i] occured
/// in the original content of the vector v will
/// be returned in counts[i] (old content of counts vector, if any, will be destroyed). 

template <class T> 
inline  void UniqueSortAndCount(vec<T> & v, vec<unsigned int> & counts) {
  UniqueSortAndCount(v,less<T>(),equal_to<T>(),counts);
}

/// When this method is applied to a vector \c v (first argument), it will be 
/// first shrunk to contain only unique elements.
/// These elements will be also sorted in 
/// ascending order and the number of 
/// times each unique element occured in the original content of the vector v 
/// will be counted. The FilterPredicate will be then applied to each such 
/// count and upon return from this method <v> and <counts> will contain only 
/// unique elements (still sorted) and counts, respectively, such that for each 
/// counts[i] the value of passes(counts[i]) is <true>

template <class T, class FilterPredicate> 
inline  void UniqueSortAndCount(vec<T> & v, vec<unsigned int> & counts, FilterPredicate pass) {
  UniqueSortAndCount(v, less<T>(), equal_to<T>(), counts, pass);
}


/////////////////////////////////////////////////////////////////////////////
//
//  SortSync - sort many vectors bases on the order of the first
//
/////////////////////////////////////////////////////////////////////////////

// forward declare for use in SortSync
template<class V, typename C, typename IDX> 
void WhatPermutation( const V & v, vec<IDX>& perm, C comparator, bool inv = true );

// SortSync( vec& v, vec& w ): sort v, moving the elements of w synchronously.

template<class S, class T, typename F, typename IDX > 
void SortSyncDispatch( vec<S>& v, vec<T>& w, F comparator ) {
  ForceAssertEq( v.size( ), w.size( ) );
  vec<IDX> perm;
  WhatPermutation<vec<S>,F, IDX>(v, perm, comparator);
  PermuteVec(v, perm);
  PermuteVec(w, perm);
}

template<class S, class T, typename F > 
void SortSync( vec<S>& v, vec<T>& w, F comparator ) {
  if (v.size() < numeric_limits<unsigned int>::max())
    SortSyncDispatch<S,T,F,unsigned int>(v, w, comparator);
  else
    SortSyncDispatch<S,T,F,typename vec<S>::size_type>(v, w, comparator);
}

template<class S, class T> 
void SortSync( vec<S>& v, vec<T>& w ) {
  SortSync( v, w, less<S>() );
}

template<class S, class T> 
void ReverseSortSync( vec<S>& v, vec<T>& w ) {
  SortSync( v, w, greater<S>() );
}

template<class S, class T> 
void UniqueSortSync( vec<S>& v, vec<T>& w ) {
  SortSync( v, w, less<S>() );
  typename vec<S>::size_type count = 0;
  for ( typename vec<S>::size_type i = 0; i < v.size( ); i++ )
  {    if ( i > 0 && v[i] == v[i-1] ) continue;
       if ( count != i )
       {    v[count] = v[i];
            w[count] = w[i];    }
       ++count;    }
  v.resize(count), w.resize(count);
}

// SortSync( v, w, x )-type functions.

template<class S, class T, class U, typename F, typename IDX > 
void SortSyncDispatch( vec<S>& v, vec<T>& w, vec<U>& x, F comparator ) {
  ForceAssertEq( v.size( ), w.size( ) );
  ForceAssertEq( v.size( ), x.size( ) );
  vec<IDX> perm;
  WhatPermutation<vec<S>,F, IDX>(v, perm, comparator);
  PermuteVec(v, perm);
  PermuteVec(w, perm);
  PermuteVec(x, perm);   
}

template<class S, class T, class U, typename F > 
void SortSync( vec<S>& v, vec<T>& w, vec<U>& x, F comparator ) {
  if (v.size() < numeric_limits<unsigned int>::max())
    SortSyncDispatch<S,T,U,F,unsigned int>(v, w, x,comparator);
  else
    SortSyncDispatch<S,T,U,F,typename vec<S>::size_type>(v, w, x, comparator);
}

template<class S, class T, class U> 
void SortSync( vec<S>& v, vec<T>& w, vec<U>& x ) {
  SortSync( v, w, x, less<S>() );
}

template<class S, class T, class U> 
void ReverseSortSync( vec<S>& v, vec<T>& w, vec<U>& x ) {
  SortSync( v, w, x, greater<S>() );
}

// SortSync( v, w, x, y )-type functions.

template<class S, class T, class U, class V, typename F, typename IDX > 
void SortSyncDispatch( vec<S>& v, vec<T>& w, vec<U>& x, vec<V>& y, F comparator ) {
  ForceAssertEq( v.size( ), w.size( ) );
  ForceAssertEq( v.size( ), x.size( ) );
  ForceAssertEq( v.size( ), y.size( ) );
  vec<IDX> perm;
  WhatPermutation<vec<S>,F, IDX>(v, perm, comparator);
  PermuteVec(v, perm);
  PermuteVec(w, perm);
  PermuteVec(x, perm);
  PermuteVec(y, perm);
}

template<class S, class T, class U, class V, typename F > 
void SortSync( vec<S>& v, vec<T>& w, vec<U>& x, vec<V>& y, F comparator ) {
  if (v.size() <= numeric_limits<unsigned int>::max())
    SortSyncDispatch<S,T,U,V,F,unsigned int>(v, w, x, y, comparator);
  else
    SortSyncDispatch<S,T,U,V,F,typename vec<S>::size_type>(v, w, x, y, comparator);
}

template<class S, class T, class U, class V> 
void SortSync( vec<S>& v, vec<T>& w, vec<U>& x, vec<V>& y ) {
  SortSync( v, w, x, y, less<S>() );
}

template<class S, class T, class U, class V> 
void ReverseSortSync( vec<S>& v, vec<T>& w, vec<U>& x, vec<V>& y ) {
  SortSync( v, w, x, y, greater<S>() );    
}

/////////////////////////////////////////////////////////////////////////////
//
//  SortIndex
//
/////////////////////////////////////////////////////////////////////////////

// forward declare for use in SortIndex
template<class V, class IDX>
void PermuteVec(V & v, const vec<IDX> & permutation);


// SortIndex( v, index ): calculates an index into v such that
// v[index[n]] returns the nth value of a sorted v.
template<class S, typename F, typename IDX > 
void SortIndex( const vec<S>& v, vec<IDX>& index, F comparator ) {
  AssertLe(v.size(), static_cast<size_t>(numeric_limits<IDX>::max()));
  index.resize(v.size());
  iota(index.begin(), index.end(), 0);
  vec<IDX> perm;
  WhatPermutation<vec<S>,F, IDX>(v, perm, comparator);
  PermuteVec(index, perm);
}

template<class S, typename IDX> 
void SortIndex( const vec<S>& v, vec<IDX>& index ) {
  SortIndex( v, index, less<S>() );
}



/////////////////////////////////////////////////////////////////////////////
//
//  Intersection & Meet - functions to find shared elements
//
/////////////////////////////////////////////////////////////////////////////

/// Intersection: make pairs of all elements in the first vector that also
/// occur in the second (which can result in duplicates).  Put the pairs
/// in the result. Both vectors must be SORTED.  This version is meant
/// for situations where two objects may compare equal but contain different
/// information, and we want to preserve both copies of that information.

template<class T> void
Intersection( const vec<T>& v1, const vec<T>& v2, vec<pair<T,T> > & result )
{    
  result.clear();
  typename vec<T>::const_iterator v1iter = v1.begin();
  typename vec<T>::const_iterator v2iter = v2.begin();
  while ( v1iter != v1.end() && v2iter != v2.end() )
  {
    if ( *v1iter < *v2iter )
      ++v1iter;
    else if ( *v2iter < *v1iter )
      ++v2iter;
    else
      result.push_back( make_pair(*v1iter++,*v2iter) );
  }
}

/// Intersection: copy all the elements in the first vector that also
/// occur in the second (which can result in duplicates).  Both
/// vectors must be SORTED.  There are two versions: the first accepts
/// as an argument a container for the result, the second simply
/// returns the result.  The elements in the result are taken from the
/// first vector.

template<class T, class R> void
Intersection( const vec<T>& v1, const vec<T>& v2, vec<R>& result )
{    
  result.clear();
  typename vec<T>::const_iterator v1iter = v1.begin();
  typename vec<T>::const_iterator v2iter = v2.begin();
  while ( v1iter != v1.end() && v2iter != v2.end() )
  {
    if ( *v1iter < *v2iter )
      ++v1iter;
    else if ( *v2iter < *v1iter )
      ++v2iter;
    else
      result.push_back( *v1iter++ );
  }
}

/// Intersection: copy all the elements in the first vector that also
/// occur in the second (which can result in duplicates).  Both
/// vectors must be SORTED.  Returns the result.  
/// The elements in the result are taken from the
/// first vector.

template<class T> vec<T> Intersection( const vec<T>& v1, const vec<T>& v2 )
{    
  vec<T> w;
  Intersection( v1, v2, w ); 
  return w;    
}

// Find the intersection of a family of sorted vectors.

template<class T> void Intersection( const vec< vec<T> >& x, vec<T>& y )
{    ForceAssert( x.nonempty( ) );
     y = x[0];
     for ( typename vec<T>::size_type j = 1; j < x.size( ); j++ )
     {    vec<T> z;
          Intersection( x[j], y, z );
          y = z;    }    }

/// Meet: determine if two SORTED vectors have an element in common.  

template<class T> Bool Meet( const vec<T>& v1, const vec<T>& v2 )
{
  typename vec<T>::const_iterator v1iter = v1.begin();
  typename vec<T>::const_iterator v2iter = v2.begin();
  while ( v1iter != v1.end() && v2iter != v2.end() )
  {
    if ( *v1iter < *v2iter )
      ++v1iter;
    else if ( *v2iter < *v1iter )
      ++v2iter;
    else
      return true;
  }
  return false;
}


/////////////////////////////////////////////////////////////////////////////
//
//  Helpers for SortSync and SortIndex
//
/////////////////////////////////////////////////////////////////////////////

/**Permute input vector v in place according to permutation.
   Preconditions:
   - v.size() == permutation.size()
   If the permutation contains a -1, the position corresponding to that
   is essentially ignored and ends up in one of the available empty spaces

   This works for all std::vectors and vecs.
*/
template<class V, typename IDX>
void PermuteVec(V & v, const vec<IDX> & permutation) {
  AssertEq(v.size(), permutation.size());
  AssertLe(v.size(),static_cast<size_t>(numeric_limits<IDX>::max()));  // make sure the index type is large enough
  IDX n = static_cast<IDX>(v.size());
  vec<IDX> o = permutation;
  const bool unsigned_index = (numeric_limits<IDX>::is_signed == false);
  const IDX minus1 = static_cast<IDX>(-1);  // -1 for a signed index, don't care for unsigned
  for (IDX i = 0; i != n; ++i) {
    while (o[i] != i && (unsigned_index || o[i] != minus1)) {  // ignore -1 check for unsigned
      std::swap(v[i], v[o[i]]);
      std::swap(o[i], o[o[i]]);
    }
  }
}

// Instantiate with a vector v, then use in place of operator<.
// It will tell you i<j if v[i]<v[j].  Helper for WhatPermutation.
template<class S, typename C> 
struct indirect_compare : public binary_function<ulonglong, ulonglong, bool> {
  const vec<S>& v;
  C comp;
  indirect_compare(const vec<S>& v ) : v(v) {}
  indirect_compare(const vec<S>& v, C comparator) : v(v), comp(comparator){}
  bool operator() ( ulonglong i, ulonglong j ) const { return comp( v[i], v[j] ); }
};

// What permutation would we apply to V to get it sorted?
// NOTE: If v is a vec<int>, this inverts it.
template<class V, typename C, typename IDX >
void WhatPermutation( const V & v, vec<IDX>& permutation, C comparator, 
		      bool inv ) {
  AssertLe(v.size(), static_cast<size_t>(numeric_limits<IDX>::max()));
  IDX n = v.size();
  vec<IDX> perm(n);
  for (IDX i = 0; i < n; i++)
    perm[i]=i;
  sort( perm.begin(), perm.end(), 
	indirect_compare<typename V::value_type,C>(v,comparator) );
  if (inv) {
    // That's the inverse of the permutation PermuteVec takes.
    permutation.resize(n);
    for (IDX i = 0; i < n; i++)
      permutation[perm[i]]=i;
  }
  else {
    swap(perm, permutation);
  }
}

template<class V, typename IDX>
void WhatPermutation( const V & v, vec<IDX>& permutation, bool inv = true ) {
  WhatPermutation< V, less<typename V::value_type> >
    (v,permutation, less<typename V::value_type>(), inv );
}

#endif
