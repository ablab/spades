///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// This file is a combination of part of Vec.h and VecUtilities.h, modified to 
// use gcc parallelized algorithm versions.  At present, only sort functions are
// included, and only the parallelized algorithms are the sorts per se.  Probably
// more could be parallelized.

#ifndef PARALLEL_VEC_UTILITIES_H
#define PARALLEL_VEC_UTILITIES_H

#include <omp.h>

#include <parallel/algorithm>
#include <functional>

#include "Vec.h"
#include "VecUtilities.h"

template<class T> void ParallelSort( vec<T>& v )
{
  TRACEVAL_STOP_TRACING_COPIES;
  __gnu_parallel::sort( v.begin( ), v.end( ), std::less<T>() );
  TRACEVAL_START_TRACING_COPIES;
}

template<class T, class StrictWeakOrdering > 
inline void ParallelSort( vec<T>& v, StrictWeakOrdering comp )
{
  TRACEVAL_STOP_TRACING_COPIES;
  __gnu_parallel::sort( v.begin( ), v.end( ), comp );
  TRACEVAL_START_TRACING_COPIES;
}

template<class T> void ParallelReverseSort( vec<T>& v )
{
  TRACEVAL_STOP_TRACING_COPIES;
  __gnu_parallel::sort( v.rbegin( ), v.rend( ), std::less<T>() );
  TRACEVAL_START_TRACING_COPIES;
}

template<class T, class StrictWeakOrdering> 
void ParallelReverseSort( vec<T>& v, StrictWeakOrdering comp )
{
  TRACEVAL_STOP_TRACING_COPIES;
  __gnu_parallel::sort( v.rbegin( ), v.rend( ), comp );
  TRACEVAL_START_TRACING_COPIES;
}

template <class T, class StrictWeakOrdering, class EqualPredicate> 
  void ParallelUniqueSort(vec<T> & v, StrictWeakOrdering comp, EqualPredicate equal) {
  if ( v.empty() ) return;

  TRACEVAL_STOP_TRACING_COPIES;
  
  ParallelSort( v, comp );
  
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
    if ( equal(*iter,*last_unique_elem_iter ) ) continue;
    *(++move_to) = *(last_unique_elem_iter = iter);
  }
  v.erase(++move_to,iter_end);

  TRACEVAL_START_TRACING_COPIES;
}

template <class T> inline void ParallelUniqueSort(vec<T> & v) {
  ParallelUniqueSort(v,less<T>(),equal_to<T>());
}

template <class T, class StrictWeakOrdering, class EqualPredicate> 
void ParallelUniqueSortAndCount(vec<T> & v, StrictWeakOrdering comp, 
     EqualPredicate equal, vec<unsigned int> & counts) {
  counts.clear();
  if ( v.empty() ) return;

  TRACEVAL_STOP_TRACING_COPIES;
  
  ParallelSort( v, comp );
  
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

template <class T, class StrictWeakOrdering, class EqualPredicate, class FilterPredicate> 
void ParallelUniqueSortAndCount(vec<T> & v, StrictWeakOrdering comp, 
     EqualPredicate equal, vec<unsigned int> & counts, FilterPredicate pass) {
  counts.clear();
  if ( v.empty() ) return;

  TRACEVAL_STOP_TRACING_COPIES;
  
  ParallelSort( v, comp );
  
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

template <class T> 
inline  void ParallelUniqueSortAndCount(vec<T> & v, vec<unsigned int> & counts) {
  ParallelUniqueSortAndCount(v,less<T>(),equal_to<T>(),counts);
}

template <class T, class FilterPredicate> 
inline  void ParallelUniqueSortAndCount(vec<T> & v, vec<unsigned int> & counts, 
     FilterPredicate pass) {
  ParallelUniqueSortAndCount(v, less<T>(), equal_to<T>(), counts, pass);
}

template<class V, typename C, typename IDX> void ParallelWhatPermutation( 
     const V & v, vec<IDX>& perm, C comparator, bool inv = true );

template<class S, class T, typename F, typename IDX > 
void ParallelSortSyncDispatch( vec<S>& v, vec<T>& w, F comparator ) {
  ForceAssertEq( v.size( ), w.size( ) );
  vec<IDX> perm;
  ParallelWhatPermutation<vec<S>,F, IDX>(v, perm, comparator);
  ParallelPermuteVec(v, perm);
  ParallelPermuteVec(w, perm);
}

template<class S, class T, typename F > 
void ParallelSortSync( vec<S>& v, vec<T>& w, F comparator ) {
  if (v.size() < numeric_limits<unsigned int>::max())
    ParallelSortSyncDispatch<S,T,F,unsigned int>(v, w, comparator);
  else
    ParallelSortSyncDispatch<S,T,F,typename vec<S>::size_type>(v, w, comparator);
}

template<class S, class T> void ParallelSortSync( vec<S>& v, vec<T>& w ) {
  ParallelSortSync( v, w, less<S>() );
}

template<class S, class T> void ParallelReverseSortSync( vec<S>& v, vec<T>& w ) {
  ParallelSortSync( v, w, greater<S>() );
}

template<class S, class T> void ParallelUniqueSortSync( vec<S>& v, vec<T>& w ) {
  ParallelSortSync( v, w, less<S>() );
  typename vec<S>::size_type count = 0;
  for ( typename vec<S>::size_type i = 0; i < v.size( ); i++ )
  {    if ( i > 0 && v[i] == v[i-1] ) continue;
       if ( count != i )
       {    v[count] = v[i];
            w[count] = w[i];    }
       ++count;    }
  v.resize(count), w.resize(count);
}

template<class S, class T, class U, typename F, typename IDX > 
void ParallelSortSyncDispatch( vec<S>& v, vec<T>& w, vec<U>& x, F comparator ) {
  ForceAssertEq( v.size( ), w.size( ) );
  ForceAssertEq( v.size( ), x.size( ) );
  vec<IDX> perm;
  ParallelWhatPermutation<vec<S>,F, IDX>(v, perm, comparator);
  ParallelPermuteVec(v, perm);
  ParallelPermuteVec(w, perm);
  ParallelPermuteVec(x, perm);   
}

template<class S, class T, class U, typename F > 
void ParallelSortSync( vec<S>& v, vec<T>& w, vec<U>& x, F comparator ) {
  if (v.size() < numeric_limits<unsigned int>::max())
    ParallelSortSyncDispatch<S,T,U,F,unsigned int>(v, w, x,comparator);
  else
    ParallelSortSyncDispatch<S,T,U,F,typename vec<S>::size_type>(v, w, x, comparator);
}

template<class S, class T, class U> 
void ParallelSortSync( vec<S>& v, vec<T>& w, vec<U>& x ) {
  ParallelSortSync( v, w, x, less<S>() );
}

template<class S, class T, class U> 
void ParallelReverseSortSync( vec<S>& v, vec<T>& w, vec<U>& x ) {
  ParallelSortSync( v, w, x, greater<S>() );
}

template<class S, class T, class U, class V, typename F, typename IDX > 
void ParallelSortSyncDispatch( vec<S>& v, vec<T>& w, vec<U>& x, vec<V>& y, F comparator ) {
  ForceAssertEq( v.size( ), w.size( ) );
  ForceAssertEq( v.size( ), x.size( ) );
  ForceAssertEq( v.size( ), y.size( ) );
  vec<IDX> perm;
  ParallelWhatPermutation<vec<S>,F, IDX>(v, perm, comparator);
  ParallelPermuteVec(v, perm);
  ParallelPermuteVec(w, perm);
  ParallelPermuteVec(x, perm);
  ParallelPermuteVec(y, perm);
}

template<class S, class T, class U, class V, typename F > 
void ParallelSortSync( vec<S>& v, vec<T>& w, vec<U>& x, vec<V>& y, F comparator ) {
  if (v.size() <= numeric_limits<unsigned int>::max())
    ParallelSortSyncDispatch<S,T,U,V,F,unsigned int>(v, w, x, y, comparator);
  else
    ParallelSortSyncDispatch<S,T,U,V,F,typename vec<S>::size_type>(v, w, x, y, comparator);
}

template<class S, class T, class U, class V> 
void ParallelSortSync( vec<S>& v, vec<T>& w, vec<U>& x, vec<V>& y ) {
  ParallelSortSync( v, w, x, y, less<S>() );
}

template<class S, class T, class U, class V> 
void ParallelReverseSortSync( vec<S>& v, vec<T>& w, vec<U>& x, vec<V>& y ) {
  ParallelSortSync( v, w, x, y, greater<S>() );    
}

template<class V, class IDX>
void ParallelPermuteVec(V & v, const vec<IDX> & permutation);

template<class S, typename F, typename IDX > 
void ParallelSortIndex( const vec<S>& v, vec<IDX>& index, F comparator ) {
  AssertLe(v.size(), static_cast<size_t>(numeric_limits<IDX>::max()));
  index.resize(v.size());
  iota(index.begin(), index.end(), 0);
  vec<IDX> perm;
  ParallelWhatPermutation<vec<S>,F, IDX>(v, perm, comparator);
  ParallelPermuteVec(index, perm);
}

template<class S, typename IDX> 
void ParallelSortIndex( const vec<S>& v, vec<IDX>& index ) {
  ParallelSortIndex( v, index, less<S>() );
}

template<class V, typename IDX>
void ParallelPermuteVec(V & v, const vec<IDX> & permutation) {
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

template<class V, typename C, typename IDX >
void ParallelWhatPermutation( const V & v, vec<IDX>& permutation, C comparator, 
		      bool inv ) {
  AssertLe(v.size(), static_cast<size_t>(numeric_limits<IDX>::max()));
  IDX n = v.size();
  vec<IDX> perm(n);
  for (IDX i = 0; i < n; i++)
    perm[i]=i;
  __gnu_parallel::sort( perm.begin(), perm.end(), 
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
void ParallelWhatPermutation( const V & v, vec<IDX>& permutation, bool inv = true ) {
  ParallelWhatPermutation< V, less<typename V::value_type> >
    (v,permutation, less<typename V::value_type>(), inv );
}

#endif
