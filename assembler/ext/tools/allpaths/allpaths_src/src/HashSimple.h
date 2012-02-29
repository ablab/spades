/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef HASH_SIMPLE_H
#define HASH_SIMPLE_H

#include "Vec.h"
#include <ext/hash_set>

using __gnu_cxx::hash;

/// Very simple linear probing hash table for positive ints only.
///
/// \class HashSimple
///
/// Primitive, possible future building block,
/// meant mainly to measure lower bound for speed of hashing because
/// hash_set was clearly too slow.
///
/// This class is faster than hash_set by about 2x

class HashSimple {
  typedef int T;
  static const T BAD=-1;
private:
  vec<T> data;
  T * start_;
  T * end_;
  int size_;
  hash<T> h_;
  int MAXSIZE_;

  void Rehash() {
    vec<T> newdata(data.size() * 2, BAD);
    data.swap(newdata);
    MAXSIZE_ = capacity() * 2 / 3;
    start_ = &data[0];
    end_ = start_ + data.size();
    const int S=newdata.size();
    for (vec<T>::iterator i=newdata.begin(); i != newdata.end(); ++i) {
      if (*i != BAD) Insert(*i);
    }
  }

  int * FindInternal(const T & t) const {
    int * current= start_ + h_(t) % capacity();
    while (*current != t && *current != BAD) {
      if (++current == end_) current = start_;
    }
    return current;
  }

public:
  class iterator;
  friend class iterator;

  HashSimple(int capacity=100): 
    data(capacity,BAD), start_(&data[0]), end_(start_ + capacity), 
    size_(0), h_(), MAXSIZE_(capacity*2/3) {}

  /// Return false if t was already in the set.
  bool Insert(const T & t) {
    AssertNe(t, BAD);
    int * here = FindInternal(t);
    if (*here != BAD) return false;
    
    *here = t; 
    if (++size_ > MAXSIZE_) { 
      Rehash();
    }
    return true;
  }

  void clear() { data.clear(); size_=0; }

  int size() const { return size_; }

  int capacity() const { return data.size(); }

  bool Has(const T & t) const {
    return (*FindInternal(t) != BAD);
  }
  
  /// Note that this iterator always ends at end(), so it 
  /// must be started at begin(): it does not know how to circle
  /// back.
  class iterator {
  private: 
    HashSimple & h; 
    int pos; 
    const int SIZE; 

  public:
    iterator(HashSimple & h, int p=0): h(h), pos(p), SIZE(h.data.size()) {}
    iterator operator++() { 
      //PRINT4(pos, SIZE, h.data[pos], h.capacity());
      while (++pos < SIZE && h.data[pos] == -1) {}
      return *this; 
    }
    T operator*() { return h.data[pos]; }
    bool operator==(const iterator & o) { return o.pos == pos; }
    bool operator!=(const iterator & o) { return !operator==(o); }
  };

  iterator begin() { return iterator(*this); }
  iterator end() { return iterator(*this, data.size()); }
};

//const int HashSimple::BAD;

#endif // HASH_SIMPLE_H
