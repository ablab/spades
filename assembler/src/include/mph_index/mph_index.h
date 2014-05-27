//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef __CXXMPH_MPH_INDEX_H__
#define __CXXMPH_MPH_INDEX_H__

// Minimal perfect hash abstraction implementing the BDZ algorithm
//
// This is a data structure that given a set of known keys S, will create a
// mapping from S to [0..|S|). The class is informed about S through the 
// method and the mapping is queried by calling index(key).
//
// This is a pretty uncommon data structure, and if you application has a real
// use case for it, chances are that it is a real win. If all you are doing is
// a straightforward implementation of an in-memory associative mapping data
// structure, then it will probably be slower. Take a look at mph_map.h
// instead.
//
// Thesis presenting this and similar algorithms:
// http://homepages.dcc.ufmg.br/~fbotelho/en/talks/thesis2008/thesis.pdf
//
// Notes:
//
// Most users can use the SimpleMPHIndex wrapper instead of the MPHIndex which
// have confusing template parameters.
// This class only implements a minimal perfect hash function, it does not
// implement an associative mapping data structure.

#include <stdint.h>

#include <cassert>
#include <climits>
#include <cmath>
#include <unordered_map>  // for std::hash
#include <vector>

#include <iostream>

using std::cerr;
using std::endl;

#include "seeded_hash.h"
#include "mph_bits.h"
#include "trigraph.h"

namespace cxxmph {

class MPHIndex {
 public:
  MPHIndex(bool square = false, double c = 1.23, uint8_t b = 7) :
      c_(c), b_(b), m_(0), n_(0), k_(0), square_(square), r_(1),
      g_(8, true), ranktable_(0)  {
    nest_displacement_[0] = 0;
    nest_displacement_[1] = r_;
    nest_displacement_[2] = (r_ << 1);
  }

  template <class SeededHashFcn, class ForwardIterator>
  bool Reset(ForwardIterator begin, ForwardIterator end, uint32_t size);
  template <class SeededHashFcn, class Key>  // must agree with Reset
  // Get a unique identifier for k, in the range [0;size()). If x wasn't part
  // of the input in the last Reset call, returns a random value.
  uint32_t index(const Key& x) const;
  size_t size() const { return m_; }
  void clear();

  // Advanced users functions. Please avoid unless you know what you are doing.
  uint32_t perfect_hash_size() const { return n_; } 
  template <class SeededHashFcn, class Key>  // must agree with Reset
  uint32_t perfect_hash(const Key& x) const;  // way faster than the minimal
  template <class SeededHashFcn, class Key>  // must agree with Reset
  uint32_t perfect_square(const Key& x) const;  // even faster but needs square=true
  uint32_t minimal_perfect_hash_size() const { return (uint32_t) size(); }
  template <class SeededHashFcn, class Key>  // must agree with Reset
  uint32_t minimal_perfect_hash(const Key& x) const;

  template<class Writer>
  void serialize(Writer &os) const;
  template<class Reader>
  void deserialize(Reader &is);
  
 private:
  template <class SeededHashFcn, class ForwardIterator>
  bool Mapping(ForwardIterator begin, ForwardIterator end,
               std::vector<TriGraph::Edge>* edges,
               std::vector<uint32_t>* queue);
  bool GenerateQueue(TriGraph* graph, std::vector<uint32_t>* queue);
  void Assigning(const std::vector<TriGraph::Edge>& edges,
                 const std::vector<uint32_t>& queue);
  void Ranking();
  uint32_t Rank(uint32_t vertex) const;

  // Algorithm parameters
  // Perfect hash function density. If this was a 2graph,
  // then probability of having an acyclic graph would be 
  // sqrt(1-(2/c)^2). See section 3 for details.
  // http://www.it-c.dk/people/pagh/papers/simpleperf.pdf
  double c_;
  uint8_t b_;  // Number of bits of the kth index in the ranktable

  // Values used during generation
  uint32_t m_;  // edges count
  uint32_t n_;  // vertex count
  uint32_t k_;  // kth index in ranktable, $k = log_2(n=3r)\varepsilon$
  bool square_;  // make bit vector size a power of 2

  // Values used during search

  // Partition vertex count, derived from c parameter.
  uint32_t r_;
  uint32_t nest_displacement_[3];  // derived from r_

  // The array containing the minimal perfect hash function graph.
  dynamic_2bitset g_;
  uint8_t threebit_mod3[10];  // speed up mod3 calculation for 3bit ints
  // The table used for the rank step of the minimal perfect hash function
  std::vector<uint32_t> ranktable_;

  // The selected hash seed triplet for finding the edges in the minimal
  // perfect hash function graph.
  uint32_t hash_seed_[3];

public:
  size_t mem_size() const { return ranktable_.capacity()*sizeof(uint32_t) + g_.size() / 4; }
};

// Template method needs to go in the header file.
template <class SeededHashFcn, class ForwardIterator>
bool MPHIndex::Reset(ForwardIterator begin, ForwardIterator end, uint32_t size) {
  if (end == begin) {
    clear();
    return true;
  }

  m_ = size;
  r_ = static_cast<uint32_t>(ceil((c_*m_)/3));
  r_ = std::max(3u, r_);
  if ((r_ % 2) == 0) r_ += 1;

  // This can be used to speed mods, but increases occupation too much. 
  // Needs to try http://gmplib.org/manual/Integer-Exponentiation.html instead
  if (square_)
    r_ = nextpoweroftwo(r_);

  nest_displacement_[0] = 0;
  nest_displacement_[1] = r_;
  nest_displacement_[2] = (r_ << 1);
  for (uint32_t i = 0; i < sizeof(threebit_mod3); ++i)
    threebit_mod3[i] = (uint8_t)(i % 3);

  n_ = 3*r_;
  k_ = 1U << b_;

  // cerr << "m " << m_ << " n " << n_ << " r " << r_ << endl;

  int iterations;
  std::vector<TriGraph::Edge> edges;
  std::vector<uint32_t> queue;
  uint32_t seed = 0;
  for (iterations = 0; iterations < 1000; ++iterations) {
    // cerr << "Iterations missing: " << iterations << endl;
    for (unsigned i = 0; i < 3; ++i)
      hash_seed_[i] = seed++;
    if (Mapping<SeededHashFcn>(begin, end, &edges, &queue))
      break;
  }
  if (iterations == 1000)
    return false;

  Assigning(edges, queue);
  std::vector<TriGraph::Edge>().swap(edges);
  std::vector<uint32_t>().swap(queue);
  Ranking();
  return true;
}

template <class SeededHashFcn, class ForwardIterator>
bool MPHIndex::Mapping(ForwardIterator begin, ForwardIterator end,
                       std::vector<TriGraph::Edge>* edges, std::vector<uint32_t>* queue) {
  TriGraph graph(n_, m_);
  for (ForwardIterator it = begin; it != end; ++it) { 
    h128 h = SeededHashFcn().hash128(*it, hash_seed_[0]);
    // for (int i = 0; i < 3; ++i) h[i] = SeededHashFcn()(*it, hash_seed_[i]);
    uint32_t v0 = h[0] % r_;
    uint32_t v1 = h[1] % r_ + r_;
    uint32_t v2 = h[2] % r_ + (r_ << 1);
    
    // cerr << "Key: " << *it << " edge " <<  it - begin << " (" << v0 << "," << v1 << "," << v2 << ")" << endl;
    graph.AddEdge(TriGraph::Edge(v0, v1, v2));
  }

  if (GenerateQueue(&graph, queue)) {
    graph.ExtractEdgesAndClear(edges);
    return true;
  }

  return false;
}

template <class SeededHashFcn, class Key>
uint32_t MPHIndex::perfect_square(const Key& key) const {
  h128 h = SeededHashFcn().hash128(key, hash_seed_[0]);
  h[0] = (h[0] & (r_-1)) + nest_displacement_[0];
  h[1] = (h[1] & (r_-1)) + nest_displacement_[1];
  h[2] = (h[2] & (r_-1)) + nest_displacement_[2];
  assert((h[0]) < g_.size());
  assert((h[1]) < g_.size());
  assert((h[2]) < g_.size());
  uint8_t nest = threebit_mod3[g_[h[0]] + g_[h[1]] + g_[h[2]]];
  uint32_t vertex = h[nest];
  return vertex;
}

template <class SeededHashFcn, class Key>
uint32_t MPHIndex::perfect_hash(const Key& key) const {
  if (!g_.size())
    return 0;

  h128 h = SeededHashFcn().hash128(key, hash_seed_[0]);
  h[0] = (h[0] % r_) + nest_displacement_[0];
  h[1] = (h[1] % r_) + nest_displacement_[1];
  h[2] = (h[2] % r_) + nest_displacement_[2];
  assert((h[0]) < g_.size());
  assert((h[1]) < g_.size());
  assert((h[2]) < g_.size());
  uint8_t nest = threebit_mod3[g_[h[0]] + g_[h[1]] + g_[h[2]]];
  uint32_t vertex = h[nest];
  return vertex;
}

template <class SeededHashFcn, class Key>
uint32_t MPHIndex::minimal_perfect_hash(const Key& key) const {
  return Rank(perfect_hash<SeededHashFcn, Key>(key));
}

template <class SeededHashFcn, class Key>
uint32_t MPHIndex::index(const Key& key) const {
  return minimal_perfect_hash<SeededHashFcn, Key>(key);
}

// Simple wrapper around MPHIndex to simplify calling code. Please refer to the
// MPHIndex class for documentation.
template <class DefaultKey, class HashFcn = typename seeded_hash<std::hash<DefaultKey>>::hash_function>
class SimpleMPHIndex : public MPHIndex {
 public:
  SimpleMPHIndex(bool advanced_usage = false) : MPHIndex(advanced_usage) {}
  template <class ForwardIterator>
  bool Reset(ForwardIterator begin, ForwardIterator end, uint32_t size) {
    return MPHIndex::Reset<HashFcn>(begin, end, size);
  }
  template<class Key = DefaultKey>
  uint32_t index(const Key& key) const { return MPHIndex::index<HashFcn>(key); }
};

// The parameters minimal and square trade memory usage for evaluation speed.
// Minimal decreases speed and memory usage, and square does the opposite.
// Using minimal=true and square=false is the same as SimpleMPHIndex.
template <bool minimal, bool square, class Key, class HashFcn>
struct FlexibleMPHIndex {};

template <class Key, class HashFcn>
struct FlexibleMPHIndex<true, false, Key, HashFcn> 
    : public SimpleMPHIndex<Key, HashFcn> {
  FlexibleMPHIndex() : SimpleMPHIndex<Key, HashFcn>(false) {}
  uint32_t index(const Key& key) const {
      return MPHIndex::minimal_perfect_hash<HashFcn>(key); }
  uint32_t size() const { return MPHIndex::minimal_perfect_hash_size(); }
};
template <class Key, class HashFcn>
struct FlexibleMPHIndex<false, true, Key, HashFcn> 
    : public SimpleMPHIndex<Key, HashFcn> {
  FlexibleMPHIndex() : SimpleMPHIndex<Key, HashFcn>(true) {}
  uint32_t index(const Key& key) const {
      return MPHIndex::perfect_square<HashFcn>(key); }
  size_t size() const { return MPHIndex::perfect_hash_size(); }
};
template <class Key, class HashFcn>
struct FlexibleMPHIndex<false, false, Key, HashFcn> 
    : public SimpleMPHIndex<Key, HashFcn> {
  FlexibleMPHIndex() : SimpleMPHIndex<Key, HashFcn>(false) {}
  uint32_t index(const Key& key) const {
      return MPHIndex::perfect_hash<HashFcn>(key); }
  size_t size() const { return MPHIndex::perfect_hash_size(); }
};
// From a trade-off perspective this case does not make much sense.
// template <class Key, class HashFcn>
// class FlexibleMPHIndex<true, true, Key, HashFcn>

#define WRITE_POD(var, os)                      \
  do {                                          \
    (os).write((char*)&(var), sizeof(var));     \
  } while(0)

#define WRITE_CLASS(var, os) \
  do {                       \
    (var).serialize(os);     \
  } while(0)

#define WRITE_VECTOR(var, os)                                         \
  do {                                                                \
    const __typeof(var) & __vect = (var);                             \
    size_t __sz = __vect.size();                                      \
    WRITE_POD(__sz, os);                                              \
    (os).write((char*)&__vect[0], __vect.size() * sizeof(__vect[0])); \
  } while(0)

#define READ_POD(var, is)                       \
  do {                                          \
    (is).read((char*)&(var), sizeof(var));      \
  } while(0)

#define READ_CLASS(var, is)    \
  do {                         \
    (var).deserialize(is);     \
  } while(0)

#define READ_VECTOR(var, is)                                          \
  do {                                                                \
    __typeof(var) & __vect = (var);                                   \
    size_t __sz = 0;                                                  \
    READ_POD(__sz, is);                                               \
    __vect.resize(__sz);                                              \
    (is).read((char*)&__vect[0], __sz * sizeof(__vect[0]));           \
  } while(0)

template<class Writer>
void MPHIndex::serialize(Writer &os) const {
  WRITE_POD(c_, os);
  WRITE_POD(b_, os);
  WRITE_POD(m_, os);
  WRITE_POD(n_, os);
  WRITE_POD(k_, os);
  WRITE_POD(square_, os);
  WRITE_POD(r_, os);
  WRITE_POD(nest_displacement_, os);
  WRITE_POD(threebit_mod3, os);
  WRITE_POD(hash_seed_, os);
  WRITE_VECTOR(ranktable_, os);
  WRITE_CLASS(g_, os);
}

template<class Reader>
void MPHIndex::deserialize(Reader &is) {
  READ_POD(c_, is);
  READ_POD(b_, is);
  READ_POD(m_, is);
  READ_POD(n_, is);
  READ_POD(k_, is);
  READ_POD(square_, is);
  READ_POD(r_, is);
  READ_POD(nest_displacement_, is);
  READ_POD(threebit_mod3, is);
  READ_POD(hash_seed_, is);
  READ_VECTOR(ranktable_, is);
  READ_CLASS(g_, is);
}

#undef WRITE_POD
#undef WRITE_CLASS
#undef WRITE_VECTOR
#undef READ_POD
#undef READ_CLASS
#undef READ_VECTOR


}  // namespace cxxmph

#endif // __CXXMPH_MPH_INDEX_H__
