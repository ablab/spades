/*
 * kmer_mapper.hpp
 *
 *  Created on: Dec 4, 2013
 *      Author: andrey
 */
//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include "omni/omni_utils.hpp"
#include "sequence/sequence_tools.hpp"
#include "omni/path_processor.hpp"

#include "runtime_k.hpp"
#include "edge_index.hpp"

#include <cstdlib>

namespace debruijn_graph {
template <class Graph, class Seq = runtime_k::RtSeq>
class KmerMapper : public omnigraph::GraphActionHandler<Graph> {
  typedef omnigraph::GraphActionHandler<Graph> base;
  typedef typename Graph::EdgeId EdgeId;
  typedef Seq Kmer;
  typedef typename runtime_k::KmerMap<Kmer, Seq> MapType;

  MapType mapping_;
  size_t k_;
  bool verification_on_;

 public:

  KmerMapper(const Graph& g, bool verification_on = true) :
      base(g, "KmerMapper"), mapping_(g.k() + 1), k_(g.k() + 1), verification_on_(verification_on) {}

  virtual ~KmerMapper() { }

  size_t get_k() const { return k_; }

  typename MapType::const_iterator begin() const {
    return mapping_.begin();
  }

  typename MapType::const_iterator end() const {
    return mapping_.end();
  }

  void Normalize() {
    std::vector<Kmer> all;
    for (auto it = begin(); it != end(); ++it) {
      all.push_back(it->first);
    }
    for (auto it = all.begin(); it != all.end(); ++it) {
      Normalize(*it);
    }
  }

  void Revert(const Kmer &kmer) {
    Kmer old_value = Substitute(kmer);
    if (old_value != kmer) {
      mapping_.erase(kmer);
      mapping_[old_value] = kmer;
    }
  }

  void Normalize(const Kmer &kmer) {
    mapping_[kmer] = Substitute(kmer);
  }

  bool CheckCanRemap(const Sequence& old_s, const Sequence& new_s) const {
    size_t old_length = old_s.size() - k_ + 1;
    size_t new_length = new_s.size() - k_ + 1;
    UniformPositionAligner aligner(old_s.size() - k_ + 1,
                                   new_s.size() - k_ + 1);
    Kmer old_kmer = old_s.start<Kmer>(k_);
    old_kmer >>= 0;

    for (size_t i = k_ - 1; i < old_s.size(); ++i) {
      old_kmer <<= old_s[i];
      size_t old_kmer_offset = i - k_ + 1;
      size_t new_kmer_offest = aligner.GetPosition(old_kmer_offset);
      if(old_kmer_offset * 2 + 1 == old_length && new_length % 2 == 0) {
        Kmer middle(k_ - 1, new_s, new_length / 2);
        if (typename Kmer::less2()(middle, !middle)) {
          new_kmer_offest = new_length - 1 - new_kmer_offest;
        }
      }
      Kmer new_kmer(k_, new_s, new_kmer_offest);
      auto it = mapping_.find(new_kmer);
      if (it != mapping_.end()) {
        if (Substitute(new_kmer) != old_kmer) {
          return false;
        }
      }
    }
    return true;
  }

  void RemapKmers(const Sequence& old_s, const Sequence& new_s) {
    VERIFY(this->IsAttached());
    size_t old_length = old_s.size() - k_ + 1;
    size_t new_length = new_s.size() - k_ + 1;
    UniformPositionAligner aligner(old_s.size() - k_ + 1,
                                   new_s.size() - k_ + 1);
    Kmer old_kmer = old_s.start<Kmer>(k_);

    for (size_t i = k_ - 1; i < old_s.size(); ++i) {
      // Instead of shifting right
      if (i != k_ - 1) {
        old_kmer <<= old_s[i];
      }

      size_t old_kmer_offset = i - k_ + 1;
      size_t new_kmer_offest = aligner.GetPosition(old_kmer_offset);
      if(old_kmer_offset * 2 + 1 == old_length && new_length % 2 == 0) {
        Kmer middle(unsigned(k_ - 1), new_s, new_length / 2);
        if(typename Kmer::less2()(middle, !middle)) {
          new_kmer_offest = new_length - 1 - new_kmer_offest;
        }
      }
      Kmer new_kmer(unsigned(k_), new_s, new_kmer_offest);
      auto it = mapping_.find(new_kmer);
      if (it != mapping_.end()) {
    	if(verification_on_)
    		VERIFY(Substitute(new_kmer) == old_kmer);
        mapping_.erase(it);
      }
      if(old_kmer != new_kmer)
            mapping_[old_kmer] = new_kmer;
    }
  }

  virtual void HandleGlue(EdgeId new_edge, EdgeId edge1, EdgeId edge2) {
    VERIFY(this->g().EdgeNucls(new_edge) == this->g().EdgeNucls(edge2));
    RemapKmers(this->g().EdgeNucls(edge1), this->g().EdgeNucls(edge2));
  }

  Kmer Substitute(const Kmer& kmer) const {
    VERIFY(this->IsAttached());
    Kmer answer = kmer;
    auto it = mapping_.find(answer);
    while (it != mapping_.end()) {
      if(verification_on_)
        VERIFY(it.first() != it.second());
      answer = it.second();
      it = mapping_.find(answer);
    }
    return answer;
  }

  void BinWrite(std::ostream& file) const {
    u_int32_t size = (u_int32_t) mapping_.size();
    file.write((const char *) &size, sizeof(u_int32_t));

    for (auto iter = mapping_.begin(); iter != mapping_.end(); ++iter) {
      Kmer::BinWrite(file, iter.first());
      Kmer::BinWrite(file, iter.second());
    }
  }

  void BinRead(std::istream& file) {
    mapping_.clear();
    u_int32_t size;
    file.read((char *) &size, sizeof(u_int32_t));

    for (u_int32_t i = 0; i < size; ++i) {
      Kmer key(k_);
      Kmer value(k_);
      Kmer::BinRead(file, &key);
      Kmer::BinRead(file, &value);
      mapping_[key] = value;
    }
  }

  bool CompareTo(KmerMapper<Graph, Kmer> const& m) {
    if (mapping_.size() != m.mapping_.size()) {
      INFO("Unequal sizes");
    }
    for (auto iter = mapping_.begin(); iter != mapping_.end(); ++iter) {
      auto cmp = m.mapping_.find(iter.first());
      if (cmp == m.mapping_.end() || cmp.second() != iter.second()) {
        return false;
      }
    }
    return true;
  }

  void clear() {
    mapping_.clear();
  }

  size_t size() const {
    return mapping_.size();
  }

  // "turn on = true" means turning of all verifies
  void SetUnsafeMode(bool turn_on){
          verification_on_ = !turn_on;
  }
};

}
