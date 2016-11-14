//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "perfect_hash_map.hpp"
#include "edge_info_updater.hpp"
#include "edge_position_index.hpp"

#include <folly/SmallLocks.h>

namespace debruijn_graph {

template<class IdType>
class EdgeInfoStorage {
public:
    typedef vector<EdgeInfo<IdType>> Content;
    typedef typename Content::iterator iterator;
    typedef typename Content::const_iterator const_iterator;
    Content content_;
    folly::MicroSpinLock lock_;

    EdgeInfoStorage(const Content &content) : content_(content) {
        lock_.init();
    }

    EdgeInfoStorage() {
        lock_.init();
    }

    EdgeInfo<IdType> &operator[](size_t i) {
        return content_[i];
    }

    iterator begin() {
        return content_.begin();
    }

    iterator end() {
        return content_.end();
    }

    const_iterator begin() const {
        return content_.cbegin();
    }

    const_iterator end() const {
        return content_.cend();
    }

    iterator find(const EdgeInfo<IdType> &info) {
        return content_.find(info);
    }

    const_iterator find(const EdgeInfo<IdType> &info) const {
        return content_.find(info);
    }

    void push_back(const EdgeInfo<IdType> &info) {
        folly::MSLGuard g(lock_);
        content_.push_back(info);
    }

    template<class... Args>
    void emplace_back(Args&&... args) {
        folly::MSLGuard g(lock_);
        content_.emplace_back(std::forward<Args>(args)...);
    }

    size_t size() const{
        return content_.size();
    }

    bool valid() const {
        //what's invalid edge info storage?
        return true;
    }

    EdgeInfoStorage conjugate(size_t k) const {
        EdgeInfoStorage result;
        for(auto it = content_.rbegin(); it != content_.rend(); ++it) {
            result.push_back(it->conjugate(k));
        }
        return result;
    }
};

//todo it is not handling graph events!!!
template<class IdType, class Seq = RtSeq,
    class traits = kmer_index_traits<Seq>,  class StoringType = SimpleStoring >
class DeBruijnEdgeMultiIndex : public KeyStoringMap<Seq, EdgeInfoStorage<IdType>, traits, StoringType > {
  typedef KeyStoringMap<Seq, EdgeInfoStorage<IdType>, traits, StoringType > base;
 public:
  typedef StoringType storing_type;
  typedef typename base::traits_t traits_t;
  typedef typename base::KMer KMer;
  typedef typename base::KMerIdx KMerIdx;
  typedef typename  base::KeyWithHash KeyWithHash;
  typedef EdgeInfoStorage<IdType> Value;

  using base::ConstructKWH;
//  typedef typename base::IdType IdType;
  //todo move this typedef up in hierarchy (need some c++ tricks)

  DeBruijnEdgeMultiIndex(unsigned k, const std::string &workdir)
      : base(k, workdir) {
      INFO("Constructing multi-kmer index");
  }

  ~DeBruijnEdgeMultiIndex() {}


  Value get(const KeyWithHash &kwh) const {
    VERIFY(contains(kwh));
    return base::get_value(kwh);
  }

  bool contains(const KeyWithHash &kwh) const {
      if (!base::valid(kwh))
          return false;
      return this->get_raw_value_reference(kwh).valid();
  }

  bool valid(const KMer &kmer) const {
      KeyWithHash kwh = base::ConstructKWH(kmer);
      return base::valid(kwh);
  }

  void PutInIndex(const KeyWithHash &kwh, IdType id, size_t offset) {
      if (!contains(kwh))
          return;

      EdgeInfoStorage<IdType> &entry = this->get_raw_value_reference(kwh);
      entry.emplace_back(id, (unsigned int)offset);
  }

  const EdgeInfoStorage<IdType> get(const KMer& kmer) const {
      auto kwh = base::ConstructKWH(kmer);
      auto entry = this->get_value(kwh);
      return entry;
  }

  //todo delete if equal seems to work improperly!!!
  bool DeleteIfEqual(const KeyWithHash &, IdType) {
      VERIFY(false);
      return false;
  }

};

}
