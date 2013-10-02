#pragma once
/*
 * edge_multi_index.hpp
 *
 *  Created on: May 24, 2013
 *      Author: anton
 */
#include "perfect_hash_map.hpp"
#include "edge_info_updater.hpp"
#include "kmer_splitters.hpp"

namespace debruijn_graph {

template<class IdType>
class EdgeInfoStorage {
public:
	typedef typename vector<EdgeInfo<IdType>> Content;
	typedef typename Content::iterator iterator;
	typedef typename Content::const_iterator const_iterator;
	Content content_;

	EdgeInfoStorage(const Content &content) : content_(content) {
	}

	EdgeInfoStorage() {
	}

	EdgeInfo<IdType> &operator[](size_t i) {
		return content_[i];
	}

	iterator begin() {
		return content_.begin();
	}

	iterator end() {
		return content_.begin();
	}

	const_iterator begin() const {
		return content_.begin();
	}

	const_iterator end() const {
		return content_.begin();
	}

	iterator find(const EdgeInfo<IdType> &info) {
		return content_.find(info);
	}

	const_iterator find(const EdgeInfo<IdType> &info) const {
		return content_.find(info);
	}

	void push_back(const EdgeInfo<IdType> &info) {
		content_.push_back(info);
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
template<class IdType, class Seq = runtime_k::RtSeq,
    class traits = kmer_index_traits<Seq> >
class DeBruijnEdgeMultiIndex : public KeyStoringMap<PerfectHashMap<Seq, EdgeInfoStorage<IdType>, traits>> {
  typedef KeyStoringMap<vector<EdgeInfo<IdType>>, traits> base;
 public:
  typedef typename base::traits_t traits_t;
  typedef typename base::KMer KMer;
  typedef typename base::KMerIdx KMerIdx;
	typedef typename  base::KeyWithHash KeyWithHash;
  using base::ConstructKWH;
//  typedef typename base::IdType IdType;
  //todo move this typedef up in hierarchy (need some c++ tricks)

  DeBruijnEdgeMultiIndex(unsigned k, const std::string &workdir)
      : base(k, workdir) {}

  ~DeBruijnEdgeMultiIndex() {}

  std::vector<EdgePosition> get(const KeyWithHash &kwh) const {
    VERIFY(contains(kwh));
    return base::operator[](kwh);
  }
//
//  std::vector<EdgePosition> get(KMerIdx idx) const {
//    VERIFY(valid_idx(idx));
//    return base::operator[](idx);
//  }

  void PutInIndex(const KMer &kmer, IdType id, size_t offset) {
    KeyWithHash kwh = base::ConstructKWH(kmer);
    if (base::contains(kwh)) {
      EdgeInfoStorage<IdType> &entry = this->get_raw_value_reference(kwh);
      EdgeInfo<IdType> new_entry;
      new_entry.edge_id = id;
      new_entry.offset = offset;
      entry.push_back(new_entry);
    }
  }

  //todo delete if equal seems to work improperly!!!
  bool DeleteIfEqual(const KeyWithHash &kwh, IdType id) {
      VERIFY(false);
      return false;
  }

};

}
