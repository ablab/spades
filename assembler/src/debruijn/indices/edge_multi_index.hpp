#pragma once
/*
 * edge_multi_index.hpp
 *
 *  Created on: May 24, 2013
 *      Author: anton
 */
#include "debruijn_kmer_index.hpp"
#include "edge_info_updater.hpp"
#include "kmer_splitters.hpp"

namespace debruijn_graph {

//todo it is not handling graph events!!!
template<class IdType, class Seq = runtime_k::RtSeq,
    class traits = kmer_index_traits<Seq> >
class DeBruijnEdgeMultiIndex : public DeBruijnKMerIndex<KmerStoringIndex<vector<EdgeInfo<IdType>>, traits>> {
  typedef DeBruijnKMerIndex<KmerStoringIndex<vector<EdgeInfo<IdType>>, traits>> base;
 public:
 public:
  typedef typename base::traits_t traits_t;
  typedef typename base::KMer KMer;
  typedef typename base::KMerIdx KMerIdx;
//  typedef typename base::IdType IdType;
  //todo move this typedef up in hierarchy (need some c++ tricks)

  DeBruijnEdgeMultiIndex(unsigned k, const std::string &workdir)
      : base(k, workdir) {
	  INFO("DeBruijnEdgeMultiIndex constructing");
  }

  ~DeBruijnEdgeMultiIndex() {}

  std::vector<EdgePosition> get(const KMer &kmer) const {
    VERIFY(contains(kmer));
    return base::operator[](kmer);
  }

  std::vector<EdgePosition> get(KMerIdx idx) const {
    VERIFY(valid_idx(idx));
    return base::operator[](idx);
  }

  void PutInIndex(const KMer &kmer, IdType id, int offset,  bool /*ignore_new_kmer*/ = false) {
    size_t idx = base::seq_idx(kmer);
    TRACE("put in index");
    if (base::valid_key(idx, kmer))  {
		std::vector<EdgeInfo<IdType> > &entry = base::operator[](idx);
		EdgeInfo<IdType> new_entry;
		new_entry.edge_id = id;
		new_entry.offset = offset;
		entry.push_back(new_entry);
    }
  }

  //todo delete if equal seems to work improperly!!!
  bool DeleteIfEqual(const KMer &kmer, IdType id) {
      VERIFY(false);
      return false;
//    typename base::KMerIdx idx = base::seq_idx(kmer);
//
//    // Early exit if kmer has not been seen at all
//    if (idx == base::InvalidKMerIdx)
//      return false;
//
//    // Now we know that idx is in range. Check the edge id.
//    typename base::KMerIndexValueType &entry = base::operator[](idx);
//
//    if (entry.edgeId_ == id) {
//      entry.offset_ = -1;
//      return true;
//    }
//
//    return false;
  }

};

}
