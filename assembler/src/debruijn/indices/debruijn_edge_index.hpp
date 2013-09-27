#pragma once
/*
 * edge_index.hpp
 *
 *  Created on: May 24, 2013
 *      Author: anton
 */

#include "debruijn_kmer_index.hpp"
#include "edge_info_updater.hpp"
#include "kmer_splitters.hpp"

namespace debruijn_graph {

template<class IdType>
struct EdgeInfo {
    IdType edge_id;
    unsigned offset;
    unsigned count;

    EdgeInfo(IdType edge_id_ = IdType(), unsigned offset_ = -1u, unsigned count_ = 0) :
            edge_id(edge_id_), offset(offset_), count(count_) { }
};

//fixme name
//todo reduce number of template parameters
template<class Graph, class Seq = runtime_k::RtSeq, class traits = kmer_index_traits<Seq>>
class KmerFreeDeBruijnEdgeIndex : public DeBruijnKMerIndex<KmerFreeIndex<EdgeInfo<typename Graph::EdgeId>, traits>> {
    typedef DeBruijnKMerIndex<KmerFreeIndex<EdgeInfo<typename Graph::EdgeId>, traits>> base;
    const Graph &graph_;

  public:
    typedef typename base::traits_t traits_t;
    typedef typename base::KMer KMer;
    typedef typename base::KMerIdx KMerIdx;
    typedef Graph GraphT;
    typedef typename Graph::EdgeId IdType;
    using base::valid_idx;
    using base::seq_idx;
    using base::BinWriteKmers;

    KmerFreeDeBruijnEdgeIndex(unsigned k, const Graph &graph, const std::string &workdir)
            : base(k, workdir), graph_(graph) {}

    /**
     * Shows if kmer has some entry associated with it
     */
    bool contains(KMerIdx idx, const KMer &k) const {
        // Sanity check
        if (!valid_idx(idx))
            return false;

        const typename base::ValueType &entry = base::operator[](idx);

        if (entry.offset == -1u)
            return false;

        return k == KMer(this->k(), graph_.EdgeNucls(entry.edge_id), entry.offset);
    }

    KMer kmer(typename base::KMerIdx idx) const {
        VERIFY(valid_idx(idx));
        const typename base::ValueType &entry = base::operator[](idx);
        VERIFY(entry.offset != -1u);
        return KMer(this->k(), graph_.EdgeNucls(entry.edge_id), entry.offset);
    }

    //todo current strategy of putting in index if slot is vacant,
    //can lead to funny behavior during gap closing regarding coverage of new k-mers!
    //Currently used both for filling and update
    void PutInIndex(const KMer &kmer, IdType e, size_t offset) {
    	TRACE("Put in KmerFreeDeBruijnEdgeIndex");
        KMerIdx idx = seq_idx(kmer);
        if (!valid_idx(idx))
            return;
        EdgeInfo<IdType>& entry = this->operator[](idx);
        // if slot is empty or it contains information about this k-mer
        // (second condition is almost always not useful)
        if (entry.offset == -1u || contains(idx, kmer)) {
            entry.edge_id = e;
            entry.offset = (unsigned) offset;
        }
    }

};

//fixme name
template<class Graph, class Seq = runtime_k::RtSeq, class traits = kmer_index_traits<Seq>>
class KmerStoringDeBruijnEdgeIndex : public DeBruijnKMerIndex<KmerStoringIndex<EdgeInfo<typename Graph::EdgeId>, traits>> {
  typedef DeBruijnKMerIndex<KmerStoringIndex<EdgeInfo<typename Graph::EdgeId>, traits>> base;

 public:
  typedef typename base::traits_t traits_t;
  typedef typename base::KMer KMer;
  typedef typename base::KMerIdx KMerIdx;
  typedef Graph GraphT;
  typedef typename Graph::EdgeId IdType;

  KmerStoringDeBruijnEdgeIndex(size_t K, const Graph& , const std::string &workdir)
          : base(K, workdir) {}

  ~KmerStoringDeBruijnEdgeIndex() {}

  /**
   * Shows if kmer has some entry associated with it
   */
  bool contains(KMerIdx idx, const KMer &k) const {
		  if (!base::valid_key(idx, k))
          return false;

      const typename base::ValueType &entry = base::operator[](idx);
      return entry.offset != -1u;
  }

  void PutInIndex(const KMer &kmer, IdType id, int offset, bool /*ignore_new_kmer*/ = false) {
	TRACE("Put in KmerStoringDeBruijnEdgeIndex");
    size_t idx = base::seq_idx(kmer);
    if (base::valid_key(idx, kmer)) {
      EdgeInfo<IdType> &entry = base::operator[](idx);
      entry.edge_id = id;
      entry.offset = offset;
    }
  }
};

//todo rename to account for the fact that it is not multi index
template<class Index>
class DeBruijnEdgeIndex : public Index {
    typedef Index base;

  public:
    typedef typename base::traits_t traits_t;
    typedef typename base::KMer KMer;
    typedef typename base::KMerIdx KMerIdx;
    typedef typename base::GraphT GraphT;
    typedef typename base::IdType IdType;

    using base::contains;

    DeBruijnEdgeIndex(unsigned K, const GraphT &graph, const std::string &workdir)
            : base(K, graph, workdir) {}

    //todo why do we need to check equality???!!!
    bool DeleteIfEqual(const KMer &kmer, EdgeId e) {
        KMerIdx idx = this->seq_idx(kmer);
        if (!contains(idx, kmer))
            return false;

        EdgeInfo<EdgeId> &entry = this->operator[](idx);
        if (entry.edge_id == e) {
            entry.offset = -1u;
            return true;
        }
        return false;
    }

    //todo change to unsigned
    std::pair<IdType, size_t> get(KMerIdx idx, const KMer &kmer) const {
        VERIFY(contains(idx, kmer));

        const EdgeInfo<IdType> &entry = base::operator[](idx);
        return std::make_pair(entry.edge_id, (size_t)entry.offset);
    }

    //todo change to unsigned
    std::pair<IdType, size_t> get(const KMer &kmer) const {
        typename base::KMerIdx idx = base::seq_idx(kmer);
        return get(idx, kmer);
    }

    /**
     * Shows if kmer has some entry associated with it
     */
    bool contains(const KMer& kmer) const {
			  KMerIdx idx = this->seq_idx(kmer);
        return contains(idx, kmer);
    }

    template<class Writer>
    void BinWrite(Writer &writer) const {
        this->index_.serialize(writer);
        size_t sz = this->data_.size();
        writer.write((char*)&sz, sizeof(sz));
        for (size_t i = 0; i < sz; ++i)
            writer.write((char*)&(this->data_[i].count), sizeof(this->data_[0].count));
        this->BinWriteKmers(writer);
    }

    template<class Reader>
    void BinRead(Reader &reader, const std::string &FileName) {
        this->clear();
        this->index_.deserialize(reader);
        size_t sz = 0;
        reader.read((char*)&sz, sizeof(sz));
        this->data_.resize(sz);
        for (size_t i = 0; i < sz; ++i)
            reader.read((char*)&(this->data_[i].count), sizeof(this->data_[0].count));
        this->BinReadKmers(reader, FileName);
    }

};

}
