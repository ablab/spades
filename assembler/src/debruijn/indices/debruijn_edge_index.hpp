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
class KmerFreeDeBruijnEdgeIndex : public DeBruijnKMerIndex<EdgeInfo<typename Graph::EdgeId>, traits> {
    typedef DeBruijnKMerIndex<EdgeInfo<typename Graph::EdgeId>, traits> base;
    const Graph &graph_;

 protected:
    template<class Writer>
    void BinWriteKmers(Writer &writer) const {
        //empty
    }

    template<class Reader>
    void BinReadKmers(Reader &reader, const std::string &FileName) {
        this->kmers = NULL;
    }

  public:
    typedef typename base::traits_t traits_t;
    typedef typename base::KMer KMer;
    typedef typename base::KMerIdx KMerIdx;
    typedef Graph GraphT;
    typedef typename Graph::EdgeId IdType;
    typedef InnerDeBruijnTotallyKMerFreeIndexBuilder<KmerFreeDeBruijnEdgeIndex> BuilderT;

    KmerFreeDeBruijnEdgeIndex(unsigned K, const Graph &graph, const std::string &workdir)
            : base(K, workdir), graph_(graph) {}

    /**
     * Doesn't work without already constructed condensed graph!!!
     */
    bool contains(KMerIdx idx, const KMer &k) const {
        // Sanity check
        if (!valid_idx(idx))
            return false;

        const typename base::ValueType &entry = base::operator[](idx);

        if (entry.offset == -1u)
            return false;

        return k == KMer(this->K(), graph_.EdgeNucls(entry.edge_id), entry.offset);
    }

    bool contains(const KMer& kmer) const {
        KMerIdx idx = seq_idx(kmer);
        return contains(idx, kmer);
    }

    KMer kmer(typename base::KMerIdx idx) const {
        VERIFY(valid_idx(idx));
        const typename base::ValueType &entry = base::operator[](idx);
        VERIFY(entry.offset != -1u);
        return KMer(this->K(), graph_.EdgeNucls(entry.edge_id), entry.offset);
    }

    //todo current strategy of putting in index if slot is vacant,
    //can lead to funny behavior during gap closing regarding coverage of new k-mers!
    //Currently used both for filling and update
    void PutInIndex(const KMer &kmer, IdType e, size_t offset) {
        KMerIdx idx = seq_idx(kmer);
        if (!valid_idx(idx))
            return;
        EdgeInfo<IdType>& entry = operator[](idx);
        // if slot is empty or it contains information about this k-mer
        // (second condition is almost always not useful)
        if (entry.offset == -1u || contains(idx, kmer)) {
            entry.edge_id = e;
            entry.offset = offset;
        }
    }

};

////fixme name
//template<class Graph, class Seq = runtime_k::RtSeq, class traits = kmer_index_traits<Seq>>
//class KmerStoringDeBruijnEdgeIndex : public DeBruijnKMerIndex<EdgeInfo<typename Graph::EdgeId>, traits> {
//  typedef DeBruijnKMerIndex<EdgeInfo<typename Graph::EdgeId>, traits> base;
//
// protected:
//  template<class Writer>
//  void BinWriteKmers(Writer &writer) const {
//      traits_t::raw_serialize(writer, this->kmers);
//  }
//
//  template<class Reader>
//  void BinReadKmers(Reader &reader, const std::string &FileName) {
//      this->kmers = traits_t::raw_deserialize(reader, FileName);
//  }
//
// public:
//  typedef typename base::traits_t traits_t;
//  typedef typename base::KMer KMer;
//  typedef typename base::KMerIdx KMerIdx;
//  typedef Graph GraphT;
//  typedef typename Graph::EdgeId IdType;
//  typedef InnerDeBruijnKMerStoringIndexBuilder<KmerStoringDeBruijnEdgeIndex> BuilderT;
//
//  KmerStoringDeBruijnEdgeIndex(unsigned K, const Graph& , const std::string &workdir)
//          : base(K, workdir) {}
//
//  ~KmerStoringDeBruijnEdgeIndex() {}
//
//  bool contains(KMerIdx idx, const KMer &k) const {
//      if (!valid_idx(idx))
//        return false;
//
//      auto it = this->kmers->begin() + idx;
//      return (typename traits::raw_equal_to()(k, *it));
//  }
//
//  KMer kmer(typename base::KMerIdx idx) const {
//      VERIFY(valid_idx(idx));
//
//      auto it = this->kmers->begin() + idx;
//      return (typename traits::raw_create()(this->K(), *it));
//  }
//
//// todo discuss with AntonK old strange version?
////  template<class Writer>
////  void BinWrite(Writer &writer) const {
////    base::index_.serialize(writer);
////    size_t sz = base::data_.size();
////    writer.write((char*)&sz, sizeof(sz));
////    for (size_t i = 0; i < sz; ++i)
////      writer.write((char*)&(base::data_[i].count_), sizeof(base::data_[0].count_));
////    sz = base::push_back_buffer_.size();
////    writer.write((char*)&sz, sizeof(sz));
////    for (size_t i = 0; i < sz; ++i)
////      writer.write((char*)&(base::push_back_buffer_[i].count_), sizeof(base::push_back_buffer_[0].count_));
////    for (auto it = base::push_back_index_.left.begin(), e = base::push_back_index_.left.end(); it != e; ++it) {
////      size_t idx = it->second;
////      KMer::BinWrite(writer, it->first);
////      writer.write((char*)&idx, sizeof(idx));
////      sz -= 1;
////    }
////    VERIFY(sz == 0);
////    traits::raw_serialize(writer, base::kmers);
////  }
////
////  template<class Reader>
////  void BinRead(Reader &reader, const std::string &FileName) {
////    base::clear();
////    base::index_.deserialize(reader);
////    size_t sz = 0;
////    reader.read((char*)&sz, sizeof(sz));
////    base::data_.resize(sz);
////    for (size_t i = 0; i < sz; ++i)
////      reader.read((char*)&(base::data_[i].count_), sizeof(base::data_[0].count_));
////    reader.read((char*)&sz, sizeof(sz));
////    base::push_back_buffer_.resize(sz);
////    for (size_t i = 0; i < sz; ++i)
////      reader.read((char*)&(base::push_back_buffer_[i].count_), sizeof(base::push_back_buffer_[0].count_));
////    for (size_t i = 0; i < sz; ++i) {
////      KMer s(base::K_);
////      size_t idx;
////
////      s.BinRead(reader);
////      reader.read((char*)&idx, sizeof(idx));
////
////      base::push_back_index_.insert(typename base::KMerPushBackIndexType::value_type(s, idx));
////    }
////    base::kmers = traits::raw_deserialize(reader, FileName);
////  }
//
//
//  void PutInIndex(const KMer &kmer, IdType id, int offset, bool ignore_new_kmer = false) {
//    size_t idx = base::seq_idx(kmer);
//    if (base::contains(idx, kmer)) {
//      EdgeInfo<IdType> &entry = base::operator[](idx);
//      entry.edge_id = id;
//      entry.offset = offset;
//    }
//  }
//};

//fixme name
template<class Graph, class Seq = runtime_k::RtSeq, class traits = kmer_index_traits<Seq>>
class KmerStoringDeBruijnEdgeIndex : public KmerStoringIndex<EdgeInfo<typename Graph::EdgeId>, traits> {
  typedef KmerStoringIndex<EdgeInfo<typename Graph::EdgeId>, traits> base;

 public:
  typedef typename base::traits_t traits_t;
  typedef typename base::KMer KMer;
  typedef typename base::KMerIdx KMerIdx;
  typedef Graph GraphT;
  typedef typename Graph::EdgeId IdType;
  typedef InnerDeBruijnKMerStoringIndexBuilder<KmerStoringDeBruijnEdgeIndex> BuilderT;

  KmerStoringDeBruijnEdgeIndex(size_t K, const Graph& , const std::string &workdir)
          : base(K, workdir) {}

  ~KmerStoringDeBruijnEdgeIndex() {}

// todo discuss with AntonK old strange version?
//  template<class Writer>
//  void BinWrite(Writer &writer) const {
//    base::index_.serialize(writer);
//    size_t sz = base::data_.size();
//    writer.write((char*)&sz, sizeof(sz));
//    for (size_t i = 0; i < sz; ++i)
//      writer.write((char*)&(base::data_[i].count_), sizeof(base::data_[0].count_));
//    sz = base::push_back_buffer_.size();
//    writer.write((char*)&sz, sizeof(sz));
//    for (size_t i = 0; i < sz; ++i)
//      writer.write((char*)&(base::push_back_buffer_[i].count_), sizeof(base::push_back_buffer_[0].count_));
//    for (auto it = base::push_back_index_.left.begin(), e = base::push_back_index_.left.end(); it != e; ++it) {
//      size_t idx = it->second;
//      KMer::BinWrite(writer, it->first);
//      writer.write((char*)&idx, sizeof(idx));
//      sz -= 1;
//    }
//    VERIFY(sz == 0);
//    traits::raw_serialize(writer, base::kmers);
//  }
//
//  template<class Reader>
//  void BinRead(Reader &reader, const std::string &FileName) {
//    base::clear();
//    base::index_.deserialize(reader);
//    size_t sz = 0;
//    reader.read((char*)&sz, sizeof(sz));
//    base::data_.resize(sz);
//    for (size_t i = 0; i < sz; ++i)
//      reader.read((char*)&(base::data_[i].count_), sizeof(base::data_[0].count_));
//    reader.read((char*)&sz, sizeof(sz));
//    base::push_back_buffer_.resize(sz);
//    for (size_t i = 0; i < sz; ++i)
//      reader.read((char*)&(base::push_back_buffer_[i].count_), sizeof(base::push_back_buffer_[0].count_));
//    for (size_t i = 0; i < sz; ++i) {
//      KMer s(base::K_);
//      size_t idx;
//
//      s.BinRead(reader);
//      reader.read((char*)&idx, sizeof(idx));
//
//      base::push_back_index_.insert(typename base::KMerPushBackIndexType::value_type(s, idx));
//    }
//    base::kmers = traits::raw_deserialize(reader, FileName);
//  }


  void PutInIndex(const KMer &kmer, IdType id, int offset, bool ignore_new_kmer = false) {
    size_t idx = base::seq_idx(kmer);
    if (base::contains(idx, kmer)) {
      EdgeInfo<IdType> &entry = base::operator[](idx);
      entry.edge_id = id;
      entry.offset = offset;
    }
  }
};

template<class Index>
class DeBruijnEdgeIndex : public Index {
    typedef Index base;

  public:
    typedef typename base::traits_t traits_t;
    typedef typename base::KMer KMer;
    typedef typename base::KMerIdx KMerIdx;
    typedef typename base::GraphT GraphT;
    typedef typename base::IdType IdType;
    typedef typename base::BuilderT BuilderT;

    using base::contains;

    DeBruijnEdgeIndex(unsigned K, const GraphT &graph, const std::string &workdir)
            : base(K, graph, workdir) {}

    //todo why do we need to check equality???!!!
    bool DeleteIfEqual(const KMer &kmer, EdgeId e) {
        KMerIdx idx = seq_idx(kmer);
        if (!contains(idx, kmer))
            return false;

        EdgeInfo<EdgeId> &entry = operator[](idx);
        if (entry.edge_id == e) {
            entry.offset = -1u;
            return true;
        }
        return false;
    }

    /**
     * Number of edges coming into param edge's end
     */
    unsigned RivalEdgeCount(const KMer &kmer) const {
      KMer kmer2 = kmer << 'A';
      unsigned res = 0;
      for (char c = 0; c < 4; ++c)
        if (contains(kmer2 >> c))
          res += 1;

      return res;
    }

    /**
     * Number of edges going out of the param edge's end
     */
    unsigned NextEdgeCount(const KMer &kmer) const {
      unsigned res = 0;
      for (char c = 0; c < 4; ++c)
        if (contains(kmer << c))
          res += 1;

      return res;
    }

    KMer NextEdge(const KMer &kmer) const { // returns any next edge
      for (char c = 0; c < 4; ++c) {
        KMer s = kmer << c;
        if (contains(s))
          return s;
      }

      VERIFY_MSG(false, "Couldn't find requested edge!");
      return KMer(base::K());
      // no next edges (we should request one here).
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

    template<class Writer>
    void BinWrite(Writer &writer) const {
        this->index_.serialize(writer);
        size_t sz = this->data_.size();
        writer.write((char*)&sz, sizeof(sz));
        for (size_t i = 0; i < sz; ++i)
            writer.write((char*)&(this->data_[i].count), sizeof(this->data_[0].count));
        BinWriteKmers(writer);
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
        BinReadKmers(reader, FileName);
    }

};

}
