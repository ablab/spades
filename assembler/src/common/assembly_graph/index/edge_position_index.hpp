//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "sequence/rtseq.hpp"
#include "utils/ph_map/perfect_hash_map.hpp"

#include <folly/SmallLocks.h>

namespace debruijn_graph {

template<class IdType>
class EdgeInfo {
    IdType edge_id_;

    typedef folly::PicoSpinLock<uint32_t> OffsetType;
    OffsetType offset_with_lock_;
    // This is a bit hacky...
    static constexpr unsigned CLEARED = -1u & ~OffsetType::kLockBitMask_;
    static constexpr unsigned TOMBSTONE = -2u & ~OffsetType::kLockBitMask_;

 public:
    EdgeInfo(IdType e = IdType(), unsigned o = CLEARED) :
        edge_id_(e) {
        offset_with_lock_.init();
        offset_with_lock_.setData(o);
        VERIFY(edge_id_ != IdType() || clean());
    }

    template<class Graph>
    EdgeInfo conjugate(const Graph &g) const {
        if (!valid())
            return EdgeInfo(IdType(), CLEARED);

        return EdgeInfo(g.conjugate(edge()), unsigned(g.length(edge()) - offset() - 1));
    }

    IdType edge() const {
        return edge_id_;
    }

    void lock() { offset_with_lock_.lock(); }
    void unlock() { offset_with_lock_.unlock(); }

    unsigned offset() const { return offset_with_lock_.getData(); }

    void clear() {
      offset_with_lock_.setData(CLEARED);
      edge_id_ = IdType();
    }
    bool clean() const { return offset() == CLEARED; }

    void remove() {
      offset_with_lock_.setData(TOMBSTONE);
      edge_id_ = IdType();
    }
    bool removed() const { return offset() == TOMBSTONE; }

    bool valid() const {
        return !clean() && !removed();
    }
};

template<class stream, class IdType>
stream &operator<<(stream &s, const EdgeInfo<IdType> &info) {
    return s << "EdgeInfo[" << info.edge() << ", " << info.offset() << "]";
}

template<class Graph,
         typename IdType = typename Graph::EdgeId>
struct GraphInverter {
    const Graph &g_;

    GraphInverter(const Graph &g)
            : g_(g) {}

    template<class K>
    EdgeInfo<IdType> operator()(const EdgeInfo<IdType>& v, const K&) const {
        return v.conjugate(g_);
    }
};


template<class Graph, class StoringType = utils::DefaultStoring>
class KmerFreeEdgeIndex : public utils::PerfectHashMap<RtSeq,
                                                       EdgeInfo<typename Graph::EdgeId>,
                                                       utils::kmer_index_traits<RtSeq>, StoringType> {
    typedef utils::PerfectHashMap<RtSeq, EdgeInfo<typename Graph::EdgeId>,
            utils::kmer_index_traits<RtSeq>, StoringType> base;
    const Graph &graph_;

public:
    typedef typename base::traits_t traits_t;
    typedef StoringType storing_type;
    typedef typename base::KeyType KMer;
    typedef typename base::IdxType KMerIdx;
    typedef Graph GraphT;
    typedef typename Graph::EdgeId IdType;
    typedef typename base::KeyWithHash KeyWithHash;
    typedef EdgeInfo<typename Graph::EdgeId> KmerPos;
    using base::valid;
    using base::ConstructKWH;

public:

    KmerFreeEdgeIndex(const Graph &graph)
            : base(unsigned(graph.k() + 1)), graph_(graph) {}


    KmerPos get_value(const KeyWithHash &kwh) const {
        return base::get_value(kwh, GraphInverter<Graph>(graph_));
    }

    void put_value(const KeyWithHash &kwh, const KmerPos &pos) {
        base::put_value(kwh, pos, GraphInverter<Graph>(graph_));
    }

    /**
     * Shows if kmer has some entry associated with it
     */
    bool contains(const KeyWithHash &kwh) const {
        // Sanity check
        if (!valid(kwh))
            return false;

        KmerPos entry = get_value(kwh);
        if (!entry.valid())
            return false;

        return graph_.EdgeNucls(entry.edge()).contains(kwh.key(), entry.offset());
    }

    void PutInIndex(KeyWithHash &kwh, IdType id, size_t offset) {
        if (!valid(kwh))
            return;

        KmerPos &entry = this->get_raw_value_reference(kwh);
        if (entry.removed())
            return;

        entry.lock();
        if (entry.clean()) {
            // Note that this releases the lock as well!
            put_value(kwh, KmerPos(id, (unsigned)offset));
        } else if (contains(kwh)) {
            entry.remove();
        }
        entry.unlock();
    }
};

template<class Graph, class StoringType = utils::DefaultStoring>
class KmerStoringEdgeIndex : public utils::KeyStoringMap<RtSeq, EdgeInfo<typename Graph::EdgeId>,
        utils::kmer_index_traits<RtSeq>, StoringType> {
  typedef utils::KeyStoringMap<RtSeq, EdgeInfo<typename Graph::EdgeId>,
          utils::kmer_index_traits<RtSeq>, StoringType> base;

public:
  typedef typename base::traits_t traits_t;
  typedef StoringType storing_type;
  typedef typename base::KMer KMer;
  typedef typename base::KMerIdx KMerIdx;
  typedef Graph GraphT;
  typedef typename Graph::EdgeId IdType;
  typedef typename base::KeyWithHash KeyWithHash;
  typedef EdgeInfo<typename Graph::EdgeId> KmerPos;
  using base::valid;
  using base::ConstructKWH;


  KmerStoringEdgeIndex(const Graph& g)
          : base(unsigned(g.k() + 1)) {}

  ~KmerStoringEdgeIndex() {}

  /**
   * Shows if kmer has some entry associated with it
   */
  bool contains(const KeyWithHash &kwh) const {
      if (!base::valid(kwh))
          return false;
      return this->get_raw_value_reference(kwh).valid();
  }

  template<class Writer>
  void BinWrite(Writer &writer) const {
      this->index_ptr_->serialize(writer);
      size_t sz = this->data_.size();
      writer.write((char*)&sz, sizeof(sz));
      this->BinWriteKmers(writer);
  }

  template<class Reader>
  void BinRead(Reader &reader, const std::string &FileName) {
      this->clear();
      this->index_ptr_->deserialize(reader);
      size_t sz = 0;
      reader.read((char*)&sz, sizeof(sz));
      this->data_.resize(sz);
      this->BinReadKmers(reader, FileName);
  }

  void PutInIndex(KeyWithHash &kwh, IdType id, size_t offset) {
      // Here valid already checks equality of query-kmer and stored-kmer sequences
      if (!base::valid(kwh))
        return;

      KmerPos &entry = this->get_raw_value_reference(kwh);
      if (entry.removed())
        return;

      entry.lock();
      if (entry.clean()) {
        put_value(kwh, KmerPos(id, (unsigned)offset));
      } else {
        entry.remove();
      }
      entry.unlock();
  }
};

}
