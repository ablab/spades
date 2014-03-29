#pragma once
/*
 * edge_index.hpp
 *
 *  Created on: May 24, 2013
 *      Author: anton
 */

#include "perfect_hash_map.hpp"
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

    template<class KWH>
    EdgeInfo conjugate(const KWH &kwh) const {
        return conjugate(kwh.key().size());
    }

    EdgeInfo conjugate(size_t k) const {
        if(!valid()) {
            return EdgeInfo(IdType(0), unsigned(-1), count);
        } else {
            return EdgeInfo(edge_id->conjugate(), (unsigned)edge_id->length(k) - offset, count);
        }
    }

    void invalidate() {
        offset = unsigned(-1);
    }

    bool valid() const {
        return offset != unsigned(-1);
    }
};

template<class stream, class IdType>
stream &operator<<(stream &s, const EdgeInfo<IdType> &info) {
    return s << "EdgeInfo[" << info.edge_id << ", " << info.offset << ", " << info.count << "]"; 
}

template<class Graph, class Seq = runtime_k::RtSeq, class traits = kmer_index_traits<Seq>, class StoringType = DefaultStoring>
class KmerFreeEdgeIndex : public KeyIteratingMap<Seq, EdgeInfo<typename Graph::EdgeId>, traits, StoringType> {
    typedef KeyIteratingMap<Seq, EdgeInfo<typename Graph::EdgeId>, traits, StoringType> base;
    const Graph &graph_;

public:
    typedef typename base::traits_t traits_t;
    typedef StoringType storing_type;
    typedef typename base::KMer KMer;
    typedef typename base::KMerIdx KMerIdx;
    typedef Graph GraphT;
    typedef typename Graph::EdgeId IdType;
    typedef typename base::KeyWithHash KeyWithHash;
    typedef EdgeInfo<typename Graph::EdgeId> Value;
    using base::valid;
    using base::ConstructKWH;

public:

    KmerFreeEdgeIndex(const Graph &graph, const std::string &workdir)
            : base(unsigned(graph.k() + 1), workdir), graph_(graph) {}

    /**
     * Shows if kmer has some entry associated with it
     */
    bool contains(const KeyWithHash &kwh) const {
        // Sanity check
        if (!valid(kwh)) {
            return false;
        }

        Value entry = base::get_value(kwh);

        if (entry.offset == -1u) {
            return false;
        }

        return kwh.key() == KMer(this->k(), graph_.EdgeNucls(entry.edge_id), entry.offset);
    }

    void PutInIndex(KeyWithHash &kwh, EdgeId id, size_t offset) {
        if (valid(kwh)) {
            auto &entry = this->get_raw_value_reference(kwh);
            if (!entry.valid() || contains(kwh)) {
                this->put_value(kwh, Value(id, (unsigned)offset, entry.count));
            }
        }
    }

    //Only coverage is loaded
    template<class Writer>
    void BinWrite(Writer &writer) const {
        this->index_.serialize(writer);
        size_t sz = this->data_.size();
        writer.write((char*)&sz, sizeof(sz));
        for (size_t i = 0; i < sz; ++i)
            writer.write((char*)&(this->data_[i].count), sizeof(this->data_[0].count));
    }

    template<class Reader>
    void BinRead(Reader &reader, const std::string/* &FileName*/) {
        this->clear();
        this->index_.deserialize(reader);
        size_t sz = 0;
        reader.read((char*)&sz, sizeof(sz));
        this->data_.resize(sz);
        for (size_t i = 0; i < sz; ++i)
            reader.read((char*)&(this->data_[i].count), sizeof(this->data_[0].count));
    }
};

template<class Graph, class Seq = runtime_k::RtSeq, class traits = kmer_index_traits<Seq>, class StoringType = DefaultStoring>
class KmerStoringEdgeIndex : public KeyStoringMap<Seq, EdgeInfo<typename Graph::EdgeId>, traits, StoringType> {
  typedef KeyStoringMap<Seq, EdgeInfo<typename Graph::EdgeId>, traits, StoringType> base;

public:
  typedef typename base::traits_t traits_t;
  typedef StoringType storing_type;
  typedef typename base::KMer KMer;
  typedef typename base::KMerIdx KMerIdx;
  typedef Graph GraphT;
  typedef typename Graph::EdgeId IdType;
  typedef typename base::KeyWithHash KeyWithHash;
  typedef EdgeInfo<typename Graph::EdgeId> Value;
  using base::valid;
  using base::ConstructKWH;


  KmerStoringEdgeIndex(const Graph& g, const std::string &workdir)
          : base(unsigned(g.k() + 1), workdir) {}

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
  void PutInIndex(KeyWithHash &kwh, EdgeId id, size_t offset) {
      if (valid(kwh)) {
          auto &entry = this->get_raw_value_reference(kwh);
          if (!entry.valid() || contains(kwh)) {
              this->put_value(kwh, Value(id, (unsigned)offset, entry.count));
          }
      }
  }
};

}
