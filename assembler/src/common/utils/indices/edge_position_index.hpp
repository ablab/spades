//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "perfect_hash_map.hpp"
#include "io/reads/single_read.hpp"

namespace debruijn_graph {

template<class IdType>
struct EdgeInfo {
    IdType edge_id;
    unsigned offset;
    unsigned count;

    EdgeInfo(IdType edge_id_ = IdType(), unsigned offset_ = unsigned(-1), unsigned count_ = 0) :
            edge_id(edge_id_), offset(offset_), count(count_) {
        VERIFY(edge_id != IdType() || clean());
    }

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

    void clear() {
        offset = unsigned(-1);
    }

    bool clean() const {
        return offset == unsigned(-1);
    }

    void remove() {
        offset = unsigned(-2);
    }

    bool removed() const {
        return offset == unsigned(-2);
    }

    bool valid() const {
        return !clean() && !removed();
    }
};

template<class stream, class IdType>
stream &operator<<(stream &s, const EdgeInfo<IdType> &info) {
    return s << "EdgeInfo[" << info.edge_id.int_id() << ", " << info.offset << ", " << info.count << "]";
}

template<class Graph, class StoringType = DefaultStoring>
class KmerFreeEdgeIndex : public KeyIteratingMap<RtSeq, EdgeInfo<typename Graph::EdgeId>,
        kmer_index_traits<RtSeq>, StoringType> {
    typedef KeyIteratingMap<RtSeq, EdgeInfo<typename Graph::EdgeId>,
            kmer_index_traits<RtSeq>, StoringType> base;
    const Graph &graph_;

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

public:

    KmerFreeEdgeIndex(const Graph &graph, const std::string &workdir)
            : base(unsigned(graph.k() + 1), workdir), graph_(graph) {}

    /**
     * Shows if kmer has some entry associated with it
     */
    bool contains(const KeyWithHash &kwh) const {
        // Sanity check
        if (!valid(kwh))
            return false;

        KmerPos entry = base::get_value(kwh);
        if (!entry.valid())
            return false;
        return graph_.EdgeNucls(entry.edge_id).contains(kwh.key(), entry.offset);
    }

    void PutInIndex(KeyWithHash &kwh, IdType id, size_t offset) {
        if (!valid(kwh))
            return;
        
        KmerPos &entry = this->get_raw_value_reference(kwh);
        if (entry.removed()) {
            //VERIFY(false);
            return;
        }
        if (entry.clean()) {
            //put verify on this conversion!
            this->put_value(kwh, KmerPos(id, (unsigned)offset, entry.count));
        } else if (contains(kwh)) {
            //VERIFY(false);
            entry.remove();
        } else {
            //VERIFY(false);
            //FIXME bad situation; some other kmer is there; think of putting verify
        }
    }

    //Only coverage is loaded
    template<class Writer>
    void BinWrite(Writer &writer) const {
        this->index_ptr_->serialize(writer);
        size_t sz = this->data_.size();
        writer.write((char*)&sz, sizeof(sz));
        for (size_t i = 0; i < sz; ++i)
            writer.write((char*)&(this->data_[i].count), sizeof(this->data_[0].count));
    }

    template<class Reader>
    void BinRead(Reader &reader, const std::string/* &FileName*/) {
        this->clear();
        this->index_ptr_->deserialize(reader);
        size_t sz = 0;
        reader.read((char*)&sz, sizeof(sz));
        this->data_.resize(sz);
        for (size_t i = 0; i < sz; ++i)
            reader.read((char*)&(this->data_[i].count), sizeof(this->data_[0].count));
    }
};

template<class Graph, class StoringType = DefaultStoring>
class KmerStoringEdgeIndex : public KeyStoringMap<RtSeq, EdgeInfo<typename Graph::EdgeId>,
        kmer_index_traits<RtSeq>, StoringType> {
  typedef KeyStoringMap<RtSeq, EdgeInfo<typename Graph::EdgeId>,
          kmer_index_traits<RtSeq>, StoringType> base;

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
      this->index_ptr_->serialize(writer);
      size_t sz = this->data_.size();
      writer.write((char*)&sz, sizeof(sz));
      for (size_t i = 0; i < sz; ++i)
          writer.write((char*)&(this->data_[i].count), sizeof(this->data_[0].count));
      this->BinWriteKmers(writer);
  }

  template<class Reader>
  void BinRead(Reader &reader, const std::string &FileName) {
      this->clear();
      this->index_ptr_->deserialize(reader);
      size_t sz = 0;
      reader.read((char*)&sz, sizeof(sz));
      this->data_.resize(sz);
      for (size_t i = 0; i < sz; ++i)
          reader.read((char*)&(this->data_[i].count), sizeof(this->data_[0].count));
      this->BinReadKmers(reader, FileName);
  }

  void PutInIndex(KeyWithHash &kwh, IdType id, size_t offset) {
      //here valid already checks equality of query-kmer and stored-kmer sequences
      if (base::valid(kwh)) {
          KmerPos &entry = this->get_raw_value_reference(kwh);
          if (entry.removed())
              return;
          if (!entry.clean()) {
              this->put_value(kwh, KmerPos(id, (unsigned)offset, entry.count));
          } else {
              entry.remove();
          }
      }
  }
};

}
