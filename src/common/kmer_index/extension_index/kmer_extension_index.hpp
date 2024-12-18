//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2015-2022 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "inout_mask.hpp"

#include "kmer_index/ph_map/perfect_hash_map.hpp"
#include "kmer_index/ph_map/kmer_maps.hpp"
#include "kmer_index/ph_map/storing_traits.hpp"

#include "sequence/rtseq.hpp"
#include "sequence/sequence.hpp"
#include "utils/stl_utils.hpp"
#include <bitset>

namespace kmers {

template<class Seq>
struct slim_kmer_index_traits : public kmers::kmer_index_traits<Seq> {
    typedef kmers::kmer_index_traits<Seq> __super;

    template<class Writer>
    static void raw_serialize(Writer&, typename __super::RawKMerStorage*) {
        VERIFY(false && "Cannot save extension index");
    }

    template<class Reader>
    static typename __super::RawKMerStorage *raw_deserialize(
            Reader&, const std::string &) {
        VERIFY(false && "Cannot load extension index");
        return NULL;
    }

};

template<typename KeyWithHash>
struct AbstractDeEdge {
    KeyWithHash start;
    KeyWithHash end;
    AbstractDeEdge(KeyWithHash s, KeyWithHash e)
            : start(std::move(s)), end(std::move(e)) {}

    bool operator==(const AbstractDeEdge &other) {
        return start == other.start && end == other.end;
    }

    bool operator!=(const AbstractDeEdge &other) {
        return !(*this == other);
    }
};

template<class stream, class KWH>
stream &operator<<(stream &s, const AbstractDeEdge<KWH> de_edge) {
    return s << "DeEdge[" << de_edge.start << ", " << de_edge.end << "]";
}

template<class traits = slim_kmer_index_traits<RtSeq>, class StoringType = DefaultStoring>
class DeBruijnExtensionIndex : public KeyIteratingMap<typename traits::SeqType, InOutMask, traits, StoringType> {
    typedef KeyIteratingMap<typename traits::SeqType, InOutMask, traits, StoringType> base;

public:
    typedef DeBruijnExtensionIndex<traits, StoringType> This;
    typedef typename base::traits_t traits_t;
    typedef StoringType storing_type;
    typedef typename base::KeyType KMer;
    typedef typename base::IdxType KMerIdx;
    typedef typename base::KeyWithHash KeyWithHash;
    typedef AbstractDeEdge<KeyWithHash> DeEdge;
    using base::ConstructKWH;
    unsigned k_size_;

    DeBruijnExtensionIndex(unsigned K)
            : base(K) {
        k_size_ = K;
    }

    using PerfectHashMap<typename traits::SeqType, InOutMask, traits, StoringType>::raw_data;
    using PerfectHashMap<typename traits::SeqType, InOutMask, traits, StoringType>::raw_size;

    This &operator|=(const char *data) {
        memor(this->raw_data(), data, this->raw_size());
        return *this;
    }

    This &operator|=(const This &rhs) {
        VERIFY(this->size() == rhs.size());
        return *this |= rhs.raw_data();
    }

    This &operator&=(const char *data) {
        memand(this->raw_data(), data, this->raw_size());
        return *this;
    }

    This &operator&=(const This &rhs) {
        VERIFY(this->size() == rhs.size());
        return *this &= rhs.raw_data();
    }

    void AddOutgoing(const KeyWithHash &kwh, char nucl) {
        TRACE("Add outgoing " << kwh << " " << ::nucl(nucl) << " " << kwh.is_minimal());
        this->get_raw_value_reference(kwh).AddOutgoing(nucl, kwh.is_minimal());
    }

    void AddIncoming(const KeyWithHash &kwh, char nucl) {
        TRACE("Add incoming " << kwh << " " << ::nucl(nucl) << " " << kwh.is_minimal());
        this->get_raw_value_reference(kwh).AddIncoming(nucl, kwh.is_minimal());
    }

    void DeleteOutgoing(const KeyWithHash &kwh, char nucl) {
        TRACE("Delete outgoing " << kwh << " " << ::nucl(nucl) << " " << kwh.is_minimal());
        this->get_raw_value_reference(kwh).DeleteOutgoing(nucl, kwh.is_minimal());
    }

    void DeleteIncoming(const KeyWithHash &kwh, char nucl) {
        TRACE("Delete incoming " << kwh << " " << ::nucl(nucl) << " " << kwh.is_minimal());
        this->get_raw_value_reference(kwh).DeleteIncoming(nucl, kwh.is_minimal());
    }

    void IsolateVertex(const KeyWithHash &kwh) {
        TRACE("Isolate vertex " << kwh);
        this->get_raw_value_reference(kwh).IsolateVertex();
    }

    void RemoveSequence(const Sequence &sequence) {
        RtSeq kmer = sequence.start<RtSeq>(k_size_);
        KeyWithHash kwh = ConstructKWH(kmer);
        IsolateVertex(kwh);
        for (size_t pos = k_size_; pos < sequence.size(); pos++) {
            kwh = kwh << sequence[pos];
            IsolateVertex(kwh);
        }
    }

    void RemoveSequences(const std::vector<Sequence> &sequences) {
#       pragma omp parallel for schedule(guided)
        for (size_t i = 0; i < sequences.size(); ++i) {
            RemoveSequence(sequences[i]);
            RemoveSequence(!sequences[i]);
        }
    }

    bool CheckOutgoing(const KeyWithHash &kwh, char nucl) const {
        return this->get_value(kwh).CheckOutgoing(nucl);
    }

    KeyWithHash GetOutgoing(const KeyWithHash &kwh, char nucl) const {
        return kwh << nucl;
    }

    bool CheckIncoming(const KeyWithHash &kwh, char nucl) const {
        return this->get_value(kwh).CheckIncoming(nucl);
    }

    KeyWithHash GetIncoming(const KeyWithHash &kwh, char nucl) const {
        return kwh >> nucl;
    }

    bool IsDeadEnd(const KeyWithHash &kwh) const {
        return this->get_value(kwh).IsDeadEnd();
    }

    bool IsDeadStart(const KeyWithHash &kwh) const {
        return this->get_value(kwh).IsDeadStart();
    }

    bool IsJunction(const KeyWithHash &kwh) const {
        return this->get_value(kwh).IsJunction();
    }

    bool CheckUniqueOutgoing(const KeyWithHash &kwh) const {
        return this->get_value(kwh).CheckUniqueOutgoing();
    }

    KeyWithHash GetUniqueOutgoing(const KeyWithHash &kwh) const {
        return GetOutgoing(kwh, this->get_value(kwh).GetUniqueOutgoing());
    }

    bool CheckUniqueIncoming(const KeyWithHash &kwh) const {
        return this->get_value(kwh).CheckUniqueIncoming();
    }

    KeyWithHash GetUniqueIncoming(const KeyWithHash &kwh) const {
        return GetIncoming(kwh, this->get_value(kwh).GetUniqueIncoming());
    }

    size_t OutgoingEdgeCount(const KeyWithHash &kwh) const {
        return this->get_value(kwh).OutgoingEdgeCount();
    }

    size_t IncomingEdgeCount(const KeyWithHash &kwh) const {
        return this->get_value(kwh).IncomingEdgeCount();
    }

    ~DeBruijnExtensionIndex() { }

private:
   DECL_LOGGER("ExtensionIndex");
};

}
