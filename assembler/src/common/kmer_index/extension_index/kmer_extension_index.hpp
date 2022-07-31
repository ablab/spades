//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "kmer_index/ph_map/perfect_hash_map.hpp"
#include "kmer_index/ph_map/kmer_maps.hpp"
#include "kmer_index/ph_map/storing_traits.hpp"
#include "sequence/rtseq.hpp"
#include "utils/stl_utils.hpp"

#include <bitset>

namespace kmers {

inline uint8_t invert_byte_slow(uint8_t a) {
    size_t res = 0;
    for (size_t i = 0; i < 8; i++) {
        res <<= 1;
        res += a & 1;
        a = uint8_t(a >> 1);
    }
    return uint8_t(res);
}

inline std::vector<uint8_t> count_invert_byte() {
    std::vector<uint8_t> result;
    for (size_t a = 0; a < 256; a++) {
        result.push_back(invert_byte_slow((uint8_t)a));
    }
    return result;
}

inline uint8_t invert_byte(uint8_t a) {
    static std::vector<uint8_t> precalc = count_invert_byte();
    return precalc[a];
}

// TODO Force data alignment and implement more efficient bitwise ops
// TODO Use OpenMP here
inline void memor(char *self, const char *rhs, size_t count) {
    for (size_t i = 0; i < count; ++i) {
        self[i] |= rhs[i];
    }
}

inline void memand(char *self, const char *rhs, size_t count) {
    for (size_t i = 0; i < count; ++i) {
        self[i] &= rhs[i];
    }
}

class InOutMask {
private:
    uint8_t mask_;

    static constexpr bool CheckUnique(uint8_t mask) {
        constexpr bool unique[] =
                { 0, 1, 1, 0, 1, 0, 0, 0,
                  1, 0, 0, 0, 0, 0, 0, 0 };
        return unique[mask];
    }

    static constexpr char GetUnique(uint8_t mask) {
        constexpr char next[] =
                { -1,  0,  1, -1,  2, -1, -1, -1,
                  3, -1, -1, -1, -1, -1, -1, -1 };
        return next[mask];
    }

    static constexpr size_t Count(uint8_t mask) {
        constexpr char count[] =
                { 0, 1, 1, 2, 1, 2, 2, 3,
                  1, 2, 2, 3, 2, 3, 3, 4 };
        return count[mask];
    }

    static constexpr char inv_position(char nucl, bool as_is) {
        return as_is ? nucl: char(7 - nucl);
    }

    uint8_t outgoing() const {
        return mask_ & 0xF;
    }

    uint8_t incoming() const {
        return (mask_ >> 4) & 0xF;
    }

public:
    explicit InOutMask(uint8_t mask = 0)
            : mask_(mask) {}

    uint8_t get_mask() const {
        return mask_;
    }

    template<class Key>
    InOutMask conjugate(const Key & /*k*/) const {
        return InOutMask(invert_byte(mask_));
    }

    void AddOutgoing(char nnucl, bool as_is) {
        unsigned nmask = (unsigned) (1 << inv_position(nnucl, as_is));
        if (!(mask_ & nmask)) {
#           pragma omp atomic
            mask_ |= (unsigned char) nmask;
        }
    }

    void AddIncoming(char pnucl, bool as_is) {
        unsigned pmask = (unsigned) (1 << inv_position(char(pnucl + 4), as_is));
        if (!(mask_ & pmask)) {
#           pragma omp atomic
            mask_ |= (unsigned char) pmask;
        }
    }

    void DeleteOutgoing(char nnucl, bool as_is) {
        unsigned nmask = (1 << inv_position(nnucl, as_is));
        if (mask_ & nmask) {
#           pragma omp atomic
            mask_ &= (unsigned char) ~nmask;
        }
    }

    void DeleteIncoming(char pnucl, bool as_is) {
        unsigned pmask = (1 << inv_position(char(pnucl + 4), as_is));
        if (mask_ & pmask) {
#           pragma omp atomic
            mask_ &= (unsigned char) ~pmask;
        }
    }

    void IsolateVertex() {
        mask_ = 0;
    }

    bool CheckOutgoing(char nucl) const {
        return mask_ & (1 << nucl);
    }

    bool CheckIncoming(char nucl) const {
        return mask_ & (1 << (4 + nucl));
    }

    bool IsDeadEnd() const {
        return outgoing() == 0;
    }

    bool IsDeadStart() const {
        return incoming() == 0;
    }

    bool IsJunction() const {
        return !CheckUniqueOutgoing() || !CheckUniqueIncoming();
    }

    bool CheckUniqueOutgoing() const {
        return CheckUnique(outgoing());
    }

    bool CheckUniqueIncoming() const {
        return CheckUnique(incoming());
    }

    char GetUniqueOutgoing() const {
        return GetUnique(outgoing());
    }

    char GetUniqueIncoming() const {
        return GetUnique(incoming());
    }

    size_t OutgoingEdgeCount() const {
        return Count(outgoing());
    }

    size_t IncomingEdgeCount() const {
        return Count(incoming());
    }

    // We pass rhs by value (not by reference) below since InOutMask is just one byte
    InOutMask &operator|=(const InOutMask rhs) {
        mask_ |= rhs.mask_;
        return *this;
    }

    InOutMask &operator&=(const InOutMask rhs) {
        mask_ &= rhs.mask_;
        return *this;
    }

    InOutMask operator|(const InOutMask rhs) const {
        InOutMask result(*this);
        result |= rhs;
        return result;
    }

    InOutMask operator&(const InOutMask rhs) const {

        InOutMask result(*this);

        result &= rhs;
        return result;
    }
};

template<class Stream>
Stream &operator<<(Stream& stream, const InOutMask &mask) {
    return stream << std::bitset<8>(mask.get_mask());
}

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

    DeBruijnExtensionIndex(unsigned K)
            : base(K) {}

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
