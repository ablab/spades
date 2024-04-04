//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2015-2022 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <bitset>
#include <vector>
#include <cstdlib>
#include <cstdint>

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

    static constexpr int8_t GetUnique(uint8_t mask) {
        constexpr int8_t next[] =
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

}
