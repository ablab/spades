//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2015-2022 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "seq.hpp"
#include "rtseq.hpp"

#include "utils/verify.hpp"

#include <llvm/ADT/IntrusiveRefCntPtr.h>
#include <llvm/Support/TrailingObjects.h>

#include <vector>
#include <string>
#include <memory>
#include <cstring>

// Silence bogus gcc warnings
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"

namespace detail {
/// A smart pointer to a reference-counted object that inherits from
/// RefCountedBase or ThreadSafeRefCountedBase.
///
/// This class increments its pointee's reference count when it is created, and
/// decrements its refcount when it's destroyed (or is changed to point to a
/// different object).
template <typename T>
class RefCntPtr {
    T *Obj = nullptr;

public:
    using element_type = T;

    explicit RefCntPtr() = default;
    RefCntPtr(T *obj) : Obj(obj) { retain(); }
    RefCntPtr(const RefCntPtr &S) : Obj(S.Obj) { retain(); }
    RefCntPtr(RefCntPtr &&S) : Obj(S.Obj) { S.Obj = nullptr; }

    template <class X>
    RefCntPtr(RefCntPtr<X> &&S) : Obj(S.get()) {
        S.Obj = nullptr;
    }

    template <class X>
    RefCntPtr(const RefCntPtr<X> &S) : Obj(S.get()) {
        retain();
    }

    ~RefCntPtr() { release(); }

    RefCntPtr &operator=(RefCntPtr S) {
        swap(S);
        return *this;
    }

    T &operator*() const { return *get(); }
    T *operator->() const { return get(); }
    T *get() const { return Obj; }
    explicit operator bool() const { return get(); }

    void swap(RefCntPtr &other) {
        T *tmp = other.Obj;
        other.Obj = Obj;
        Obj = tmp;
    }

    void reset() {
        release();
        Obj = nullptr;
    }

    void resetWithoutRelease() { Obj = nullptr; }

private:
    void retain() {
        if (auto obj = get()) obj->Retain();
    }
    void release() {
        if (auto obj = get()) obj->Release();
    }

    template <typename X>
    friend class RefCntPtr;
};

template <class T, class U>
inline bool operator==(const RefCntPtr<T> &A, const RefCntPtr<U> &B) {
    return A.get() == B.get();
}

template <class T, class U>
inline bool operator!=(const RefCntPtr<T> &A, const RefCntPtr<U> &B) {
    return A.get() != B.get();
}

template <class T, class U>
inline bool operator==(const RefCntPtr<T> &A, U *B) {
    return A.get() == B;
}

template <class T, class U>
inline bool operator!=(const RefCntPtr<T> &A, U *B) {
    return A.get() != B;
}

template <class T, class U>
inline bool operator==(T *A, const RefCntPtr<U> &B) {
    return A == B.get();
}

template <class T, class U>
inline bool operator!=(T *A, const RefCntPtr<U> &B) {
    return A != B.get();
}

template <class T>
bool operator==(std::nullptr_t, const RefCntPtr<T> &B) {
    return !B;
}

template <class T>
bool operator==(const RefCntPtr<T> &A, std::nullptr_t B) {
    return B == A;
}

template <class T>
bool operator!=(std::nullptr_t A, const RefCntPtr<T> &B) {
    return !(A == B);
}

template <class T>
bool operator!=(const RefCntPtr<T> &A, std::nullptr_t B) {
    return !(A == B);
}

} // namespace detail


class Sequence {
    // Type to store Seq in Sequences
    typedef seq::seq_element_type ST;
    // Number of bits in ST
    static constexpr size_t STBits = sizeof(ST) << 3;
    // Number of nucleotides in ST
    static constexpr size_t STN = (STBits >> 1);
    // Number of bits in STN (for faster div and mod)
    static constexpr size_t STNBits = log_<STN, 2>::value;

    class ManagedNuclBuffer final : public llvm::ThreadSafeRefCountedBase<ManagedNuclBuffer>,
                                    protected llvm::TrailingObjects<ManagedNuclBuffer, ST> {
        friend TrailingObjects;

        ManagedNuclBuffer() {}

        ManagedNuclBuffer(size_t nucls, ST *buf) {
            std::uninitialized_copy(buf, buf + Sequence::data_size(nucls), data());
        }

      public:
        void operator delete(void *p) { ::operator delete(p); }

        static ManagedNuclBuffer *create(size_t nucls) {
            void *mem = ::operator new(totalSizeToAlloc<ST>(Sequence::data_size(nucls)));
            return new (mem) ManagedNuclBuffer();
        }

        static ManagedNuclBuffer *create(size_t nucls, ST *data) {
            void *mem = ::operator new(totalSizeToAlloc<ST>(Sequence::data_size(nucls)));
            return new (mem) ManagedNuclBuffer(nucls, data);
        }

        const ST *data() const { return getTrailingObjects<ST>(); }
        ST *data() { return getTrailingObjects<ST>(); }
    };

    // Historically we're having very unfortunate "mixed-endian" internal represenation of the nucl buffer.
    // Essentially, everything works in words (64-bit by default) and for the uneven length, the top
    // bits of the last word are zeroed. As a result, the buffer layout on litle-endian platform is as follows:
    //
    // |0000111122223333|000011112222....|
    //
    // Sequence itsels is essentially a pair of 64-bit values: first represents packed tuple of size, start offset
    // and reverse-complementary indicator, plus a pointer to a reference-counted nucl buffer.
    // In order to make short sequence optimization and reuse these 128-bits for the nucl buffer while keepeing
    // compatible buffer layout (lots of code expects word layout of the buffer) we essentially need to:
    //  - Reserve upper 8 bits of the buffer pointer for metadata: "is_short" indicator plus short size.
    //    Fortunately, the majority of the platforms operate in 48 or 56-bit address space, so upper byte could
    //    be (somehow) used
    //  - Ensure metadata is properly restored after potentially disrupting buffer copies
    //  - Short sequences are never shared and always copied as-is
    //  - Short sequences could only be constructed, we do not create them for substrings, etc. to keep
    //    these constant-time
    // The final layout is then as follows:
    //  - Long sequence:
    //     | size: 32, from: 31, rtl: 1 | refptr : 56, metadata: 8 |
    //                                       |
    //                                       +---> | refcounter |3333222211110000|...3222211110000|
    //  - Short sequence:
    //     |0000111122223333| 000011112222, metadata:8|
    // This way we can use buffer up to 120 bits for short strings, so we can store inline up to 60 nucls
    // The metadata is as follows:
    //  - [ size: 6, short_rtl : 1, is_short : 1 |
    // When is_short == 1, then it is expected that all other bits are zero, so the "is short" check is simply
    // check that the upper byte of the pointer is zero w/o any bitfield layout implications.

    static_assert(CHAR_BIT == 8, "This implementation assumes that one byte contains 8 bits");

    struct LongRep {
        size_t size_ : 32;
        size_t from_ : 31;
        bool   rtl_  : 1;  // Right to left + complimentary
        detail::RefCntPtr<ManagedNuclBuffer> data_;

        LongRep()
                : size_(0), from_(0), rtl_(0), data_(nullptr) {}

        LongRep(size_t size, size_t from, bool rtl, ManagedNuclBuffer *data)
                : size_(size), from_(from), rtl_(rtl), data_(data) {}

        bool rtl() const { return rtl_; }
        size_t from() const { return from_; }
        size_t size() const { return size_; }
    };

    enum { ShortSize = sizeof(LongRep) / sizeof(ST),
           MetadataOffset = (sizeof(ST) - 1)*8,
           MetadataMask = 0xFF,
           DataMask = (ST(1) << MetadataOffset) - 1,
           RTLOffset = 1, RTLMask = 0x1,
           SizeOffset = 2, SizeMask = 0x3F};

    static_assert(sizeof(ST) == 8, "This implementation assumes 64-bit underlying words");

    struct ShortRep {
        ST data[ShortSize];

        ShortRep() = default;

        uint8_t metadata() const {
            return (data[ShortSize - 1] >> MetadataOffset) & MetadataMask;
        }
        void set_metadata(size_t size, bool rtl) {
            uint8_t metadata = (size & SizeMask) << SizeOffset |
                               (rtl << RTLOffset) |
                               1;
            data[ShortSize - 1] = (data[ShortSize - 1] & DataMask) |
                                  ST(metadata) << MetadataOffset;
        }
        void set_metadata(uint8_t metadata) {
            data[ShortSize - 1] = (data[ShortSize - 1] & DataMask) |
                                  ST(metadata) << MetadataOffset;
        }

        bool rtl() const {
            return (metadata() >> RTLOffset) & RTLMask;
        }
        void set_rtl(bool rtl) {
            ST RTLFullMask = ST(RTLMask) << (MetadataOffset + RTLOffset);
            if (rtl)
                data[ShortSize - 1] |= RTLFullMask;
            else
                data[ShortSize - 1] &= ~RTLFullMask;
        }
        size_t size() const {
            return (metadata() >> SizeOffset) & SizeMask;
        }
        size_t from() const { return 0; }
    };

    // We borrow one byte for size, etc.
    static constexpr size_t ShortBits = (sizeof(ShortRep) - 1) << 3;
    static constexpr size_t ShortN = ShortBits >> 1;
    static constexpr size_t ShortNBits = log_<ShortN, 2>::value;

    static_assert(SizeMask >= ShortN, "inconsistent internal bitfield widths");

    #define ALWAYS_INLINE __attribute__((always_inline))
    union Rep {
        LongRep l;
        ShortRep s;

        // We need to be extremely careful and not invoke UB via
        // type-punning. Unfortunately, due to reasons explained above we cannot
        // access inactive union member, nor we might have a common leading
        // field sequence. As a result, the only escape hatch for us is the explicit
        // byte-based representation. This assumes little endian platform (!)
        ALWAYS_INLINE uint8_t metadata() const {
            const uint8_t *bytes = reinterpret_cast<const uint8_t*>(this);
            size_t offs = offsetof(ShortRep, data[ShortSize]);
            return bytes[offs - 1];
        }

        ALWAYS_INLINE bool is_short() const {
            return metadata() != 0;
        }

        ALWAYS_INLINE Rep() noexcept : l() {  }
        ALWAYS_INLINE ~Rep() noexcept { if (!is_short()) (&l)->LongRep::~LongRep(); }
        ALWAYS_INLINE Rep(const Rep &rhs) noexcept : l() {
            if (rhs.is_short())
                s = rhs.s;
            else
                l = rhs.l;
        }

        ALWAYS_INLINE Rep(Rep &&rhs) noexcept : l() {
            if (rhs.is_short())
                s = rhs.s;
            else
                l = std::move(rhs.l);
        }

        ALWAYS_INLINE Rep &operator=(const Rep &rhs) noexcept {
            if (&rhs == this)
                return *this;

            bool this_short =  is_short();
            if (rhs.is_short()) {
                if (!this_short) (&l)->LongRep::~LongRep();
                s = rhs.s;
            } else if (this_short)
                new(&l) LongRep(rhs.l);
            else // both this and rhs are long
                l = rhs.l;

            return *this;
        }
        ALWAYS_INLINE Rep &operator=(Rep &&rhs) noexcept {
            if (&rhs == this)
                return *this;

            bool this_short =  is_short();
            if (rhs.is_short()) {
                if (!this_short) (&l)->LongRep::~LongRep();
                s = rhs.s;
            } else if (this_short)
                new(&l) LongRep(std::move(rhs.l));
            else // both this and rhs are long
                l = std::move(rhs.l);

            return *this;
        }
    };

    Rep rep_;

    static constexpr bool is_short_rep(size_t size) { return size < ShortN; }

    static constexpr size_t data_size(size_t size) {
        return (size + STN - 1) >> STNBits;
    }
    size_t data_size() const {
        return (size() + STN - 1) >> STNBits;
    }

    ALWAYS_INLINE ST *data() {
        return is_short() ? rep_.s.data : rep_.l.data_->data();
    }

    ALWAYS_INLINE size_t from() const {
        return is_short() ? rep_.s.from() : rep_.l.from();
    }

    ALWAYS_INLINE bool rtl() const {
        return is_short() ? rep_.s.rtl() : rep_.l.rtl();
    }

    // NOTE: This function would destroy medata byte for short sequences
    // It should be "restored" by hands. Also, it is only allowed to look
    // over data as other fields are not expected to be populated yet.
    template<typename S>
    void InitFromNucls(const S &s, size_t size, bool rc) {
        size_t word_size = data_size(size);
        ST *words = data();

        VERIFY(is_dignucl(s[0]) || is_nucl(s[0]));

        // Which symbols does our string contain : 0123 or ACGT?
        bool digit_str = is_dignucl(s[0]);

        // data -- one temporary variable corresponding to the i-th array element
        // and some counters
        ST data = 0;
        size_t cnt = 0;
        size_t cur = 0;

        if (rc) {
            for (int i = (int) size - 1; i >= 0; --i) {
                //VERIFY(is_dignucl(s[i]) || is_nucl(s[i]));
                char c = complement(digit_str ? s[(unsigned) i] : dignucl(s[(unsigned) i]));

                data = data | (ST(c) << cnt);
                cnt += 2;

                if (cnt == STBits) {
                    words[cur++] = data;
                    cnt = 0;
                    data = 0;
                }
            }
        } else {
            for (size_t i = 0; i < size; ++i) {
                //VERIFY(is_dignucl(s[i]) || is_nucl(s[i]));
                char c = digit_str ? s[i] : dignucl(s[i]);

                data = data | (ST(c) << cnt);
                cnt += 2;

                if (cnt == STBits) {
                    words[cur++] = data;
                    cnt = 0;
                    data = 0;
                }
            }
        }

        if (cnt != 0)
            words[cur++] = data;

        for (; cur < word_size; ++cur)
            words[cur] = 0;
    }

    inline bool ReadHeader(std::istream &file);
    inline bool WriteHeader(std::ostream &file) const;

    uint8_t metadata() const {
        return rep_.metadata();
    }

    void set_metadata(uint8_t val) {
        return rep_.s.set_metadata(val);
    }

    Sequence(size_t size) {
        if (is_short_rep(size)) {
            rep_.s.set_metadata(size, /* rtl */ false);
        } else {
            new (&rep_.l) LongRep(size, 0, false, ManagedNuclBuffer::create(size));
        }
    }

    // Low level constructor. Handle with care.
    Sequence(const Sequence &seq, size_t from, size_t size, bool rtl) {
        // We are having few importants cases here:
        //   - If seq is a long sequence, we're just reusing buffer
        //   - Otherwise, we promote to long sequence if seq.from != from
        bool is_long = !seq.is_short();
        size_t seq_size = seq.size();
        if (is_long || seq.from() != from) {
            if (is_long) {
                rep_.l.size_ = size;
                rep_.l.from_ = from;
                rep_.l.rtl_  = rtl;
                rep_.l.data_ = seq.rep_.l.data_;
            } else {
                VERIFY(seq.from() == 0);
                new (&rep_.l) LongRep(size, from, rtl, ManagedNuclBuffer::create(seq_size));
                size_t dest_words = data_size(seq_size);
                ST *dest = rep_.l.data_->data();
                memcpy(dest, seq.data(), dest_words * sizeof(ST));
                // Ensure the metadata byte is cleared if we copied the last
                // word
                if (dest_words == ShortSize) dest[ShortSize - 1] &= DataMask;
            }
        } else {
            // Now we know that seq is short with same from / size
            // Copy the data and set rtl flag
            rep_.s = seq.rep_.s;
            rep_.s.set_metadata(size, rtl);
        }
    }

public:
    /**
     * Sequence initialization (arbitrary size string)
     *
     * @param s ACGT or 0123-string
     */
    explicit Sequence(const char *s, bool rc = false)
            : Sequence(strlen(s)) {
        uint8_t m = metadata();
        InitFromNucls(s, size(), rc);
        if (m) set_metadata(m);
    }

    explicit Sequence(char *s, bool rc = false)
            : Sequence(strlen(s)) {
        uint8_t m = metadata();
        InitFromNucls(s, size(), rc);
        if (m) set_metadata(m);
    }

    template<typename S>
    explicit Sequence(const S &s, bool rc = false)
            : Sequence(s.size()) {
        uint8_t m = metadata();
        InitFromNucls(s, size(), rc);
        if (m) set_metadata(m);
    }

    Sequence()
            : Sequence(size_t(0)) {
        uint8_t m = metadata();
        memset(data(), 0, data_size() * sizeof(ST));
        if (m) set_metadata(m);
    }

    template<size_t size2_>
    explicit Sequence(const Seq<size2_> &kmer, size_t)
            : Sequence(kmer.size()) {
        uint8_t m = metadata();
        kmer.copy_data(data());
        if (m) set_metadata(m);
    }

    template<size_t size2_>
    explicit Sequence(const RuntimeSeq<size2_> &kmer, size_t)
            : Sequence(kmer.size()) {
        uint8_t m = metadata();
        kmer.copy_data(data());
        if (m) set_metadata(m);
    }

    Sequence(const Sequence &s)
            : Sequence(s, s.from(), s.size(), s.rtl()) { }

    Sequence(Sequence &&) noexcept = default;
    Sequence &operator=(const Sequence &rhs) noexcept = default;
    Sequence &operator=(Sequence &&) noexcept = default;

    ALWAYS_INLINE bool is_short() const { return rep_.is_short(); }

    #undef ALWAYS_INLINE

    Sequence operator!() const {
        return Sequence(*this, from(), size(), !rtl());
    }

    char operator[](const size_t index) const {
        VERIFY_DEV(index < size());
        const ST *words = data();
        if (rtl()) {
            size_t i = from() + size() - 1 - index;
            return complement((words[i >> STNBits] >> ((i & (STN - 1)) << 1)) & 3);
        } else {
            size_t i = from() + index;
            return (words[i >> STNBits] >> ((i & (STN - 1)) << 1)) & 3;
        }
    }

    const ST *data() const {
        return is_short() ? rep_.s.data : rep_.l.data_->data();
    }

    bool operator==(const Sequence &that) const {
        if (size() != that.size())
            return false;

        // If both sequences are long we need to check if they share the same
        // buffer and have same parameters
        if (!is_short() && !that.is_short()) {
            if (rep_.l.data_ == that.rep_.l.data_ &&
                rep_.l.from_ == that.rep_.l.from_ &&
                rep_.l.rtl_  == that.rep_.l.rtl_)
                return true;
        }

        for (size_t i = 0; i < size(); ++i) {
            if (this->operator[](i) != that[i]) {
                return false;
            }
        }
        return true;
    }

    bool operator!=(const Sequence &that) const {
        return !(operator==(that));
    }

    /**
     * @todo Might be optimized via int comparison (not so easy)
     */
    bool operator<(const Sequence &that) const {
        size_t s = std::min(size(), that.size());
        for (size_t i = 0; i < s; ++i) {
            if (this->operator[](i) != that[i]) {
                return (this->operator[](i) < that[i]);
            }
        }
        return (size() < that.size());
    }

    inline Sequence operator<<(char c) const;

    /**
     * @param from inclusive
     * @param to exclusive;
     */
    inline Sequence Subseq(size_t from, size_t to) const;
    inline Sequence Subseq(size_t from) const; // up to size_ by default

    inline Sequence operator+(const Sequence &s) const;

    inline Sequence First(size_t count) const;
    inline Sequence Last(size_t count) const;

    inline size_t find(const Sequence &t, size_t from = 0) const;

    template<size_t size2_>
    Seq<size2_> start() const;
    template<size_t size2_>
    Seq<size2_> end() const;

    template<class Seq>
    Seq start(size_t k) const;
    template<class Seq>
    Seq end(size_t k) const;

    inline std::string str() const;
    inline std::string err() const;

    size_t size() const {
        return is_short() ? rep_.s.size() : rep_.l.size();
    }

    bool empty() const {
        return size() == 0;
    }

    template<class Seq>
    bool contains(const Seq& s, size_t offset = 0) const {
        VERIFY_DEV(offset + s.size() <= size());

        for (size_t i = 0, e = s.size(); i != e; ++i)
            if (operator[](offset + i) != s[i])
                return false;

        return true;
    }

    static bool RawCompare(const Sequence &s1, const Sequence &s2) {
        VERIFY_DEV(s1.from() == 0);
        VERIFY_DEV(s2.from() == 0);
        VERIFY_DEV(!s1.rtl());
        VERIFY_DEV(!s2.rtl());

        if (s1.size() != s2.size()) {
            return s1.size() < s2.size();
        }

        size_t szwords = s1.data_size();
        const ST *d1 = s1.data();
        const ST *d2 = s2.data();
        for (size_t i = 0; i < szwords; ++i) {
            if (d1[i] != d2[i]) {
                return d1[i] < d2[i];
            }
        }
        return false;  // They are equal
    }

    inline bool BinRead(std::istream &file);
    inline bool BinWrite(std::ostream &file) const;
};

inline std::ostream &operator<<(std::ostream &os, const Sequence &s);

/**
 * start of Sequence is Seq with preferred size
 */
template<size_t size2_>
Seq<size2_> Sequence::start() const {
    //VERIFY(size2_ <= size_);
    return Seq<size2_>(*this);
}

template<size_t size2_>
Seq<size2_> Sequence::end() const {
    return Seq<size2_>(*this, size() - size2_);
}


template<class Seq>
Seq Sequence::start(size_t k) const {
    return Seq(unsigned(k), *this);
}

template<class Seq>
Seq Sequence::end(size_t k) const {
    return Seq(unsigned(k), *this, size() - k);
}

// O(1)
//including from, excluding to
//safe if not #DEFINE NDEBUG
Sequence Sequence::Subseq(size_t from, size_t to) const {
    VERIFY(to >= from);
    VERIFY(to <= size());
    //VERIFY(to - from <= size_);
    if (rtl()) {
        return Sequence(*this, this->from() + size() - to, to - from, true);
    } else {
        return Sequence(*this, this->from() + from, to - from, false);
    }
}

//including from, excluding to
Sequence Sequence::Subseq(size_t from) const {
    return Subseq(from, size());
}

Sequence Sequence::First(size_t count) const {
    return Subseq(0, count);
}

Sequence Sequence::Last(size_t count) const {
    return Subseq(size() - count);
}

/**
 * @todo : must be KMP or hashing instead of this
 */
size_t Sequence::find(const Sequence &t, size_t from) const {
    for (size_t i = from; i <= size() - t.size(); i++) {
        if (Subseq(i, i + t.size()) == t) {
            return i;
        }
    }
    return -1ULL;
}

/**
 * @todo optimize
 */
Sequence Sequence::operator+(const Sequence &s) const {
    return Sequence(str() + s.str());
    // TODO might be opposite to correct
    //    int total = size_ + s.size_;
    //    std::vector<Seq<4> > bytes((total + 3) >> 2);
    //    for (size_t i = 0; i < size_; ++i) {
    //        bytes[i / 4] = (bytes[i / 4] << operator [](i)); // TODO :-) use <<=
    //    }
    //    for (size_t i = 0, j = size_; i < s.size_; ++i, ++j) {
    //        bytes[j / 4] = (bytes[j / 4]) << s[i];
    //    }
    //    return Sequence(new Data(bytes), 0, total, false);
}

std::string Sequence::str() const {
    std::string res(size(), '-');
    for (size_t i = 0; i < size(); ++i) {
        res[i] = nucl2(this->operator[](i));
    }
    return res;
}

std::string Sequence::err() const {
    std::ostringstream oss;
    oss << (is_short() ? "{ short seq " : "{ long seq ");
    oss << "*data=" << data() <<
            ", from_=" << from() <<
            ", size_=" << size() <<
            ", rtl_=" << int(rtl()) << " }";
    return oss.str();
}

std::ostream &operator<<(std::ostream &os, const Sequence &s) {
    os << s.str();
    return os;
}

bool Sequence::ReadHeader(std::istream &file) {
    size_t size;
    file.read((char *) &size, sizeof(size));

    if (!is_short())
        rep_.l.data_.reset();

    if (is_short_rep(size)) {
        rep_.s.set_metadata(size, false);
    } else {
        new (&rep_.l) LongRep(size, 0, false, ManagedNuclBuffer::create(size));
    }

    return !file.fail();
}

bool Sequence::WriteHeader(std::ostream &file) const {
    VERIFY(from() == 0);
    VERIFY(!rtl());

    size_t size = this->size();
    file.write((const char *)&size, sizeof(size));

    return !file.fail();
}


bool Sequence::BinRead(std::istream &file) {
    ReadHeader(file);

    uint8_t m = metadata();
    file.read((char *) data(), data_size() * sizeof(ST));
    if (m) set_metadata(m);

    return !file.fail();
}


bool Sequence::BinWrite(std::ostream &file) const {
    if (from() != 0 || rtl()) {
        Sequence clear(this->str());
        return clear.BinWrite(file);
    }

    WriteHeader(file);

    file.write((const char *) data(), data_size() * sizeof(ST));

    return !file.fail();
}

/**
 * @class SequenceBuilder
 * @section DESCRIPTION
 *
 * Class was created for build sequence. It is included method: size(), append()
 */

class SequenceBuilder {
    std::vector<char> buf_;
public:
    template<typename S>
    SequenceBuilder &append(const S &s) {
        for (size_t i = 0; i < s.size(); ++i) {
            buf_.push_back(s[i]);
        }
        return *this;
    }

    SequenceBuilder &append(char c) {
        buf_.push_back(c);
        return *this;
    }

    Sequence BuildSequence() {
        return Sequence(buf_);
    }

    size_t size() const {
        return buf_.size();
    }

    void clear() {
        return buf_.clear();
    }

    char operator[](const size_t index) const {
        VERIFY_DEV(index < buf_.size());
        return buf_[index];
    }

    std::string str() const {
        std::string s(buf_.size(), '-');
        for (size_t i = 0; i < s.size(); ++i) {
            s[i] = nucl(buf_[i]);
        }
        return s;
    }
};

#pragma GCC diagnostic pop
