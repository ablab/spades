//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <vector>
#include <string>
#include <memory>
#include <cstring>

#include "seq.hpp"
#include "rtseq.hpp"

#include <llvm/ADT/IntrusiveRefCntPtr.h>
#include <llvm/Support/TrailingObjects.h>

class Sequence {
    // Type to store Seq in Sequences
    typedef seq_element_type ST;
    // Number of bits in ST
    const static size_t STBits = sizeof(ST) << 3;
    // Number of nucleotides in ST
    const static size_t STN = (STBits >> 1);
    // Number of bits in STN (for faster div and mod)
    const static size_t STNBits = log_<STN, 2>::value;

    class ManagedNuclBuffer final : public llvm::ThreadSafeRefCountedBase<ManagedNuclBuffer>,
                                    protected llvm::TrailingObjects<ManagedNuclBuffer, ST> {
        friend TrailingObjects;

        ManagedNuclBuffer() {}

        ManagedNuclBuffer(size_t nucls, ST *buf) {
            std::uninitialized_copy(buf, buf + Sequence::DataSize(nucls), data());
        }

      public:
        static ManagedNuclBuffer *create(size_t nucls) {
            void *mem = ::operator new(totalSizeToAlloc<ST>(Sequence::DataSize(nucls)));
            return new (mem) ManagedNuclBuffer();
        }

        static ManagedNuclBuffer *create(size_t nucls, ST *data) {
            void *mem = ::operator new(totalSizeToAlloc<ST>(Sequence::DataSize(nucls)));
            return new (mem) ManagedNuclBuffer(nucls, data);
        }

        const ST *data() const { return getTrailingObjects<ST>(); }
        ST *data() { return getTrailingObjects<ST>(); }
    };

    size_t from_;
    size_t size_;
    bool rtl_; // Right to left + complimentary (?)
    llvm::IntrusiveRefCntPtr<ManagedNuclBuffer> data_;

    static size_t DataSize(size_t size) {
        return (size + STN - 1) >> STNBits;
    }

    template<typename S>
    void InitFromNucls(const S &s, bool rc = false) {
        size_t bytes_size = DataSize(size_);
        ST *bytes = data_->data();

        VERIFY(is_dignucl(s[0]) || is_nucl(s[0]));

        // Which symbols does our string contain : 0123 or ACGT?
        bool digit_str = is_dignucl(s[0]);

        // data -- one temporary variable corresponding to the i-th array element
        // and some counters
        ST data = 0;
        size_t cnt = 0;
        size_t cur = 0;

        if (rc) {
            for (int i = (int) size_ - 1; i >= 0; --i) {
                //VERIFY(is_dignucl(s[i]) || is_nucl(s[i]));
                char c = complement(digit_str ? s[(unsigned) i] : dignucl(s[(unsigned) i]));

                data = data | (ST(c) << cnt);
                cnt += 2;

                if (cnt == STBits) {
                    bytes[cur++] = data;
                    cnt = 0;
                    data = 0;
                }
            }
        } else {
            for (size_t i = 0; i < size_; ++i) {
                //VERIFY(is_dignucl(s[i]) || is_nucl(s[i]));
                char c = digit_str ? s[i] : dignucl(s[i]);

                data = data | (ST(c) << cnt);
                cnt += 2;

                if (cnt == STBits) {
                    bytes[cur++] = data;
                    cnt = 0;
                    data = 0;
                }
            }
        }

        if (cnt != 0)
            bytes[cur++] = data;

        for (; cur < bytes_size; ++cur)
            bytes[cur] = 0;
    }

    inline bool ReadHeader(std::istream &file);
    inline bool WriteHeader(std::ostream &file) const;

    Sequence(size_t size, int)
            : from_(0), size_(size), rtl_(false), data_(ManagedNuclBuffer::create(size_)) {}

    //Low level constructor. Handle with care.
    Sequence(const Sequence &seq, size_t from, size_t size, bool rtl)
            : from_(from), size_(size), rtl_(rtl), data_(seq.data_) {}

public:
    /**
     * Sequence initialization (arbitrary size string)
     *
     * @param s ACGT or 0123-string
     */
    explicit Sequence(const char *s, bool rc = false)
            : Sequence(strlen(s), 0) {
        InitFromNucls(s, rc);
    }

    explicit Sequence(char *s, bool rc = false)
            : Sequence(strlen(s), 0) {
        InitFromNucls(s, rc);
    }

    template<typename S>
    explicit Sequence(const S &s, bool rc = false)
            : Sequence(s.size(), 0) {
        InitFromNucls(s, rc);
    }

    Sequence()
            : Sequence(size_t(0), 0) {
        memset(data_->data(), 0, DataSize(size_));
    }

    template<size_t size2_>
    explicit Sequence(const Seq<size2_> &kmer, size_t)
            : Sequence(kmer.size(), 0) {
        kmer.copy_data(data_->data());
    }

    template<size_t size2_>
    explicit Sequence(const RuntimeSeq<size2_> &kmer, size_t)
            : Sequence(kmer.size(), 0) {
        kmer.copy_data(data_->data());
    }

    Sequence(const Sequence &s)
            : Sequence(s, s.from_, s.size_, s.rtl_) {}

    const Sequence &operator=(const Sequence &rhs) {
        if (&rhs == this)
            return *this;

        from_ = rhs.from_;
        size_ = rhs.size_;
        rtl_ = rhs.rtl_;
        data_ = rhs.data_;

        return *this;
    }

    char operator[](const size_t index) const {
        VERIFY_DEV(index < size_);
        const ST *bytes = data_->data();
        if (rtl_) {
            size_t i = from_ + size_ - 1 - index;
            return complement((bytes[i >> STNBits] >> ((i & (STN - 1)) << 1)) & 3);
        } else {
            size_t i = from_ + index;
            return (bytes[i >> STNBits] >> ((i & (STN - 1)) << 1)) & 3;
        }
    }

    bool operator==(const Sequence &that) const {
        if (size_ != that.size_)
            return false;

        if (data_ == that.data_ && from_ == that.from_ && rtl_ == that.rtl_)
            return true;

        for (size_t i = 0; i < size_; ++i) {
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
        size_t s = std::min(size_, that.size_);
        for (size_t i = 0; i < s; ++i) {
            if (this->operator[](i) != that[i]) {
                return (this->operator[](i) < that[i]);
            }
        }
        return (size_ < that.size_);
    }

    Sequence operator!() const {
        return Sequence(*this, from_, size_, !rtl_);
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
        return size_;
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

public:
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
    return Seq<size2_>(*this, size_ - size2_);
}


template<class Seq>
Seq Sequence::start(size_t k) const {
    return Seq(unsigned(k), *this);
}

template<class Seq>
Seq Sequence::end(size_t k) const {
    return Seq(unsigned(k), *this, size_ - k);
}

// O(1)
//including from, excluding to
//safe if not #DEFINE NDEBUG
Sequence Sequence::Subseq(size_t from, size_t to) const {
    VERIFY(to >= from);
    VERIFY(to <= size_);
    //VERIFY(to - from <= size_);
    if (rtl_) {
        return Sequence(*this, from_ + size_ - to, to - from, true);
    } else {
        return Sequence(*this, from_ + from, to - from, false);
    }
}

//including from, excluding to
Sequence Sequence::Subseq(size_t from) const {
    return Subseq(from, size_);
}

Sequence Sequence::First(size_t count) const {
    return Subseq(0, count);
}

Sequence Sequence::Last(size_t count) const {
    return Subseq(size_ - count);
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
    std::string res(size_, '-');
    for (size_t i = 0; i < size_; ++i) {
        res[i] = nucl(this->operator[](i));
    }
    return res;
}

std::string Sequence::err() const {
    std::ostringstream oss;
    oss << "{ *data=" << data_->data() <<
            ", from_=" << from_ <<
            ", size_=" << size_ <<
            ", rtl_=" << int(rtl_) << " }";
    return oss.str();
}

std::ostream &operator<<(std::ostream &os, const Sequence &s) {
    os << s.str();
    return os;
}

bool Sequence::ReadHeader(std::istream &file) {
    file.read((char *) &size_, sizeof(size_));

    from_ = 0;
    rtl_ = false;

    return !file.fail();
}

bool Sequence::WriteHeader(std::ostream &file) const {
    VERIFY(from_ == 0);
    VERIFY(!rtl_);

    file.write((const char *) &size_, sizeof(size_));

    return !file.fail();
}


bool Sequence::BinRead(std::istream &file) {
    ReadHeader(file);

    data_ = llvm::IntrusiveRefCntPtr<ManagedNuclBuffer>(ManagedNuclBuffer::create(size_));
    file.read((char *) data_->data(), DataSize(size_) * sizeof(ST));

    return !file.fail();
}


bool Sequence::BinWrite(std::ostream &file) const {
    if (from_ != 0 || rtl_) {
        Sequence clear(this->str());
        return clear.BinWrite(file);
    }

    WriteHeader(file);

    file.write((const char *) data_->data(), DataSize(size_) * sizeof(ST));

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
