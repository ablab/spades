//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

/*
 * rtseq.hpp
 *
 *  Created on: Jun 28, 2012
 *      Author: andrey
 */

#ifndef RTSEQ_HPP_
#define RTSEQ_HPP_

#include <string>
#include "utils/verify.hpp"
#include <array>
#include <algorithm>
#include "nucl.hpp"
#include "math/log.hpp"
#include "seq_common.hpp"
#include "seq.hpp"
#include "simple_seq.hpp"

#include <cstring>
#include <iostream>

template<size_t max_size_, typename T = seq_element_type>
class RuntimeSeq {
public:
    /**
     * @variable Number of bits in type T (e.g. 8 for char)
     * @example 8: 2^8 = 256 or 16
     */
    const static size_t TBits = sizeof(T) << 3;

    /**
     * @variable Number of nucleotides that can be stored in one type T (e.g. 4 for char)
     * TNucl MUST be a power of two
     * @example 4: 8/2 = 4 or 16/2 = 8
     */
    const static size_t TNucl = TBits >> 1;

    /**
     * @variable Number of bits in TNucl (e.g. 2 for char). Useful for shifts instead of divisions.
     */
    const static size_t TNuclBits = log_<TNucl, 2>::value;

    const static size_t Iterations = log_<TBits, 2>::value;

    static const std::array<T, Iterations> ConstructLeftMasks() {
        std::array<T, Iterations> result;
        for (size_t i = 0; i < Iterations; i++) {
            size_t shift = 1 << i;
            T mask = T(T(1) << shift) - T(1);
            result[i] = T(mask << shift);
            for (size_t j = 0; j < i; j++) {
                result[j] += T(result[j] << shift);
            }
        }
        return result;
    }

    static const std::array<T, Iterations> ConstructRightMasks() {
        std::array<T, Iterations> result(ConstructLeftMasks());
        for (size_t i = 0; i < Iterations; i++) {
            result[i] = T(~result[i]);
        }
        return result;
    }


    RuntimeSeq<max_size_, T> FastRC() const {
        const static std::array<T, Iterations> LeftMasks(ConstructLeftMasks());
        const static std::array<T, Iterations> RightMasks(ConstructRightMasks());
        const static size_t LogTSize = log_<sizeof(T), 2>::value + 3;

        RuntimeSeq<max_size_, T> res(this->size());

        const size_t bit_size = size_ << 1;
        const size_t extra = bit_size & ((1 << LogTSize) - 1);
        const size_t to_extra = TBits - extra;
        const size_t filled = bit_size >> LogTSize;
        size_t real_length = filled;
        if (extra == 0) {
            for (size_t i = 0, j = filled - 1; i < filled; i++, j--) {
                res.data_[i] = data_[j];
            }
        } else {
            for (size_t i = 0, j = filled; i < filled && j > 0; i++, j--) {
                res.data_[i] = (data_[j] << to_extra) + (data_[j - 1] >> extra);
            }
            res.data_[filled] = (data_[0] << to_extra);
            real_length++;
        }

        for (size_t i = 0; i < real_length; i++) {
            res.data_[i] = res.data_[i] ^ T(-1);
            for (size_t it = 1; it < Iterations; it++) {
                size_t shift = 1 << it;
                res.data_[i] = T((res.data_[i] & LeftMasks[it]) >> shift) ^ T((res.data_[i] & RightMasks[it]) << shift);
            }
        }

        if (extra != 0) {
            res.data_[real_length - 1] = (res.data_[real_length - 1] & ((T(1) << extra) - 1));
        }
        return res;
    }

    /**
     * @variable Number of Ts which required to store all sequence.
     */
    const static size_t DataSize = (max_size_ + TNucl - 1) >> TNuclBits;

    /**
     * @variable Number of meaningful bytes in whick seq is stored
     */
    const static size_t TotalBytes = sizeof(T) * DataSize;

    typedef T DataType;

    static size_t GetDataSize(size_t size) {
        return (size + TNucl - 1) >> TNuclBits;
    }

private:
    /* *
     * @variable Just some prime number to count the hash function of the kmer
     * */
    const static size_t PrimeNum = 239;


    // number of nucleotides in the last data_ bucket
    static size_t NuclsRemain(size_t size) {
        return size & (TNucl - 1);
    }

    // useful mask to fill the last element of the data_ array
    static size_t MaskForLastBucket(size_t size) {
        size_t nr = NuclsRemain(size);
        return nr != 0 ? (((T) 1) << (nr << 1)) - 1 : -1ul;
    }


    /**
     * @variable Inner representation of sequence: array of Ts with length = DataSize.
     *
     * @invariant Invariant: all nucleotides >= size_ are 'A's (useful for comparison)
     */
    std::array<T, DataSize> data_;

    size_t size_;

    /**
     * Initialize data_ array of this object with C-string
     *
     * @param s C-string (ACGT chars only), strlen(s) = size_
     */
    void init(const char *s) {
        T data = 0;
        size_t cnt = 0;
        size_t cur = 0;
        for (size_t pos = 0; pos < size_; ++pos, ++s) { // unsafe!
            // VERIFY(is_nucl(*s)); // for performance
            data = data | ((T) dignucl(*s) << cnt);
            cnt += 2;
            if (cnt == TBits) {
                this->data_[cur++] = data;
                cnt = 0;
                data = 0;
            }
        }
        if (cnt != 0) {
            this->data_[cur++] = data;
        }

        for (; cur < DataSize; ++cur)
            this->data_[cur] = 0;

        VERIFY(*s == 0); // C-string always ends on 0
    }

    /**
     * Sets i-th symbol of Seq with 0123-char
     */
    inline void set(const size_t i, char c) {
        data_[i >> TNuclBits] =
                (data_[i >> TNuclBits] & ~((T) 3 << ((i & (TNucl - 1)) << 1))) | ((T) c << ((i & (TNucl - 1)) << 1));
    }

    // Template voodoo to calculate the length of the string regardless whether it is std::string or const char*
    template<class S>
    size_t size(const S &t,
                typename std::enable_if<std::is_class<S>::value, T>::type * = 0) {
        return t.size();
    }

    template<class S>
    size_t size(const S &t,
                typename std::enable_if<std::is_same<S, const char *>::value, T>::type * = 0) {
        return strlen(t);
    }


public:

    const static size_t max_size = max_size_;

    RuntimeSeq() : size_(0) {
        std::fill(data_.begin(), data_.end(), 0);
    }

    /**
     * Default constructor, fills Seq with A's
     */

    explicit RuntimeSeq(size_t k) : size_(k) {
        VERIFY_DEV(k <= max_size_);
        //VERIFY((T)(-1) >= (T)0);//be sure to use unsigned types
        std::fill(data_.begin(), data_.end(), 0);
    }

    RuntimeSeq(size_t k, const char *s) : size_(k) {
        VERIFY_DEV(k <= max_size_);
        //VERIFY((T)(-1) >= (T)0);//be sure to use unsigned types
        init(s);
    }


    explicit RuntimeSeq(size_t k, const T *data_array) : size_(k) {
        VERIFY_DEV(k <= max_size_);
        std::fill(data_.begin(), data_.end(), 0);

        size_t data_size = GetDataSize(size_);
        memcpy(data_.data(), data_array, data_size * sizeof(T));

        if (NuclsRemain(size_)) {
            data_[data_size - 1] = data_[data_size - 1] & MaskForLastBucket(size_);
        }
    }

    explicit RuntimeSeq(size_t k, T *data_array) : size_(k) {
        VERIFY_DEV(k <= max_size_);
        std::fill(data_.begin(), data_.end(), 0);

        size_t data_size = GetDataSize(size_);
        memcpy(data_.data(), data_array, data_size * sizeof(T));

        if (NuclsRemain(size_)) {
            data_[data_size - 1] = data_[data_size - 1] & MaskForLastBucket(size_);
        }
    }

    explicit RuntimeSeq(size_t k, const RuntimeSeq &seq)
            : RuntimeSeq(k, seq.data_.data()) {}

    template<size_t size2_, typename T2 = T>
    explicit RuntimeSeq(const Seq<size2_, T2> &seq, bool) : size_(size2_) {
        VERIFY_DEV(size_ <= max_size_);
        std::fill(data_.begin(), data_.end(), 0);
        seq.copy_data(data_.data());
    }

    template<size_t size2_, typename T2 = T>
    explicit RuntimeSeq(const SimpleSeq<size2_, T2> &seq, size_t k) : size_(k) {
        VERIFY_DEV(size_ <= max_size_);
        VERIFY_DEV(size2_ <= max_size_);
        std::fill(data_.begin(), data_.end(), 0);
        seq.copy_data(data_.data());
    }


    /**
     * Ultimate constructor from ACGT0123-string.
     *
     * @param s Any object with operator[], which returns 0123 chars
     * @param offset Offset when this sequence starts
     * @number_to_read A number of nucleotides, we want to fetch from this string
     * @warning assuming that s is a correct string, filled with ACGT _OR_ 0123
     * no init method, filling right here
     */
    template<typename S>
    explicit RuntimeSeq(size_t k, const S &s, size_t offset = 0) : size_(k) {
        VERIFY(size_ <= max_size_);
        //TRACE("New Constructor for seq " << s[0] << " is first symbol");
        VERIFY(size_ == 0 || is_dignucl(s[0]) || is_nucl(s[0]));
        VERIFY(offset + size_ <= this->size(s));

        // which symbols does our string contain : 0123 or ACGT?
        bool digit_str = size_ == 0 || is_dignucl(s[0]);

        // we fill everything with zeros (As) by default.
        std::fill(data_.begin(), data_.end(), 0);

        // data -- one temporary variable corresponding to the i-th array element
        // and some counters
        T data = 0;
        size_t cnt = 0;
        size_t cur = 0;

        for (size_t i = 0; i < size_; ++i) {
            //VERIFY(is_dignucl(s[i]) || is_nucl(s[i]));

            // we fill everything with zeros (As) by default.
            char c = (char) (digit_str ? s[offset + i] : dignucl(s[offset + i]));

            data = data | (T(c) << cnt);
            cnt += 2;

            if (cnt == TBits) {
                this->data_[cur++] = data;
                cnt = 0;
                data = 0;
            }
        }

        if (cnt != 0) {
            this->data_[cur++] = data;
        }

        for (; cur < DataSize; ++cur)
            this->data_[cur] = 0;
    }

    RuntimeSeq start(size_t K) const {
        return RuntimeSeq(K, data_.data());
    }
    
    /**
     *  Reads sequence from the file (in the same format as BinWrite writes it)
     *  and returns false if error occured, true otherwise.
     */
    bool BinRead(std::istream &file) {
        file.read((char *) data_.data(), sizeof(T) * GetDataSize(size_));
        return !file.fail();
    }

    /**
     *  Writes sequence to the file (in the same format as BinRead reads it)
     *  and returns false if error occured, true otherwise.
     */
    bool BinWrite(std::ostream &file) const {
        file.write((const char *) data_.data(), sizeof(T) * GetDataSize(size_));
        return !file.fail();
    }

    /**
     *  Reads sequence from the file (in the same format as BinWrite writes it)
     *  and returns false if error occured, true otherwise.
     */
    static bool BinRead(std::istream &file, RuntimeSeq<max_size_, T> *seq) {
        return seq->BinRead(file);
    }

    /**
     *  Writes sequence to the file (in the same format as BinRead reads it)
     *  and returns false if error occured, true otherwise.
     */
    static bool BinWrite(std::ostream &file, const RuntimeSeq<max_size_, T> &seq) {
        return seq.BinWrite(file);
    }


    /**
     * Get i-th symbol of Seq.
     *
     * @param i Index of the symbol (0 <= i < size_)
     * @return 0123-char on position i
     */
    char operator[](const size_t i) const {
        VERIFY_DEV(i < size_);
        return (data_[i >> TNuclBits] >> ((i & (TNucl - 1)) << 1)) & 3;
    }

    /**::
     * Reverse complement.
     *
     * @return Reverse complement Seq.
     */
    RuntimeSeq<max_size_, T> operator!() const {
//    RuntimeSeq<max_size_, T> res(*this);
//    for (size_t i = 0; i < (size_ >> 1); ++i) {
//      auto front = complement(res[i]);
//      auto end = complement(res[size_ - 1 - i]);
//      res.set(i, end);
//      res.set(size_ - 1 - i, front);
//    }
//    if ((size_ & 1) == 1) {
//      res.set(size_ >> 1, complement(res[size_ >> 1]));
//    }
        return FastRC();
//    return res;
    }

    /**
     * Is the kmer minimal among this and !this.
     *
     * @return True if kmer < !kmer and false otherwise.
     */
    bool IsMinimal() const {
        for (size_t i = 0; (i << 1) + 1 <= size_; ++i) {
            auto front = this->operator[](i);
            auto end = complement(this->operator[](size_ - 1 - i));
            if (front != end)
                return front < end;
        }
        return true;
    }

    /**
     * Shift left
     *
     * @param c New 0123 char which should be added to the right.
     * @return Shifted (to the left) sequence with 'c' char on the right.
     */
    RuntimeSeq<max_size_, T> operator<<(char c) const {
        if (is_nucl(c)) {
            c = dignucl(c);
        }

        RuntimeSeq<max_size_, T> res(*this);
        std::array<T, DataSize> &data = res.data_;

        size_t data_size = GetDataSize(size_);

        if (data_size != 0) { // unless empty sequence
            T rm = data[data_size - 1] & 3;
            T lastnuclshift_ = ((size_ + TNucl - 1) & (TNucl - 1)) << 1;
            data[data_size - 1] = (data[data_size - 1] >> 2) | ((T) c << lastnuclshift_);

            if (data_size >= 2) { // if we have at least 2 elements in data
                for (int i = (int) data_size - 2; i >= 0; --i) {
                    T new_rm = data[i] & 3;
                    data[i] = (data[i] >> 2) |
                              (rm << (TBits - 2)); // we need & here because if we shift negative, it fill with ones :(
                    rm = new_rm;
                }
            }
        }
        return res;
    }

    void operator<<=(char c) {
        if (is_nucl(c)) {
            c = dignucl(c);
        }

        size_t data_size = GetDataSize(size_);

        if (data_size == 0) {
            return;
        }

        for (size_t i = 0; i < data_size - 1; ++i) {
            data_[i] = (data_[i] >> 2) | (((T) data_[i + 1] & 3) << (TBits - 2));
        }

        T lastnuclshift_ = ((size_ + TNucl - 1) & (TNucl - 1)) << 1;
        data_[data_size - 1] = (data_[data_size - 1] >> 2) | ((T) c << lastnuclshift_);
    }

//todo naming convention violation!
    RuntimeSeq<max_size_, T> pushBack(char c) const {
        //VERIFY(size_ + 1 <= max_size_);

        if (is_nucl(c)) {
            c = dignucl(c);
        }
        //VERIFY(is_dignucl(c));
        RuntimeSeq<max_size_, T> s(size_ + 1);
        std::copy(this->data_.begin(), this->data_.end(), s.data_.begin());

        size_t data_size = GetDataSize(size_ + 1);

        s.data_[data_size - 1] |= ((T) c << ((size_ & (TNucl - 1)) << 1));

        return s; //was: Seq<size_ + 1, T>(str() + nucl(c));
    }


//todo naming convention violation!
    void pushBackThis(char c) {
        VERIFY(size_ + 1 <= max_size_);

        if (is_nucl(c)) {
            c = dignucl(c);
        }

        size_ += 1;
        size_t data_size = GetDataSize(size_);

        data_[data_size - 1] |= ((T) c << (((size_ - 1) & (TNucl - 1)) << 1));
    }

    //    /**
    //     * @todo optimize!!!
    //     */
    //    RuntimeSeq<max_size_, T> pushFront(char c) const {
    //        VERIFY(size_ + 1 < max_size_);
    //        if (is_nucl(c)) {
    //            c = dignucl(c);
    //        }
    //        VERIFY(is_dignucl(c));
    //        return RuntimeSeq<max_size_, T> (size_ + 1, nucl(c) + str());
    //    }

    //todo naming convention violation!
    RuntimeSeq<max_size_, T> pushFront(char c) const {
        VERIFY(size_ + 1 <= max_size_);
        if (is_nucl(c)) {
            c = dignucl(c);
        }
        VERIFY(is_dignucl(c));
        RuntimeSeq<max_size_, T> res(size_ + 1);

        size_t data_size = GetDataSize(size_ + 1);

        T rm = c;
        for (size_t i = 0; i < data_size; ++i) {
            T new_rm = (data_[i] >> (TBits - 2)) & 3;
            res.data_[i] = (data_[i] << 2) | rm;
            rm = new_rm;
        }

        return res;
    }

//todo naming convention violation!
    void pushFrontThis(char c) {
        VERIFY(size_ + 1 <= max_size_);

        if (is_nucl(c)) {
            c = dignucl(c);
        }

        size_ += 1;
        size_t data_size = GetDataSize(size_);

        T rm = c;
        for (size_t i = 0; i < data_size; ++i) {
            T new_rm = (data_[i] >> (TBits - 2)) & 3;
            data_[i] = (data_[i] << 2) | rm;
            rm = new_rm;
        }
    }

    /**
     * Shift right
     *
     * @param c New 0123 char which should be added to the left.
     * @return Shifted (to the right) sequence with 'c' char on the left.
     */
    RuntimeSeq<max_size_, T> operator>>(char c) const {
        if (is_nucl(c)) {
            c = dignucl(c);
        }
        VERIFY_DEV(is_dignucl(c));

        RuntimeSeq<max_size_, T> res(*this);
        size_t data_size = GetDataSize(size_);

        T rm = c;
        for (size_t i = 0; i < data_size; ++i) {
            T new_rm = (res.data_[i] >> (TBits - 2)) & 3;
            res.data_[i] = (res.data_[i] << 2) | rm;
            rm = new_rm;
        }

        res.data_[data_size - 1] &= MaskForLastBucket(size_);

        return res;
    }

    //todo remove code duplication!
    void operator>>=(char c) {
        if (is_nucl(c)) {
            c = dignucl(c);
        }
        VERIFY_DEV(is_dignucl(c));

        size_t data_size = GetDataSize(size_);

        T rm = (T) c;
        for (size_t i = 0; i < data_size; ++i) {
            T new_rm = (data_[i] >> (TBits - 2)) & 3;
            data_[i] = (data_[i] << 2) | rm;
            rm = new_rm;
        }

        data_[data_size - 1] &= MaskForLastBucket(size_);
    }

    bool operator==(const RuntimeSeq<max_size_, T> &s) const {
        VERIFY_DEV(size_ == s.size_);

        size_t data_size = GetDataSize(size_);
        for (size_t i = 0; i < data_size; ++i)
            if (data_[i] != s.data_[i])
                return false;

        return true;
    }

    /**
     * @see operator ==()
     */
    bool operator!=(const RuntimeSeq<max_size_, T> &s) const {
        return !operator==(s);
    }

    /**
     * String representation of this Seq
     *
     * @return ACGT-string of length size_
     * @see nucl()
     */
    std::string str() const {
        std::string res(size_, '-');
        for (size_t i = 0; i < size_; ++i) {
            res[i] = nucl(operator[](i));
        }
        return res;
    }

    std::string err() const {
        return "";
    }


    std::string full_str() const {
        std::string res(max_size, '-');
        for (size_t i = 0; i < max_size; ++i) {
            res[i] = nucl(operator[](i));
        }
        return res;
    }

    size_t size() const {
        return size_;
    }

    size_t data_size() const {
        return GetDataSize(size_);
    }

    const T *data() const {
        return data_.data();
    }

    template<size_t size2_, typename T2 = T>
    Seq<size2_, T2> get_seq() const {
        VERIFY(size2_ == size_);
        return Seq<size2_, T2>((T2 *) data_.data());
    }

    template<size_t size2_, typename T2 = T>
    SimpleSeq<size2_, T2> get_sseq() const {
        VERIFY(size2_ <= max_size_);
        return SimpleSeq<size2_, T2>((T2 *) data_.data());
    }

    void copy_data(void *dst) const {
        memcpy(dst, (const void *) data_.data(), GetDataSize(size_) * sizeof(T));
    }

    char last() const {
        return operator[](size_ - 1);
    }

    char first() const {
        return operator[](0);
    }

    static size_t GetHash(const DataType *data, size_t sz, uint64_t seed = 0) {
        return CityHash64WithSeed((const char *) data, sz * sizeof(DataType), 0x9E3779B9 ^ seed);
    }

    size_t GetHash(uint64_t seed = 0) const {
        return GetHash(data_.data(), GetDataSize(size_), seed);
    }

    struct hash {
        size_t operator()(const RuntimeSeq<max_size_, T> &seq, uint64_t seed = 0) const {
            return seq.GetHash(seed);
        }

        size_t operator()(const DataType *data, size_t sz, uint64_t seed = 0) {
            return GetHash(data, sz, seed);
        }
    };

    struct less2 {
        int operator()(const RuntimeSeq<max_size_, T> &l, const RuntimeSeq<max_size_, T> &r) const {
            for (size_t i = 0; i < l.size(); ++i) {
                if (l[i] != r[i]) {
                    return (l[i] < r[i]);
                }
            }
            return l.size() < r.size();
        }
    };

    /**
     * Denotes some (weird) order on k-mers. Works fast.
     */
    struct less2_fast {
        bool operator()(const RuntimeSeq<max_size_, T> &l, const RuntimeSeq<max_size_, T> &r) const {
            return 0 > memcmp(l.data(), r.data(), sizeof(T) * l.data_size());
        }
    };

    struct less3 {
        bool operator()(const RuntimeSeq<max_size_, T> &l, const RuntimeSeq<max_size_, T> &r) const {
            VERIFY_DEV(l.size() == r.size());
            const T* l_data = l.data();
            const T* r_data = r.data();
            for (size_t i = 0; i < l.data_size(); ++i)
                if (l_data[i] != r_data[i])
                    return l_data[i] < r_data[i];
            return false;
        }
    };
};

template<size_t max_size_, typename T = seq_element_type>
bool operator<(const RuntimeSeq<max_size_, T> &l, const RuntimeSeq<max_size_, T> &r) {
    for (size_t i = 0; i < l.size(); ++i) {
        if (l[i] != r[i]) {
            return (l[i] < r[i]);
        }
    }

    return l.size() < r.size();
}

template<size_t max_size_, typename T>
std::ostream &operator<<(std::ostream &os, RuntimeSeq<max_size_, T> seq) {
    os << seq.str();
    return os;
}

namespace std {
template<size_t max_size, typename T>
struct hash<RuntimeSeq<max_size, T>> {
    size_t operator()(const RuntimeSeq<max_size, T> &seq) const {
        return seq.GetHash();
    }
};

}

typedef RuntimeSeq<UPPER_BOUND> RtSeq;

#endif /* RTSEQ_HPP_ */
