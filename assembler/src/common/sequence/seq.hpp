//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

/**
 * @file    seq.hpp
 * @author  vyahhi
 * @version 1.0
 *
 * @section LICENSE
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of
 * the License, or (at your option) any later version.
 *
 * @section DESCRIPTION
 *
 * Immutable ACGT-sequence with compile-time size.
 * It compress sequence to array of Ts (default: char).
 */

#ifndef SEQ_HPP_
#define SEQ_HPP_

#include <string>
#include <array>
#include <algorithm>
#include <cstring>
#include <iostream>

#include <city/city.h>

#include "utils/verify.hpp"
#include "nucl.hpp"
#include "utils/log.hpp"
#include "seq_common.hpp"


/**
 * @param T is max number of nucleotides, type for storage
 */
template<size_t size_, typename T = seq_element_type>
class Seq {
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

    /**
     * @variable Number of Ts which required to store all sequence.
     */
    const static size_t DataSize = (size_ + TNucl - 1) >> TNuclBits;

    typedef T DataType;

    /**
     * @variable Number of meaningful bytes in whick seq is stored
     */
    const static size_t TotalBytes = sizeof(T) * DataSize;

    static size_t GetDataSize(size_t size) {
        VERIFY(size == size_);
        return (size_ + TNucl - 1) >> TNuclBits;
    }

private:
    /* *
     * @variable Just some prime number to count the hash function of the kmer
     * */
    const static size_t PrimeNum = 239;

    // number of nucleotides in the last data_ bucket
    const static size_t NuclsRemain = size_ & (TNucl - 1);

    // useful mask to fill the last element of the data_ array
    const static size_t MaskForLastBucket = (((T) 1) << (NuclsRemain << 1)) - 1;


    /**
     * @variable Inner representation of sequence: array of Ts with length = DataSize.
     *
     * @invariant Invariant: all nucleotides >= size_ are 'A's (useful for comparison)
     */
    std::array<T, DataSize> data_;

    friend class Seq<size_ - 1, T>;

    /**
     * Initialize data_ array of this object with C-string
     *
     * @param s C-string (ACGT chars only), strlen(s) = size_
     */
    void init(const char *s) {
        T data = 0;
        size_t cnt = 0;
        int cur = 0;
        for (size_t pos = 0; pos != size_; ++pos, ++s) { // unsafe!
            // VERIFY(is_nucl(*s)); // for performance
            data = data | (T) ((T) dignucl(*s) << cnt);
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
        VERIFY(*s == 0); // C-string always ends on 0
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
    /**
     * Default constructor, fills Seq with A's
     */
    Seq() {
        std::fill(data_.begin(), data_.end(), 0);
    }

    Seq(const char *s) {
        init(s);
    }

    explicit Seq(T *data_array) {
        memcpy(data_.data(), data_array, TotalBytes);
    }

    explicit Seq(unsigned, const T *data_array) {
        memcpy(data_.data(), data_array, TotalBytes);
    }


    /**
     * Ultimate constructor from ACGT0123-string.
     *
     * @param s Any object with operator[], which returns 0123 chars
     * @param offset Offset when this sequence starts
     * @number_to_read A number of nucleotides, we want to fetch from this string
     * @raw Flag whether to check for string length (e.g. via strlen, or not)
     * @warning assuming that s is a correct string, filled with ACGT _OR_ 0123
     * no init method, filling right here
     */
    template<typename S>
    explicit Seq(const S &s, size_t offset = 0, size_t number_to_read = size_,
                 bool raw = false) {
        if (this->size(s) == 0) {
            return;
        }
        VERIFY(offset < this->size(s));
        VERIFY(is_dignucl(s[offset]) || is_nucl(s[offset]));
        if (!raw)
            VERIFY(offset + number_to_read <= this->size(s));

        // which symbols does our string contain : 0123 or ACGT?
        bool digit_str = is_dignucl(s[offset]);

        // data -- one temporary variable corresponding to the i-th array element
        // and some counters
        T data = 0;
        size_t cnt = 0;
        size_t cur = 0;

        for (size_t i = 0; i < number_to_read; ++i) {
            //VERIFY(is_dignucl(s[i]) || is_nucl(s[i]));

            // we fill everything with zeros (As) by default.
            char c = digit_str ? s[offset + i] : (char) dignucl(s[offset + i]);

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

        for (; cur != DataSize; ++cur)
            this->data_[cur] = 0;
    }


    /**
     * Get i-th symbol of Seq.
     *
     * @param i Index of the symbol (0 <= i < size_)
     * @return 0123-char on position i
     */
    char operator[](const size_t i) const {
        return (data_[i >> TNuclBits] >> ((i & (TNucl - 1)) << 1)) & 3;
    }

    /**
     * Reverse complement.
     *
     * @return Reverse complement Seq.
     */
    Seq<size_, T> operator!() const {
        Seq<size_, T> res(*this);
        for (size_t i = 0; i < (size_ >> 1); ++i) {
            T front = complement(res[i]);
            T end = complement(res[size_ - 1 - i]);
            res.set(i, (char) end);
            res.set(size_ - 1 - i, (char) front);
        }
        if ((size_ & 1) == 1) {
            res.set(size_ >> 1, complement(res[size_ >> 1]));
        }
        // can be made without complement calls, but with xor on all bytes afterwards.
        return res;
    }

    /**
     * Shift left
     *
     * @param c New 0123 char which should be added to the right.
     * @return Shifted (to the left) sequence with 'c' char on the right.
     */
    Seq<size_, T> operator<<(char c) const {
        if (is_nucl(c)) {
            c = dignucl(c);
        }
        Seq<size_, T> res(*this);
        std::array<T, DataSize> &data = res.data_;
        if (DataSize != 0) { // unless empty sequence
            T rm = data[DataSize - 1] & 3;
            T lastnuclshift_ = ((size_ + TNucl - 1) & (TNucl - 1)) << 1;
            data[DataSize - 1] = (data[DataSize - 1] >> 2) | ((T) c << lastnuclshift_);

            if (DataSize >= 2) { // if we have at least 2 elements in data
                int data_size = DataSize;
                for (int i = data_size - 2; i >= 0; --i) {
                    T new_rm = data[i] & 3;
                    data[i] = (data[i] >> 2) |
                              (rm << (TBits - 2)); // we need & here because if we shift negative, it fill with ones :(
                    rm = new_rm;
                }
            }
        }
        return res;
    }

    Seq<size_ + 1, T> pushBack(char c) const {
        if (is_nucl(c)) {
            c = dignucl(c);
        }
        //VERIFY(is_dignucl(c));
        Seq<size_ + 1, T> s;
        copy(this->data_.begin(), this->data_.end(), s.data_.begin());
        s.data_[s.DataSize - 1] = s.data_[s.DataSize - 1] | ((T) c << ((size_ & (TNucl - 1)) << 1));

        return s; //was: Seq<size_ + 1, T>(str() + nucl(c));

    }

    //    /**
    //   * @todo optimize!!!
    //   */
    //    Seq<size_ + 1, T> pushFront(char c) const {
    //        if (is_nucl(c)) {
    //            c = dignucl(c);
    //        }
    //        VERIFY(is_dignucl(c));
    //        return Seq<size_ + 1, T> (nucl(c) + str());
    //    }

    Seq<size_ + 1, T> pushFront(char c) const {
        if (is_nucl(c)) {
            c = dignucl(c);
        }
        VERIFY(is_dignucl(c));
        Seq<size_ + 1, T> res;

        //if new kmer has more Ts
        if (Seq<size_ + 1, T>::DataSize > DataSize) {
            res.data_[DataSize] = (data_[DataSize - 1] >> (TBits - 2)) & 3;
        }

        T rm = c;
        for (size_t i = 0; i < DataSize; ++i) {
            T new_rm = (data_[i] >> (TBits - 2)) & 3;
            res.data_[i] = (data_[i] << 2) | rm;
            rm = new_rm;
        }

        return res;
    }

    /**
     * Shift right
     *
     * @param c New 0123 char which should be added to the left.
     * @return Shifted (to the right) sequence with 'c' char on the left.
     */
    Seq<size_, T> operator>>(char c) const {
        if (is_nucl(c)) {
            c = dignucl(c);
        }
        VERIFY(is_dignucl(c));
        Seq<size_, T> res(*this);
        T rm = c;
        for (size_t i = 0; i < DataSize; ++i) {
            T new_rm = (res.data_[i] >> (TBits - 2)) & 3;
            res.data_[i] = (res.data_[i] << 2) | rm;
            rm = new_rm;
        }
        if ((size_ & (TNucl - 1)) != 0) {
            T lastnuclshift_ = (size_ & (TNucl - 1)) << 1;
            res.data_[DataSize - 1] = res.data_[DataSize - 1] & (((T) 1
                                                                  << lastnuclshift_) - 1);
        }
        return res;
    }

    /**
     * Sets i-th symbol of Seq with 0123-char
     */
    inline void set(const size_t i, char c) {
        data_[i >> TNuclBits] =
                (data_[i >> TNuclBits] & ~((T) 3 << ((i & (TNucl - 1)) << 1))) | ((T) c << ((i & (TNucl - 1)) << 1));
    }

    bool operator==(const Seq<size_, T> &s) const {
        for (size_t i = 0; i < DataSize; ++i)
            if (data_[i] != s.data_[i])
                return false;
        return true;
    }

    /**
     * @see operator ==()
     */

    bool operator!=(const Seq<size_, T> &s) const {
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
        for (size_t i = 0; i != size_; ++i) {
            res[i] = nucl(operator[](i));
        }
        return res;
    }

    static size_t size() {
        return size_;
    }


    void copy_data(void *dst) const {
        memcpy(dst, (const void *) data_.data(), TotalBytes);
    }

    /**
     *  Reads sequence from the file (in the same format as BinWrite writes it)
     *  and returns false if error occured, true otherwise.
     */
    static bool BinRead(std::istream &file, Seq<size_> *seq) {
        file.read((char *) seq->data_.data(), sizeof(T) * DataSize);
        return !file.fail();
    }

    /**
     *  Writes sequence to the file (in the same format as BinRead reads it)
     *  and returns false if error occured, true otherwise.
     */
    static bool BinWrite(std::ostream &file, const Seq<size_> &seq) {
        file.write((const char *) seq.data_.data(), sizeof(T) * DataSize);
        return !file.fail();
    }

    /**
     *  Reads sequence from the file (in the same format as BinWrite writes it)
     *  and returns false if error occured, true otherwise.
     */
    bool BinRead(std::istream &file) {
        return BinRead(file, this);
    }

    /**
     *  Writes sequence to the file (in the same format as BinRead reads it)
     *  and returns false if error occured, true otherwise.
     */
    bool BinWrite(std::ostream &file) const {
        return BinWrite(file, *this);
    }

    /**
     * @see Seq
     */
    template<size_t size2_, typename T2 = T>
    Seq<size2_, T2> start() const {
        VERIFY(size2_ <= size_);
        return Seq<size2_, T2>(*this);
    }

    template<size_t size2_/* = size_ - 1*/, typename T2 = T>
    Seq<size2_, T2> end() const {
        VERIFY(size2_ <= size_);
        return Seq<size2_, T2>(*this, size_ - size2_);
    }

    const T *data() const {
        return data_.data();
    }

    size_t data_size() const {
        return DataSize;
    }


    char last() const {
        return operator[](size_ - 1);
    }

    char first() const {
        return operator[](0);
    }

    static size_t GetHash(const DataType *data, size_t sz = DataSize, uint32_t seed = 0) {
        return CityHash64WithSeed((const char *) data, sz * sizeof(DataType), 0x9E3779B9 ^ seed);
    }

    size_t GetHash(uint32_t seed = 0) const {
        return GetHash(data_.data(), DataSize, seed);
    }

    struct hash {
        size_t operator()(const Seq<size_, T> &seq, uint32_t seed = 0) const {
            return seq.GetHash(seed);
        }

        size_t operator()(const DataType *data, size_t sz = DataSize, uint32_t seed = 0) {
            return GetHash(data, sz, seed);
        }
    };

    struct equal_to {
        bool operator()(const Seq<size_, T> &l, const Seq<size_, T> &r) const {
            return r == l;
        }
    };

    struct less2 {
        bool operator()(const Seq<size_, T> &l, const Seq<size_, T> &r) const {
            for (size_t i = 0; i < size_; ++i) {
                if (l[i] != r[i]) {
                    return (l[i] < r[i]);
                }
            }
            return false;
        }
    };

    /**
     * Denotes some (weird) order on k-mers. Works fast.
     */
    struct less2_fast {
        bool operator()(const Seq<size_, T> &l, const Seq<size_, T> &r) const {
            return 0 > memcmp(l.data_.data(), r.data_.data(), sizeof(T) * DataSize);
        }
    };
};

template<size_t size_, typename T>
std::ostream &operator<<(std::ostream &os, Seq<size_, T> seq) {
    os << seq.str();
    return os;
}

//namespace std {
//
//template<size_t size_, typename T = seq_element_type>
//struct hash<Seq<size_, T> {
//    typedef size_t result_type;
//    typedef Seq<size_, T> argument_type;
//
//    result_type operator() (const argument_type& arg) {
//        return Seq<size_, T>::hash()(arg);
//    }
//};
//
//}

#endif /* SEQ_HPP_ */
