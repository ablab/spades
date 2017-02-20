//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

/*
 * simple_seq.hpp
 *
 *  Created on: Jul 23, 2012
 *      Author: andrey
 */

#ifndef SIMPLE_SEQ_HPP_
#define SIMPLE_SEQ_HPP_

#include <string>
#include <array>
#include <algorithm>
#include <cstring>
#include <iostream>

#include "utils/verify.hpp"
#include "nucl.hpp"
#include "math/log.hpp"
#include "seq_common.hpp"
/**
 * @param T is max number of nucleotides, type for storage
 */
template<size_t size_, typename T = seq_element_type>
class SimpleSeq {
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

private:
    // number of nucleotides in the last data_ bucket
    const static size_t NuclsRemain = size_ & (TNucl - 1);

    // useful mask to fill the last element of the data_ array
    const static size_t MaskForLastBucket = (((T) 1) << (NuclsRemain << 1) ) - 1;


    /**
     * @variable Inner representation of sequence: array of Ts with length = DataSize.
     *
     * @invariant Invariant: all nucleotides >= size_ are 'A's (useful for comparison)
     */
    std::array<T, DataSize> data_;


public:

    SimpleSeq() {
        //VERIFY((T)(-1) >= (T)0);//be sure to use unsigned types
        std::fill(data_.begin(), data_.end(), 0);
    }

    explicit SimpleSeq(T * data_array) {
        memcpy(data_.data(), data_array, TotalBytes);
    }


    char operator[](const size_t i) const {
        //VERIFY(i >= 0);
        //VERIFY(i < size_);
        return (data_[i >> TNuclBits] >> ((i & (TNucl - 1)) << 1)) & 3;
    }

    std::string str() const {
        std::string res(size_, '-');
        for (size_t i = 0; i < size_; ++i) {
            res[i] = nucl(operator[](i));
        }
        return res;
    }

    void copy_data(void * dst) const {
        memcpy(dst, (const void *) data_.data(), TotalBytes);
    }

    static size_t GetHash(const DataType *data, size_t sz, uint32_t seed = 0) {
        return CityHash64WithSeed((const char*)data, sz * sizeof(DataType), 0x9E3779B9 ^ seed);
    }

    size_t GetHash(uint32_t seed = 0) const {
        return GetHash(data_.data(), DataSize, seed);
    }

    struct hash {
        size_t operator()(const SimpleSeq<size_, T>& seq, uint32_t seed = 0) const {
            return seq.GetHash(seed);
        }

        size_t operator()(const DataType *data, size_t sz, unsigned seed = 0) {
            return GetHash(data, sz, seed);
        }
    };

    struct equal_to {
        bool operator()(const SimpleSeq<size_, T>& l, const SimpleSeq<size_, T>& r) const {
            for (size_t i = 0; i < DataSize; ++i)
                if (l.data_[i] != r.data_[i])
                    return false;
            return true;
        }
    };

    struct less2 {
        int operator()(const SimpleSeq<size_, T> &l, const SimpleSeq<size_, T> &r) const {
            for (size_t i = 0; i < size_; ++i) {
                if (l[i] != r[i]) {
                    return (l[i] < r[i]);
                }
            }
            return false;
        }
    };

};

template<size_t size_, typename T>
std::ostream& operator<<(std::ostream& os, SimpleSeq<size_, T> seq) {
    os << seq.str();
    return os;
}


#endif /* SIMPLE_SEQ_HPP_ */
