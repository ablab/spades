//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/**
 * @file    sequence_data.hpp
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
 * Sequence Data with ref counter
 */

#ifndef SEQUENCE_DATA_HPP_
#define SEQUENCE_DATA_HPP_

#include "log.hpp"
#include "seq_common.hpp"

#include <vector>
#include <string>

class SequenceData {
private:
    friend class Sequence;
    // type to store Seq in Sequences
    typedef seq_element_type ST;
    // number of bits in ST
    const static size_t STBits = sizeof(ST) << 3;
    // number of nucleotides in ST
    const static size_t STN = (STBits >> 1);
    // number of bits in STN (for faster div and mod)
    const static size_t STNBits = log_<STN, 2>::value;
    // ref counter
    size_t kCount;
    // sequence (actual data for what it's for)
    ST *bytes_;

    // methods:
    SequenceData(const SequenceData &sd); // forbidden
    SequenceData& operator=(const SequenceData&); // forbidden

    void Grab() {
#     pragma omp atomic
      kCount += 1;
    }
    void Release() {
      // FIXME: This is really not correct. Here we have race condition between Grab() and Release()
#     pragma omp atomic
      kCount -= 1;

#     pragma omp flush(kCount)
#     pragma omp critical
      if (kCount == 0) {
          delete this;
      }
    }

public:

    template<size_t size_>
    SequenceData(const Seq<size_> &kmer): kCount(0) {
        size_t size = (size_ + STN - 1) >> STNBits;

        bytes_ = (ST*) malloc(size * sizeof(ST));

        kmer.copy_data((void *) bytes_);
    }

    template<size_t size_>
    SequenceData(const RuntimeSeq<size_> &kmer): kCount(0) {
        size_t size = (kmer.size() + STN - 1) >> STNBits;

        bytes_ = (ST*) malloc(size * sizeof(ST));

        kmer.copy_data((void *) bytes_);
    }
    /**
     * Sequence initialization (arbitrary size string)
     * copypaste from the seq constructor
     *
     * @param s ACGT or 0123-string
    */
    template<typename S>
    SequenceData(const S &s, size_t size) : kCount(0)
    {
        size_t bytes_size = (size + STN - 1) >> STNBits;

        VERIFY(is_dignucl(s[0]) || is_nucl(s[0]));

        bytes_ = (ST*) malloc(bytes_size * sizeof(ST));

        // which symbols does our string contain : 0123 or ACGT?
        bool digit_str = is_dignucl(s[0]);

        // data -- one temporary variable corresponding to the i-th array element
        // and some counters
        ST data = 0;
        size_t cnt = 0;
        size_t cur = 0;

        for (size_t i = 0; i < size; ++i) {
            //VERIFY(is_dignucl(s[i]) || is_nucl(s[i]));
            char c = digit_str ? s[i] : dignucl(s[i]);

            data = data | (ST(c) << cnt);
            cnt += 2;

            if (cnt == STBits) {
                bytes_[cur++] = data;
                cnt = 0;
                data = 0;
            }
        }


        if (cnt != 0) {
            bytes_[cur++] = data;
        }

        for (; cur < bytes_size; ++cur)
            bytes_[cur] = 0;
    }


    SequenceData(size_t size_): kCount(0) {
        size_t size = size_;
        size_t bytes_size = (size + STN - 1) >> STNBits;

        bytes_ = (ST*) malloc(bytes_size * sizeof(ST));

        for (size_t i = 0; i < bytes_size; ++i) {
            bytes_[i] = 0;
        }
    }

    bool BinRead(std::istream& file, size_t size) {
        size_t bytes_size = (size + STN - 1) >> STNBits;
        for (size_t i = 0; i < bytes_size; ++i) {
            file.read((char *) &bytes_[i], sizeof(ST));
        }
        return !file.fail();
    }

    bool BinWrite(std::ostream& file, size_t size) {
        size_t bytes_size = (size + STN - 1) >> STNBits;
        for (size_t i = 0; i < bytes_size; ++i) {
            file.write((const char *) &bytes_[i], sizeof(ST));
        }
        return !file.fail();
    }

    ~SequenceData() {
        VERIFY(bytes_);
        free(bytes_);
        bytes_ = NULL;
    }

    // for a faster start method of Sequence
    template<size_t size>
    Seq<size> GetStartSeq() {
        return Seq<size>(bytes_);
    }

    void print(ST * s, size_t size) const {
        for (size_t i = 0; i < size; ++i) {
            std::cerr << nucl_map[ (s[i >> STNBits] >> ((i & (STN - 1)) << 1)) & 3 ];
        }
        std::cerr << std::endl;
    }

    template<size_t size>
    Seq<size> FastStartSeq(size_t from) const {
        ST result[(size + STN - 1) >> STNBits] = {0};

        size_t start = from >> STNBits;
        size_t end = (from + size - 1) >> STNBits;
        size_t shift = (from & (STN - 1)) << 1;

        for (size_t i = start; i <= end; ++i) {
            result[i - start] = bytes_[i] >> shift;
        }

        if (shift != 0) {
            shift = STBits - shift;

            for (size_t i = start + 1; i <= end; ++i) {
                result[i - start - 1] |= bytes_[i] << shift;
            }
        }

        return Seq<size>(result);
    }

    char operator[](const size_t i) const {
        return (bytes_[i >> STNBits] >> ((i & (STN - 1)) << 1)) & 3; // btw (i % Tnucl) <=> (i & (Tnucl-1))
    }
};

#endif /* REFCOUNT_HPP_ */
