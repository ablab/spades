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

#include <vector>
#include <string>
#include "log.hpp"

class SequenceData {
private:
    friend class Sequence;
    // type to store Seq in Sequences
    typedef u_int64_t ST;
    // number of bits in ST
	const static size_t STBits = sizeof(ST) << 3;
    // number of nucleotides in ST
    const static size_t STN = (STBits >> 1);
    // number of bits in STN (for faster div and mod)
    const static size_t STNBits = log_<STN, 2>::value;
    // ref counter
    volatile size_t kCount; 
    // sequence (actual data for what it's for)
    ST *bytes_;



    // methods:
    SequenceData(const SequenceData &sd); // forbidden
    SequenceData& operator=(const SequenceData&); // forbidden


    void Grab() {
        ++kCount;
    }
    void Release() {
        --kCount;
        if (kCount == 0) {
            delete this;
        }
    }

public:
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

    SequenceData(size_t size_): count(0) {
        size_t size = size_;
        size_t bytes_size = (size + STN - 1) >> STNbits;
        bytes_ = (Seq<STN, ST>*) malloc(bytes_size * sizeof(Seq<STN, ST> )); // it's a bit faster than new
        for (size_t i = 0; i < bytes_size; ++i) {
            bytes_[i] = Seq<STN, ST>();
        }
    }

    bool BinRead(std::istream& file, size_t size) {
        size_t bytes_size = (size + STN - 1) >> STNbits;
        for (size_t i = 0; i < bytes_size; ++i) {
            bytes_[i].BinRead(file);
        }
        return !file.fail();
    }

    bool BinWrite(std::ostream& file, size_t size) {
        size_t bytes_size = (size + STN - 1) >> STNbits;
        for (size_t i = 0; i < bytes_size; ++i) {
            bytes_[i].BinWrite(file);
        }
        return !file.fail();
    }

    ~SequenceData() {
        free(bytes_);
    }

    // for a faster start method of Sequence
    template<size_t size>
    Seq<size> GetStartSeq() {
        return Seq<size>(bytes_);
    }

    char operator[](const size_t i) const {
		return (bytes_[i >> STNBits] >> ((i & (STN - 1)) << 1)) & 3; // btw (i % Tnucl) <=> (i & (Tnucl-1))
    }
};

#endif /* REFCOUNT_HPP_ */
