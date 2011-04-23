/*
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
using namespace std;

class SequenceData {
private:
	friend class Sequence;
	// type to store Seq in Sequences
	typedef unsigned int ST;
	// number of nucleotides in ST
	static const size_t STN = (sizeof(ST) * 4);
	// number of bits in STN (for faster div and mod)
	static const size_t STNbits = log_<STN, 2>::value;
	// ref counter
	size_t count; // should be volatile some day for multi-threading
	// sequence (actual data for what it's for)
	Seq<STN, ST> *bytes_;
	// methods:
	SequenceData(const SequenceData &sd); // forbidden
	SequenceData& operator=(const SequenceData&); // forbidden
	void Grab() {
		++count;
	}
	void Release() {
		--count;
		if (count == 0) {
			delete this;
		}
	}

        template<typename S>
	string SubString(const S &s, size_t offset, size_t size) {
		string str;
		for (size_t i = offset; i < size; ++i) {
			str += s[i];
		}
		return str;
	}
public:
	/**
	 * Sequence initialization (arbitrary size string)
	 *
	 * @param s ACGT or 0123-string
	 */
	template<typename S>
	SequenceData(const S &s, size_t size_) :
		count(0) {
		size_t size = size_;
		size_t bytes_size = (size + STN - 1) >> STNbits;
		//bytes_ = NULL;
		bytes_ = (Seq<STN, ST>*) malloc(bytes_size * sizeof(Seq<STN, ST> )); // it's a bit faster than new
		//bytes_ = new Seq<STN,ST>[bytes_size];
		for (size_t i = 0; i < (size >> STNbits); ++i) {
		  bytes_[i] = Seq<STN, ST> (s, i << STNbits);
		}
		if (size & (STN - 1)) {
			// fill with As for not breaking contract of Seq
			//todo think of something better than using string
		  string s2 = SubString(s, size & ~(STN - 1), size);
			size_t count = STN - (size & (STN - 1));
			s2.append(count, 'A');
			bytes_[size >> STNbits] = Seq<STN, ST> (s2);
		}
	}

	~SequenceData() {
	        free(bytes_);
        }

	char operator[](const size_t i) const {
		return bytes_[i >> STNbits][i & (STN - 1)];
	}

};

#endif /* REFCOUNT_HPP_ */
