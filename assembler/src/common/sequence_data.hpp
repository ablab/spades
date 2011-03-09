/*
 * Sequence Data with ref counter
 *
 * Some of the source code was taken from the book "C++ for Real Programmers"
 *
 *  Created on: 09.03.2011
 *      Author: vyahhi
 */

#ifndef REFCOUNT_HPP_
#define REFCOUNT_HPP_

#include <vector>
#include <string>
#include "log.hpp"
using namespace std;

class SequenceData {
private:
	// type to store Seq in Sequences
	typedef int ST;
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
public:
	void Grab() {
		count++;
	}
	void Release() {
		if (count > 0) {
			count--;
		}
		if (count == 0) {
			delete this;
		}
	}

	~SequenceData() {
		free(bytes_);
		//delete[] bytes_;
	}
	/**
	 * Sequence initialization (arbitrary size string)
	 *
	 * @param s ACGT-string
	 */
	SequenceData(const std::string &s) :
		count(0) {
		size_t size_ = s.size();
		size_t bytes_size = (size_ + STN - 1) >> STNbits;
		//bytes_ = NULL;
		bytes_ = (Seq<STN, ST>*) malloc(bytes_size * sizeof(Seq<STN, ST> )); // it's a bit faster than new
		//bytes_ = new Seq<STN,ST>[bytes_size];
		for (size_t i = 0; i < (size_ >> STNbits); ++i) {
			bytes_[i] = Seq<STN, ST> (s, i << STNbits);
		}
		if (size_ & (STN - 1)) {
			// fill with As for not breaking contract of Seq
			string s2 = s.substr(size_ & ~(STN - 1), size_ & (STN - 1));
			size_t count = STN - (size_ & (STN - 1));
			s2.append(count, 'A');
			bytes_[size_ >> STNbits] = Seq<STN, ST> (s2);
		}
	}

	inline char operator[](const size_t i) const {
		return bytes_[i >> STNbits][i & (STN - 1)];
	}

};

#endif /* REFCOUNT_HPP_ */
