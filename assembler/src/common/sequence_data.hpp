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
using namespace std;

class SequenceData {
private:
	// type to store Seq in Sequences
	typedef short ST;
	// number of nucleotides in ST
	static const int STN = (sizeof(ST) * 4);
	// ref counter
	size_t count; // should be volatile some day for multi-threading
	// sequence (actual data for what it's for)
	Seq<STN,ST> *bytes_;
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
		//free(bytes_);
		delete[] bytes_;
	}
	/**
	 * Sequence initialization (arbitrary size string)
	 *
	 * @param s ACGT-string
	 */
	SequenceData(const std::string &s) :
		count(0) {
		size_t size_ = s.size();
		size_t bytes_size = (size_ + STN - 1) / STN;
		bytes_ = NULL;
		//void* x = malloc(bytes_size * sizeof(Seq<STN, ST>));
		//bytes_ = (Seq<STN,ST>*) x; // it's a bit faster than ne
		bytes_ = new Seq<STN,ST>[bytes_size];
		for (size_t i = 0; i < size_ / STN; ++i) {
			bytes_[i] = Seq<STN, ST>(s.substr(i * STN, STN));
		}
		if (size_ % STN != 0) {
			// fill with As for not breaking contract of Seq
			string s2 = s.substr((size_ / STN) * STN, size_ % STN);
			size_t count = STN - (size_ % STN);
			s2.append(count, 'A');
			bytes_[size_ / STN] = Seq<STN, ST>(s2);
		}

	}

	char operator[](const size_t i) const {
		return bytes_[i / STN][i % STN];
	}
};

#endif /* REFCOUNT_HPP_ */
