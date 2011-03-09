/*
 * sequence.hpp
 *
 *  Created on: 01.03.2011
 *      Author: vyahhi
 */

#ifndef SEQUENCE_HPP_
#define SEQUENCE_HPP_

#include "seq.hpp"
#include "sequence_data.hpp"
#include <vector>
#include <string>
#include <iostream>
using namespace std;

/**
 * Immutable runtime length sequence (someway slow)
 */

class Sequence {
private:
	SequenceData *data_;
	const size_t from_;
	const size_t size_;
	const bool rtl_; // Right to left + complimentary (?)
	Sequence(const Sequence &seq, size_t from, size_t size, bool rtl);
	Sequence& operator=(const Sequence &); // forbidden
public:
	// constructors:
	template<size_t _size>
	Sequence(const Seq<_size> seq) :
		from_(0), size_(seq.size()), rtl_(false) { // TODO: optimize
		data_ = new SequenceData(seq.str());
		data_->Grab();
	}
	Sequence(const Sequence &s);
	Sequence(const string &s);
	~Sequence();

	// other methods:
	char operator[](const size_t index) const;
	bool operator==(const Sequence &that) const;
	Sequence operator!() const;
	/**
	 * @param from inclusive
	 * @param to exclusive;
	 */
	Sequence Subseq(size_t from, size_t to) const;
	Sequence Subseq(size_t from) const; // up to size_ by default
	Sequence operator+(const Sequence &s) const;
	int find(const Sequence &t, int from = 0) const;
	int similar(const Sequence &t, int k, char directed = 0) const;
	int leftSimilar(const Sequence &t, int k) const;
	int rightSimilar(const Sequence &t, int k) const;
	string str() const;
	size_t size() const;
};

#endif /* SEQUENCE_HPP_ */
