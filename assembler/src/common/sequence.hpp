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
#include <string.h>
using namespace std;

/**
 * Immutable runtime length sequence (someway slow)
 */

class Sequence {
private:
	SequenceData *data_;
	size_t from_;
	size_t size_;
	bool rtl_; // Right to left + complimentary (?)
	Sequence(const Sequence &seq, size_t from, size_t size, bool rtl);
public:
	// constructors:
	//	template<size_t _size>
	//	Sequence(const Seq<_size> seq) :
	//		from_(0), size_(seq.size()), rtl_(false) { // TODO: optimize
	//		data_ = new SequenceData(seq);
	//		data_->Grab();
	//	}
	const Sequence& operator=(const Sequence &rhs) {
		if (data_ == rhs.data_) {
			return *this;
		}
		data_->Release();
		data_ = rhs.data_;
		data_->Grab();
		from_ = rhs.from_;
		size_ = rhs.size_;
		rtl_ = rhs.rtl_;
		return *this;
	}

	Sequence(char* s) :
		from_(0), size_(strlen(s)), rtl_(false) {
		data_ = new SequenceData(s, size_);
		data_->Grab();
	}

	Sequence(const char* s) :
		from_(0), size_(strlen(s)), rtl_(false) {
		data_ = new SequenceData(s, size_);
		data_->Grab();
	}

	template<typename S>
	explicit Sequence(const S &s) :
		from_(0), size_(s.size()), rtl_(false) {
		data_ = new SequenceData(s, size_);
		data_->Grab();
	}

	Sequence(const Sequence &s);
	~Sequence();

	// other methods:
	char operator[](const size_t index) const;
	bool operator==(const Sequence &that) const;
	bool operator!=(const Sequence &that) const;
	bool operator<(const Sequence &that) const;
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

	/// returns true if two sequences intersect
	bool intersects(const Sequence &t) const;

	//	template<size_t size2_>
	//	Seq<size2_> start() const;
	//
	//	template<size_t size2_>
	//	Seq<size2_> end() const;

	template<size_t size2_>
	Seq<size2_> start() const;
	template<size_t size2_>
	Seq<size2_> end() const;
	string str() const;
	size_t size() const;
};

ostream& operator<<(ostream& os, const Sequence& s);

// TODO: optimize
template<size_t size2_>
Seq<size2_> Sequence::start() const {
	assert(size2_ <= size_);
	return Seq<size2_> (*this);
}

// TODO: optimize
template<size_t size2_>
Seq<size2_> Sequence::end() const {
	assert(size2_ <= size_);
	return Seq<size2_> (*this, size_ - size2_);
}

class SequenceBuilder {
	vector<char> buf_;
public:
	template<typename S>
	SequenceBuilder& append(const S &s) {
		for (size_t i = 0; i < s.size(); ++i) {
			buf_.push_back(s[i]);
		}
		return *this;
	}

	SequenceBuilder& append(char c) {
		buf_.push_back(c);
		return *this;
	}

	Sequence BuildSequence() {
		return Sequence(buf_);
	}

	size_t size() const {
		return buf_.size();
	}

	char operator[](const size_t index) const {
		assert(index < buf_.size());
		return buf_[index];
	}
};

#endif /* SEQUENCE_HPP_ */
