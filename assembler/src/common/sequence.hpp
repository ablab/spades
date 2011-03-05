/*
 * sequence.hpp
 *
 *  Created on: 01.03.2011
 *      Author: vyahhi
 */

#ifndef SEQUENCE_HPP_
#define SEQUENCE_HPP_

#include "seq.hpp"
#include <vector>
#include <string>
#include <iostream>
using namespace std;

typedef char SequenceT;
#define SequenceTnucl 4

//SEQUENCE IS IMMUTABLE!!!
class Sequence { // immutable runtime length sequence (slow!!!)

private:
	/**
	 * Data class for reference counting, contains sequence and ref count.
	 */
	class Data {
	public:
		Seq<4,SequenceT> *bytes_;
		size_t ref_;
		inline void IncRef() {++ref_;}
		inline void DecRef() {--ref_;}
		Data(const std::vector<Seq<4,SequenceT> > &bytes) : ref_(1) {
			bytes_ = new Seq<4,SequenceT>[bytes.size()];
			copy(bytes.begin(), bytes.end(), bytes_);
		}
		~Data() {
			delete[] bytes_;
		}
	};

	void init(const std::string &s);

	Data* data_;
	const size_t from_; // should be const
	const size_t size_; // should be const
	/**
	 * Right to left + complimentary?
	 */
	const bool rtl_; // should be const
	Sequence(Data* data, size_t from, size_t size, bool rtl);
	Sequence& operator=(const Sequence& that);

public:
	template <size_t _size>
	Sequence(const Seq<_size> seq) : from_(0), size_(seq.size()), rtl_(false) {
		init(seq.str());
	}

	Sequence(const std::string &s);
	~Sequence();
	char operator[](const size_t index) const;
	bool operator==(const Sequence& that) const;
	Sequence operator!() const;

	/**
	 * @param from inclusive
	 * @param to exclusive;
	 */
	Sequence Subseq(size_t from, size_t to) const;

	Sequence Subseq(size_t from) const; // up to size_ by default
	Sequence operator+ (const Sequence &s) const;
	Sequence(const Sequence& s);
	int find(const Sequence& t, int from = 0) const;
	int similar(const Sequence &t, int k, char directed = 0) const;
	std::string str() const;
	int leftSimilar(const Sequence &t, int k) const;
	int rightSimilar(const Sequence &t, int k) const;
	size_t size() const;
};

#endif /* SEQUENCE_HPP_ */
