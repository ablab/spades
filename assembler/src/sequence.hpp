/*
 * sequence.hpp
 *
 *  Created on: 01.03.2011
 *      Author: vyahhi
 */

#ifndef SEQUENCE_HPP_
#define SEQUENCE_HPP_

#include "seq.hpp"

//SEQUENCE IS IMMUTABLE!!!
class Sequence { // immutable runtime length sequence (slow!!!)
	class Data {
	public:
		vector < Seq<4> > bytes_;
		int ref_;
		void IncRef() {++ref_;}
		void DecRef() {--ref_;}
		Data(vector< Seq<4> > bytes) : bytes_(bytes) , ref_(1) {}
	};

public:
//	template <size_t _size> static Sequence constr(const Seq<_size> seq);

	void init(const std::string &s);

	Sequence(const std::string &s);

	template <size_t _size>
	Sequence(const Seq<_size> seq) {
		init(seq.str());
	}

	~Sequence();
	char operator[](const size_t index) const;
	bool operator==(const Sequence& that) const;
	const Sequence operator!() const;

	/**
	 * @param from inclusive
	 * @param to exclusive;
	 */
	Sequence Subseq(size_t from, size_t to) const;
	Sequence operator+ (const Sequence &s) const;
	Sequence(const Sequence& s);
	int find(const Sequence& t, int from = 0) const;
	int similar(const Sequence &t, int k, char directed = 0) const;
	std::string Str() const;
	size_t size() const;
private:
	Data* data_;
	size_t from_; // should be const
	size_t size_; // should be const
	/**
	 * Right to left + complimentary?
	 */
	bool rtl_; // should be const
//	Sequence(const Sequence *svl, bool reverse); // reverse
	Sequence(Data* data, size_t from, size_t size, bool rtl);
	Sequence& operator=(const Sequence& that) {cerr << "Don't call operator= for Sequence"; return *this;};
};

#endif /* SEQUENCE_HPP_ */
