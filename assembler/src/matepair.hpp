/*
 * matepair.hpp
 *
 *  Created on: 03.03.2011
 *      Author: vyahhi
 */

#ifndef MATEPAIR_HPP_
#define MATEPAIR_HPP_

#include "seq.hpp"

template <int size_, typename T = char> // max number of nucleotides in each read
class MatePair {
public:
	MatePair(const char *s1, const char *s2, const int id) : id_(id), seq1_(s1), seq2_(s2) {
		// nothing
	}

	bool isNull() {
		return this == null;
	}

	int id() const {
		return id_;
	}

	Seq<size_,T> seq1() const {
		return seq1_;
	}

	Seq<size_,T> seq2() const {
		return seq2_;
	}

	const static MatePair<size_, T> null;

private: // make private!
	int id_; // consecutive number from input file :)
	Seq<size_,T> seq1_;
	Seq<size_,T> seq2_;
};

template <int size_, typename T>
const MatePair<size_, T> MatePair<size_, T>::null = MatePair<size_,T>("", "", -1);

#endif /* MATEPAIR_HPP_ */
