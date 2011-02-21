/*
 * seq.hpp
 *
 *  Created on: 21.02.2011
 *      Author: vyahhi
 */

#ifndef SEQ_HPP_
#define SEQ_HPP_

#include <string>

template <int size> // max number of nucleotides
class Seq {
public:
	Seq(const std::string &s);
	char operator[] (const int &index) const;
	std::string str() const;
public:
	char bytes[(size << 2) + ((size & 3) != 0)]; // little-endian
};

template <int size> // max number of nucleotides in each read
class MatePair {
public:
	MatePair(const std::string &s1, const std::string &s2, const int id_);
	MatePair(const MatePair &mp);
// private:
	int id; // consecutive number from input file :)
	Seq<size> seq1;
	Seq<size> seq2;
};

// ******************** //
// * TEMPLATE METHODS * //
// ******************** //

template <int size>
Seq<size>::Seq (const std::string &s) {
	char byte = 0;
	int cnt = 6;
	int cur = 0;
	for (std::string::const_iterator si = s.begin(); si != s.end(); si++) {
		switch (*si) {
			case 'C': byte |= (1 << cnt); break;
			case 'G': byte |= (2 << cnt); break;
			case 'T': byte |= (3 << cnt); break;
		}
		cnt -= 2;
		if (cnt < 0) {
			this->bytes[cur++] = byte;
			cnt = 6;
			byte = 0;
		}
	}
	if (cnt != 6) {
		this->bytes[cur++] = byte;
	}
}

template <int size>
char Seq<size>::operator[] (const int &index) const {
	switch ( ( this->bytes[index >> 2] >> ((3-(index%4))*2) ) & 3) { // little endian!
		case 0: return 'A'; break;
		case 1: return 'C'; break;
		case 2: return 'G'; break;
		case 3: return 'T'; break;
		default: return 'N';
	}
}

template <int size>
std::string Seq<size>::str() const {
	std::string res = "";
	for (int i = 0; i < size; ++i) {
		res += this->operator[](i);
	}
	return res;
 }

template <int size>
MatePair<size>::MatePair(const std::string &s1, const std::string &s2, const int id_) : id(id_), seq1(s1), seq2(s2) {
	//
}

template <int size>
MatePair<size>::MatePair(const MatePair &mp) : id(mp.id), seq1(mp.seq1), seq2(mp.seq2) {
	//
}

#endif /* SEQ_HPP_ */
