/*
 * seq.hpp
 *
 *  Created on: 21.02.2011
 *      Author: vyahhi
 */

#ifndef SEQ_HPP_
#define SEQ_HPP_

#include <string>

enum Nucl {
	A, T, G, C
};

char complement(char c);

template <int size> // max number of nucleotides
class Seq {
public:
	//Seq(const std::string &s);
	Seq(const char* s);
	char operator[] (const int index) const;
	std::string str() const;
private:
	char _bytes[(size >> 2) + ((size & 3) != 0)]; // little-endian
};

//template <>
//class Seq<0> {};

/*struct A {
	int _cnt;
	int _len;
	Seq<4>* _bytes;
};

struct B {
	int from, size;
	A* a;
};*/


class Sequence { // runtime length sequence (slow!!!)
public:
	Sequence(const std::string &s);
	~Sequence();
	char operator[](int index) const;
	Sequence& operator!() const;
//	SeqVarLen operator+ (const SeqVarLen &svl1, const SeqVarLen &svl2) const;
	std::string str() const;
	int len() const;
public:
	Seq<4>* _bytes;
	int _len;
	bool _reverse;
	Sequence(const Sequence *svl, bool reverse); // reverse
};

template <int size> // max number of nucleotides in each read
class MatePair {
public:
	MatePair(const char *s1, const char *s2, const int id_);
	MatePair(const MatePair &mp);
	const static MatePair<size> null;
//private:
	int id; // consecutive number from input file :)
	Seq<size> seq1;
	Seq<size> seq2;
};

template <int size>
const MatePair<size> MatePair<size>::null = MatePair<size>("", "", -1);

// ******************** //
// * TEMPLATE METHODS * //
// ******************** //

/*template <int size>
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
			this->_bytes[cur++] = byte;
			cnt = 6;
			byte = 0;
		}
	}
	if (cnt != 6) {
		this->_bytes[cur++] = byte;
	}
}*/

template <int size>
Seq<size>::Seq (const char* s) {
	char byte = 0;
	int cnt = 6;
	int cur = 0;
	for (const char* si = s; *si != 0; si++) { // unsafe!
		switch (*si) {
			case 'C': byte |= (1 << cnt); break;
			case 'G': byte |= (2 << cnt); break;
			case 'T': byte |= (3 << cnt); break;
		}
		cnt -= 2;
		if (cnt < 0) {
			this->_bytes[cur++] = byte;
			cnt = 6;
			byte = 0;
		}
	}
	if (cnt != 6) {
		this->_bytes[cur++] = byte;
	}
}

template <int size>
char Seq<size>::operator[] (const int index) const {
	switch ( ( this->_bytes[index >> 2] >> ((3-(index%4))*2) ) & 3) { // little endian!
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
MatePair<size>::MatePair(const char *s1, const char *s2, const int id_) : id(id_), seq1(s1), seq2(s2) {
	//
}

template <int size>
MatePair<size>::MatePair(const MatePair &mp) : id(mp.id), seq1(mp.seq1), seq2(mp.seq2) {
	//
}

template <int size> // max number of nucleotides
class Kmer {
public:
	Kmer(const char* s);
	char operator[] (const int index) const;
	std::string str() const;
	void PushBack(Nucl c);
	void PushFront(Nucl c);
	bool Follows(Kmer prev);
	Kmer RevCompl();
	bool Canonic();
private:

};
#endif /* SEQ_HPP_ */
