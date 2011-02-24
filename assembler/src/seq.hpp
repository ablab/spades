/*
 * seq.hpp
 *
 *  Created on: 21.02.2011
 *      Author: vyahhi
 */

#ifndef SEQ_HPP_
#define SEQ_HPP_

#include <string>

char complement(char c); // 0123 -> 3210

char nucl(char c); // 0123 -> ACGT

class Sequence;

template <int size> // max number of nucleotides
class Seq {
public:
	Seq(const char* s);
	Seq(const Seq<size> &s);
	char operator[] (const int index) const;
	Seq<size> operator!() const;
	Seq<size> shift_right(char c) const; // char should be 0123
	Seq<size> shift_left(char c) const; // char should be 0123
	Sequence substring(int from, int to) const;
	std::string str() const;
	static int len();
private:
	const static int _bytelen = (size >> 2) + ((size & 3) != 0);
	char _bytes[_bytelen]; // little-endian
};

class Sequence { // runtime length sequence (slow!!!)
public:
	Sequence(const std::string &s);
	~Sequence();
	char operator[](int index) const;
	Sequence& operator!() const;
//	SeqVarLen operator+ (const SeqVarLen &svl1, const SeqVarLen &svl2) const;
	std::string str() const;
	int len() const;
private:
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

template <int size>
Seq<size>::Seq (const Seq<size> &s) {
	memcpy(_bytes, s._bytes, (size >> 2) + ((size & 3) != 0));
}

template <int size>
Seq<size>::Seq (const char* s) {
	char byte = 0;
	int cnt = 6;
	int cur = 0;
	for (const char* si = s; *si != 0; si++) { // unsafe!
		switch (*si) {
			case 'C': case '1': byte |= (1 << cnt); break;
			case 'G': case '2': byte |= (2 << cnt); break;
			case 'T': case '3': byte |= (3 << cnt); break;
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
	return ((_bytes[index >> 2] >> ((3-(index%4))*2) ) & 3);
}

template <int size>
std::string Seq<size>::str() const {
	std::string res = "";
	for (int i = 0; i < size; ++i) {
		res += nucl((*this)[i]);
	}
	return res;
}

template <int size>
int Seq<size>::len() {
	return size;
}


template <int size>
Seq<size> Seq<size>::operator! () const {
	Sequence s = Sequence(this->str());
	return Seq<size>((!s).str());
}

template <int size>
Seq<size> Seq<size>::shift_right(char c) const {
	Seq<size> res = *this;
	c <<= (((4-(size%4))%4)*2); // omg >.<
	for (int i = _bytelen - 1; i >= 0; --i) {
		char rm = (_bytes[i] & 192) >> 6;
		res._bytes[i] <<= 2;
		res._bytes[i] &= 252;
		res._bytes[i] |= c;
		c = rm;
	}
	return res;
}

template <int size>
Seq<size> Seq<size>::shift_left(char c) const {
	Seq<size> res = *this;
	for (int i = 0; i < _bytelen; ++i) {
		char rm = _bytes[i] & 3;
		res._bytes[i] >>= 2;
		res._bytes[i] &= 63;
		res._bytes[i] |= (c << 6);
		c = rm;
	}
	return res;
}


template <int size>
Sequence Seq<size>::substring(int from, int to) const {
	std::string s = str();
	s = s.substr(from, to);
	return Sequence(s);
}

template <int size>
MatePair<size>::MatePair(const char *s1, const char *s2, const int id_) : id(id_), seq1(s1), seq2(s2) {
	//
}

template <int size>
MatePair<size>::MatePair(const MatePair &mp) : id(mp.id), seq1(mp.seq1), seq2(mp.seq2) {
	//
}


/*
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
	bool IsCanonic();
private:

};*/

#endif /* SEQ_HPP_ */
