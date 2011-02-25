/*
 * seq.hpp
 *
 *  Created on: 21.02.2011
 *      Author: vyahhi
 */

#ifndef SEQ_HPP_
#define SEQ_HPP_

#include <string>
#include <functional>
#include <memory>

char complement(char c); // 0123 -> 3210
char nucl(char c); // 0123 -> ACGT

template <size_t _size> // max number of nucleotides
class Seq;

// *****************************************

class Sequence { // runtime length sequence (slow!!!)
public:
	Sequence(const std::string &s);
	~Sequence();
	char operator[](int index) const;
	bool operator==(const Sequence &that) const;
	Sequence& operator!() const;
	std::string str() const;
	int size() const;
private:
	Seq<4>* _bytes;
	int _size;
	bool _reverse;
	Sequence(const Sequence *svl, bool reverse); // reverse
};

// *****************************************

template <size_t _size> // max number of nucleotides
class Seq {
private:
	const static size_t _byteslen = (_size >> 2) + ((_size & 3) != 0);
	char _bytes[_byteslen]; // little-endian
public:
	Seq() {}; // random Seq, use with care!

	Seq(const char* s) {
		char byte = 0;
		int cnt = 6;
		int cur = 0;
		for (size_t pos = 0; pos < _size; ++pos, ++s) { // unsafe!
			switch (*s) {
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

	Seq(const Seq<_size> &s) {
		memcpy(_bytes, s._bytes, (_size >> 2) + ((_size & 3) != 0));
	}

	char operator[] (const size_t index) const { // 0123
		return ((_bytes[index >> 2] >> ((3-(index%4))*2) ) & 3);
	}

	Seq<_size> operator!() const { // TODO: optimize
		Sequence s = Sequence(this->str());
		return Seq<_size>((!s).str());
	}

	Seq<_size> shift_right(char c) const { // char should be 0123
		Seq<_size> res = *this;
		c <<= (((4-(_size%4))%4)*2); // omg >.<
		for (int i = _byteslen - 1; i >= 0; --i) {
			char rm = (_bytes[i] & 192) >> 6;
			res._bytes[i] <<= 2;
			res._bytes[i] &= 252;
			res._bytes[i] |= c;
			c = rm;
		}
		return res;
	}

	Seq<_size> shift_left(char c) const { // char should be 0123
		Seq<_size> res = *this;
		for (size_t i = 0; i < _byteslen; ++i) {
			char rm = _bytes[i] & 3;
			res._bytes[i] >>= 2;
			res._bytes[i] &= 63;
			res._bytes[i] |= (c << 6);
			c = rm;
		}
		return res;
	}

	Sequence substring(int from, int to) const { // TODO: optimize
		std::string s = str();
		s = s.substr(from, to);
		return Sequence(s);
	}

	// string representation of Seq - only for debug and output purposes
	std::string str() const {
		std::string res = "";
		for (int i = 0; i < _size; ++i) {
			res += nucl((*this)[i]);
		}
		return res;
	}

	static int size() {
		return _size;
	}

	struct hash {
		 size_t operator() (const Seq<_size> &seq) const {
			size_t h = 0;
			for (size_t i = 0; i < _byteslen; ++i) {
				h += seq._bytes[i];
			}
			return h;
		}
	};

	struct equal_to {
		bool operator() (const Seq<_size> &l, const Seq<_size> &r) const {
			 return 0 == memcmp(l._bytes, r._bytes, _byteslen);
		}
	};

	struct less {
		int operator() (const Seq<_size> &l, const Seq<_size> &r) const {
			 return 0 > memcmp(l._bytes, r._bytes, _byteslen);
		}
	};

};

/*namespace std { // hacking standard hash function for Seq<_size>
	namespace tr1 { // hacking standard hash function for Seq<_size>
		template <>
		template <int _size>
		struct hash<Seq<_size> > {
			size_t operator()(const Seq<_size> &seq) const{
				return 1;
			}
		};
	}
}*/

// *****************************************


template <int _size> // max number of nucleotides in each read
class MatePair {
public:
	MatePair(const char *s1, const char *s2, const int id_) : id(id_), seq1(s1), seq2(s2) {};
	MatePair(const MatePair &mp) : id(mp.id), seq1(mp.seq1), seq2(mp.seq2) {};
	const static MatePair<_size> null = MatePair<_size>("", "", -1);
private:
	int id; // consecutive number from input file :)
	Seq<_size> seq1;
	Seq<_size> seq2;
};

// *****************************************

#endif /* SEQ_HPP_ */
