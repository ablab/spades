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
#include <iostream> // for debug
#include <cstring>
#include <cassert>
#include <cstring>
#include <array>
#include <vector>
#include "nucl.hpp"

typedef long long word;

using namespace std;

template <size_t size_, typename T = char> // max number of nucleotides, type for storage
class Seq {
private:
	// compile-time static constants.
	const static size_t Tbits = sizeof(T) << 3; // ex. 8: 2^8 = 256 or 16
	const static size_t Tnucl = Tbits >> 1; // ex. 4: 8/2 = 4 or 16/2 = 8
	const static size_t Tnucl_bits = 1 + sizeof(T); // ex. 2: 2^2 = 4 or 3: 2^3 = 8
	const static size_t data_size_ = (size_ + Tnucl - 1) >> Tnucl_bits;
	std::array<T,data_size_> data_; // 0 bits overhead

	//const static char _lastnuclshift = ((size_ + 3) & 3) << 1;

	void init(const char* s) {
		T data = 0;
		size_t cnt = 0;
		int cur = 0;
		for (size_t pos = 0; pos < size_; ++pos, ++s) { // unsafe!
			switch (*s) {
				case 'C': case '1': case 1: data |= (1 << cnt); break;
				case 'G': case '2': case 2: data |= (2 << cnt); break;
				case 'T': case '3': case 3: data |= (3 << cnt); break;
			}
			cnt += 2;
			if (cnt == Tbits) {
				this->data_[cur++] = data;
				cnt = 0;
				data = 0;
			}
		}
		if (cnt != 0) {
			this->data_[cur++] = data;
		}
	}

	Seq(std::array<T,data_size_> bytes): data_(bytes) {};

public:
	Seq() {}; // random Seq, use with care!

	Seq(const char* s) {
		init(s);
	}

	Seq(const Seq<size_> &seq) : data_(seq.data_) {
		//memcpy?
	}

	template <size_t _bigger_size>
	Seq(const Seq<_bigger_size>& seq) {
		assert(_bigger_size > size_);
		std::copy(data_, seq.data, size_);
		//memcpy(data_, seq._bytes, data_size_); ?
	}

	template <typename T2>
	Seq(const T2& t, size_t offset = 0) {
		char a[size_];
		for (size_t i = 0; i < size_; ++i) {
			a[i] = t[offset + i];
		}
		init(a);
	}

	char operator[] (const size_t index) const { // 0123
		int ind = index >> Tnucl_bits;
		return (data_[ind] >> ((index % Tnucl)*2)) & 3;
	}

	Seq<size_> operator!() const { // TODO: optimize
		// TODO!!!
		return *this;
		//Sequence s = Sequence(this->str());
		//return Seq<_size>((!s).Str());
	}

//	// add nucleotide to the right
//	Seq<_size> shift_right(char c) const { // char should be 0123
//		assert(c <= 3);
//		Seq<_size> res = *this; // copy constructor
//		c <<= (((4-(_size%4))%4)*2); // omg >.<
//		for (int i = _byteslen - 1; i >= 0; --i) { // don't make it size_t :)
//			char rm = (res._bytes[i] >> 6) & 3;
//			res._bytes[i] <<= 2;
//			//res._bytes[i] &= 252;
//			res._bytes[i] |= c;
//			c = rm;
//		}
//		return res;
//	}

	//todo optimize via machine words;
	/**
	 * add one nucl to the right, shifting seq
	 */
	Seq<size_> operator<<(char c) const {
		/*std::array<T,data_size_> new_a(data_);
		char buf = new_a[data_size_ - 1] & 3;
		new_a[data_size_ - 1] >>= 2;
		new_a[data_size_ - 1] |= c << _lastnuclshift;
		for (size_t i = data_size_ - 2; i >= 0; --i) { // bug!
			char new_buf = new_a[i] & 3;
			new_a[i] >>= 2;
			new_a[i] |= buf << 6;
			buf = new_buf;
		}
		return Seq<size_>(new_a);*/
		// TODO!
		return *this;
	}

//	// add nucleotide to the left
//	Seq<_size> shift_left(char c) const { // char should be 0123
//		Seq<_size> res = *this; // copy constructor
//		// TODO: clear last nucleotide
//		for (size_t i = 0; i < _byteslen; ++i) {
//			char rm = res._bytes[i] & 3;
//			res._bytes[i] >>= 2;
//			//res._bytes[i] &= 63;
//			res._bytes[i] |= (c << 6);
//			c = rm;
//		}
//		return res;
//	}

	// string representation of Seq - only for debug and output purposes
	std::string str() const {
		std::string res(size_, '-');
		for (size_t i = 0; i < size_; ++i) {
			res[i] = nucl(this->operator[](i));
		}
		return res;
	}

	static int size() {
		return size_;
	}

	struct hash {
		 size_t operator() (const Seq<size_> &seq) const {
			size_t h = 0;
			for (size_t i = 0; i < data_size_; ++i) {
				h += seq.data_[i];
			}
			return h;
		}
	};

	struct equal_to {
		bool operator() (const Seq<size_> &l, const Seq<size_> &r) const {
			return l.data_ == r.data_;
			//return 0 == memcmp(l._bytes.data(), r._bytes.data(), _byteslen);
		}
	};

	struct less {
		int operator() (const Seq<size_> &l, const Seq<size_> &r) const {
			return l.data_ < r.data_;
			//return 0 > memcmp(l._bytes.data(), r._bytes.data(), _byteslen);
		}
	};

	template <int _size2>
	Seq<_size2> head() { // TODO: optimize (Kolya)
		std::string s = str();
		return Seq<_size2>(s.substr(0, _size2).c_str());
	}

	template <int _size2>
	Seq<_size2> tail() const { // TODO: optimize (Kolya)
		std::string s = str();
		return Seq<_size2>(s.substr(size_ - _size2, _size2).c_str());
	}

};

// *****************************************


template <int _size> // max number of nucleotides in each read
class MatePair {
public:
	MatePair(const char *s1, const char *s2, const int id_) : id(id_), seq1(s1), seq2(s2) {};
	MatePair(const MatePair &mp) : id(mp.id), seq1(mp.seq1), seq2(mp.seq2) {};
	const static MatePair<_size> null;
public: // make private!
	int id; // consecutive number from input file :)
	Seq<_size> seq1;
	Seq<_size> seq2;
};

template <int _size>
const MatePair<_size> MatePair<_size>::null = MatePair<_size>("", "", -1);

// *****************************************

#endif /* SEQ_HPP_ */
