/*
 * Sequence class with compile-time size (immutable)
 *
 *  Created on: 20.02.2011
 *      Author: vyahhi
 */

#ifndef SEQ_HPP_
#define SEQ_HPP_

#include <string>
#include <cassert>
#include <array>
#include <algorithm>
#include "nucl.hpp"
#include "log.hpp"

template<size_t size_, typename T = char> // max number of nucleotides, type for storage
class Seq {
private:
	// compile-time static constants.
	const static size_t Tbits = sizeof(T) << 3; // ex. 8: 2^8 = 256 or 16
	const static size_t Tnucl = Tbits >> 1; // ex. 4: 8/2 = 4 or 16/2 = 8
	const static size_t Tnucl_bits = log_<Tnucl,2>::value; // ex. 2: 2^2 = 4 or 3: 2^3 = 8
	const static size_t data_size_ = (size_ + Tnucl - 1) >> Tnucl_bits;
	std::array<T,data_size_> data_; // invariant: all nucleotides >= size_ are 'A's

	/*
	 * gets ACGT-string and initialize data_ array on the object
	 */
	void init(const char* s) {
		T data = 0;
		size_t cnt = 0;
		int cur = 0;
		for (size_t pos = 0; pos < size_; ++pos, ++s) { // unsafe!
			assert(is_nucl(*s));
			data = data | ((T)denucl(*s) << cnt);
			cnt += 2;
			if (cnt == Tbits) {
				this->data_[cur++] = data;
				cnt = 0;
				data = 0;
			}
		}
		assert(*s == 0);//c string always ends on 0
		if (cnt != 0) {
			this->data_[cur++] = data;
		}
	}

	Seq(std::array<T, data_size_> data) :
		data_(data) {
	}
	;

	void set(const size_t index, char c) {
		data_[index >> Tnucl_bits] = ( data_[index >> Tnucl_bits]
								   & ~((T)3 << ((index % Tnucl) << 1)) )
								   | ((T)c << ((index % Tnucl) << 1));
	}

public:
	Seq() {
		std::fill(data_.begin(), data_.end(), 0); // fill with A-s
	}

	Seq(const char* s) {
		init(s);
	}

	Seq(const Seq<size_,T> &seq): data_(seq.data_) {
	 //does memcpy faster?
	}

	template <typename S> Seq(const S& s, size_t offset = 0) { // TODO: optimize
		char a[size_ + 1];
		for (size_t i = 0; i < size_; ++i) {
			a[i] = nucl(s[offset + i]);
		}
		a[size_] = 0;
		init(a);
	}

	template<size_t _bigger_size, T>
	Seq(const Seq<_bigger_size, T>& seq) {// TODO: optimize (Kolya)
		assert(_bigger_size > size_);
		init(seq.str().substr(0, size_).c_str());
	}

	char operator[](const size_t index) const { // 0123
		assert(index >= 0);
		assert(index < size_);
		T ind = index;
		return (data_[ind >> Tnucl_bits] >> ((ind % Tnucl) << 1)) & 3;
	}

	/*
	 * reverse complement from the Seq
	 */
	Seq<size_,T> operator!() const { // TODO: optimize
		Seq<size_, T> res(data_);
		for(size_t i = 0; i < (size_ >> 1); ++i) {
			T front = complement(res[i]);
			T end = complement(res[size_ - 1 - i]);
			res.set(i, end);
			res.set(size_ - 1 - i, front);
		}
		if ((size_ & 1) == 1) {
			res.set(size_ >> 1, complement(res[size_ >> 1]));
		}
		// can be made without complement calls, but with xor on all bytes afterwards.
		return res;
	}

	/**
	 * add one nucl to the right, shifting seq to the left
	 */
	Seq<size_, T> operator<<(char c) const {
		//todo talk with vyahhi about this
		assert(is_nucl(c) || c < 4);
		Seq<size_, T> res(data_);
		if (data_size_ != 0) { // unless empty sequence
			T rm = res.data_[data_size_ - 1] & 3;
			T lastnuclshift_ = ((size_ + Tnucl - 1) % Tnucl) << 1;
			res.data_[data_size_ - 1] = (res.data_[data_size_ - 1] >> 2) | ((T)(is_nucl(c) ? denucl(c) : c) << lastnuclshift_);
			if (data_size_ >= 2) { // if we have at least 2 elements in data
				size_t i = data_size_ - 1;
				do {
					--i;
					T new_rm = res.data_[i] & 3;
					res.data_[i] = ( (res.data_[i] >> 2) & (((T)1 << (Tbits - 2)) - 1) ) | (rm << (Tbits - 2));	// we need & here because if we shift negative, it fill with ones :(
					rm = new_rm;
				} while (i != 0);
			}
		}
		return res;
	}

	/**
	 * add one nucl to the left, shifting seq to the right
	 */
	Seq<size_,T> operator>>(char c) {	// TODO: better name
		assert(is_nucl(c));
        Seq<size_, T> res(data_);
		T rm = denucl(c);
		for (size_t i = 0; i < data_size_; ++i) {
			T new_rm = (res.data_[i] >> (Tbits - 2)) & 3;
			res.data_[i] = (res.data_[i] << 2) | rm;
			rm = new_rm;
		}
		if (size_ % Tnucl != 0) {
			T lastnuclshift_ = ((size_  % Tnucl) + 1) << 1;
			res.data_[data_size_ - 1] = res.data_[data_size_ - 1] & (((T)1 << lastnuclshift_) - 1);
		}
        return res;
	}

	bool operator==(const Seq<size_, T> s) const {	// TODO: optimize
		return s.data_ == data_;
		//return this->equal_to()(s);
	}

	// string representation of Seq - only for debug and output purposes
	std::string str() const {
		std::string res(size_, '-');
		for (size_t i = 0; i < size_; ++i) {
			res[i] = nucl(operator[](i));
		}
		return res;
	}

	static size_t size() {
		return size_;
	}

	struct hash {
		size_t operator()(const Seq<size_> &seq) const {
			size_t h = 0;
			for (size_t i = 0; i < data_size_; ++i) {
				h += seq.data_[i];
			}
			return h;
		}
	};

	struct equal_to {
		bool operator()(const Seq<size_> &l, const Seq<size_> &r) const {
			return l.data_ == r.data_;
			//return 0 == memcmp(l._bytes.data(), r._bytes.data(), _byteslen);
		}
	};

	struct less {
		int operator()(const Seq<size_> &l, const Seq<size_> &r) const {
			return l.data_ < r.data_;
			//return 0 > memcmp(l._bytes.data(), r._bytes.data(), _byteslen);
		}
	};

	template<int size2, typename T2 = char>
	Seq<size2, T2> head() { // TODO: optimize (Kolya)
		std::string s = str();
		return Seq<size2, T2> (s.substr(0, size2).c_str());
	}

	template<int size2, typename T2 = char>
	Seq<size2, T2> tail() const { // TODO: optimize (Kolya)
		std::string s = str();
		return Seq<size2, T2> (s.substr(size_ - size2, size2).c_str());
	}

};

// *****************************************
// LEGACY CODE

/*template <typename T2>
 Seq(const T2& t, size_t offset = 0) {
 char a[size_];
 for (size_t i = 0; i < size_; ++i) {
 a[i] = t[offset + i];
 }
 init(a);
 }*/

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
//
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

#endif /* SEQ_HPP_ */
