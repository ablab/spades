/*
 *
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

using namespace std;

#define HASH_SEED 239

/**
 * Immutable ACGT-sequence with compile-time size.
 * It compress sequence to array of Ts (default: char).
 */
template<size_t size_, typename T = int> // max number of nucleotides, type for storage
class Seq {
private:
	/**
	 * Number of bits in type T (e.g. 8 for char)
	 */
	const static size_t Tbits = sizeof(T) << 3; // ex. 8: 2^8 = 256 or 16

	/**
	 * Number of nucleotides that can be stored in one type T (e.g. 4 for char)
	 */
	const static size_t Tnucl = Tbits >> 1; // ex. 4: 8/2 = 4 or 16/2 = 8

	/**
	 * Number of bits in Tnucl (e.g. 2 for char). Useful for shifts instead of divisions.
	 */
	const static size_t Tnucl_bits = log_<Tnucl, 2>::value;

	/**
	 * Number of Ts which required to store all sequence.
	 */
	const static size_t data_size_ = (size_ + Tnucl - 1) >> Tnucl_bits;

	/**
	 * Inner representation of sequence: array of Ts with length = data_size_.
	 *
	 * Invariant: all nucleotides >= size_ are 'A's (useful for comparison)
	 */
	std::array<T, data_size_> data_;

	/**
	 * Initialize data_ array of this object with C-string
	 *
	 * @param s C-string (ACGT chars only), strlen(s) = size_
	 */
	void init(const char* s) {
		T data = 0;
		size_t cnt = 0;
		int cur = 0;
		for (size_t pos = 0; pos < size_; ++pos, ++s) { // unsafe!
			assert(is_nucl(*s));
			data = data | ((T) dignucl(*s) << cnt);
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
		assert(*s == 0); // C-string always ends on 0
	}

	/**
	 * Constructor from std::array
	 */
	Seq(std::array<T, data_size_> data) :
		data_(data) {
		;
	}

	/**
	 * Sets i-th symbol of Seq with 0123-char
	 */
	inline void set(const size_t i, char c) {
		data_[i >> Tnucl_bits] = (data_[i >> Tnucl_bits] & ~((T) 3 << ((i
				% Tnucl) << 1))) | ((T) c << ((i % Tnucl) << 1));
	}

public:

	/*
	 * Default constructor, fills Seq with A's
	 */
	Seq() {
		std::fill(data_.begin(), data_.end(), 0);
	}

	/*
	 * Copy constructor
	 */
	Seq(const Seq<size_, T> &seq) :
		data_(seq.data_) {
	}

	Seq(const char* s) {
		init(s);
	}

	/**
	 * Ultimate constructor from ACGT0123-string.
	 *
	 * @param s Any object with operator[], which returns 0123 chars
	 * @param offset Offset when this sequence starts
	 */
	template<typename S>
	Seq(const S &s, size_t offset = 0) {
	  //		assert(size_ + offset <= s.size());
		char a[size_ + 1];
		for (size_t i = 0; i < size_; ++i) {
			char c = s[offset + i];
			assert(is_nucl(c) || is_dignucl(c));
			if (is_dignucl(c)) {
				c = nucl(c);
			}
			a[i] = c;
		}
		a[size_] = 0;
		init(a);
	}

	/**
	 * Get i-th symbol of Seq.
	 *
	 * @param i Index of the symbol (0 <= i < size_)
	 * @return 0123-char on position i
	 */
	char operator[](const size_t i) const {
		assert(i >= 0);
		assert(i < size_);
		return (data_[i >> Tnucl_bits] >> ((i & (Tnucl - 1)) << 1)) & 3; // btw (i % Tnucl) <=> (i & (Tnucl-1))
	}

	/**
	 * Reverse complement.
	 *
	 * @return Reverse complement Seq.
	 */
	Seq<size_, T> operator!() const {
		Seq<size_, T> res(data_);
		for (size_t i = 0; i < (size_ >> 1); ++i) {
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
	 * Shift left
	 *
	 * @param c New 0123 char which should be added to the right.
	 * @return Shifted (to the left) sequence with 'c' char on the right.
	 */
	Seq<size_, T> operator<<(char c) const {
		if (is_nucl(c)) {
			c = dignucl(c);
		}
		assert(is_dignucl(c));
		Seq<size_, T> res(data_);
		if (data_size_ != 0) { // unless empty sequence
			T rm = res.data_[data_size_ - 1] & 3;
			T lastnuclshift_ = ((size_ + Tnucl - 1) % Tnucl) << 1;
			res.data_[data_size_ - 1] = (res.data_[data_size_ - 1] >> 2)
					| ((T) (c) << lastnuclshift_);
			if (data_size_ >= 2) { // if we have at least 2 elements in data
				size_t i = data_size_ - 1;
				do {
					--i;
					T new_rm = res.data_[i] & 3;
					res.data_[i] = ((res.data_[i] >> 2) & (((T) 1
							<< (Tbits - 2)) - 1)) | (rm << (Tbits - 2)); // we need & here because if we shift negative, it fill with ones :(
					rm = new_rm;
				} while (i != 0);
			}
		}
		return res;
	}

	/**
	 * Shift right
	 *
	 * @param c New 0123 char which should be added to the left.
	 * @return Shifted (to the right) sequence with 'c' char on the left.
	 */
	Seq<size_, T> operator>>(char c) const {
		if (is_nucl(c)) {
			c = dignucl(c);
		}
		assert(is_dignucl(c));
		Seq<size_, T> res(data_);
		T rm = c;
		for (size_t i = 0; i < data_size_; ++i) {
			T new_rm = (res.data_[i] >> (Tbits - 2)) & 3;
			res.data_[i] = (res.data_[i] << 2) | rm;
			rm = new_rm;
		}
		if (size_ % Tnucl != 0) {
			T lastnuclshift_ = (size_ % Tnucl) << 1;
			res.data_[data_size_ - 1] = res.data_[data_size_ - 1] & (((T) 1
					<< lastnuclshift_) - 1);
		}
		return res;
	}

	bool operator==(const Seq<size_, T> s) const {
		return s.data_ == data_;
		//return this->equal_to()(s);
	}

	bool operator!=(const Seq<size_, T> s) const {
		return s.data_ != data_;
		//return this->equal_to()(s);
	}

	/**
	 * String representation of this Seq
	 *
	 * @return ACGT-string of length size_
	 */
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

	template<size_t size2_, typename T2 = T>
	Seq<size2_,T2> start() const {
		assert(size2_ <= size_);
		return Seq<size2_,T2> (*this);
	}

	template<size_t size2_, typename T2 = T>
	Seq<size2_,T2> end() const {
		assert(size2_ <= size_);
		return Seq<size2_,T2> (*this, size_ - size2_);
	}

	//	template<size_t HASH_SEED>
	struct hash {
		size_t operator()(const Seq<size_> &seq) const {
			size_t h = HASH_SEED;
			for (size_t i = 0; i < seq.data_size_; i++) {
				h = ((h << 5) - h) + seq.data_[i];
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

};

template<size_t size_, typename T = int>
ostream& operator<<(ostream& os, Seq<size_, T> seq) {
	os << seq.str();
	return os;
}

// *****************************************
// LEGACY CODE

/*
 template<size_t _bigger_size, T>
 Seq(const Seq<_bigger_size, T>& seq) {
 assert(_bigger_size > size_);
 init(seq.str().substr(0, size_).c_str());
 }

 template<int size2, typename T2 = char>
 Seq<size2, T2> head() {
 std::string s = str();
 return Seq<size2, T2> (s.substr(0, size2).c_str());
 }

 template<int size2, typename T2 = char>
 Seq<size2, T2> tail() const {
 std::string s = str();
 return Seq<size2, T2> (s.substr(size_ - size2, size2).c_str());
 }

 */

/*
 * Constructor from ACGT C-string
 */
/*Seq(const char* s) {
 init(s);
 }*/

/*
 * Constructor from ACGT std::string
 */
/*Seq(std::string s) {
 init(s.c_str());
 }*/

#endif /* SEQ_HPP_ */
