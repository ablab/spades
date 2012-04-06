//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/**
 * @file    seq.hpp
 * @author  vyahhi
 * @version 1.0
 *
 * @section LICENSE
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of
 * the License, or (at your option) any later version.
 *
 * @section DESCRIPTION
 *
 * Immutable ACGT-sequence with compile-time size.
 * It compress sequence to array of Ts (default: char).
 */

#ifndef SEQ_HPP_
#define SEQ_HPP_

#include <string>
#include "verify.hpp"
#include <array>
#include <algorithm>
#include "sequence/nucl.hpp"
#include "log.hpp"
#include <cstring>
#include <iostream>

/**
 * @param T is max number of nucleotides, type for storage
 */
template<size_t size_, typename T = u_int64_t>
class Seq {
private:
	/**
	 * @variable Number of bits in type T (e.g. 8 for char)
	 * @example 8: 2^8 = 256 or 16
	 */
	const static size_t Tbits = sizeof(T) << 3;

	/**
	 * @variable Number of nucleotides that can be stored in one type T (e.g. 4 for char)
     * Tnucl MUST be a power of two
	 * @example 4: 8/2 = 4 or 16/2 = 8
	 */
	const static size_t Tnucl = Tbits >> 1;

	/**
	 * @variable Number of bits in Tnucl (e.g. 2 for char). Useful for shifts instead of divisions.
	 */
	const static size_t Tnucl_bits = log_<Tnucl, 2>::value;

	/**
	 * @variable Number of Ts which required to store all sequence.
	 */
	const static size_t data_size_ = (size_ + Tnucl - 1) >> Tnucl_bits;

    
    /* *
     * @variable Just some prime number to count the hash function of the kmer
     * */    
    const static size_t prime_num = 239;
	

    /**
	 * @variable Inner representation of sequence: array of Ts with length = data_size_.
	 *
	 * @invariant Invariant: all nucleotides >= size_ are 'A's (useful for comparison)
	 */
	std::array<T, data_size_> data_;

    //size_t seq_hash_;

	friend class Seq<size_ - 1, T> ;

	/**
	 * Initialize data_ array of this object with C-string
	 *
	 * @param s C-string (ACGT chars only), strlen(s) = size_
	 */
	void init(const char* s) {
        //seq_hash_ = prime_num;
		T data = 0;
		size_t cnt = 0;
		int cur = 0;
		for (size_t pos = 0; pos < size_; ++pos, ++s) { // unsafe!
			// VERIFY(is_nucl(*s)); // for performance
			data = data | ((T) dignucl(*s) << cnt);
			cnt += 2;
			if (cnt == Tbits) {
				//seq_hash_ = ((seq_hash_ << 5) - seq_hash_) + data;
                this->data_[cur++] = data;
				cnt = 0;
				data = 0;
			}
		}
		if (cnt != 0) {
			this->data_[cur++] = data;
			//seq_hash_ = ((seq_hash_ << 5) - seq_hash_) + data;
		}
		VERIFY(*s == 0); // C-string always ends on 0
	}

	/**
	 * Initialize hash function for the Seq with the data_ given
	 *
	 */
    //void init_hash() {
        ////Assuming the data is okay and already filled
        //size_t seq_hash = prime_num;
        //for (size_t i = 0; i < data_size_; ++i) {
            //seq_hash = ((seq_hash << 5) - seq_hash) + data_[i];
        //}
        //seq_hash_ = seq_hash;
    //}

	/**
	 * Constructor from std::array
	 */
    //Seq(std::array<T, data_size_> data) : data_(data) {
        //WARN("constructor from array");
        //init_hash();
    //}


	/**
	 * Sets i-th symbol of Seq with 0123-char
	 */
	inline void set(const size_t i, char c) {
		data_[i >> Tnucl_bits] = (data_[i >> Tnucl_bits] & ~((T) 3 << ((i & (Tnucl - 1)) << 1))) | ((T) c << ((i & (Tnucl - 1)) << 1));
	}

public:

    /**
     *  Reads sequence from the file (in the same format as BinWrite writes it)
     *  and returns false if error occured, true otherwise.
     */
	static bool BinRead(std::istream& file, Seq<size_> *seq) {
	    //std::cerr << "SEQ: " << sizeof(T) << " " << data_size_ << " " << sizeof(T) * data_size_ << std::endl;
		file.read((char *) seq->data_.data(), sizeof(T) * data_size_);
		return !file.fail();
	}

	/**
     *  Writes sequence to the file (in the same format as BinRead reads it)
     *  and returns false if error occured, true otherwise.
     */
	static bool BinWrite(std::ostream& file, const Seq<size_> &seq) {
	    //std::cerr << "SEQ: " << sizeof(T) << " " << data_size_ << " " << sizeof(T) * data_size_ << std::endl;
		file.write((const char *) seq.data_.data(), sizeof(T) * data_size_);
		return !file.fail();
	}

    /**
     *  Reads sequence from the file (in the same format as BinWrite writes it)
     *  and returns false if error occured, true otherwise.
     */
    bool BinRead(std::istream& file) {
        return BinRead(file, this);
    }

    /**
     *  Writes sequence to the file (in the same format as BinRead reads it)
     *  and returns false if error occured, true otherwise.
     */
    bool BinWrite(std::ostream& file) {
        return BinWrite(file, *this);
    }

	/**
	 * Default constructor, fills Seq with A's
	 */
	Seq() {
		//VERIFY((T)(-1) >= (T)0);//be sure to use unsigned types
		std::fill(data_.begin(), data_.end(), 0);
        //init_hash();
	}

	Seq(const char* s) {
		//VERIFY((T)(-1) >= (T)0);//be sure to use unsigned types
		init(s);
	}

	/**
	 * Ultimate constructor from ACGT0123-string.
	 *
	 * @param s Any object with operator[], which returns 0123 chars
	 * @param offset Offset when this sequence starts
     * @number_to_read A number of nucleotides, we want to fetch from this string
     * @warning assuming that s is a correct string, filled with ACGT _OR_ 0123 
     * no init method, filling right here
	 */
	template<typename S>
	explicit Seq(const S &s, size_t offset = 0, size_t number_to_read = size_) : data_() {
		
        TRACE("New Constructor for seq " << s[0] << " is first symbol");
        VERIFY(is_dignucl(s[0]) || is_nucl(s[0]));

        // which symbols does our string contain : 0123 or ACGT?
        bool digit_str = is_dignucl(s[0]); 

        // data -- one temporary variable corresponding to the i-th array element
        // and some counters
        T data = 0;
		size_t cnt = 0;
		size_t cur = 0;

		for (size_t i = 0; i < number_to_read; ++i) {
            VERIFY(is_dignucl(s[i]) || is_nucl(s[i]));
            
            // we fill everything with zeros (As) by default. 
		    char c = digit_str ? s[offset + i] : dignucl(s[offset + i]);
            
	        data = data | (T(c) << cnt);
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

        for (; cur < data_size_; ++cur)
            this->data_[cur] = 0;

        //init_hash();
		
        ////init
		//for (size_t pos = 0; pos < size_; ++pos, ++s) { // unsafe!
			//data = data | ((T) dignucl(*s) << cnt);
			//cnt += 2;
			//if (cnt == Tbits) {
				//this->data_[cur++] = data;
				//cnt = 0;
				//data = 0;
			//}
		//}
		//if (cnt != 0) {
			//this->data_[cur++] = data;
		//}
        
        // init_hash
	}

	/**
	 * Get i-th symbol of Seq.
	 *
	 * @param i Index of the symbol (0 <= i < size_)
	 * @return 0123-char on position i
	 */
	char operator[](const size_t i) const {
		//VERIFY(i >= 0);
		//VERIFY(i < size_);
		return (data_[i >> Tnucl_bits] >> ((i & (Tnucl - 1)) << 1)) & 3; // btw (i % Tnucl) <=> (i & (Tnucl-1))
	}

	/**
	 * Reverse complement.
	 *
	 * @return Reverse complement Seq.
	 */
	Seq<size_, T> operator!() const {
		Seq<size_, T> res(*this);
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
		//res.init_hash();
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
		//VERIFY(is_dignucl(c));
		Seq<size_, T> res(*this);
        std::array<T, data_size_>& data = res.data_;
		if (data_size_ != 0) { // unless empty sequence
			T rm = data[data_size_ - 1] & 3;
			T lastnuclshift_ = ((size_ + Tnucl - 1) & (Tnucl - 1)) << 1;
			data[data_size_ - 1] = (data[data_size_ - 1] >> 2) | ((T) c << lastnuclshift_);
			if (data_size_ >= 2) { // if we have at least 2 elements in data
			    for (int i = data_size_ - 2; i >= 0; --i){
					T new_rm = data[i] & 3;
					data[i] = (data[i] >> 2) | (rm << (Tbits - 2)); // we need & here because if we shift negative, it fill with ones :(
					rm = new_rm;
			    }
			}
		}
        //res.init_hash();
		return res;
	}

	Seq<size_ + 1, T> pushBack(char c) const {
		if (is_nucl(c)) {
			c = dignucl(c);
		}
		//VERIFY(is_dignucl(c));
		Seq<size_ + 1, T> s;
		copy(this->data_.begin(), this->data_.end(), s.data_.begin());
		s.data_[s.data_size_ - 1] = s.data_[s.data_size_ - 1] | ((T) c << ((size_ & (Tnucl - 1)) << 1));
        //s.init_hash();

        return s; //was: Seq<size_ + 1, T>(str() + nucl(c));

	}

	/**
	 * @todo optimize!!!
	 */
	Seq<size_ + 1, T> pushFront(char c) const {
		if (is_nucl(c)) {
			c = dignucl(c);
		}
		VERIFY(is_dignucl(c));
        return Seq<size_ + 1, T> (nucl(c) + str());
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
		VERIFY(is_dignucl(c));
		Seq<size_, T> res(*this);
		T rm = c;
		for (size_t i = 0; i < data_size_; ++i) {
			T new_rm = (res.data_[i] >> (Tbits - 2)) & 3;
			res.data_[i] = (res.data_[i] << 2) | rm;
			rm = new_rm;
		}
		if ((size_ & (Tnucl - 1)) != 0) {
			T lastnuclshift_ = (size_ & (Tnucl - 1)) << 1;
			res.data_[data_size_ - 1] = res.data_[data_size_ - 1] & (((T) 1
					<< lastnuclshift_) - 1);
		}
		//res.init_hash();
        return res;
	}

	bool operator==(const Seq<size_, T>& s) const {
		return 0 == memcmp(data_.data(), s.data_.data(), sizeof(T) * data_size_);
	}

	/**
	 * @see operator ==()
	 */

	bool operator!=(const Seq<size_, T>& s) const {
		return 0 != memcmp(data_.data(), s.data_.data(), sizeof(T) * data_size_);
	}

	//	/*
	//	 * now usual order, but some linear order on Seq which works fast
	//	 */
	//	bool operator<(const Seq<size_, T> that) const {
	//		return 0 > memcmp(data_.data(), that.data_.data(), sizeof(T) * data_size_);
	//	}


	/**
	 * String representation of this Seq
	 *
	 * @return ACGT-string of length size_
	 * @see nucl()
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

	/**
	 * @see Seq
	 */
	template<size_t size2_, typename T2 = T>
	Seq<size2_, T2> start() const {
		VERIFY(size2_ <= size_);
		return Seq<size2_, T2> (*this);
	}

	template<size_t size2_/* = size_ - 1*/, typename T2 = T>
	Seq<size2_, T2> end() const {
		VERIFY(size2_ <= size_);
		return Seq<size2_, T2> (*this, size_ - size2_);
	}

	char last() const {
		return operator[](size_ - 1);
	}

	char first() const {
		return operator[](0);
	}

    size_t GetHash() const {
        size_t hash = prime_num;
        for (size_t i = 0; i < data_size_; i++) {
            hash = ((hash << 5) - hash) + data_[i];
        }
        return hash;
    }

    struct hash {
        //size_t operator()(const Seq<size_, T>& seq) const {
            //return seq.seq_hash_; 
        //}
        size_t operator()(const Seq<size_, T>& seq) const {
            size_t hash = prime_num;
            for (size_t i = 0; i < seq.data_size_; i++) {
                hash = ((hash << 5) - hash) + seq.data_[i];
            }
            return hash;
        }
	};

	struct multiple_hash {
		size_t operator()(const Seq<size_, T>& seq, size_t hash_num, 
				size_t h) const {
            WARN("using multiple hash");
			++hash_num;
			for (size_t i = 0; i < seq.data_size_; i++) {
				h = (h << hash_num) + seq.data_[i];
			}
			return h;
		}
	};

	struct equal_to {
		bool operator()(const Seq<size_, T>& l, const Seq<size_, T>& r) const {
			return memcmp(l.data_.data(), r.data_.data(), sizeof(T) * data_size_) == 0;
		}
	};

	struct less2 {
		int operator()(const Seq<size_, T> &l, const Seq<size_, T> &r) const {
			for (size_t i = 0; i < size_; ++i) {
				if (l[i] != r[i]) {
					return (l[i] < r[i]);
				}
			}
			return false;
		}
	};

	//	/**
	//	 * Denotes some (weird) order on k-mers. Works fast.
	//	 */
	//	struct less {
	//		int operator()(const Seq<size_> &l, const Seq<size_> &r) const {
	//			return 0 > memcmp(data_.data(), that.data_.data(), sizeof(T) * data_size_);
	//		}
	//	};

};

template<size_t size_, typename T = int>
ostream& operator<<(ostream& os, Seq<size_, T> seq) {
	os << seq.str();
	return os;
}

#endif /* SEQ_HPP_ */
