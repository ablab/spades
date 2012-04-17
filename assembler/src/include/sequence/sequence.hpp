//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/**
 * @file    sequence.hpp
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
 *		data_(seq.data_), from_(from), size_(size), rtl_(rtl) {
			data_->Grab();
		}
 * Immutable runtime length sequence (someway slow)		data_(seq.data_), from_(from), size_(size), rtl_(rtl) {
			data_->Grab();
		}
 */


#ifndef SEQUENCE_HPP_
#define SEQUENCE_HPP_

#include "sequence/seq.hpp"
#include "sequence/sequence_data.hpp"
#include <vector>
#include <string>
#include <cstring>
using namespace std;

class Sequence {
private:
    SequenceData* data_;
    size_t from_;
    size_t size_;
    bool rtl_; // Right to left + complimentary (?)
    inline Sequence(const Sequence &seq, size_t from, size_t size, bool rtl);
public:
    const Sequence& operator=(const Sequence &rhs) {
        from_ = rhs.from_;
        size_ = rhs.size_;
        rtl_ = rhs.rtl_;
        if (data_ != rhs.data_) {
            data_->Release();
            data_ = rhs.data_;
            data_->Grab();
        }
        return *this;
    }

    /**
         * Sequence initialization (arbitrary size string)
         *
         * @param s ACGT or 0123-string
         */
    explicit Sequence(char* s) :
            from_(0), size_(strlen(s)), rtl_(false) {
        data_ = new SequenceData(s, size_);
        data_->Grab();
    }

    Sequence() :
             from_(0), size_(0), rtl_(false) {
         data_ = new SequenceData("", size_);
         data_->Grab();
    }

    explicit Sequence(const char* s) :
            from_(0), size_(strlen(s)), rtl_(false) {
        data_ = new SequenceData(s, size_);
        data_->Grab();
    }

    template<typename S>
    explicit Sequence(const S &s) :
            from_(0), size_(s.size()), rtl_(false) {
        data_ = new SequenceData(s, size_);
        data_->Grab();
    }

    inline Sequence(const Sequence &s);
    inline ~Sequence();

    inline char operator[](const size_t index) const;
    inline bool operator==(const Sequence &that) const;
    inline bool operator!=(const Sequence &that) const;
    inline bool operator<(const Sequence &that) const;
    inline Sequence operator!() const;

    /**
	 * @param from inclusive
	 * @param to exclusive;
	 */
    inline Sequence Subseq(size_t from, size_t to) const;
    inline Sequence Subseq(size_t from) const; // up to size_ by default
    inline Sequence operator+(const Sequence &s) const;

    /////todo what are these methods???
    inline int find(const Sequence &t, int from = 0) const;
    inline int similar(const Sequence &t, int k, char directed = 0) const;
    inline int leftSimilar(const Sequence &t, int k) const;
    inline int rightSimilar(const Sequence &t, int k) const;

    /**
         * @param from inclusive
         * @param to exclusive;
         * @return true if two sequences intersect
         */
    inline bool intersects(const Sequence &t) const;
    /////

    //	template<size_t size2_>
    //	Seq<size2_> start() const;
    //
    //	template<size_t size2_>
    //	Seq<size2_> end() const;

    template<size_t size2_>
    Seq<size2_> start() const;
    template<size_t size2_>
    Seq<size2_> end() const;
    
    inline string str() const;
    inline size_t size() const;
};

inline ostream& operator<<(ostream& os, const Sequence& s);

/**
 * @todo optimize
 */
template<size_t size2_>
Seq<size2_> Sequence::start() const {
    VERIFY(size2_ <= size_);
    return Seq<size2_> (*this);
}

/**
 * @todo optimize
 */
template<size_t size2_>
Seq<size2_> Sequence::end() const {
    VERIFY(size2_ <= size_);
    return Seq<size2_> (*this, size_ - size2_);
}

/**
 * @class SequenceBuilder
 * @section DESCRIPTION
 *
 * Class was created for build sequence. It is included method: size(), append()
 */


class SequenceBuilder {
    vector<char> buf_;
public:
    template<typename S>
    SequenceBuilder& append(const S &s) {
        for (size_t i = 0; i < s.size(); ++i) {
            buf_.push_back(s[i]);
        }
        return *this;
    }

    SequenceBuilder& append(char c) {
        buf_.push_back(c);
        return *this;
    }

    Sequence BuildSequence() {
        return Sequence(buf_);
    }

    size_t size() const {
        return buf_.size();
    }

    char operator[](const size_t index) const {
        VERIFY(index < buf_.size());
        return buf_[index];
	}

	string str() const {
		string s(buf_.size(), '-');
		for (size_t i = 0; i < s.size(); ++i) {
			s[i] = nucl(buf_[i]);
		}
		return s;
    }
};

/////////////////////////////////////////////////////////
// impl

Sequence::Sequence(const Sequence &seq, size_t from, size_t size, bool rtl) :
	data_(seq.data_), from_(from), size_(size), rtl_(rtl) {
	data_->Grab();
}

Sequence::Sequence(const Sequence &s) :
	data_(s.data_), from_(s.from_), size_(s.size_), rtl_(s.rtl_) {
	data_->Grab();
}

Sequence::~Sequence() {
	data_->Release();
}

char Sequence::operator[](const size_t index) const {
	VERIFY(index >= 0);
	VERIFY(index < size_);
	if (rtl_) {
		int i = from_ + size_ - 1 - index;
		return complement(data_->operator[](i));
	} else {
		int i = from_ + index;
		return data_->operator[](i);
	}
}

bool Sequence::operator==(const Sequence &that) const {
	if (size_ != that.size_) {
		return false;
	}
	if (data_ == that.data_ && from_ == that.from_ && rtl_ == that.rtl_) {
		return true;
	}
	for (size_t i = 0; i < size_; ++i) {
		if (this->operator[](i) != that[i]) {
			return false;
		}
	}
	return true;
}

bool Sequence::operator!=(const Sequence &that) const {
	return !(*this == that);
}

bool Sequence::intersects(const Sequence &t) const {
	for (size_t i = 0; i < min(size_, t.size_); ++i) {
		if (this->operator[](i) == t[i]) {
			return true;
		}
	}
	return false;
}

/**
  * @todo Might be optimized via int comparison (not so easy)
  */
bool Sequence::operator<(const Sequence &that) const {
	size_t s = min(size_, that.size_);
	for (size_t i = 0; i < s; ++i) {
		if (this->operator[](i) != that[i]) {
			return (this->operator[](i) < that[i]);
		}
	}
	return (size_ < that.size_);
}

Sequence Sequence::operator!() const {
	return Sequence(*this, from_, size_, !rtl_);
}

// O(1)
//including from, excluding to
//safe if not #DEFINE NDEBUG
Sequence Sequence::Subseq(size_t from, size_t to) const {
//	cerr << endl<<"subseq:" <<   from <<" " << to << " " <<  this->str() << endl;
	VERIFY(to >= from);
	VERIFY(from >= 0);
	VERIFY(to <= size_);
	//VERIFY(to - from <= size_);
	if (rtl_) {
		return Sequence(*this, from_ + size_ - to, to - from, true);
	} else {
		return Sequence(*this, from_ + from, to - from, false);
	}
}

//including from, excluding to
Sequence Sequence::Subseq(size_t from) const {
	return Subseq(from, size_);
}

/**
* @todo : must be KMP or hashing instead of this shit
*/
int Sequence::find(const Sequence &t, int from) const {
	for (size_t i = from; i <= size() - t.size(); i++) {
		if (Subseq(i, i + t.size()) == t) {
			return i;
		}
	}
	return -1;
}
/**
 *
 *@param k  minimal intersection of sequences
 *@param directed  LEFT means that after intersection t continues to left over _this and matches perfectly with _this on overlaping
 *@return 0 - undirected similarity, 1: t extends this to right, -1: this extends t
 *
 */
int Sequence::similar(const Sequence &t, int k, char directed) const {
//	cerr << endl << t.str()<< "similar started" <<k << endl;
	int result = 0;
	if (directed != -1)
		result |= rightSimilar(t, k);
	if (directed != 1)
		result |= leftSimilar(t, k);
	return result;
}

int Sequence::leftSimilar(const Sequence &t, int k) const {
	return t.rightSimilar(*this, k);
}

int Sequence::rightSimilar(const Sequence &t, int k) const {
	int tsz = t.size();
	int sz = size();
	Sequence d(t.Subseq(0, k));
	for (int res = find(d, 0); res != -1; res = find(d, res + 1)) {
		if (res + tsz < sz)
			continue;
		int i;
		for (i = k; i + res < sz; i++) {
			if (t[i] != this->operator[](i + res)) {
				break;
			};
		}
		if (i == sz - res)
			return 1;
	}
	return 0;
}

/**
* @todo optimize
  */
Sequence Sequence::operator+(const Sequence &s) const {
	return Sequence(str() + s.str());
	// TODO might be opposite to correct
	//	int total = size_ + s.size_;
	//	std::vector<Seq<4> > bytes((total + 3) >> 2);
	//	for (size_t i = 0; i < size_; ++i) {
	//		bytes[i / 4] = (bytes[i / 4] << operator [](i)); // TODO :-) use <<=
	//	}
	//	for (size_t i = 0, j = size_; i < s.size_; ++i, ++j) {
	//		bytes[j / 4] = (bytes[j / 4]) << s[i];
	//	}
	//	return Sequence(new Data(bytes), 0, total, false);
}

std::string Sequence::str() const {
	std::string res(size_, '-');
	for (size_t i = 0; i < size_; ++i) {
		res[i] = nucl(this->operator[](i));
	}
	return res;
}

ostream& operator<<(ostream& os, const Sequence& s) {
	os << s.str();
	return os;
}

size_t Sequence::size() const {
	return size_;
}

#endif /* SEQUENCE_HPP_ */
