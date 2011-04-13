/*
 * sequence.cpp
 *
 *  Created on: 01.03.2011
 *      Author: vyahhi
 */

#include "sequence.hpp"
#include "nucl.hpp"

#include <ostream>

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
	assert(index >= 0);
	assert(index < size_);
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

bool Sequence::operator<(const Sequence &that) const {
	for (size_t i = 0; i < size_; ++i) {
		if (i > that.size_) return true;
		if (this->operator[](i) < that[i]) return true;
	}
	if (size_ == that.size_) return false;
	return true;
}

Sequence Sequence::operator!() const {
	return Sequence(*this, from_, size_, !rtl_);
}

// O(1)
//including from, excluding to
//safe if not #DEFINE NDEBUG
Sequence Sequence::Subseq(size_t from, size_t to) const {
//	cerr << endl<<"subseq:" <<   from <<" " << to << " " <<  this->str() << endl;
	assert(to >= from);
	assert(from >= 0);
	assert(to <= size_);
	//assert(to - from <= size_);
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

//TODO: must be KMP or hashing instead of this shit
int Sequence::find(const Sequence &t, int from) const {
	for (size_t i = from; i <= size() - t.size(); i++) {
		if (Subseq(i, i + t.size()) == t) {
			return i;
		}
	}
	return -1;
}
/*
 *@param k  minimal intersection of sequences
 *@param directed  LEFT means that after intersection t continues to left over _this and matches perfectly with _this on overlaping
 *
 *
 */
// 0 - undirected similarity, 1: t extends this to right, -1: this extends t
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

// TODO optimize
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

