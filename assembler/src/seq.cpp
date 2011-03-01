/*
 * seq.cpp
 *
 *  Created on: 21.02.2011
 *      Author: vyahhi
 */

#include <iostream>
#include <cstdlib>
#include "seq.hpp"
#include <string>

using namespace std;

Sequence::Sequence(Data* data, size_t from, size_t size, bool rtl) :
	data_ (data),
	from_ (from),
	size_ (size),
	rtl_ (rtl)
{
	data_->IncRef();
}

Sequence::Sequence(const Sequence &s) :
	data_ (s.data_),
	from_ (s.from_),
	size_ (s.size_),
	rtl_ (s.rtl_)
{
	data_->IncRef();
}

Sequence::Sequence(const std::string &s): from_(0), size_(s.size()), rtl_(false) {  //accepts both 0123 and ACGT
	vector< Seq<4> > inner_data((size_ + 3) >> 2);
	for (size_t i = 0; i < size_ / 4; ++i) {
		inner_data[i] = Seq<4>(s.substr(i * 4, 4).c_str());
	}
	if (size_ & 3) {
		inner_data[size_ / 4] = Seq<4>(s.substr((size_ / 4) * 4, size_ & 3).c_str());
	}
	data_ = new Data(inner_data);
}

Sequence::~Sequence() {
	if (--data_->ref_ == 0) {
		delete data_;
	}
}

char Sequence::operator[] (const size_t index) const {
	if (rtl_) {
		int i = from_ + size_ - 1 - index;
		return complement(data_->bytes_[i / 4][i % 4]);
	}
	else {
		int i = from_ + index;
		return data_->bytes_[i / 4][i % 4];
	}
}

bool Sequence::operator== (const Sequence &that) const {
	if (size_ != that.size_) {
		return false;
	}
	for (size_t i = 0; i < size_; ++i) {
		if (operator [](i) != that[i]) {
			return false;
		}
	}
	return true;
}

const Sequence Sequence::operator! () const {
	return Sequence(data_, from_, size_, !rtl_);
}

Sequence Sequence::Subseq(size_t from, size_t to) const {
	if (rtl_) {
		return Sequence(data_, from_ + size_ - to, to - from, true);
	} else {
		return Sequence(data_, from_ + from, to - from, false);
	}
}
//TODO: must be KMP or hashing instead of this shit
int Sequence::find (const Sequence &t) {
	for(int i = 0; i < size()- t.size() + 1; i++) {
		if (Subseq(i, i + t.size()) == t) {
			return 1;
		}
	}
	return 0;
}
int similar(const Sequence &a, const Sequence &b, int k){
	Sequence c(a.Subseq(0, k));
//	if (b.find(c)) {
		return 1;
//	}
	Sequence d(b.Subseq(0, k));
//	if (a.find(d)) {
		return 1;
//	}
	return 0;
}

// TODO optimize
// TODO might be opposite to correct
Sequence Sequence::operator+ (const Sequence &s) const {
	int total = size_ + s.size_;
	vector< Seq<4> > bytes((total + 3) >> 2);
	for (size_t i = 0; i < size_; ++i) {
		bytes[i / 4] = bytes[i / 4].shift_left((operator [](i))); // TODO :-) use <<=
	}
	for (size_t i = 0, j = size_; i < s.size_; ++i, ++j) {
		bytes[j / 4] = (bytes[j / 4]) << s[i];
	}
	return Sequence(new Data(bytes), 0, total, false);
}

std::string Sequence::Str() const {
	std::string res = "";
	for (size_t i = 0; i < size_; ++i) {
		res += nucl(operator [](i));
	}
	return res;
}

size_t Sequence::size() const {
	return size_;
}

char complement(char c) {
	return c ^ 3;
}

char nucl(char c) {
	switch(c) {
		case 0: return 'A';
		case 1: return 'C';
		case 2: return 'G';
		case 3: return 'T';
		default: return 'N';
	}
}

