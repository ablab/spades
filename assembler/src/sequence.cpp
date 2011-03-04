/*
 * sequence.cpp
 *
 *  Created on: 01.03.2011
 *      Author: vyahhi
 */

#include "sequence.hpp"
#include "nucl.hpp"

Sequence::Sequence(Data* data, size_t from, size_t size, bool rtl) :
	data_(data), from_(from), size_(size), rtl_(rtl) {
	data_->IncRef();
}

Sequence::Sequence(const Sequence &s) :
	data_(s.data_), from_(s.from_), size_(s.size_), rtl_(s.rtl_) {
	data_->IncRef();
}

void Sequence::init(const std::string &s) { //accepts both 0123 and ACGT
	from_ = 0;
	size_ = s.size();
	rtl_ = false;
	std::vector<Seq<4> > inner_data((size_ + 3) >> 2);
	for (size_t i = 0; i < size_ / 4; ++i) {
		inner_data[i] = Seq<4> (s.substr(i * 4, 4).c_str());
	}
	if (size_ & 3) {
		inner_data[size_ / 4] = Seq<4> (
				s.substr((size_ / 4) * 4, size_ & 3).c_str());
	}
	data_ = new Data(inner_data);
}

Sequence::Sequence(const std::string &s) {
	init(s);
}

Sequence::~Sequence() {
	if (--data_->ref_ == 0) {
		delete data_;
	}
}

char Sequence::operator[](const size_t index) const {
	if (rtl_) {
		int i = from_ + size_ - 1 - index;
		return complement(data_->bytes_[i / 4][i % 4]);
	} else {
		int i = from_ + index;
		return data_->bytes_[i / 4][i % 4];
	}
}

bool Sequence::operator==(const Sequence &that) const {
	if (size_ != that.size_) {
		return false;
	}

	for (size_t i = 0; i < size_; ++i) {
		if (operator[](i) != that[i]) {
			return false;

		}
	}
	return true;
}

const Sequence Sequence::operator!() const {
	return Sequence(data_, from_, size_, !rtl_);
}

//including from, excluding to
//safe if not #DEFINE NDEBUG
Sequence Sequence::Subseq(size_t from, size_t to) const {
	assert(to >= from);
	assert(from >= 0);
	assert(to - from <= size_);
	if (rtl_) {
		return Sequence(data_, from_ + size_ - to, to - from, true);
	} else {
		return Sequence(data_, from_ + from, to - from, false);
	}
}
//TODO: must be KMP or hashing instead of this shit
int Sequence::find (const Sequence &t, int from) const{
	for(size_t i = from; i <= size()- t.size(); i++) {
		if (Subseq(i, i + t.size()) == t) {
			return i;
		}
	}
	return -1;
}
/*int Sequence::findAfter(const Sequence &t, int pos){


 return -1;
 }*/
//0 - undirected similarity, 1: t extends this to right, -1: this extends t
int Sequence::similar(const Sequence &t, int k, char directed) const {
	Sequence c(Subseq(0, k));

	//	cerr << endl <<"trying 1" << endl;
	int tsz = t.size();
	int sz = size();
	int res;
	int i;
	int from = 0;
	if (directed != -1) {
		Sequence d(t.Subseq(0, k));
		//	cerr << endl <<"trying 2" << endl;
		while ((res = find(d, from)) != -1) {
			//		cerr <<"Res=" << res <<endl;
			from = res + 1;
			for (i = k; i < sz - res; i++) {
				if (i >= tsz)
					break;
				//			cerr << t[i] + '0'<< " ";
				//			cerr << this->operator[](i + res) + '0' << endl;
				if (t[i] != this->operator[](i + res)) {
					break;
				};
			}
			if (i == sz - res)
				return 1;
		}
	}
	from = 0;

	if (directed != 1) {
		while ((res = t.find(c, from)) != -1) {
			//		cerr <<"Res=" << res <<endl;
			from = res + 1;
			for (i = k; i < tsz - res; i++) {
				if (i >= sz)
					break;
				//			cerr << t[i +res]+'0'<< " ";
				//			cerr << this->operator[](i) +'0'<< endl;
				if (t[i + res] != this->operator[](i)) {
					break;
				}
			}
			if (i == tsz - res)
				return 1;

		}
	}
	return 0;
}

// TODO optimize
// TODO might be opposite to correct
Sequence Sequence::operator+(const Sequence &s) const {
	int total = size_ + s.size_;
	std::vector<Seq<4> > bytes((total + 3) >> 2);
	for (size_t i = 0; i < size_; ++i) {
		bytes[i / 4] = (bytes[i / 4] << operator [](i)); // TODO :-) use <<=
	}
	for (size_t i = 0, j = size_; i < s.size_; ++i, ++j) {
		bytes[j / 4] = (bytes[j / 4]) << s[i];
	}
	return Sequence(new Data(bytes), 0, total, false);
}

std::string Sequence::Str() const {
	std::string res(size_, '-');
	for (size_t i = 0; i < size_; ++i) {
		res[i] = nucl(operator[](i));
	}
	return res;
}

size_t Sequence::size() const {
	return size_;
}
