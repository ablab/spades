/*
 * strobe_read.hpp
 *
 *  Created on: 03.03.2011
 *      Author: vyahhi
 */

#ifndef STROBE_READ_HPP_
#define STROBE_READ_HPP_

#include "seq.hpp"
using namespace std;

template <int size, int cnt = 1, typename T = char>
class strobe_read {
public:
	strobe_read(array<Seq<size,T>,cnt> data, const int id = -1) : id_(id), data_(data) {
		;
	}
	strobe_read() : id_(-1) {
		;
	}
	Seq<size,T> get(int i) {
		return data_[i];
	}
	void put(int i, const Seq<size,T> &s) {
		assert(i < cnt);
		data_[i] = s;
	}
	bool valid() {
		return id != -1;
	}
	int id() {
		assert(valid());
		return id;
	}
	void invalidate() {
		id = -1;
	}
private:
	int id_;
	array<Seq<size,T>,cnt> data_;
};

#endif /* STROBE_READ_HPP_ */
