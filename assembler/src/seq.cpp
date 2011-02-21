/*
 * seq.cpp
 *
 *  Created on: 21.02.2011
 *      Author: vyahhi
 */

#include <cstdlib>
#include "seq.hpp"

SeqVarLen::SeqVarLen(const std::string &s): _len(s.size()) {
	_bytes = (Seq<4>*) malloc(this->_len >> 2); // sizeof(Seq<4>()) == 1;
	for (int i = 0; i < _len / 4; ++i) {
		_bytes[i] = Seq<4>(s.substr(i*4, 4));
	}
}

SeqVarLen::~SeqVarLen() {
	free(_bytes);
}

char SeqVarLen::operator[] (const int &index) const {
	return _bytes[index / 4][index % 4];
}

std::string SeqVarLen::str() const {
	std::string res = "";
	for (int i = 0; i < this->_len; ++i) {
		res += this->operator[](i);
	}
	return res;
}

int SeqVarLen::len() const {
	return _len;
}
