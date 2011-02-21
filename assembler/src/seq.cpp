/*
 * seq.cpp
 *
 *  Created on: 21.02.2011
 *      Author: vyahhi
 */

#include <iostream>
#include <cstdlib>
#include "seq.hpp"

SeqVarLen::SeqVarLen(const std::string &s): _len(s.size()), _reverse(false) {
	_bytes = (Seq<4>*) malloc(this->_len >> 2); // sizeof(Seq<4>()) == 1;
	for (int i = 0; i < _len / 4; ++i) {
		_bytes[i] = Seq<4>(s.substr(i*4, 4));
	}
}

SeqVarLen::~SeqVarLen() {
	if (!_reverse) { // cheat, free only one memory! should be implemented with pointer counters
		free(_bytes);
	}
}

char SeqVarLen::operator[] (int index) const {
	if (_reverse) {
		index = _len - index - 1;
		return complement(_bytes[index / 4][index % 4]);
	}
	else {
		return _bytes[index / 4][index % 4];
	}
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

SeqVarLen& SeqVarLen::operator! () const {
	SeqVarLen* res = new SeqVarLen(this, true);
	return *res;
}

SeqVarLen::SeqVarLen(const SeqVarLen *svl, bool reverse = false): _bytes(svl->_bytes), _len(svl->_len), _reverse(svl->_reverse) {
	if (reverse) {
		this->_reverse = !this->_reverse;
	}
}
