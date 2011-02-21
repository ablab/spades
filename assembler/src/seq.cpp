/*
 * seq.cpp
 *
 *  Created on: 21.02.2011
 *      Author: vyahhi
 */

#include "seq.hpp"

SeqVarLen::SeqVarLen(const std::string &s, int len) {
	// TODO
}

char SeqVarLen::operator[] (const int &index) const {
	// TODO
	return 'A';
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
