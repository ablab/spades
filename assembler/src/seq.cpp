/*
 * seq.cpp
 *
 *  Created on: 21.02.2011
 *      Author: vyahhi
 */

#include <iostream>
#include <cstdlib>
#include "seq.hpp"

Sequence::Sequence(const std::string &s): _size(s.size()), _reverse(false) {
	_bytes = (Seq<4>*) malloc(this->_size >> 2); // sizeof(Seq<4>()) == 1;
	for (int i = 0; i < _size / 4; ++i) {
		_bytes[i] = Seq<4>(s.substr(i*4, 4).c_str());
	}
}

Sequence::~Sequence() {
	if (!_reverse) { // cheat, free only one memory! should be implemented with pointer counters
		free(_bytes);
	}
}

char Sequence::operator[] (int index) const {
	if (_reverse) {
		index = _size - index - 1;
		return complement(_bytes[index / 4][index % 4]);
	}
	else {
		return _bytes[index / 4][index % 4];
	}
}

// TODO rewrite this!
bool Sequence::operator== (const Sequence &that) const {
	return this->str() == that.str();
}

std::string Sequence::str() const {
	std::string res = "";
	for (int i = 0; i < this->_size; ++i) {
		res += nucl((*this)[i]);
	}
	return res;
}

int Sequence::size() const {
	return _size;
}

Sequence& Sequence::operator! () const {
	Sequence* res = new Sequence(this, true);
	return *res;
}

Sequence::Sequence(const Sequence *svl, bool reverse = false): _bytes(svl->_bytes), _size(svl->_size), _reverse(svl->_reverse) {
	if (reverse) {
		_reverse = !_reverse;
	}
}

char complement(char c) {
	return 3 - c;
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

