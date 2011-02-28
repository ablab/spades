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


Sequence::Sequence(const std::string &s): _size(s.size()), _reverse(false) {  //accepts both 0123 and ACGT
	_bytes = (Seq<4>*) malloc((this->_size >> 2) + 1); // sizeof(Seq<4>()) == 1;
	for (int i = 0; i < _size / 4; ++i) {
		_bytes[i] = Seq<4>(s.substr(i*4, 4).c_str());
	}
	if (_size & 3) {
		_bytes[_size/4] = Seq<4>(s.substr((_size/4)*4, _size & 3).c_str());
	}
}
Sequence Sequence::shift_right() const{
	return *this;
}

Sequence Sequence::shift_left() const{
	return *this;
}
Sequence::~Sequence() {
	if (!_reverse) { // cheat, free only one memory! should be implemented with pointer counters
		free(_bytes);
	}
}

char Sequence::operator[] (const size_t index) const {
	int i = index;
	if (_reverse) {
		i = _size - i - 1;
		return complement(_bytes[i / 4][i % 4]);
	}
	else {
		return _bytes[i / 4][i % 4];
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

Sequence Sequence::substr(int start, int end) const {
	//todo rewrite when other constructors will exist
	string s = "";
	for (int i = start; i < end; ++i) {
		s += operator[](i);
	}
	return Sequence(s);
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

