//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include "single_read.hpp"

namespace io {

class ofastastream {
private:
	ofstream ofstream_;
public:
	ofastastream(const string& filename): ofstream_(filename.c_str()) {
	}

	ofastastream& operator<<(const SingleRead& read) {
		ofstream_ << ">" << read.name() << endl;
		size_t cur = 0;
		string s = read.GetSequenceString();
		while (cur < s.size()) {
			ofstream_ << s.substr(cur, 60) << endl;
			cur += 60;
		}
		return *this;
	}
};

}
