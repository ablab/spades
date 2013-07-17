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

	virtual ofastastream& operator<<(const SingleRead& read) {
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

class RCRemovingOFastaStream: public ofastastream {
    typedef ofastastream base;
    size_t cnt_;
public:
    RCRemovingOFastaStream(const string& filename): base(filename), cnt_(0) {
	}

	/*virtual*/ RCRemovingOFastaStream& operator<<(const SingleRead& read) {
	    if (++cnt_ % 2 == 1 /*cnt ^= 1*/) {
	        base::operator<<(read);
	    }
	    return *this;
	}
};

}
