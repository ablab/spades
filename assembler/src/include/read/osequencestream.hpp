//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/*
 * oreadstream.hpp
 *
 *  Created on: 23.06.2011
 *      Author: vyahhi
 */

#ifndef OSEQUENCESTREAM_HPP_
#define OSEQUENCESTREAM_HPP_

#include <fstream>

using std::string;
using std::endl;


class osequencestream {

protected:
	std::ofstream ofstream_;

	int id_;


	void write_str(const string& s) {
        size_t cur = 0;
        while (cur < s.size()) {
            ofstream_ << s.substr(cur, 60) << endl;
            cur += 60;
        }
	}

	virtual void write_header(const string& s) {
        // Velvet format: NODE_1_length_24705_cov_358.255249
	    ofstream_ << ">NODE_" << id_++ << "_length_" << s.size() << endl;
	}

public:
	osequencestream(const string& filename): id_(0) {
		ofstream_.open(filename.c_str());
	}

	virtual ~osequencestream() {
		ofstream_.close();
	}

    osequencestream& operator<<(const string& s) {
        write_header(s);
        write_str(s);
        return *this;
    }

	osequencestream& operator<<(const Sequence& seq) {
	    std::string s = seq.str();
		return operator <<(s);
	}
};


class osequencestream_cov: public osequencestream {
protected:

	double coverage_;


    virtual void write_header(const string& s) {
        // Velvet format: NODE_1_length_24705_cov_358.255249
        ofstream_ << ">NODE_" << id_++ << "_length_" << s.size() << "_cov_" << coverage_ << endl;
    }


public:
	osequencestream_cov(const string& filename): osequencestream(filename), coverage_(0.) {
	}

    virtual ~osequencestream_cov() {
        ofstream_.close();
    }

	osequencestream_cov& operator<<(double coverage) {
		coverage_ = coverage;
		return *this;
	}

	osequencestream_cov& operator<<(const string& s) {
        write_header(s);
        write_str(s);
        return *this;
    }

	osequencestream_cov& operator<<(const Sequence& seq) {
        std::string s = seq.str();
        return operator <<(s);
    }

};


class osequencestream_with_id: public osequencestream {
protected:
	int uid_;

	double cov_;

    virtual void write_header(const string& s) {
        ofstream_ << ">NODE_" << id_++ << "_length_" << s.size() << "_cov_" << cov_ << "_ID_" << uid_ << endl;
    }

public:
	osequencestream_with_id(const string& filename): osequencestream(filename), uid_(0), cov_(0.0) {
	}

    virtual ~osequencestream_with_id() {
        ofstream_.close();
    }


	void setCoverage(double c) {
		cov_ = c;
	}

	void setID(int uid) {
		uid_ = uid;
	}

	osequencestream_with_id& operator<<(const string& s) {
        write_header(s);
        write_str(s);
        return *this;
    }

	osequencestream_with_id& operator<<(const Sequence& seq) {
        std::string s = seq.str();
        return operator <<(s);
    }

};


class osequencestream_with_data_for_scaffold: public osequencestream_with_id  {
protected:
    std::ofstream scstream_;

    virtual void write_header(const string& s) {
        scstream_ << id_ << "\tNODE_" << id_ << "\t" << s.size() << "\t" << (int) round(cov_) << endl;
        ofstream_ << ">NODE_" << id_++ << "_length_" << s.size() << "_cov_" << cov_ << "_ID_" << uid_ << endl;
    }

public:
	osequencestream_with_data_for_scaffold(const string& filename): osequencestream_with_id(filename) {
	    id_ = 1;
		std::string sc_filename = filename + ".info";
		scstream_.open(sc_filename.c_str());
	}

	virtual ~osequencestream_with_data_for_scaffold() {
	    ofstream_.close();
		scstream_.close();
	}

	osequencestream_with_data_for_scaffold& operator<<(const string& s) {
        write_header(s);
        write_str(s);
        return *this;
    }

	osequencestream_with_data_for_scaffold& operator<<(const Sequence& seq) {
        std::string s = seq.str();
        return operator <<(s);
    }

};

#endif /* OSEQUENCESTREAM_HPP_ */
