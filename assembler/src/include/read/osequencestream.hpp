//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
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

class osequencestream {
private:
	ofstream ofstream_;
	int id_;
public:
	osequencestream(const string& filename): id_(0) {
		ofstream_.open(filename.c_str());
	}

	virtual ~osequencestream() {
		ofstream_.close();
	}

	osequencestream& operator<<(const Sequence& seq) {
//		DEBUG("outputting");
		string s = seq.str();
		ofstream_ << ">NODE_" << id_++ << "_length_" << s.size() << endl;
		// Velvet format: NODE_1_length_24705_cov_358.255249
		size_t cur = 0;
		while (cur < s.size()) {
			ofstream_ << s.substr(cur, 60) << endl;
			cur += 60;
		}
		return *this;
	}
};

class osequencestream_cov {
private:
	ofstream ofstream_;
	int id_;
	double coverage_;
public:
	osequencestream_cov(const string& filename): id_(0), coverage_(0.) {
		ofstream_.open(filename.c_str());
	}

	virtual ~osequencestream_cov() {
		ofstream_.close();
	}

	osequencestream_cov& operator<<(double coverage) {
		coverage_ = coverage;
		return *this;
	}

	osequencestream_cov& operator<<(const Sequence& seq) {
		string s = seq.str();
		ofstream_ << ">NODE_" << id_++ << "_length_" << s.size() << "_cov_" << coverage_ << endl;
		// Velvet format: NODE_1_length_24705_cov_358.255249
		size_t cur = 0;
		while (cur < s.size()) {
			ofstream_ << s.substr(cur, 60) << endl;
			cur += 60;
		}
		return *this;
	}
};


class osequencestream_with_id {
private:
	ofstream ofstream_;
	int id_;

	void* uid_;
	double cov_;

public:
	osequencestream_with_id(const string& filename): id_(0), uid_(0), cov_(0.0) {
		ofstream_.open(filename.c_str());
	}

	virtual ~osequencestream_with_id() {
		ofstream_.close();
	}

	void setCoverage(double c) {
		cov_ = c;
	}

	void setID(void* uid) {
		uid_ = uid;
	}

	osequencestream_with_id& operator<<(const Sequence& seq) {
		string s = seq.str();
		ofstream_ << ">NODE_" << id_++ << "_length_" << s.size() << "_cov_" << cov_ << "_ID_" << uid_ << endl;
		size_t cur = 0;
		while (cur < s.size()) {
			ofstream_ << s.substr(cur, 60) << endl;
			cur += 60;
		}
		return *this;
	}
};


class osequencestream_with_data_for_scaffold {
private:
	ofstream ofstream_;
	ofstream scstream_;
	int id_;

	int uid_;
	double cov_;

public:
	osequencestream_with_data_for_scaffold(const string& filename): id_(1), uid_(0), cov_(0.0) {
		ofstream_.open(filename.c_str());
		string sc_filename = filename + ".info";
		scstream_.open(sc_filename);
	}

	virtual ~osequencestream_with_data_for_scaffold() {
		ofstream_.close();
		scstream_.close();
	}

	void setCoverage(double c) {
		cov_ = c;
	}

	void setID(size_t uid) {
		uid_ = uid;
	}

	osequencestream_with_data_for_scaffold& operator<<(const Sequence& seq) {
		string s = seq.str();
		scstream_ << id_ << "\tNODE_" << id_ << "\t" << s.size() << "\t" << (int) round(cov_) << endl;
		ofstream_ << ">NODE_" << id_++ << "_length_" << s.size() << "_cov_" << cov_ << "_ID_" << uid_ << endl;

		size_t cur = 0;
		while (cur < s.size()) {
			ofstream_ << s.substr(cur, 60) << endl;
			cur += 60;
		}
		return *this;
	}
};

#endif /* OSEQUENCESTREAM_HPP_ */
