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
	void * ptr;

	osequencestream(const string& filename): id_(0), ptr(0) {
		ofstream_.open(filename.c_str());
	}

	virtual ~osequencestream() {
		ofstream_.close();
	}

	osequencestream& operator<<(const Sequence& seq) {
//		DEBUG("outputting");
		string s = seq.str();
		ofstream_ << ">NODE_" << id_++ << "_length_" << s.size() << "_ID_" << ptr << endl;
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
	osequencestream_cov(const string& filename): id_(0) {
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

#endif /* OSEQUENCESTREAM_HPP_ */
