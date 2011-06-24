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
		ofstream_.open(filename);
	}

	virtual ~osequencestream() {
		ofstream_.close();
	}

	osequencestream& operator<<(const Sequence& seq) {
		ofstream_ << "> contig #" << id_++ << endl;
		string s = seq.str();
		size_t cur = 0;
		while (cur < s.size()) {
			ofstream_ << s.substr(cur, 70) << endl;
			cur += 70;
		}
		return *this;
	}


};

#endif /* OSEQUENCESTREAM_HPP_ */
