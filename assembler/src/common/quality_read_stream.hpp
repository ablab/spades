/*
 * quality_read_stream.hpp
 *
 *  Created on: Mar 22, 2011
 *      Author: sergey
 */

#ifndef QUALITY_READ_STREAM_HPP_
#define QUALITY_READ_STREAM_HPP_

#define QUAL_SHIFT 33

#include <algorithm>
#include <string>
#include <zlib.h>
#include <cstdlib>
#include <cstdarg>
#include <iostream>
#include <vector>
#include "ifaststream.hpp"
#include "strobe_read.hpp"

using namespace std;

/*
 * reads 'read's from multiple FASTQ files (simultaneously)
 *
 * skips any reads with Ns!
 */

class QualityReadStream {
private:
	//	vector<ifaststream*> ifs_;
	ifaststream* ifs_;
	bool eof_;
	bool is_open_;

	string read_;
	string qual_;

	bool ValidNucls(const string& s) {
		for (size_t j = 0; j < s.size(); ++j) { // if at least one letter isn't ACGT
			if (!is_nucl(s[j])) {
				return false;
			}
		}
		return true;
	}

	bool ValidRead() {
		if (!is_open() || eof()) {
			return false;
		}
		string name, seq, qual;
		*ifs_ >> name >> seq >> qual;

		if (ifs_->eof()) {
			eof_ = true;
		}

		bool valid = ValidNucls(seq);

		if (valid) {
			read_ = seq;
			qual_ = qual;
		}
		return valid;
	}

	void read_ahead() {
		if (!is_open()) {
			return;
		}
		while (!eof() && !ValidRead()) {
			;
		}
	}

	int ToPhred(char c) {
		return (int)c - QUAL_SHIFT;
	}

	vector<int> ConvertToPhred(const string& qual) {
		vector<int> answer;
		for (size_t i = 0; i < qual.size(); ++i) {
			answer.push_back(ToPhred(qual[i]));
		}
		return answer;
	}

public:

	QualityReadStream(const char *filename) {
		ifs_ = new ifaststream(filename);
		is_open_ = true;
		eof_ = false;
		read_ahead();
	}

	bool eof() const {
		return eof_;
	}

	bool is_open() const {
		return is_open_;
	}

	void reset() {
		ifs_->reset();
		is_open_ = true;
		eof_ = false;
		read_ahead();
	}

	virtual ~QualityReadStream() {
		close();
	}

	void close() {
		if (is_open()) {
			delete ifs_;
			is_open_ = false;
		}
	}

	pair<Sequence, vector<int> > Next() {
		pair<Sequence, vector<int> > answer = make_pair(Sequence(read_), ConvertToPhred(qual_));
		read_ahead();
		return answer;
	}

/*
	vector<pair<Sequence, vector<char> > >* readAll(int number = -1) { // read count `reads`, default: 2^32 - 1
		vector<pair<Sequence, vector<char> > > *v = new vector<pair<Sequence, vector<char> > > ();
		while (!eof_ && number--) {
			v->push_back(Next());
		}
		return v;
	}
*/
};
#endif /* QUALITY_READ_STREAM_HPP_ */
