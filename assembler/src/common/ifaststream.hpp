/*
 * ifastqstream.hpp
 *
 *  Created on: 03.03.2011
 *      Author: vyahhi
 */

#ifndef IFASTSTREAM_HPP_
#define IFASTSTREAM_HPP_

#include "libs/kseq/kseq.h"
#include <zlib.h>
#include <cassert>

using namespace std;

// STEP 1: declare the type of file handler and the read() function
// KSEQ_INIT(gzFile, gzread)

/*
 * Read name, seq and qual strings from FASTQ data (one by one)
 */
class ifaststream {

public:
	ifaststream(const string& filename) {
		filename_ = filename;
		is_open_ = open(filename);
		assert(is_open_); // Fails if there is no such file
	}

	virtual ~ifaststream() {
		close();
	}

	bool is_open() const {
		return is_open_;
	}

	bool eof() const {
		return eof_;
	}

	ifaststream& operator>>(string &s) {
		if (!is_open() || eof()) {
			return *this;
		}
		// if there is 'N' in sequence, then throw out this mate read
		switch (state_) {
		case 0:
			s = seq_->name.s;
			break;
		case 1:
			s = seq_->seq.s;
			break;
		case 2:
			s = seq_->qual.s;
			read_ahead(); // make actual read for the next result
			break;
		}
		state_ = (state_ + 1) % 3; // next state
		return *this;
	}

	void close() {
		if (is_open()) {
			kseq_destroy(seq_); // STEP 5: destroy seq
			gzclose(fp_); // STEP 6: close the file handler
			is_open_ = false;
		}
	}

	void reset() {
		close();
		open(filename_);
	}

private:
	std::string filename_;
	gzFile fp_;
	kseq_t* seq_;
	bool is_open_;
	bool eof_;
	char state_; // 0 - name, 1 - seq, 2 - qual

	/*
	 * open i's file with FASTQ reads,
	 * return true if it opened file, false otherwise
	 */
	bool open(string filename) {
		fp_ = gzopen(filename.c_str(), "r"); // STEP 2: open the file handler
		if (!fp_) {
			return false;
		}
		is_open_ = true;
		seq_ = kseq_init(fp_); // STEP 3: initialize seq
		eof_ = false;
		state_ = 2;
		read_ahead();
		state_ = 0;
		return true;
	}

	void read_ahead() {
		assert(is_open());
		assert(!eof());
		assert(state_ == 2);
		if (kseq_read(seq_) < 0) {
			eof_ = true;
		}
	}
};

#endif /* IFASTSTREAM_HPP_ */
