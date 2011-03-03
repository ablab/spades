/*
 * ifastqstream.hpp
 *
 *  Created on: 03.03.2011
 *      Author: vyahhi
 */

#ifndef IFASTSTREAM_HPP_
#define IFASTSTREAM_HPP_

#include "libs/kseq/kseq.h"

using namespace std;

// STEP 1: declare the type of file handler and the read() function
KSEQ_INIT(gzFile, gzread)

class ifaststream {
public:
	ifaststream(const char* filename) {
		is_open_ = open(filename);
		if (is_open_){
			eof_ = false;
			do_read();
		}
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
			do_read(); // make actual read for the next result
			break;
		}
		state_++;
		return *this;
	}

private:
	gzFile fp_;
	kseq_t* seq_;
	bool is_open_;
	bool eof_;
	char state_; // 0 - name, 1 - seq, 2 - qual

	/*
	 * open i's file with FASTQ reads,
	 * return true if it opened file, false otherwise
	 */
	bool open(const char *filename) {
		fp_ = gzopen(filename, "r"); // STEP 2: open the file handler
		if (!fp_) {
			return false;
		}
		seq_ = kseq_init(fp_); // STEP 3: initialize seq
		return true;
	}

	void close() {
		if (is_open()) {
			kseq_destroy(seq_); // STEP 5: destroy seq
			gzclose(fp_); // STEP 6: close the file handler
		}
	}

	void do_read() {
		assert(is_open());
		assert(!eof());
		if (kseq_read(seq_) >= 0) {
			state_ = 0;
		}
		else {
			eof_ = true;
		}
	}
};

#endif /* IFASTSTREAM_HPP_ */
