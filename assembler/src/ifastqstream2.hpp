/*
 * ifastqstream.hpp
 *
 *  Created on: 03.03.2011
 *      Author: vyahhi
 */

#ifndef IFASTQSTREAM_HPP_
#define IFASTQSTREAM_HPP_

#include <string>
#include <cstdarg>
#include <cstdio>
#include <iostream>
#include <vector>
#include "libs/kseq/kseq.h"
#include "strobe_read.hpp"

using namespace std;

// STEP 1: declare the type of file handler and the read() function
KSEQ_INIT(gzFile, gzread)

template <size_t size, size_t cnt = 1, typename T = char>
class ifastqstream : basic_istream<Seq<size,T> > { // TODO: inheritance from basic_istream
private:
	array<gzFile,cnt> fp;
	array<kseq_t*,cnt> seq;
	bool eofbit;
	bool failbit;
	int id;

	/*
	 * open i's file with FASTQ reads,
	 * return true if it opened file, false otherwise
	 */
	bool open(const string &filename, int i) {
		fp[i] = gzopen(filename.c_str(), "r"); // STEP 2: open the file handler
		if (!fp[i]) {
			//cerr << "File " << filename << " not found!" << endl;
			return false;
		}
		seq[i] = kseq_init(fp[i]); // STEP 3: initialize seq
		return true;
	}

	/*
	 * close all FASTQ files
	 */
	void close() {
		assert(!fail());
		for (size_t i = 0; i < cnt; ++i) {
			kseq_destroy(seq[i]); // STEP 5: destroy seq
			gzclose(fp[i]); // STEP 6: close the file handler
		}
		failbit = true;
	}

	/*
	 * do actual read, it's one sequence ahead of read()
	 */
	void do_read() {
		assert(good());
		for (size_t i = 0; i < cnt; ++i) {
			if (kseq_read(seq[i]) < 0) {
				eofbit = true;
			}
		}
	}
public:
	ifastqstream(const char *filename, ...) {
		va_list ap;
		va_start(ap, filename);
		for (size_t i = 0; i < cnt; ++i) {
			if (!open(filename, i)) {
				failbit = true;
				return;
			}
			filename = va_arg(ap, const char *);
		}
		va_end(ap);
		failbit = false;
		eofbit = false;
		id = 0;
		do_read();
	}

	~ifastqstream() {
		close();
	}

	bool fail() const {
		return failbit;
	}

	bool eof() const {
		return eofbit;
	}

	bool good() const {
		return !fail() && !eof();
	}

	ifastqstream<size,cnt>& operator<<(strobe_read<size,cnt,T> &sr) {
		assert(good());
		// invariant: actual read from file was already done by do_read();
		// assert that they are mate pairs (have: name/1 and name/2)!
		// assert(strncmp(_seq1->name.s, _seq2->name.s, _seq1->name.l - 1) == 0);
		// assume that all paired reads have the same length
		// assert(_seq1->seq.l == _seq2->seq.l);

		// fill result
		for (size_t i = 0; i < cnt; ++i) {
			// if there is 'N' in sequence, then throw out this mate read
			for (char *si = seq[i]->seq.s; *si != 0; ++si) {
				if (*si != 'A' && *si != 'C' && *si != 'G' && *si != 'T') {
					sr.invalidate();
					break;
				}
			}
			if (!sr.valid()) {
				break;
			}
			sr.put(i, Seq<size,T>(seq[i]->seq.s));
		}
		// make actual read for the next result
		do_read();
		return *this;
	}
};


#endif /* IFASTQSTREAM_HPP_ */
