#ifndef PARSER_HPP
#define PARSER_HPP

#include <algorithm>
#include <string>
#include <zlib.h>
#include <cstdlib>
#include <iostream>
#include "seq.hpp"
#include "libs/kseq/kseq.h"

using namespace std;

// STEP 1: declare the type of file handler and the read() function 
KSEQ_INIT(gzFile, gzread)

template <int size> // size of reads in base pairs
class FASTQParser {
public:

	FASTQParser() {
		_opened = false;
	}

	static vector<MatePair<size> >* readAll(const string &filename1, const string &filename2) {
		FASTQParser<size>* fqp = new FASTQParser<size>();
		vector<MatePair<size> >* res = new vector<MatePair<size> >;
		fqp->open(filename1, filename2);
		while (!fqp->eof()) {
			MatePair<size> mp = fqp->read(); // is it copy? :)
			if (mp.id != -1) { // don't have 'N' in reads
				res->push_back(mp);
			}
		}
		fqp->close();
		return res;
	}

public: // make it private!

	void open(string filename1, string filename2) { // open two files with mate paired reads
		_fp1 = gzopen(filename1.c_str(), "r"); // STEP 2: open the file handler
		_fp2 = gzopen(filename2.c_str(), "r"); // STEP 2: open the file handler
		bool fail = false;
		if (!_fp1) { cerr << "File " << filename1 << " not found!" << endl; fail = true; }
		if (!_fp2) { cerr << "File " << filename2 << " not found!" << endl; fail = true; }
		if (fail) {
			_opened = false;
			return;
		} // behave like it's empty file
		_seq1 = kseq_init(_fp1); // STEP 3: initialize seq
		_seq2 = kseq_init(_fp2); // STEP 3: initialize seq
		_opened = true;
		_eof = false;
		_cnt = 0;
		do_read();
	}

	MatePair<size> read() {
		assert(_opened);
		// invariant: actual read from file was already done by do_read();
		// assert that they are mate pairs (have: name/1 and name/2)!
		assert(strncmp(_seq1->name.s, _seq2->name.s, _seq1->name.l - 1) == 0);
		// assume that all paired reads have the same length
		assert(_seq1->seq.l == _seq2->seq.l);
		// if there is 'N' in sequence, then throw out this mate read
		for (char *si1 = _seq1->seq.s, *si2 = _seq2->seq.s; *si1 != 0; ++si1, ++si2) {
			if (*si1 == 'N' || *si2 == 'N') {
				this->do_read();
				return MatePair<size>::null;
			}
		}
		// fill result
		MatePair<size> res(_seq1->seq.s, _seq2->seq.s, _cnt++);
		// make actual read for the next result
		do_read();
		return res;
	}

	bool eof() {
		return !_opened || _eof;
	}

	void close() {
		if (!_opened) {
			return;
		}
		kseq_destroy(_seq1); // STEP 5: destroy seq
		kseq_destroy(_seq2); // STEP 5: destroy seq
		gzclose(_fp1); // STEP 6: close the file handler
		gzclose(_fp2); // STEP 6: close the file handler
		_opened = false;
	}

private:
	gzFile _fp1, _fp2;
	kseq_t *_seq1, *_seq2;
	bool _eof;
	bool _opened; // ready for read
	int _cnt;
	void do_read() { // do actual read, it's one sequence ahead of read()
		if (kseq_read(_seq1) < 0 || kseq_read(_seq2) < 0) { // STEP 4: read sequence
			_eof = true;
		}
	}
};

#endif
