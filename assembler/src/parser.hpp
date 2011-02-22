#ifndef PARSER_HPP
#define PARSER_HPP

#include <algorithm>
#include <string>
#include <zlib.h>
#include <cstdlib>
#include <iostream>
#include "seq.hpp"
#include "libs/kseq/kseq.h"

// STEP 1: declare the type of file handler and the read() function 
KSEQ_INIT(gzFile, gzread)

template <int size> // size of reads in base pairs
class FASTQParser {
public:
	FASTQParser();
	void open(std::string filename1, std::string filename2); // open two files with mate paired reads
	MatePair<size> read();
	bool eof();
	void close();
private:
	gzFile fp1, fp2;
	kseq_t *seq1, *seq2;
	bool _eof_bit;
	bool _opened; // ready for read
	int _cnt;
	void do_read(); // do actual read, it's one sequence ahead of read()
};

// ******************** //
// * TEMPLATE METHODS * //
// ******************** //

template <int size>
MatePair<size> FASTQParser<size>::read() {
	assert(this->_opened);
	// invariant: actual read from file was already done by do_read();
	// assert that they are mate pairs (have: name/1 and name/2)!
	assert(strncmp(this->seq1->name.s, this->seq2->name.s, this->seq1->name.l - 1) == 0);
	// assume that all paired reads have the same length
	assert(this->seq1->seq.l == this->seq2->seq.l);
	// if there is 'N' in sequence, then throw out this mate read
	for (char *si1 = this->seq1->seq.s, *si2 = this->seq2->seq.s; *si1 != 0; ++si1, ++si2) {
		if (*si1 == 'N' || *si2 == 'N') {
			this->do_read();
			return MatePair<size>::null;
		}
	}
	// fill result
	MatePair<size> res(this->seq1->seq.s, this->seq2->seq.s, this->_cnt++);
	// make actual read for the next result
	this->do_read();
	return res;
}

template <int size>
FASTQParser<size>::FASTQParser() {
	this->_opened = false;
}

template <int size>
void FASTQParser<size>::open(std::string filename1, std::string filename2) {
	this->fp1 = gzopen(filename1.c_str(), "r"); // STEP 2: open the file handler
	this->fp2 = gzopen(filename2.c_str(), "r"); // STEP 2: open the file handler
	bool fail = false;
	if (!this->fp1) { std::cerr << "File " << filename1 << " not found!" << std::endl; fail = true; }
	if (!this->fp2) { std::cerr << "File " << filename2 << " not found!" << std::endl; fail = true; }
	if (fail) { this->_eof_bit = true; this->_opened = true; return; } // behave like it's empty file
	this->seq1 = kseq_init(fp1); // STEP 3: initialize seq
	this->seq2 = kseq_init(fp2); // STEP 3: initialize seq
	this->_eof_bit = false;
	this->_cnt = 0;
	this->do_read();
	this->_opened = true;
}

template <int size>
bool FASTQParser<size>::eof() {
	return this->_eof_bit;
}

template <int size>
void FASTQParser<size>::close() {
	if (!this->_opened) {
		return;
	}
	kseq_destroy(this->seq1); // STEP 5: destroy seq
	kseq_destroy(this->seq2); // STEP 5: destroy seq
	gzclose(this->fp1); // STEP 6: close the file handler
	gzclose(this->fp2); // STEP 6: close the file handler
	this->_opened = false;
}

template <int size>
void FASTQParser<size>::do_read() {
	if (kseq_read(this->seq1) < 0 || kseq_read(this->seq2) < 0) { // STEP 4: read sequence
		this->_eof_bit = true;
	}
}
#endif
