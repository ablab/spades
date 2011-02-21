#ifndef PARSER_HPP
#define PARSER_HPP

#include <string>
#include <zlib.h>
#include <iostream>
#include "seq.hpp"
#include "libs/kseq/kseq.h"

// STEP 1: declare the type of file handler and the read() function 
KSEQ_INIT(gzFile, gzread)

template <int size> // size of reads
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
	bool eof_bit;
	bool opened; // ready for read
	int cnt;
	void do_read(); // do actual read, it's one sequence ahead of read()
};

// ******************** //
// * TEMPLATE METHODS * //
// ******************** //

template <int size>
MatePair<size> FASTQParser<size>::read() {
	assert(this->opened);
	// invariant: actual read from file was already done by do_read();
	// assert that they are mate pairs (have: name/1 and name/2)!
	assert(strncmp(this->seq1->name.s, this->seq2->name.s, this->seq1->name.l - 1) == 0);
	// if there is 'N' in sequence, then throw out this mate read
	for (unsigned int i = 0; i < this->seq1->seq.l; ++i) {
		if (this->seq1->seq.s[i] == 'N') {
			this->do_read();
			return MatePair<size>("", "", -1);
		}
	}
	for (unsigned int i = 0; i < this->seq2->seq.l; ++i) {
		if (this->seq2->seq.s[i] == 'N') {
			this->do_read();
			return MatePair<size>("", "", -1);
		}
	}
	// fill result
	MatePair<size> res(this->seq1->seq.s, this->seq2->seq.s, this->cnt++);
	// make actual read for the next result
	this->do_read();
	return res;
}

template <int size>
FASTQParser<size>::FASTQParser() {
	this->opened = false;
}

template <int size>
void FASTQParser<size>::open(std::string filename1, std::string filename2) {
	this->fp1 = gzopen(filename1.c_str(), "r"); // STEP 2: open the file handler
	this->fp2 = gzopen(filename2.c_str(), "r"); // STEP 2: open the file handler
	bool fail = false;
	if (!this->fp1) { std::cerr << "File " << filename1 << " not found!" << std::endl; fail = true; }
	if (!this->fp2) { std::cerr << "File " << filename2 << " not found!" << std::endl; fail = true; }
	if (fail) { this->eof_bit = true; this->opened = true; return; } // behave like it's empty file
	this->seq1 = kseq_init(fp1); // STEP 3: initialize seq
	this->seq2 = kseq_init(fp2); // STEP 3: initialize seq
	this->eof_bit = false;
	this->cnt = 0;
	this->do_read();
	this->opened = true;
}

template <int size>
bool FASTQParser<size>::eof() {
	return this->eof_bit;
}

template <int size>
void FASTQParser<size>::close() {
	kseq_destroy(this->seq1); // STEP 5: destroy seq
	kseq_destroy(this->seq2); // STEP 5: destroy seq
	gzclose(this->fp1); // STEP 6: close the file handler
	gzclose(this->fp2); // STEP 6: close the file handler
	this->opened = false;
}

template <int size>
void FASTQParser<size>::do_read() {
	if (kseq_read(this->seq1) < 0 || kseq_read(this->seq2) < 0) { // STEP 4: read sequence
		eof_bit = true;
	}
}
#endif
