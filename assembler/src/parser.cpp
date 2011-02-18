#include <zlib.h>
#include <stdio.h>
#include <cassert>
#include <iostream> // for debug
#include "parser.hpp"

FASTQParser::FASTQParser() {
	this->opened = false;
}

void FASTQParser::open(std::string filename1, std::string filename2) {
	this->fp1 = gzopen(filename1.c_str(), "r"); // STEP 2: open the file handler 
	this->fp2 = gzopen(filename2.c_str(), "r"); // STEP 2: open the file handler 
	this->seq1 = kseq_init(fp1); // STEP 3: initialize seq 
	this->seq2 = kseq_init(fp2); // STEP 3: initialize seq 
	this->eof_bit = false;
	this->cnt = 0;
	this->do_read();
	this->opened = true;
}

MatePair FASTQParser::read() {
	assert(this->opened);
	// invariant: actual read from file was already done by do_read();	
	
	// assert that they are mate pairs (have: name/1 and name/2)!
	assert(strncmp(this->seq1->name.s, this->seq2->name.s, this->seq1->name.l - 1) == 0);
	
	// fill result
	MatePair res;
	res.id = this->cnt++;
	res.seq1 = this->seq1->seq.s;
	res.seq2 = this->seq2->seq.s;
	
	// make actual read for the next result
	this->do_read();
	
	return res;
}

bool FASTQParser::eof() {
	return this->eof_bit;
}

void FASTQParser::close() {
	kseq_destroy(this->seq1); // STEP 5: destroy seq
	kseq_destroy(this->seq2); // STEP 5: destroy seq
	gzclose(this->fp1); // STEP 6: close the file handler
	gzclose(this->fp2); // STEP 6: close the file handler
	this->opened = false;
}

void FASTQParser::do_read() {
	if (kseq_read(this->seq1) < 0 || kseq_read(this->seq2) < 0) { // STEP 4: read sequence 
		eof_bit = true;
	}
}
