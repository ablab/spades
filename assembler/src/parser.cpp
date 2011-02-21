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

	// if there is 'N' in sequence, then throw out this mate read
	//std::cerr << this->seq1->seq.s << std::endl;
	for (int i = 0; i < this->seq1->seq.l; ++i) {
		if (this->seq1->seq.s[i] == 'N') {
			this->do_read();
			return MatePair("", "", -1);
		}
	}
	for (int i = 0; i < this->seq2->seq.l; ++i) {
		if (this->seq2->seq.s[i] == 'N') {
			this->do_read();
			return MatePair("", "", -1);
		}
	}
	
	// fill result
	MatePair res(this->seq1->seq.s, this->seq2->seq.s, this->cnt++);
	
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

MatePair::MatePair(const std::string &s1, const std::string &s2, const int id_) : id(id_), seq1(s1), seq2(s2) {
}

MatePair::MatePair(const MatePair &mp) : id(mp.id), seq1(mp.seq1), seq2(mp.seq2) {
}

Read::Read (const std::string &s) {
	char byte = 0;
	int cnt = 6;
	int cur = 0;
	for (std::string::const_iterator si = s.begin(); si != s.end(); si++) {
		switch (*si) {
			case 'C': byte |= (1 << cnt); break;
			case 'G': byte |= (2 << cnt); break;
			case 'T': byte |= (3 << cnt); break;
		}
		cnt -= 2;
		if (cnt < 0) {
			this->bytes[cur++] = byte;
			cnt = 6;
			byte = 0;
		}
	}
	if (cnt != 6) {
		this->bytes[cur++] = byte;
	}
}

char Read::operator[] (const int &index) const {
	switch ( ( this->bytes[index/4] >> ((3-(index%4))*2) ) & 3) { // little endian!
		case 0: return 'A'; break;
		case 1: return 'C'; break;
		case 2: return 'G'; break;
		case 3: return 'T'; break;
		default: return 'N';
	}
}

std::string Read::str() const {
	std::string res = "";
	for (int i = 0; i < 100; ++i) {
		res += this->operator[](i);
	}
	return res;
 }
