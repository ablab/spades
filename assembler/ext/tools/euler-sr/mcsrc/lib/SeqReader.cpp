/***************************************************************************
 * Title:          SeqReader.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  10/31/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <sstream>
#include <iostream>
#include "SeqReader.h"
#include "utils.h"

ssize_t  SeqReader::_maskRepeats   = 0;
ssize_t  SeqReader::_removeRepeats = 1;
ssize_t  SeqReader::noConvert = 0;
ssize_t  SeqReader::lineNumber = 0;
#define BUF_LEN 1024

#ifdef _INTEL_
#define CROPE std::crope
#else
#define CROPE __gnu_cxx::crope
#endif


ssize_t SeqReader::GetRead(DNASequence &sequence) {
	return GetRead(*_in, sequence);
}

ssize_t SeqReader::GetRead(std::istream &in, DNASequence &sequence) {
	std::string title;
	if (!in || !in.good())
		return 0;
	while (in and in.peek() != '>') in.get();

	std::getline(in, title);
	title = title.substr(1);
	std::string line;
	std::string read;
	char ws;
	while(in and in.peek() != '>' and !in.eof()) {
		std::getline(in, line);
		lineNumber++;
		read.append(line);
		while (in and ((ws = in.peek()) == ' ' or
									 ws == '\n' or
									 ws == '\t')) {
			in.get();
		}
		//		in.get(); // discard '\n'
	}
	
	//UNUSED// ssize_t pos;
	//UNUSED// ssize_t noGapPos;
	// Look for gaps in the reads.
	/*	noGapPos = 0;
	for (pos = 0; pos < read.size(); pos++) {
		if (read[pos] != ' ' and 
				read[pos] != '\t' and
				read[pos] != '\r' and
				read[pos] != '\n') {
			read[noGapPos] = read[pos];
			noGapPos++;
		}
		}*/
	sequence.Reset(read.size());
	ssize_t p;

	memcpy(sequence.seq, read.c_str(), read.size());
	sequence.length = read.size();
	for (p = 0; p < sequence.length; p++) {
		sequence.seq[p] = toupper(sequence.seq[p]);
	}
	sequence.namestr = title;
	while (in and in.peek() != '>' and !in.eof())
		in.get();
	return 1;
}

ssize_t SeqReader::GetSeq(DNASequence &sequence, ssize_t convert) {
  return GetSeq(*_in, sequence, convert);
}
ssize_t SeqReader::GetSeq(std::string seqName, DNASequence &sequence, ssize_t convert,
											std::ostream &report) {
  std::ifstream in;
  // in.open(seqName.c_str(), std::ios::binary | std::ios::in);
	//  if (!in.good()){

	if (!openck(seqName, in, std::ios::binary | std::ios::in, report, 0)) {
    std::cout << "could not open sequence file " << seqName << std::endl;
    return 0;
  }
  ssize_t retval = GetSeq(in, sequence, convert);
  in.close();
  in.clear();
  return retval;
}

ssize_t SeqReader::GetSeq(std::ifstream &in, DNASequence *&seqPtr, ssize_t convert) {
  seqPtr = NULL;
  if (!in or !in.good() or in.peek() == EOF) {
    return 0;
  }
  else {
    seqPtr = new DNASequence;
    return GetSeq(in, *seqPtr, convert);
  }
}
    

ssize_t SeqReader::GetSeq(std::ifstream &in, DNASequence &sequence, ssize_t convert){
  // early exit if this is the end of the file.
  if (!in or !in.good() or in.peek() == EOF) return 0;

  // Read the '>' part of the fasta file.
  if (!GetSeqName(in, sequence.namestr)) {
		return 0;
	}

  // Read in the fasta file. 
  //  std::string *line;
  char *line;
  //UNUSED// ssize_t tlc = 0;
  //UNUSED// ssize_t l,l1;
  //UNUSED// char c;
  std::streampos strmPos, endOfTitlePos;
  endOfTitlePos = strmPos = in.tellg();
  char p;
  while ( (p = in.peek()) != EOF and p != '>') in.get();
  if (p == EOF) {
    in.clear(); // want to continue reading although
    // fail bit is set on eof
  }

  ssize_t length;
  length = in.tellg() - endOfTitlePos;
  
  line = new char[length];
  in.seekg(endOfTitlePos, std::ios::beg);
  in.read(line, length);
    
  ssize_t curp, wsp;
  curp = wsp = 0;
  while (wsp < length) {
    if (line[wsp] == '\n' or 
				line[wsp] == ' ' or 
				line[wsp] == '\0') {
			if (line[wsp] == '\n')
				lineNumber++;
      wsp++;
		}
    else {
      line[curp] = line[wsp];
			// a hack for solexa sequences.
			if (line[curp] == '.')
				line[curp] = 'N';
      curp++;
      wsp++;
    }
  }
  sequence.Reset(curp);
  memcpy(&sequence.seq[0], line, curp);
  sequence.length = curp;
  delete[] line;

  ssize_t i;
  // It is possible to mask different things. 
  // use repeat masking.
  if (_maskRepeats)
    sequence.StoreAsMasked(REPEAT_MASKED);
	else 
		sequence.RemoveRepeats();
  
  // Store sequence as binary representation rather than character.
  if (convert) {
    for (i = 0; i < sequence.length; i++) {
      sequence.seq[i] = sequence.CharToNumeric(sequence.seq[i]);
    }
    sequence._ascii = 0;
  }
  else
    sequence._ascii = 1;
  return 1;
}

ssize_t SeqReader::GetSeq(DNASequence* &sequence, ssize_t convert ) {
  sequence = new DNASequence;
  ssize_t worked;
  worked =  GetSeq(*_in, *sequence, convert);
  if ( ! worked ) {
    delete sequence;
    sequence = NULL;
  }
  return worked;
}

ssize_t SeqReader::GetSeqName(std::ifstream &in, std::string &name) {
  //UNUSED// char *s;
  std::stringbuf strbuf;
	if (in.peek() != '>') {
		std::cout << "ERROR, titles in FASTA format must begin with a \">\"."<<std::endl;
		std::cout << "at line: " << lineNumber << std::endl;
		std::string line;
		std::getline(in, line);
		std::cout << "The offending line is: " << std::endl;
		std::cout << line << std::endl;
		exit(1);
	}
  if (in.peek() == '>') in.get();
	std::getline(in, name);
	lineNumber++;
  return in.peek() != EOF;
}

ssize_t ReadDNASequences(std::string &inFileName, DNASequenceList &sequences, ssize_t numSeq,
										 std::ostream &report) {
  DNASequence seq;
  std::ifstream seqIn;
  openck(inFileName, seqIn, std::ios::in, report);
	ssize_t count = 0;
  while(SeqReader::GetSeq(seqIn, seq, SeqReader::noConvert) and 
				(numSeq == -1 or count < numSeq) ) {
    sequences.push_back(seq);
		count++;
  }
  return sequences.size();
}

void ReadNamedSequences(std::string &inFileName, NamedSequenceList &sequences,
												std::ostream &report) {
  std::ifstream seqIn;
  openck(inFileName, seqIn, std::ios::in, report);
  DNASequence seq;
	NamedSequence namedSeq;
  while(SeqReader::GetSeq(seqIn, seq, SeqReader::noConvert)) {
		namedSeq.seq = new unsigned char[seq.length];
		memcpy(namedSeq.seq, seq.seq, seq.length);
		namedSeq.length= seq.length;
		namedSeq.namestr = seq.namestr;
		sequences.push_back(namedSeq);
	}
}

ssize_t ReadSimpleSequences(std::string &inFileName, std::vector<SimpleSequence> &sequences,
												std::ostream &report) {
  DNASequence seq;
  SimpleSequence simpleSeq;
  std::ifstream seqIn;
  openck(inFileName, seqIn, std::ios::in, report);
  while(SeqReader::GetSeq(seqIn, seq, SeqReader::noConvert)) {
    simpleSeq.seq = new unsigned char[seq.length];
    memcpy(simpleSeq.seq, seq.seq, seq.length);
		ssize_t s;
		for (s = 0; s < seq.length; s++) {
			simpleSeq.seq[s] = toupper(simpleSeq.seq[s]);
		}
    simpleSeq.length = seq.length;
    sequences.push_back(simpleSeq);
  }
  return sequences.size();
}
