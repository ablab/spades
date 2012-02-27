/***************************************************************************
 * Title:          SeqReader.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  10/31/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef SEQ_READER
#define SEQ_READER

#include <fstream>
#include <vector>
#include "DNASequence.h"

class SeqReader {
private:
  char *_buffer;
  std::ifstream *_in;
  ssize_t _buffFull;
  ssize_t _buffLen;
  static ssize_t GetSeqName(std::ifstream &in, std::string &name);
  static ssize_t _maskRepeats;
	static ssize_t lineNumber;
  ssize_t _lineNumber;
	static ssize_t _removeRepeats;
public: 
  static ssize_t noConvert;
 SeqReader(std::ifstream *in) :  _in(in), _buffFull(0), _buffLen(0), _lineNumber(0) {};
  
  // Public utility functions
	 static ssize_t GetSeq(std::string seqName, DNASequence &sequence, ssize_t c=0, std::ostream &report=std::cout);
  static ssize_t GetSeq(std::ifstream &in, DNASequence &sequence, ssize_t c=1);
  static ssize_t GetSeq(std::ifstream &in, DNASequence *&seqPtr, ssize_t c=1);
  // Functions for working on one stream per seqreader
  void SetSource(std::ifstream *source) { _in = source; _lineNumber = 0; _buffLen = 0; }
  //  int GetSeq(char* & seq, char* &name, long & len, int convert=1);
  ssize_t GetSeq(DNASequence &s, ssize_t convert=1);
  ssize_t GetSeq(DNASequence* &s, ssize_t convert=1);
	static ssize_t GetRead(std::istream &in, DNASequence &sequence);
	ssize_t GetRead(DNASequence &sequence);
	
  ssize_t GetLine(char* &line);
  static void MaskRepeats() { SeqReader::_maskRepeats = 1;}
	static void UnmaskRepeats() { SeqReader::_maskRepeats = 0;}
	void KeepRepeats() {_removeRepeats = 0;}
};

ssize_t ReadSequences(std::string &inFileName, DNASequenceList &sequences);
ssize_t ReadDNASequences(std::string &inFileName, DNASequenceList &sequences, ssize_t numToRead = -1, std::ostream &report = std::cout );
void ReadNamedSequences(std::string &inFileName, NamedSequenceList &sequences, std::ostream &report = std::cout);
ssize_t ReadSimpleSequences(std::string &inFileName, std::vector<SimpleSequence> &sequences, std::ostream &report=std::cout);
#endif
