/***************************************************************************
 * Title:          IntegralCombineTupleLists.cpp 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  02/28/2010
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "SimpleSequence.h"
#include "SeqReader.h"
#include "DNASequence.h"
#include "NumericHashedSpectrum.h"
#include "utils.h"
#include "IntegralTupleStatic.h"



class TupleBuffer {
	ssize_t curTuple;
	ssize_t bufferSize;
	ssize_t curBufferSize;
public:
	ssize_t totalFlushed;
	ssize_t index;
	std::ifstream *file;
	std::ofstream *out;
	std::vector<std::string> tupleBuf;
	std::vector<ssize_t> multBuf;
	ssize_t bufIndex;
	ssize_t readyBuf;
	ssize_t FetchTupleMult(std::string &tuple, ssize_t &mult);
	ssize_t ReloadBuffer();
	void Init(std::ifstream *tupleFile, ssize_t bufSize) {
		file = tupleFile;
		bufferSize = bufSize;
		curBufferSize = 0;
		curTuple = 0;
	}
	void OutputInit(std::ofstream *outputFile, ssize_t bufSize) {
		out = outputFile;
		bufferSize = bufSize;
		tupleBuf.resize(bufSize);
		multBuf.resize(bufSize);
		curBufferSize = 0;
		index  = 0;
		totalFlushed = 0;
	}

	TupleBuffer(std::ifstream *tupleFile, ssize_t bufSize) {
		Init(tupleFile, bufSize);
	}
	TupleBuffer() {
		curBufferSize = 0;
	}
	void Flush() {
		ssize_t i;
		totalFlushed += curBufferSize;
		for (i = 0 ; i < curBufferSize; i++ ){ 
			*out << tupleBuf[i] << " " << multBuf[i] << std::endl;
		}
		curBufferSize = 0;
		index = 0;
	}
	void Append(std::string tuple, ssize_t mult) {
		tupleBuf[index] = tuple;
		multBuf[index] = mult;
		index++;
		curBufferSize++; // I don't remember the difference between curbuffersize and index.
		if (index == bufferSize) {
			Flush();
		}
	}
};
