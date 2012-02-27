/***************************************************************************
 * Title:          CombineTupleLists.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  10/24/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <iostream> 
#include <fstream>
#include <string>
#include <assert.h>
#include "utils.h"
#include <unistd.h>
#include <string.h>

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


ssize_t TupleBuffer::FetchTupleMult(std::string &tuple, ssize_t &mult) {
	if (curTuple >= curBufferSize) {
		if (!ReloadBuffer()) {
			return 0;
		}
		else {
			curTuple = 0;
		}
	}
	tuple = tupleBuf[curTuple];
	mult  = multBuf[curTuple];
	curTuple++;
	return 1;
}

ssize_t TupleBuffer::ReloadBuffer() {
	//	std::cout << "reloading buffer: " << index << std::endl;
	if ( !(*file).good()) return 0;
	ssize_t i;
	if (tupleBuf.size() != bufferSize) {
		tupleBuf.resize(bufferSize);
		multBuf.resize(bufferSize);
	}
	std::string tuple;
	ssize_t mult;
	i = 0;
	while (i < bufferSize and (*file).good()) {
		if (! ( (*file) >> tuple >> mult)) 
			break;
		tupleBuf[i] = tuple;
		multBuf[i]  = mult;
		i++;
	}
	curBufferSize = i;
	if (curBufferSize > 0)
		return 1;
	else
		return 0;
}
		



int main(int argc, char* argv[]) {

	std::string fileA, fileB, outName;
	if (argc < 4) {
		std::cout << "usage: combinetupleLists combined list1 list2 ... [-minMult m] [-listOut file]" << std::endl;
		std::cout << " -minMult m    Only store tuples with frequency greater or equal to m." << std::endl;
		std::cout << " -bufferOutput Periodically write output to a file." << std::endl;
		exit(1);
		// TODO: "-listOut" is in usage but isn't implemented
	}
	int argi = 1;
	outName = argv[argi++];
	std::vector<std::string> files;
	ssize_t minMult = 0;
	ssize_t useBuffering = 0;
	//UNUSED// ssize_t i;
	while (argi < argc) {
		if (strcmp(argv[argi],"-minMult") == 0) {
			minMult = atoi(argv[++argi]);
		}
		else if (strcmp(argv[argi], "-bufferOutput") == 0) {
			useBuffering = 1;
		}
		else {
			files.push_back(argv[argi]);
		}
		++argi;
	}
	if (files.size() < 2) {
		std::cout << "must merge at least two files" << std::endl;
		exit(1);
	}
	std::vector<std::ifstream*> listFiles;
	std::vector<char> listStat;
	std::vector<TupleBuffer> tupleBuffers;
	std::vector<std::string> curTuples;
	std::vector<ssize_t> curMults;
	
	listFiles.resize(files.size());
	listStat.resize(files.size());
	tupleBuffers.resize(files.size());
	curTuples.resize(files.size());
	curMults.resize(files.size());
	ssize_t f;
	ssize_t size;
	std::string outBufName;
	std::ofstream out, outBuf;
	outBufName = outName + ".buffer";

	TupleBuffer outputBuffer;
	openck(outName, out, std::ios::out);

	if (useBuffering) {
		std::cout << "opening buffer " << outBufName << std::endl;
		openck(outBufName, outBuf, std::ios::out);
		outputBuffer.OutputInit(&outBuf, 1000);
	}
		
	for (f = 0; f < files.size(); f++) {
		listFiles[f] = new std::ifstream;
		openck(files[f], *listFiles[f], std::ios::in);
		*(listFiles[f]) >> size;
		tupleBuffers[f].Init(listFiles[f], 1000);
		tupleBuffers[f].index = f;
		listStat[f] = 1;
		curTuples[f] = "";
		curMults[f]  = 0;
	}

	std::vector<std::string> mergedTuples;
	std::vector<ssize_t> mergedMults;
	std::string curTuple;
	ssize_t unemptyFiles = files.size();
	// Initialize the current tuples.
	for (f = 0; f < files.size(); f++ ) {
		if (!( tupleBuffers[f].FetchTupleMult(curTuples[f],curMults[f]))) {
			listStat[f] = 0;
		}
	}
	
	ssize_t curTupleMult;
	ssize_t numIts = 0;
	ssize_t numSkipped = 0;
	while (unemptyFiles > 0) {
		// Find which tuple to start on;
	  if (numIts % 1000 == 999){
	    std::cout << "iter: " <<numIts+1 <<" " << numSkipped << " / " << mergedTuples.size() << std::endl;
	  }
	  numIts ++;
		curTuple = "";
		for (f = 0; f < files.size(); f++) {
			if (listStat[f] and (curTuple == "" or curTuples[f] < curTuple)) {
				curTuple = curTuples[f];
			}
		}
		
		// Find the multiplicity of the current tuple
		// Also, any file that contains a tuple that is equivalent 
		// to the current tuple should be advanced to its next tuple 
		curTupleMult = 0;
		for (f = 0; f < files.size(); f++ ) {
			if (listStat[f] and curTuples[f] == curTuple) {
				curTupleMult += curMults[f];
				if (! (tupleBuffers[f].FetchTupleMult(curTuples[f], curMults[f]))) 
					listStat[f] = 0;
			}
		}
		
		if (curTupleMult >= minMult) {
			if (!useBuffering) {
				mergedTuples.push_back(curTuple);
				mergedMults.push_back(curTupleMult);
			}
			else {
				outputBuffer.Append(curTuple, curTupleMult);
			}
		}
		else {
		  numSkipped++;
		}


		// Now do a check to see if we are done
		unemptyFiles = 0;
		for (f =0 ; f < files.size(); f++ ) {
			unemptyFiles += listStat[f];
		}
	}
	
	if (!useBuffering) {
		ssize_t i;
		out << mergedTuples.size() << std::endl;
		for (i = 0; i < mergedTuples.size(); i++ ) {
			out << mergedTuples[i] << " " << mergedMults[i] << std::endl;
		}
	}
	else {
		outputBuffer.Flush();
		std::cout << "flushed: " << outputBuffer.totalFlushed << " tuples" << std::endl;
		// Now append the output 
		out << outputBuffer.totalFlushed << std::endl;
		std::ifstream bufIn;
		outBuf.close();
		openck(outBufName, bufIn, std::ios::in);
		out << bufIn.rdbuf();
		std::string removeOutCmd;
		removeOutCmd = "rm " + outBufName;
		system(removeOutCmd.c_str());
	}
	
	return 0;

}
