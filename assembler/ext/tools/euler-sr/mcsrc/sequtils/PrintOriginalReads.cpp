/***************************************************************************
 * Title:          PrintOriginalReads.cpp 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  12/10/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "DNASequence.h"
#include "SeqReader.h"
#include "SeqUtils.h"
#include "utils.h"
#include "ParseTitle.h"
#include <string>

void PrintUsage() {
		std::cout << "usage: printOriginalReads refFile readFile origReadFile"<<std::endl;
		std::cout << " [-minLength len]   Print at least 'len' from every read." << std::endl;
}

int main(int argc, char* argv[]) {
	std::string refFileName, readFileName, origReadFileName;
	if (argc < 3) {
		PrintUsage();
		exit(1);
	}
	
	refFileName      = argv[1];
	readFileName     = argv[2];
	origReadFileName = argv[3];
	int argi = 4;
	ssize_t minLength = 0;
	while (argi < argc) {
		if (strcmp("-minLength", argv[argi]) == 0) {
			minLength = atoi(argv[++argi]);
		}
		else {
			PrintUsage();
			std::cout << "bad option: " << argv[argi] << std::endl;
			exit(0);
		}
		++argi;
	}

	DNASequence genome;
	SeqReader::GetSeq(refFileName, genome, SeqReader::noConvert);
	std::ifstream readsIn;
	std::ofstream readsOut;
	openck(readFileName, readsIn, std::ios::in);
	openck(origReadFileName, readsOut, std::ios::out);

	DNASequence read, origRead, origReadRC;
	ssize_t pos;
	ssize_t strand;
	origRead._ascii = 1;
	origReadRC._ascii = 1;
	while(SeqReader::GetSeq(readsIn, read, SeqReader::noConvert)) {
		if (ParseKeyword(read.namestr, "pos", pos) and
				ParseKeyword(read.namestr, "strand", strand)) {
			if (read.length > minLength) 
				origRead.length = read.length;
			else {
				if (pos + minLength > genome.length)
					origRead.length = genome.length - pos;
				else
					origRead.length = minLength;
			}
			origRead.seq    = &genome.seq[pos];
			origRead.namestr = read.namestr;
			if (strand == 0) {
				origRead.PrintlnSeq(readsOut);
			}
			else {
				MakeRC(origRead, origReadRC);
				origReadRC.namestr = origRead.namestr;
				origReadRC.PrintlnSeq(readsOut);
			}
		}
	}
}
