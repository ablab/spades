/***************************************************************************
 * Title:          PrintReadKmerCoverage.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  02/27/2010
 *
 * Copyright (c) 2007-2010 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "Spectrum.h"
#include "SimpleSequence.h"
#include "SeqReader.h"
#include "ListSpectrum.h"
#include "Tuple.h"
#include <iostream>
#include "utils.h"
#include "StringMultTuple.h"
#include "IntegralTupleStatic.h"
using namespace std;

int main(int argc, char* argv[]) {

  std::string spectrumFileName, readFileName, coverageFileName;

  if (argc < 5) {
    std::cout << "usage: iPrintReadKmerCov spectrumFile readsFile vertexSize " <<std::endl
							<< "  [-printCoverage file] print covereage of each nucleotide to 'file'" << std::endl
							<< "  [-printContigs] print covered contigs to coverageOut " << std::endl
							<< "  [-minMult m]" << std::endl
							<< "  [-printGaps file] print gaps to 'file'" << std::endl
							<< "  [-coveredReads file] print covered reads to 'file'" << std::endl;
		
    return 1;
  }

  int tupleSize;

  spectrumFileName = argv[1];
  readFileName     = argv[2];
  tupleSize        = atoi(argv[3]);

	
	std::string gapsFileName;
	std::string coveredReadsName;
	ssize_t printGaps = 0;
	int argi = 4;
	ssize_t printContigs = 0;
	ssize_t minMult = 0;
	ssize_t printCoverage = 0;
	ssize_t printCovered = 0;
	while (argi < argc) {
		if (strcmp("-printCoverage", argv[argi]) == 0) {
			printCoverage = 1;
			coverageFileName = argv[++argi];
		}
		else if (strcmp("-coveredReads", argv[argi]) == 0) {
			printCovered = 1;
			coveredReadsName = argv[++argi];
		}
		else if (strcmp("-printContigs", argv[argi]) == 0) {
			printContigs = 1;
		}
		else if (strcmp("-minMult", argv[argi]) == 0) {
			minMult = atoi(argv[++argi]);
		}
		else if (strcmp("-printGaps", argv[argi]) == 0) {
			printGaps = 1;
			gapsFileName = argv[++argi];
		}
		else {
			std::cout << "bad option: " << argv[argi] << std::endl;
			exit(1);
		}
		++argi;
	}
	//	IntegralTuple::tupleSize = tupleSize;
	IntegralTuple::SetTupleSize(tupleSize);
	vector<CountedIntegralTuple> spectrum;
	ssize_t spectrumSize;
	ReadMultBoundedBinaryTupleList(spectrumFileName, minMult, spectrum);
	spectrumSize = spectrum.size();
  std::ofstream covOut, gapsOut;
	std::ofstream readsOut;

	if (printCoverage) 
		openck(coverageFileName, covOut, std::ios::out);
	
	if (printCovered) 
		openck(coveredReadsName, readsOut, std::ios::out);

	if (printGaps) {
		openck(gapsFileName, gapsOut, std::ios::out);
	}
	//UNUSED// int p;
  ssize_t s=-1; 	// TODO: verify fix: s wasn't initialized but was used.  Added this & s++ in loop.
	CountedIntegralTuple tuple;
  ssize_t index;
	ssize_t mult, prevMult;
	ssize_t start=0, end=-1;
	prevMult = -1;

	DNASequence read;
	std::ifstream readFile;
	openck(readFileName, readFile, std::ios::in);
	while (SeqReader::GetSeq(readFile, read)) {
		s++;
		if (printCoverage)
			covOut << read.namestr << std::endl;
		ssize_t allCovered = 1;
		ssize_t p;
		ssize_t nextInvalidNuc = -1;
		for (p = 0; p < read.length - IntegralTuple::tupleSize + 1 and p < IntegralTuple::tupleSize; p++) {
			if (nucToIndex[read.seq[p]] >= 4)
				nextInvalidNuc = p;
		}
					 
		ssize_t valid = 0;
    for (p = 0; p < read.length - IntegralTuple::tupleSize + 1; p++ ) {
			if (nucToIndex[read.seq[p+IntegralTuple::tupleSize - 1]] >= 4)
				nextInvalidNuc = p;

			if (p <= nextInvalidNuc)
				valid = 0;
			else 
				valid = 1;

			tuple.StringToTuple((unsigned char*)&read.seq[p]);
			//				std::cout << p;
			index = -1;
			if (valid) 
				index = LookupBinaryTuple(spectrum, tuple);
				
			if (index >= 0) {
				mult = spectrum[index].GetMult();;
				/*
					if (mult > 1000) {
					string tuplestr;
					tuplestr.assign((char*) read.seq[p], IntegralTuple::tupleSize);
					cout << "mult of " << tuplestr << " " << tuple.tuple << " " << tuple.count << endl;
				}
				*/
			}
			else {
				mult = 0;
				allCovered = 0;
			}
			if (printCoverage) {
				covOut << mult << " ";
			}
			else {
				if (mult >= minMult and prevMult < minMult) {
					start = p;
					ssize_t gapEnd = start;
					ssize_t gapStart = end - IntegralTuple::tupleSize - 20;
					if (printGaps and start > 0) {
						DNASequence gap;
						assert(start - end + 1 > 0);
						if (gapStart < 0) 
							gapStart = 0;
						gap.seq = &read.seq[gapStart];
						gap.length = gapEnd - gapStart + 1;
						std::stringstream namestrm;
						namestrm << s << "_" << gapStart << "_" << gapEnd << " pos=" << gapStart << " strand=0 len="<< gap.length;
						gap.namestr = namestrm.str();
						gap.PrintlnSeq(gapsOut);
					}							
				}
				else if ((mult < minMult and minMult <= prevMult) or
								 (p == read.length - IntegralTuple::tupleSize) ) {
					end = p;
					DNASequence contig;
					contig.seq = &read.seq[start];
					contig.length = end - start + 1;
					std::stringstream namestrm;
					namestrm << s << "_" << start << "_" << end << " pos=" << start << " strand=0 len="<< end - start + 1;
					contig.namestr = namestrm.str();
					if (printCoverage) 
						contig.PrintlnSeq(covOut);
				}
				prevMult = mult;
			}

			if (p > 0 and p % 50 == 0 and !printContigs)
				covOut << std::endl;
    }
		if (printCoverage)
			covOut << std::endl;
		if (allCovered and printCovered) {
			read.PrintlnSeq(readsOut);
		}
  }
}

  
