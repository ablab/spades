/***************************************************************************
 * Title:          PrintReadKmerCoverage.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  10/24/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
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


int main(int argc, char* argv[]) {

  std::string spectrumFileName, readFileName, coverageFileName;

  if (argc < 5) {
    std::cout << "usage: printReadKmerCoverage spectrumFile readsFile vertexSize " <<std::endl
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

  DNASequenceList seqList;
  ReadDNASequences(readFileName, seqList);

  ListSpectrum<StringMultTuple> spectrum;
  spectrum.tupleSize = tupleSize;
  spectrum.Read(spectrumFileName, 0);

  std::ofstream covOut, gapsOut;
	std::ofstream readsOut;

	if (printCoverage) 
		openck(coverageFileName, covOut, std::ios::out);
	
	if (printCovered) 
		openck(coveredReadsName, readsOut, std::ios::out);

	if (printGaps) {
		openck(gapsFileName, gapsOut, std::ios::out);
	}
  ssize_t s, p;
  StringMultTuple tuple;
  ssize_t index;
	ssize_t mult, prevMult;
	ssize_t start = 0;
	ssize_t end = -1;
	prevMult = -1;
  for (s =0 ; s < seqList.size(); s++ ) {
		if (printCoverage)
			covOut << seqList[s].namestr << std::endl;
		ssize_t allCovered = 1;
    for (p = 0; p < seqList[s].length - spectrum.tupleSize + 1; p++ ) {
			tuple.assign((char*)&seqList[s].seq[p]);
      if (tuple.Valid()) { //GetHashValue(seqList[s], p, spectrum.tupleSize, tuple);
				//				std::cout << p;
				index = spectrum.FindTuple(tuple);

				if (index >= 0) {
					mult = spectrum[index].GetMult();;
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
						ssize_t gapStart = end - spectrum.tupleSize - 20;
						if (printGaps and start > 0) {
							DNASequence gap;
							assert(start - end + 1 > 0);
							if (gapStart < 0) 
								gapStart = 0;
							gap.seq = &seqList[s].seq[gapStart];
							gap.length = gapEnd - gapStart + 1;
							std::stringstream namestrm;
							namestrm << s << "_" << gapStart << "_" << gapEnd << " pos=" << gapStart << " strand=0 len="<< gap.length;
							gap.namestr = namestrm.str();
							gap.PrintlnSeq(gapsOut);
						}							
					}
					else if ((mult < minMult and minMult <= prevMult) or
									 (p == seqList[s].length - spectrum.tupleSize) ) {
						end = p;
						DNASequence contig;
						contig.seq = &seqList[s].seq[start];
						contig.length = end - start + 1;
						std::stringstream namestrm;
						namestrm << s << "_" << start << "_" << end << " pos=" << start << " strand=0 len="<< end - start + 1;
						contig.namestr = namestrm.str();
						if (printCoverage) 
							contig.PrintlnSeq(covOut);
					}
					prevMult = mult;
				}
			}
			if (p > 0 and p % 50 == 0 and !printContigs)
				covOut << std::endl;
    }
		if (printCoverage)
			covOut << std::endl;
		if (allCovered and printCovered) {
			seqList[s].PrintlnSeq(readsOut);
		}
  }
}

  
