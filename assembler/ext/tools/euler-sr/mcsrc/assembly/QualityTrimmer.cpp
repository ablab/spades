/***************************************************************************
 * Title:          QualityTrimmer.cpp 
 * Author:         Mark Chaisson, Glenn Tesler
 * Created:        2007
 * Last modified:  01/28/2010
 *
 * Copyright (c) 2007-2010 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "SeqReader.h"
#include "DNASequence.h"
#include "utils.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <cmath>


using namespace std;
void PrintUsage() {
  cout << "usages:" << std::endl;
	cout << "  qualityTrimmer -fasta inFastA -qual inQual -outFasta out" << std::endl;
	cout << "  qualityTrimmer -fastq inFastQ -outFasta out" << std::endl;
  cout << "    -minQual qual   Trim off leading and trailing positions below 'qual'" << std::endl;
  cout << "                    (default 12)." << std::endl;
  cout << "    -span span      Trim off any region that does not average at least 'qual'" << std::endl;
	cout << "                    over a span of 'span' (default 10)." << std::endl; 
	cout << "    -maxTrim trim   If at least 'trim' bases are removed, discard the read." << std::endl;
  cout << "                    Set to 0 to disable (default)." << std::endl;
	cout << "    -type           illumina|sanger   (default sanger)." << std::endl;
 	cout << "                    Use either Illumina or Sanger scaling of quality values." << std::endl;
	cout << "    -trimFront N    Automatically trim N bases from the front (default 0)." << std::endl;
	cout << "    -trimEnd   N    Automatically trim N bases from the end (default 0)." << std::endl;
}


ssize_t ReadQualityValues(std::istream &in, std::string &qualTitle, std::vector<ssize_t> &qualValues) {
  std::string line;
  if (!in ) {
    return 0;
  }
  if (in.peek() == '>') {
    in.get();
    in >> qualTitle;
    std::getline(in, line);
  }
  qualValues.clear();
  ssize_t qualValue;
  
  while (in and in.peek() != '>') {
    std::getline(in, line);
    std::stringstream lineStrm(line);
    // Parse the line and append values in it
    while (lineStrm >> qualValue) { 
      qualValues.push_back(qualValue);
      qualValue = -1;
    }
    if (in.eof()) 
      return 1;
    line = "";
  }
  return 1;
}

ssize_t TrimRead(DNASequence &origRead, std::vector<ssize_t> &qualValues,
						 ssize_t minQual, ssize_t span, ssize_t minFrontTrim, ssize_t minEndTrim,
						 DNASequence &newRead, ssize_t &trimFront, ssize_t &trimEnd) {
  // copy the title for completeness 
	//  newRead.namestr = origRead.namestr;
	//  int s;

	//  double totQual;
	//  DNASequence maskedSeq;
	//  maskedSeq = origRead;

	//  int qualSpan;

  // Find the first valid start position (this excludes parts
  // that just have 1 or 2 masked nucleotides)
	ssize_t start = minFrontTrim;
	ssize_t end = origRead.length - 1 - minEndTrim;

  for ( ; start <= end; start++) {
				if (qualValues[start] >=minQual) break;
  }

	// find last valid position 
  for ( ; end >= start; end--) {
      if (qualValues[end] >= minQual)
				break;
  }

	if (start > end) {
		// read is bad, end now
		return 0;
	}

	// Compute averages on windows of length = span;
	// in the tail, the window length is smaller.
	// Find a window with avg QV >= minQual
  ssize_t firstSolidPos, lastSolidPos;
  firstSolidPos = -1;
  lastSolidPos  = -1;
	ssize_t p,p1,nterms;
	ssize_t sum, avg;

	for (p = start; p <= end; p++) {
		nterms = origRead.length - p;
		if (nterms > span)
			nterms = span;

		sum = 0;

		for (p1 = 0; p1 < nterms; p1++) {
			sum += qualValues[p+p1];
		}
		avg = sum / nterms;

		if (avg >= minQual) {
			if (firstSolidPos == -1) {
				firstSolidPos = p;
			}
			lastSolidPos = p;
		} else {
			if (firstSolidPos != -1)
				break;
		}
	}

  if (firstSolidPos == -1	or lastSolidPos == -1	or firstSolidPos >= lastSolidPos) {
    // The read is just plain bad, exit now
		return 0;
  }

  // Now produce new sequence that has the bad ends of the read trimmed off
	//    int trimmedLength;
	//    trimmedLength = lastSolidPos - firstSolidPos;
	trimFront = firstSolidPos;
	trimEnd   = lastSolidPos;
	newRead._ascii = 1;
	newRead.seq    = origRead.seq + trimFront;
	newRead.length = trimEnd - trimFront + 1;
	return 1;
}

int main(int argc, char* argv[]) {
  
  std::string seqName, qualName, outName, qualType;
  ssize_t minQual, span;
  minQual = 12, span=10;

  if (argc < 4) {
    PrintUsage();
    exit(0);
  }
  int argi;
  argi = 1;
  ssize_t maxTrim = 0;
	ssize_t fileType = -1; // type not known
	ssize_t fastA = 1;
	ssize_t fastQ = 2;
	seqName  = "";
	qualName = "";
	outName  = "";
	qualType = "";
	ssize_t minTrimFront = 0;
	ssize_t minTrimEnd   = 0;
  while (argi < argc) {
		if (strcmp(argv[argi], "-fasta") == 0) {
			seqName = argv[++argi];
			fileType = fastA;
		}
		else if (strcmp(argv[argi], "-qual") == 0) {
			qualName = argv[++argi];
		}
		else if (strcmp(argv[argi], "-outFasta") == 0) {
			outName = argv[++argi];
		}
		else if (strcmp(argv[argi], "-fastq") == 0) {
			seqName = argv[++argi];
			fileType = fastQ;
		}
    else if (strcmp(argv[argi], "-minQual") == 0) {
      minQual = atoi(argv[++argi]);
    }
    else if (strcmp(argv[argi], "-span") == 0) {
      span = atoi(argv[++argi]);
    }
		else if (strcmp(argv[argi], "-maxTrim") == 0) {
			maxTrim = atoi(argv[++argi]);
		}
		else if (strcmp(argv[argi], "-type") == 0) {
			qualType = argv[++argi];
			if (qualType != "illumina" and
					qualType != "sanger") {
				PrintUsage();
				cout << "Type '"<< qualType << "' must be either 'illumina' or 'sanger'" << endl;
				exit(1);
			}
		}
		else if (strcmp(argv[argi], "-trimFront") == 0) {
			minTrimFront = atoi(argv[++argi]);
		}
		else if (strcmp(argv[argi], "-trimEnd") == 0) {
			minTrimEnd = -atoi(argv[++argi]);
		}
    else {
      PrintUsage();
			std::cout << "bad option: " << argv[argi] << std::endl;
      exit(1);
    }
    ++argi;
  }
	if (seqName == outName) {
		std::cout << "Error! input and output file names are the same." << std::endl;
		exit(0);
	}
	if (seqName == ""|| outName == "") {
		PrintUsage();
		std::cout << "Missing input or output name." << std::endl;
		exit(0);
	}
	if (qualName == "" and fileType == fastA) {
		std::cout << "You must specify a quality file if the input is FASTA" << std::endl;
		exit(0);
	}

  std::ofstream seqOut;
  openck(outName, seqOut, std::ios::out);

	
  std::ifstream qualIn;
	if (fileType == fastA)
		openck(qualName, qualIn, std::ios::in);

  std::ifstream seqIn;
  openck(seqName, seqIn, std::ios::in);
	std::vector<ssize_t> qualValues;
	//	int *qualValuesPtr;
	//	int qualValuesLength, qualValuesBufLength;
	//	qualValuesBufLength = 0;
	//	qualValuesLength = 0;
	//  std::vector<int> qualValues;
  std::string qualTitle;
  DNASequence seq, trimmedSeq;

  ssize_t totalReads = 0;       // # of reads in input
  ssize_t nDiscards = 0;        // # of reads discarded completely
	ssize_t nTrimmed = 0;         // # reads trimmed

	ssize_t totalBases = 0;       // # of bases in input
  ssize_t trimmedBases = 0;     // # bases trimmed from reads that are retained
  ssize_t discardedBases = 0;   // # bases in reads that are discarded

	seq._masked = 0;
	trimmedSeq._masked = 0;

	ssize_t illmnQualScale[256];
	ssize_t i;
	for (i = 0; i < 256; i++) {
		illmnQualScale[i] = (ssize_t) (10*log(1+ std::pow(10, ((double)i-64)/10.0))/log(10) + 0.499);
	}

  while (true) {
		if (fileType == fastA) {
			if (!SeqReader::GetSeq(seqIn, seq, SeqReader::noConvert))
				break;
			if (!ReadQualityValues(qualIn, qualTitle, qualValues)) {
				std::cout << "error reading quality values for sequence " << seq.namestr << std::endl;
				return 1;
			}
			/*
			if (qualValues.size() != seq.length) {
				std::cout << "error, there are " << qualValues.size() << " values for " << qualTitle
									<< " but the read " << seq.namestr << " is of length " << seq.length << std::endl;
				return 1;
			}
			*/
		}
		else if (fileType == fastQ) {
			std::string title, seqString, qualString, qualTitle;
			// Get rid of '@'
			if (seqIn.get() == EOF)
				break;
			// Read the title.
			if (!std::getline(seqIn, title)) break;

			// discard newline
			//			if (seqIn.get() == EOF) break;

			// Read the sequence.
			if (!std::getline(seqIn, seqString)) break;
			
			//			seq.Reset(seqString.size());
			if (seq.length < seqString.size()) {
				seq.Reset(seqString.size());
			}
			seq.length = seqString.size();

			memcpy(seq.seq, seqString.c_str(), seqString.size());
			//			seq.namestr = title.substr(1);    // not needed, we already deleted '@'
			seq.namestr = title;

			// Discard the quality title.
			if (!std::getline(seqIn, qualTitle)) break;
			// Read the quality string.
			if (!std::getline(seqIn, qualString)) break;
			
      if (qualValues.size() < qualString.size()) {
				qualValues.resize(qualString.size());
			}
			ssize_t qualValuesLength = qualString.size();
			ssize_t i;
			//			int qualLength = qualValues.size();
			char * qualPtr = (char*) qualString.c_str();

			// Guess the fastq type: illumina or sanger
			ssize_t sanger = 1;
			ssize_t illumina = 2;
			ssize_t type = sanger;
			
			// Guess the quality type if necessary
			if (qualType == "") {
				for (i = 0; i < qualValuesLength; i++ ){
					if (qualPtr[i] > 74) {
						type = illumina;
						break;
					}
				}
			}
			else if (qualType == "illumina") {
				type = illumina;
			}
			for (i = 0; i < qualValuesLength; i++ ){
				if (type == illumina) {
					qualValues[i] =  illmnQualScale[(unsigned char) qualPtr[i]]; //10*log(1+ std::pow(10, (qualPtr[i]-64)/10.0))/log(10) + 0.499;
				}
				else {
					qualValues[i] = qualPtr[i] - 33;
				}
			}

			//			std::cout << std::endl;
		}
		PrintStatus(totalReads, 100000);
    totalReads++;
		totalBases += seq.length;
		//    std::cout << "trimming " << seq.namestr << std::endl;

		ssize_t trimFront, trimEnd;
		trimFront = 0;
		trimEnd = 0;
		
    if (TrimRead(seq, qualValues, minQual, span, minTrimFront, minTrimEnd, 
								 trimmedSeq, trimFront, trimEnd)) {
			// std::stringstream trimKW;

			if (maxTrim == 0 or (seq.length - trimmedSeq.length) < maxTrim) {
				/*if (trimFront >= 0 or trimEnd >= 0) {
					trimKW << " qualFront=" << trimFront << " qualEnd=" << trimEnd;
					trimmedSeq.namestr += trimKW.str();
				}
				*/
				//seq.PrintSeq(seqOut);
				seqOut << ">" << seq.namestr << endl;
				trimmedSeq.PrintSeq(seqOut);
				seqOut << std::endl;
				trimmedBases += seq.length - trimmedSeq.length;
				if (seq.length != trimmedSeq.length) {
					nTrimmed++;
				}
			}
			else {
				//if (maxTrim > 0 and (seq.length-trimEnd) + trimFront > maxTrim) {
				++nDiscards;
				discardedBases += seq.length;
			}
    }
    else {
      ++nDiscards;
			discardedBases += seq.length;
    }
  }

	ssize_t nOutput = totalReads - nDiscards;
  std::cout << std::endl;
	std::cout << "Input " << totalReads << " reads, "
						<< "comprising " << totalBases << " nucleotides" << std::endl;
	std::cout << "Discarded " << nDiscards << " reads, "
						<< "comprising " << discardedBases << " nucleotides" << std::endl;
	std::cout << "Trimmed " << nTrimmed << " reads; "
						<< trimmedBases << " nucleotides trimmed" << std::endl;
	std::cout << "Output " << nOutput << " reads; "
						<< nTrimmed << " trimmed reads, "
						<< (nOutput - nTrimmed) << " whole reads" << std::endl;

	//  std::cout << "discarded a total of " << nDiscards << " / " << totalReads << " reads" << std::endl;
	//  std::cout << "trimmed a total of " << trimmedBases << " nucleotides" << std::endl;
  return 0;
}
