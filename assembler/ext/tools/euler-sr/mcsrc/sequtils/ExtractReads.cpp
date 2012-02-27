/***************************************************************************
 * Title:          ExtractReads.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  01/22/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <set>
#include <map>
#include <string>
#include "SeqReader.h"
#include "DNASequence.h"
#include "utils.h"
#include "ParseTitle.h"

int main(int argc, char *argv[]) {
  std::string seqFileName, seqNameFileName;
  std::string havesFileName, haveNotsFileName;

  if (argc < 4) {
    std::cout << "usage: extractReads seqFile seqNamesFile havesFile haveNotsFile " << std::endl
							<< " [-ordered] seqNames file is in the same order as seqFile " << std::endl
							<< " [-orderByNames] order the output file according to seqNamesFile, " << std::endl
							<< "                 which is not in the same order as seqFile." << std::endl;
    exit(1);
  }
	ssize_t printExcluded = 0;
  seqFileName = argv[1];
  seqNameFileName = argv[2];
  havesFileName  = argv[3];
	int argi = 4;
	if (argc >= 5 && argv[4][0] != '-') {
		printExcluded = 1;
		haveNotsFileName = argv[4];
		argi = 5;
	}
	ssize_t doOrderedSearch = 0;
	ssize_t doOrderByNames  = 0;
	while (argi < argc) {
		if (strcmp(argv[argi], "-ordered") == 0) {
			doOrderedSearch = 1;
		}
		else if (strcmp(argv[argi], "-orderByNames") == 0) {
			doOrderByNames = 1;
		}
		argi++;
	}

  std::set<std::string> names;
	std::vector<std::string> orderedNames;
  std::ifstream seqIn, namesIn;
  std::ofstream haveOut, haveNotOut;
	DNASequenceList sequences;
	std::map<std::string, ssize_t> index;
	std::string name, fastaName;
	if (doOrderByNames) {
		ReadDNASequences(seqFileName, sequences);
		ssize_t i;
		for (i = 0; i < sequences.size(); i++) {
			name = sequences[i].namestr;
			if (name.c_str()[0] == '>')
				name.replace(0,1, "");
			ParseTitle(name, fastaName);
			index[name] = i;
		}
	}

  openck(seqFileName, seqIn, std::ios::in);
  openck(seqNameFileName, namesIn, std::ios::in);
  openck(havesFileName, haveOut, std::ios::out);
	if (printExcluded)
		openck(haveNotsFileName, haveNotOut, std::ios::out);

  // Read all the names
  std::string line;
	ssize_t seqIndex;
  while (namesIn) {
    namesIn >> name;
    std::getline(namesIn, line);
    if (name.c_str()[0] == '>') {
      name.replace(0, 1, "");
			ParseTitle(name, fastaName);
			if (doOrderedSearch) {
				orderedNames.push_back(fastaName);
			}
			else if (doOrderByNames) {
				if (index.find(fastaName) != index.end()) {
					seqIndex = index[fastaName];
					sequences[seqIndex].PrintlnSeq(haveOut);
				}
			} else {
				names.insert(fastaName);
			}
		}
  }

	if (doOrderByNames) {
		// Already done with the printing. 
		haveOut.close();
		return 0;
	}
	
  DNASequence seq;
	ssize_t curName = 0;
	// In the case of an ordered search the names file
	// is a strict subset of the reads file.
	// Search 
  while (SeqReader::GetSeq(seqIn, seq)) {
    ParseTitle(seq.namestr, fastaName);
		if (doOrderedSearch) {
			//			std::cout << curName << " " << orderedNames[curName] << " " << fastaName << std::endl;
			if (fastaName == orderedNames[curName]) {
				seq.PrintlnSeq(haveOut);
				curName++;
			}
			else {
				if (printExcluded) {
					seq.PrintlnSeq(haveNotOut);
				}
			}
		}
		else {
			if (names.find(fastaName) != names.end()) {
				seq.PrintlnSeq(haveOut);
			}
			else {
				if (printExcluded) {
					seq.PrintlnSeq(haveNotOut);
				}
			}
		}
  }
  
  return 1;
}
