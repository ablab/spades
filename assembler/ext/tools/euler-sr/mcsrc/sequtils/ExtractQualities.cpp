/***************************************************************************
 * Title:          ExtractQualities.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  10/24/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <set>
#include <string>
#include "SeqReader.h"
#include "DNASequence.h"
#include "utils.h"
#include "ParseTitle.h"


int main(int argc, char *argv[]) {
  std::string seqFileName, seqNameFileName;
  std::string havesFileName, haveNotsFileName;

  if (argc < 4) {
    std::cout << "usage: extractQuals qualFile seqNamesFile havesFile, haveNotsFile" << std::endl;
    exit(1);
  }

  seqFileName = argv[1];
  seqNameFileName = argv[2];
  havesFileName  = argv[3];
  haveNotsFileName = argv[4];

  std::set<std::string> names;

  std::ifstream seqIn, namesIn;
  std::ofstream haveOut, haveNotOut;
  openck(seqFileName, seqIn, std::ios::in);
  openck(seqNameFileName, namesIn, std::ios::in);

  // Read all the names
  std::string line, name;
  while (namesIn) {
    namesIn >> line;
		if (line.c_str()[0] == '>') {
			ParseTitle(line, name);
		}
		names.insert(name);
  }

	openck(havesFileName, haveOut, std::ios::out);
	openck(haveNotsFileName, haveNotOut, std::ios::out);
	std::string  readQual, readTitle;
	readTitle = "";
	while (seqIn) {
		line = "";
		// Fetch the next line
		while (line.size() == 0 and seqIn.good() and std::getline(seqIn, line)) ;
		if (line.size() == 0) break;

		if (line.c_str()[0] == '>') {
			if (readTitle != "") {
				// check to see if the last qual should be printed
				if (names.find(readTitle) != names.end()) {
					haveOut << readQual;
				}
				else {
					haveNotOut << readQual;
				}
			}
			readQual = line + "\n";
			ParseTitle(line, readTitle);
		}
		else {
			readQual += line + "\n";
		}
	}
	if (readTitle != "") {
		if (names.find(readTitle) != names.end())
			haveOut << readQual;
		else
			haveNotOut << readQual;
	}
 
  return 1;
}
