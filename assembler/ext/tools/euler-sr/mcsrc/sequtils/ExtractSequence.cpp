/***************************************************************************
 * Title:          ExtractSequence.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  02/04/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "DNASequence.cpp"
#include "SeqReader.h"
#include "utils.h"
#include <string>
#include <iostream>
#include <fstream>
#include <vector>

int main(int argc, char* argv[]) {
  if (argc < 4) {
    std::cout << "usage: extractseq sequencein start end sequenceout" 
								<< std::endl;
		std::cout << "or: extractseq sequence coordsFile out" << std::endl;
    exit(1);
  }
	//UNUSED// ssize_t readName = 0;
  std::string infile = argv[1];
  std::vector<ssize_t> start, end;
  std::string coordsFileName;
  std::string outFileName = "";
  std::ofstream outFile;
	std::vector<std::string> titles;
  if (argc == 5) {
    start.push_back(atoi(argv[2]));
    end.push_back(atoi(argv[3]));
    outFileName = argv[4];
    openck(outFileName, outFile, std::ios::out);
  }
  else {
    std::ifstream coordsFile;
    coordsFileName = argv[2];
		outFileName    = argv[3];
    openck(coordsFileName, coordsFile, std::ios::in);
    openck(outFileName, outFile, std::ios::out);
    ssize_t s,e;
		std::string line, title;
    while(!coordsFile.eof() and coordsFile) {
      coordsFile >> s >> e;
			std::getline(coordsFile, line); // discard the newline.
			if (line != "") {
				std::stringstream linestr(line);
				linestr >> title;
			}
      if (coordsFile) {
				/*				if (title != "") {
					std::cout << " title: "<< title<<std::endl;
				}
				else {
					std::cout<< std::endl;
				}
				*/
				start.push_back(s);
				end.push_back(e);
				//				titles.push_back(title);
      }
    }
  }
  DNASequence seq, out;

  std::ifstream inFile;
  openck(infile, inFile);

  SeqReader::MaskRepeats();

  SeqReader::GetSeq(inFile, seq);

  ssize_t i;
  ssize_t s, e;

	if (outFileName == "" and start.size() > 100){ 
		std::cout << "WARNING! You are about to generate: " << start.size() << "files!" << std::endl;
		std::string value;
		std::cout << "Type something to continue." << std::endl;
		std::cin >> value;
	}
			
  for (i = 0; i < start.size(); i++ ) {
    std::stringstream outFileStream;
    s = start[i];
    e = end[i];
		std::stringstream namestrm;
    if (outFileName == "") {
      std::string coordsToOutName;
      outFileStream << coordsFileName << "." << s << "." << e << ".fasta";
      coordsToOutName = outFileStream.str();
      openck(coordsToOutName, outFile, std::ios::out);
    }
    if (s >= 0 and e < seq.length) {
			if (titles.size() == 0) {
				namestrm.str("");
				namestrm << seq.namestr << "_" << start[i] << "_" << end[i] 
								 << " pos=" << start[i] << "  strand=0";
				out.namestr = namestrm.str();
			}
			else {
				out.namestr = titles[i];
			}
//				std::cout << "printing from " << s << " to " << e << std::endl;
      out._ascii = seq._ascii;
      out.seq    = &seq.seq[s];
      out.length = e - s + 1;
      out.PrintSeq(outFile);
			outFile << std::endl;
			if (outFileName == "" ) {
				outFile.close();
			}

    }
		/*
			Case when the bounds of the sequence aren't quite right. Just skip
			those for now.

    else {
      std::cout << "error: could not open " << coordsToOutName << std::endl;
      return 0;
    }
		*/
  }
	if (outFileName != "") {
		outFile.close();
	}
  return 0;

}
