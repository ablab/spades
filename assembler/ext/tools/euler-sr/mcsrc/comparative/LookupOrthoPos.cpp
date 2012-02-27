/***************************************************************************
 * Title:          LookupOrthoPos.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <iostream>
#include <ostream>
#include <stdio.h>
#include <unistd.h>
#include <cmath>
 
#include <algorithm>

#include "SeqReader.h"
#include "alignutils.h"
#include "lav/LAVFile.h"
#include "lav/LAVReader.h"


void InitEnv(int argc, char* argv[], 
	     std::string &alignFile, 
	     std::string &posFile,
	     std::vector<ssize_t> &pos, ssize_t &useRef, ssize_t &mapOnly);


void PrintUsage();


ssize_t GetQryPos(ssize_t pos, ssize_t length, ssize_t strand) {
  if (strand == 0) return pos;
  else 
    return length - pos;
}
void GetNextPos(std::ifstream &in, ssize_t &pos, std::string &preFormatting);

int main(int argc, char* argv[]) {

  std::string alignFile, posFileName;
  ssize_t useRef;
  ssize_t mapOnly;
  useRef = 0;
  mapOnly = 0;
  alignFile = "";
  std::vector<ssize_t> positions;
  InitEnv(argc, argv, alignFile, posFileName, positions, useRef, mapOnly);
  // Get input.
  std::ifstream in;
  if (posFileName != "") {
    in.open(posFileName.c_str());
    if (!in.good() ) {
      std::cout << "could not open " << posFileName << std::endl;
      exit(0);
    }
  }
  else {
    std::cout << "must enter a pos file " << std::endl;
    exit(0);
  }
  LAVFile lavAlignment;
  LAVReader::ReadLAVFile(alignFile, lavAlignment);
  // Get the query sequence length from the align file.
  ssize_t qryLength, qryEnd, qryBegin;

  if (lavAlignment.alignments.size() <= 0) {
    std::cout << "Need at least one alignment " << std::endl;
    exit(0);
  }
 


  // Look for them in the alignment.

  ssize_t i;
  ssize_t a, b,p;
  LAVAlignedContig *alignedContig;
  LAVBlock *lavBlock;
  
  std::vector<ssize_t> orthoPositions;
  ssize_t dir;
  std::vector<ssize_t> notFound;
  ssize_t searchedAll;
  ssize_t orthoPos;
  i = 0;
  ssize_t firstpass;
  std::string format;
  ssize_t pos;
  GetNextPos(in, pos, format);
  while (!in.eof() ) {
    searchedAll = 1;
    orthoPos = -1;
    firstpass = 1;
    for (a = 0; a < lavAlignment.alignments.size() && !in.eof(); a++) {
      alignedContig = lavAlignment.alignments[a];
      dir = alignedContig->qryContig.strand;
      qryLength = alignedContig->qryContig.end - alignedContig->qryContig.start + 1;
      qryEnd   = alignedContig->qryContig.end; 
      qryBegin = alignedContig->qryContig.start;
      ssize_t qryLen, refLen; // different from qryLength
      for (b = 0; b < alignedContig->alignments.size() && !in.eof(); b++) {
	lavBlock = alignedContig->alignments[b];
	for (p = 0; p < lavBlock->size() && !in.eof(); p++) {
	  if (useRef) {
	    // position is in an aligned block
	    if (lavBlock->refALBegin[p] <= pos &&
		lavBlock->refALEnd[p] >= pos) {
	      // store the position
	      orthoPositions.push_back(pos);
	      // store the orthologous position.
	      refLen = lavBlock->refALEnd[p] - lavBlock->refALBegin[p];
	      qryLen = lavBlock->qryALEnd[p] - lavBlock->qryALBegin[p];
	      
	      orthoPos = (ssize_t) std::floor((double(pos - lavBlock->refALBegin[p]) / refLen) * (qryLen))
					+ lavBlock->qryALBegin[p];
	      orthoPositions.push_back(orthoPos);
	      std::cout << format << orthoPos;
	      GetNextPos(in, pos, format);
	      searchedAll = 0;
	    }
	    else if (p < lavBlock->size() - 1 &&
		     lavBlock->refALEnd[p] <= pos &&
		     lavBlock->refALBegin[p+1] >= pos) {

	      // position is inbetween aligned blocks
	      orthoPositions.push_back(pos);
	      // store the orthologous position.
	      refLen = lavBlock->refALBegin[p+1] - lavBlock->refALEnd[p];
	      qryLen = lavBlock->qryALEnd[p+1] - lavBlock->qryALEnd[p];
	      
	      orthoPos = (ssize_t) std::floor((double(pos - lavBlock->refALEnd[p]) / refLen) * (qryLen)) 
		+ lavBlock->qryALEnd[p];
	      orthoPositions.push_back(orthoPos);
	      std::cout << format << orthoPos;
	      GetNextPos(in, pos, format);
	      searchedAll = 0;
	    }
	  }
	  else { 
	    // use query
	    if (dir == 0) {
	      if (lavBlock->qryALBegin[p] <= pos &&
		  lavBlock->qryALEnd[p] >= pos) {
		// store the position
		orthoPositions.push_back(pos);
		// store the orthologous position.
		refLen = lavBlock->refALEnd[p] - lavBlock->refALBegin[p];
		qryLen = lavBlock->qryALEnd[p] - lavBlock->qryALBegin[p];
		orthoPos = (ssize_t) std::floor((double(pos - lavBlock->qryALBegin[p]) / qryLen) * (refLen)) + lavBlock->refALBegin[p];
		orthoPositions.push_back(orthoPos);
		std::cout << format << orthoPos;
		GetNextPos(in, pos, format);
		searchedAll = 0;
	      }
	      else if (p < lavBlock->size() - 1 &&
		       lavBlock->qryALEnd[p] <= pos &&
		       lavBlock->qryALBegin[p+1] >= pos) {
		// position is inbetween aligned blocks
		orthoPositions.push_back(pos);
	      // store the orthologous position.
		refLen = lavBlock->refALBegin[p+1] - lavBlock->refALEnd[p];
		qryLen = lavBlock->qryALEnd[p+1] - lavBlock->qryALEnd[p];
		
		orthoPos = (ssize_t) std::floor((double(pos - lavBlock->qryALEnd[p]) / qryLen) * (refLen)) + lavBlock->refALEnd[p];
		orthoPositions.push_back(orthoPos);
		std::cout << format << orthoPos;
		GetNextPos(in, pos, format);
		searchedAll = 0;
	      }
	    }
	    else {
	      /*	      std::cout << GetQryPos(lavBlock->qryALEnd[p], qryLength, dir) << " " << pos << " "
			<< GetQryPos(lavBlock->qryALBegin[p], qryLength, dir) << " " 
			<< GetQryPos(lavBlock->qryALBegin[p+1], qryLength, dir) << std::endl;
	      */
	      if (GetQryPos(lavBlock->qryALEnd[p], qryLength, dir) <= pos &&
		  GetQryPos(lavBlock->qryALBegin[p], qryLength, dir) >= pos) {
		// store the position
		orthoPositions.push_back(pos);
		// store the orthologous position.
		refLen = lavBlock->refALEnd[p] - lavBlock->refALBegin[p];
		qryLen = lavBlock->qryALEnd[p] - lavBlock->qryALBegin[p];
		orthoPos = (ssize_t) std::floor(((pos - double(GetQryPos(lavBlock->qryALEnd[p], qryLength, dir))) / qryLen) * (refLen)) + lavBlock->refALBegin[p];
		orthoPositions.push_back(orthoPos);
		std::cout << format << orthoPos;
		GetNextPos(in, pos, format);
		searchedAll = 0;
	      }
	      else if (p < lavBlock->size() - 1 &&
		       GetQryPos(lavBlock->qryALBegin[p], qryLength, dir) <= pos &&
		       GetQryPos(lavBlock->qryALBegin[p+1], qryLength, dir) >= pos) {
		// position is inbetween aligned blocks
		orthoPositions.push_back(pos);
	      // store the orthologous position.
		refLen = lavBlock->refALBegin[p+1] - lavBlock->refALEnd[p];
		qryLen = lavBlock->qryALEnd[p+1] - lavBlock->qryALEnd[p];
		
		orthoPos = (ssize_t) std::floor(((pos - double(GetQryPos(lavBlock->qryALEnd[p], qryLength, dir))) / qryLen) * (refLen)) + lavBlock->refALEnd[p];
		orthoPositions.push_back(orthoPos);
		std::cout << format << orthoPos;
		GetNextPos(in, pos, format);
		searchedAll = 0;
	      }
	    } // end reverse complement search
	  } // end looking in qry alignment
	} // endlookign through blocks
      } // end looking through contigs
    } // end looking through alignments (should be 1 or 2).
    
    if (searchedAll) {
      notFound.push_back(pos);
      std::cout << format << "NE" ;
      GetNextPos(in, pos, format);
    }
  }
  return 0;
  /*
  i = 0;
  // print the result
  while (i < (ssize_t)orthoPositions.size() - 1 ) {
    if (mapOnly)
      std::cout << orthoPositions[i+1] << " ";
    else
      std::cout << orthoPositions[i] << "\t" << orthoPositions[i+1] << std::endl;
    i+= 2;
  }

  if (mapOnly)
    std::cout << std::endl;
  */
}



void InitEnv(int argc, char* argv[], 
	     std::string &alignFile,
	     std::string &posFileName,
	     std::vector<ssize_t> &pos, 
	     ssize_t &useRef, ssize_t &mapOnly) {

  ssize_t copt;
  ssize_t i;
  while ( (copt=getopt(argc, argv, "rmp:")) != EOF){
    switch(copt) {
    case 'r':
      useRef = 1;
      continue;
    case 'p':
      posFileName = optarg;
      continue;
    case 'm':
      mapOnly = 1;
      continue;
    default:
      PrintUsage();
      exit(1);
    }
  }
  i = optind;
  if (i >= argc) {
    std::cout << "You must specify a reference seq and query seq. " << std::endl;
    PrintUsage();
    exit(1);
  }
  alignFile = argv[i];
  i++;
  while (i < argc) {
    pos.push_back(atoi(argv[i]));
    i++;
  }
}


void PrintUsage() {
  std::cout << "loop\tLook up orthologous positions not very efficiently.\n\tGet in the loop. " 
	    << std::endl; 
  std::cout << "usage: luopos [rmp] alignfile pos1 pos2 ... " << std::endl;
  std::cout << " -r    use the reference position " << std::endl;
  std::cout << " -m    (deprecated) " << std::endl;
  std::cout << " -p posfile   read positions from posfile" << std::endl;
}


 void GetNextPos(std::ifstream &in, ssize_t &pos, std::string &preFormatting) {
   preFormatting = "";
   char c;
   while (!in.eof()) {
     c = in.peek();
     if ((ssize_t)c < (ssize_t) '0' || 
	 (ssize_t)c > (ssize_t) '9') {
       preFormatting += c;
       in.get();
     }
     else {
       in >> pos;
       return;
     }
   }
 }
