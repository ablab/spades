/***************************************************************************
 * Title:          TransCoords.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  01/08/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <iostream>
#include <ostream>
#include <vector>
#include <map>
#include <string>
#include <unistd.h>

#include "mysql/mysql.h"
#include "mysql++/mysql++.h"

#include "DNASequence.h"
#include "utils.h"

typedef std::map<std::string, ssize_t> OffsetMap;

void ReadCMLFile(std::string cmlFileName, 
		 std::vector<std::string> &contigs, 
		 std::vector<ssize_t> &cmlSize);

void InitEnv(int argc, char* argv[], 
	     std::string &contigMap,
	     std::string &coordinates);

void PrintUsage();

int main(int argc, char* argv[]) {
  OffsetMap offsets;
  std::string qrySeqName, refSeqName;
  DNASequence refSeq, qrySeq; 
  std::string coordsFileName, contigOffsetsName;
  qrySeqName = "";
  refSeqName = "";

  InitEnv(argc, argv, contigOffsetsName, coordsFileName );
  // Get input.

  std::vector<std::string> contigs;
  std::vector<ssize_t> cmlContigs;
  ReadCMLFile(contigOffsetsName, contigs, cmlContigs);
  mysqlpp::Connection sqlCon("hsapiens");

  mysqlpp::Query query = sqlCon.query();

  std::ifstream coordsFile;
  openck(coordsFileName, coordsFile);

  ssize_t coord;
  ssize_t sign, chr;
  std::string start, signSt, chromosome;
  mysqlpp::Row row;
  mysqlpp::Result res;
  ssize_t i;

  ssize_t contigPos, contigStart, offset;
  while (!coordsFile.eof()) {
    if (coordsFile.peek() == EOF) 
      break;
    coordsFile >> chromosome >> chr >> signSt >> sign >> start >> coord;

    if (coord > cmlContigs[cmlContigs.size()-1]) {
      std::cout << "coordinate " << coord 
		<< " exceeds size of chromosome " << cmlContigs[cmlContigs.size()-1] << std::endl;
    }
    // find the contig.
    for (i = 0; i < contigs.size() && coord > cmlContigs[i]; i++);

    assert(i < cmlContigs.size());
    //  build the query
    query << "select chromStart from ctgpos where contig=\"" << contigs[i] << "\"";
    // executee the query
    res = query.store();
    
    if (res.size() > 1) {
      std::cout << "multiple matches for contig: " << contigs[i] << std::endl;
    }
    else if (res.size() == 0) {
      std::cout << "did not find match for contig: " << contigs[i] 
		<< std::endl;
    }
    else {
      row = res[0];
      if (row.size() != 1) {
	std::cout << "invalid row format  for " << contigs[i] << std::endl;
	exit(0);
      }
      contigStart = atoi(row[0]);
      if (i == 0) 
	offset = coord;
      else
	offset = coord - cmlContigs[i-1];
      contigPos = contigStart + offset;
      
      std::cout << chromosome << " " << chr << " " << signSt << " " << sign << " " << start << " " << contigPos << std::endl;
    }
  }
  coordsFile.close();
  return 0;
}

void ReadCMLFile(std::string cmlFileName, 
		 std::vector<std::string> &contigs, 
		 std::vector<ssize_t> &cmlSize) {

  std::ifstream cmlFile;
  openck(cmlFileName, cmlFile);
  std::string contig;
  ssize_t size;
  while (!cmlFile.eof()) {
    if (cmlFile.peek() == EOF)
      break;
    cmlFile >> contig >> size;
    contigs.push_back(contig);
    cmlSize.push_back(size);
  }
  cmlFile.close();
}

void InitEnv(int argc, char* argv[], 
	     std::string &offsetFilename,
	     std::string &coordsFile) {

  ssize_t copt;
  ssize_t i;
  while ( (copt=getopt(argc, argv, "")) != EOF){
    switch(copt) {
    default:
      PrintUsage();
      exit(1);
    }
  }
  i = optind;
  if (i >= argc) {
    std::cout << "You must specify an offset file and coordinate file. " << std::endl;
    PrintUsage();
    exit(1);
  }
  offsetFilename = argv[i];
  i++;
  if (i >= argc) {
    std::cout << "You must specify a coordinate file. " << std::endl;
    PrintUsage();
    exit(1);
  }
  coordsFile = argv[i];
}


void PrintUsage() {
  std::cout << "stub.  your description here. " << std::endl;
  std::cout << "usage: your usage here " << std::endl;
}


