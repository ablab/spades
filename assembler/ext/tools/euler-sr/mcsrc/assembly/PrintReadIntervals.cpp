/***************************************************************************
 * Title:          PrintReadIntervals.cpp
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  10/31/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <iostream>
#include <string>
#include <vector>
#include "DeBruijnGraph.h"
#include "ReadIntervals.h"
#include "SeqReader.h"
#include "SeqUtils.h"
#include "SimpleSequence.h"
#include "IntegralTupleStatic.h"




void PrintUsage() {
  std::cout << "usage: printReadIntervals overlapFileName edgeFile sequences "
	    << "vertexSize intervalFileName pathFileName" << std::endl;
}

int main(int argc, char* argv[] ) {

  std::string overlapFileName;
  std::string edgeFileName;
  std::string sequenceListFileName;
  std::string readIntervalFileName;
  std::string intervalFileName;
	std::string pathFileName;
  int vertexSize, overlapSize;

  if (argc < 5) {
    PrintUsage();
    exit(1);
  }

  int argi = 1;
  overlapFileName = argv[argi++];
  edgeFileName       = argv[argi++];
  sequenceListFileName = argv[argi++];
  vertexSize = atoi(argv[argi++]);
  readIntervalFileName = argv[argi++];
	pathFileName = argv[argi++];
  overlapSize = vertexSize + 1;
	
	std::string reportFileName = FormReportName(sequenceListFileName);
	std::ofstream reportFile;
	openck(reportFileName, reportFile, std::ios::app, reportFile);
	BeginReport(argc, argv, reportFile);
  SimpleSequenceList sequences, edges;
  // Get the edges (don't care about their reverse complement, since it should be here
  ReadSimpleSequences(edgeFileName, edges, reportFile);
	
	// Try skipping these for now... it may be faster
	//  ReadSimpleSequences(sequenceListFileName, sequences, reportFile);
	//  AppendReverseComplements(sequences);
  
  // Read in the vertex list
  std::ifstream overlapIn;
  ReadPositions overlaps;
  openck(overlapFileName, overlapIn, std::ios::in, reportFile);
  overlapIn >> overlaps;
  overlapIn.close();

	//  std::vector<int> mult;
	//  EdgeIntervalListList edgeIntervals;
	//  EdgeIntervalList::tupleSize = overlapSize;
	PathIntervalList paths;
	PathLengthList pathLengths;
	std::vector<BareReadIntervalList> edgeReadIntervals;
	//	std::vector<BareReadIntervalIndexList> edgeReadIntervalIndices;
	edgeReadIntervals.resize(edges.size());

	std::cout << "storing read intervals. " << std::endl;
	ssize_t skipGapped = 1;
	StoreAllReadIntervals(overlaps, overlapSize, edges, sequenceListFileName, 
												edgeReadIntervals, paths, pathLengths, skipGapped);
	
	std::cout << "sorting edge read interval lists." << edgeReadIntervals.size() << endl;
	ssize_t e, i;
	for (e = 0; e < edgeReadIntervals.size(); e++) {
		SortBareReadIntervalsByReadPos(edgeReadIntervals[e]);
	}
	std::cout << "done, printing\n";
  std::ofstream readIntervalOut;
  openck(readIntervalFileName, readIntervalOut, std::ios::out, reportFile);

  ssize_t edgeIndex;
  ssize_t intervalIndex;
  intervalIndex = 0;
	ssize_t ri;
  for (edgeIndex = 0; edgeIndex < edges.size();edgeIndex++) {
    readIntervalOut << "EDGE " << edgeIndex 
										<< " Length " << edges[edgeIndex].length
										<< " Multiplicity " << edgeReadIntervals[edgeIndex].size()
										<< std::endl;
		cout << edgeReadIntervals[edgeIndex].size() << endl;
		for (ri = 0; ri < edgeReadIntervals[edgeIndex].size(); ri++) {
      readIntervalOut << "INTV " << edgeReadIntervals[edgeIndex][ri].read 
											<< " " << edgeReadIntervals[edgeIndex][ri].readPos 
											<< " " << edgeReadIntervals[edgeIndex][ri].length
											<< " " << edgeReadIntervals[edgeIndex][ri].edgePos << std::endl;
    }
  }
  readIntervalOut.close();
	std::vector<std::vector<ssize_t> > edgeIntervalIndices;
	edgeIntervalIndices.resize(edgeReadIntervals.size());
	for (e = 0; e < edgeReadIntervals.size(); e++) {
		SortBareReadIntervalsByEdgePos(edgeReadIntervals[e]);
		edgeIntervalIndices[e].resize(edgeReadIntervals[e].size());
		for (i = 0; i < edgeIntervalIndices[e].size(); i++ ) {
			edgeIntervalIndices[e][i] = i;
		}
		SortReadIntervalIndicesByRead(edgeReadIntervals[e], edgeIntervalIndices[e]);
	}

	PathReadPosToIntervalIndex(edgeReadIntervals, edgeIntervalIndices, paths, pathLengths);
	WriteReadPaths(pathFileName, paths, pathLengths, reportFile);
	EndReport(reportFile);
	return 0;
}


