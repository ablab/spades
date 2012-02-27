/***************************************************************************
 * Title:          EdgesToOverlapList.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  10/31/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <string>
#include "DeBruijnGraph.h"
#include "SeqReader.h"

#include "IntervalGraph.h"
#include "IntegralTupleStatic.h"

int main(int argc, char* argv[]) {
  std::string edgeFileName;
  ssize_t overlapLength;
  std::string overlapPosFileName;
 
  if (argc < 4) {
    std::cout << "usage: edgesToOverlapList edgeFileName vertexSize "
	      << "overlapListFile " << std::endl;
    exit(1);
  }

  int argi = 1;
  edgeFileName = argv[argi++];
  overlapLength = atoi(argv[argi++]) + 1;
  overlapPosFileName = argv[argi++];

	std::string reportFileName = FormReportName(edgeFileName);
	ofstream report;
	openck(reportFileName, report, std::ios::app, report);
	BeginReport(argc, argv, report);


  // every overlap should be unique, so there is no reason to have 
  // any special data structure to find unique overlaps.

  ssize_t numOverlaps = 0;
  //UNUSED+// ssize_t  s;
  ssize_t e, p;

  SimpleSequenceList edges;
  ReadSimpleSequences(edgeFileName, edges, report);

  for (e = 0; e < edges.size(); e++ ){ 
    numOverlaps += edges[e].length - overlapLength + 1;
  }
  
  ReadPositions readPositions;

  readPositions.resize(numOverlaps);
  ssize_t ovp = 0;
  for (e = 0; e < edges.size(); e++ ) {
    for (p = 0; p < edges[e].length - overlapLength + 1; p++) {
      readPositions[ovp].read = e;
      readPositions[ovp].pos  = p;
      ovp++;
    }
  }

  CompareTuples<SimpleSequenceList> comp;
  comp.sequencesPtr = &edges;
  comp.length = overlapLength;
  std::sort(readPositions.begin(), readPositions.end(), comp);
  
  std::ofstream overlapOut;
  openck(overlapPosFileName, overlapOut, std::ios::out, report);

  overlapOut << readPositions;
  overlapOut.close();

	EndReport(report);
	report.close();

  return 0;
}
