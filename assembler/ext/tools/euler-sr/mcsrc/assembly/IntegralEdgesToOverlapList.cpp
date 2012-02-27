/***************************************************************************
 * Title:          IntegralEdgesToOverlapList.cpp 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  02/28/2010
 *
 * Copyright (c) 2007-2010 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/

#include <string>
#include "compatibility.h"
#include "DeBruijnGraph.h"
#include "SeqReader.h"
#include "IntegralTupleStatic.h"
#include "IntervalGraph.h"

int main(int argc, char* argv[]) {
  std::string edgeFileName;
  _INT_ overlapLength;
  std::string overlapPosFileName;
 
  if (argc < 4) {
    std::cout << "usage: integralEdgesToOverlapList edgeFileName vertexSize "
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
	
  
  EdgePosIntegralTupleList readPositions;

  readPositions.resize(numOverlaps);
	//	EdgePosIntegralTuple::tupleSize = overlapLength;
	EdgePosIntegralTuple::SetTupleSize(overlapLength);
  ssize_t ovp = 0;
  for (e = 0; e < edges.size(); e++ ) {
    for (p = 0; p < edges[e].length - overlapLength + 1; p++) {
			readPositions[ovp].StringToTuple(&edges[e].seq[p]);
      readPositions[ovp].edge = e;
      readPositions[ovp].pos  = p;
      ovp++;
    }
  }
	std::sort(readPositions.begin(), readPositions.end());

  
  std::ofstream overlapOut;
  openck(overlapPosFileName, overlapOut, std::ios::out | std::ios::binary, report);
	//int numOverlaps = readPositions.size();
	// TODO: should we add tupleSize to output, for consistency with .spect files?
	overlapOut.write((char*) &numOverlaps, sizeof(ssize_t));
	overlapOut.write((char*) &readPositions[0], sizeof(EdgePosIntegralTuple) * numOverlaps);
  overlapOut.close();

	EndReport(report);
	report.close();

  return 0;
}
