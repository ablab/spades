/***************************************************************************
 * Title:          SortVertexList.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  10/31/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <vector>
#include "utils.h"
#include "DeBruijnGraph.h"
#include "DNASequence.h"
#include "SeqReader.h"
#include "SeqUtils.h"
#include "IntervalGraph.h"
#include "IntegralTupleStatic.h"


void PrintUsage() {
  std::cout << "usage: sortVertexList readPos.in sequences vertexSize readPos.out " << std::endl;
  std::cout
		<< "  readPos.in    The list of read positions.  The first line has the number of " << std::endl
		<< "                read positions, and each line contains one read position, " << std::endl
		<< "                (read, pos). " << std::endl
		<< "  sequences     The list of sequences (without their reverse complement " << std::endl
		<< "  vertexSize    The length of the tuple to sort on " << std::endl
		<< "  readPos.out   Sort, and store in this file. " << std::endl;
}

int main(int argc, char* argv[]) {

  if (argc < 5) {
    PrintUsage();
    exit(1);
  }
  std::string inFileName, outFileName, readFileName;
  int tupleSize;
  inFileName = argv[1];
  readFileName = argv[2];
  tupleSize = atoi(argv[3]);
  outFileName = argv[4];

	std::string reportFileName = FormReportName(inFileName);
	ofstream report;
	openck(reportFileName, report, std::ios::app, report);
	BeginReport(argc, argv, report);



  ReadPositions readPositions;
  std::ifstream in;
  
  openck(inFileName, in, std::ios::in, report);

  // Read in the read positions
  in >> readPositions;

  // Read in all reads
  SimpleSequenceList sequences;
  ReadSimpleSequences(readFileName, sequences, report);
  AppendReverseComplements(sequences);

  // Set up the functor to compare the tuples
  CompareTuples<SimpleSequenceList> comp;
  comp.sequencesPtr = &sequences;
  comp.length = tupleSize;

  // Do quicksort on the read positions
  std::sort(readPositions.begin(), readPositions.end(), comp);

  // Done, print them
  std::ofstream out;
  openck(outFileName, out, std::ios::out, report);
  out << readPositions;
  out.close();

	EndReport(report);
	report.close();

  return 0;
}
