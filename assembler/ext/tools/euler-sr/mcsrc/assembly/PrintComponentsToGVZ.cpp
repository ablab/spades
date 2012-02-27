/***************************************************************************
 * Title:          PrintComponentsToGVZ.cpp
 * Author:         Boyko Kakaradov, Mark Chaisson, Glenn Tesler
 * Created:        2009
 * Last modified:  12/15/2009
 *
 * Copyright (c) 2007-2009 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "DeBruijnGraph.h"
#include "IntervalGraph.h"
#include <string>


int main(int argc, char* argv[]) {
  std::string graphFile;

  if (argc < 2) {
    std::cout << "usage: printComponentsToGVZ reads_file" << std::endl;
    exit(1);
  }
  graphFile = argv[1];
  graphFile += ".bgraph";

	std::string reportFileName = FormReportName(argv[1]);
	ofstream report;
	openck(reportFileName, report, std::ios::app, report);
	BeginReport(argc, argv, report);


  TVertexList vertices;
  TEdgeList edges;
  ReadBGraph(graphFile, vertices, edges, report);

  std::cout << "Done reading " << graphFile << std::endl;

  graphFile = argv[1];
  graphFile += ".dot";
  GVZPrintBGraph(vertices, edges, graphFile, report);

	EndReport(report);
	report.close();

  return 0;
}
