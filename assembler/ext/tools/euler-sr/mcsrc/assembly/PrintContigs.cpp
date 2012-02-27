/***************************************************************************
 * Title:          PrintContigs.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  10/31/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "DeBruijnGraph.h"
#include "IntervalGraph.h"
#include "SeqReader.h"
#include "IntegralTupleStatic.h"



int main(int argc, char* argv[]) {

	BVertexList vertices;
	BEdgeList edges;
  std::string baseInName, bGraphName, contigOutName, edgeFileName;

	int argi = 1;
	if (argc <= 1) {
		std::cout << "usage: printContigs graphBase" << std::endl;
		std::cout << "prints the contigs out to graphBase.contig" << std::endl;
		exit(0);
	}
	baseInName       = argv[argi++];
  
	bGraphName   = baseInName  + ".bgraph";
  edgeFileName = baseInName  + ".edge";
	contigOutName = baseInName + ".contig";

	std::string reportFileName = FormReportName(baseInName);
	ofstream report;
	openck(reportFileName, report, std::ios::app, report);
	BeginReport(argc, argv, report);

	ReadBGraph(bGraphName, vertices, edges, report);
  ReadSequences(edgeFileName, edges, report);


	std::ofstream contigOut;
	openck(contigOutName, contigOut, std::ios::out, report);

	ssize_t e;
	Unmark(edges);
	std::stringstream titleStrm;
	for (e = 0; e < edges.size(); e++ ) {
		if (edges[e].marked != GraphEdge::Marked) {
			titleStrm.str("");
			titleStrm << e << " " << edges[e].length << " " << edges[e].multiplicity;
			edges[e].seq.PrintSeq(contigOut, titleStrm.str());
			contigOut << std::endl;
			edges[e].marked = GraphEdge::Marked;
			edges[edges[e].balancedEdge].marked = GraphEdge::Marked;
		}
	}

	EndReport(report);
	report.close();

	return 0;
}
