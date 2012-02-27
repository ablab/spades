/***************************************************************************
 * Title:          MergeAdjacentEdges.cpp 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  10/31/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "IntervalGraph.h"
#include <stdlib.h>
#include "compatibility.h"
#include "IntegralTupleStatic.h"

using namespace std;

void PrintUsage() {
	std::cout << "usage: mergeAdjacentEdges inGraph outGraph [-maxLength l] [-skipIntervals]" 
						<< std::endl;
}

int main(int argc, char* argv[]) {
	if (argc < 3) {
		PrintUsage();
		exit(0);
	}
	std::string inGraphName, outGraphName;
	int vertexSize;
	inGraphName = argv[1];
	outGraphName= argv[2];

	int argi = 3;
	ssize_t maxLength = 0;
	_INT_ skipIntervals = 0;
	while (argi < argc) {
		if (strcmp(argv[argi], "-maxLength") == 0) {
			maxLength = atoi(argv[++argi]);
		}
		else if (strcmp(argv[argi], "-skipIntervals") == 0) {
			skipIntervals = 1;
		}
		else {
			PrintUsage();
			std::cout << "bad option: " << argv[argi] << std::endl;
			exit(1);
		}
		++argi;
	}
	std::string reportFileName = FormReportName(inGraphName);
	std::ofstream report;
	openck(reportFileName, report, std::ios::app, report);
	BeginReport(argc, argv, report);

	IntervalGraph graph;
	ReadIntervalGraph(inGraphName, graph, vertexSize, skipIntervals, report);
	
	//UNUSED// ssize_t edgesMerged = 0;

	TEdgeList        &edges       = graph.edges;
	TVertexList      &vertices    = graph.vertices;
	std::vector<ssize_t> removedVertices, removedEdges;
	std::cout <<" max length: " << maxLength<< std::endl;
	do {
		//		graph.CalcDMST();
		ssize_t v;
		removedVertices.clear();
		removedEdges.clear();
		for (v =0 ; v < vertices.size(); v++) {

			if (v == graph.LookupBalVertex(v))
				continue;

			// TODO: move the test for balVertex == v somewhere up around here

			ssize_t outEdgeIndex = vertices[v].EndOut();
			// 
			ssize_t nextOutEdgeIndex = vertices[v].NextOut(outEdgeIndex);
			ssize_t outEdge, nextOutEdge;

			for (outEdgeIndex = vertices[v].FirstOut();
					 outEdgeIndex != vertices[v].EndOut();
					 outEdgeIndex = vertices[v].NextOut(outEdgeIndex)) {
				for (nextOutEdgeIndex = vertices[v].NextOut(outEdgeIndex);
						 nextOutEdgeIndex < vertices[v].EndOut();
						 nextOutEdgeIndex = vertices[v].NextOut(nextOutEdgeIndex)) {

					outEdge = vertices[v].out[outEdgeIndex];
					nextOutEdge = vertices[v].out[nextOutEdgeIndex];

					if (edges[outEdge].src == edges[outEdge].dest or
							edges[nextOutEdge].src == edges[nextOutEdge].dest)
						continue;

					std::cout << edges[outEdge].length << " " << edges[nextOutEdge].length << std::endl;
					DNASequence out, next;
					out.namestr = "out";
					next.namestr= "next";
					out.seq     = edges[outEdge].seq.seq;
					out.length  = edges[outEdge].length;
					next.seq    = edges[nextOutEdge].seq.seq;
					next.length = edges[nextOutEdge].length;

					if (!maxLength or 
							(edges[outEdge].length < maxLength and
							 edges[nextOutEdge].length < maxLength)) {

						// Make sure we aren't joining edges with 
						// their balanced edge for now, since I don't
						// know how to do that.

						ssize_t outBal, nextOutBal, balVertex;
						nextOutBal = edges[nextOutEdge].balancedEdge;
						outBal     = edges[outEdge].balancedEdge;

						balVertex = edges[outBal].dest;

						// TODO: could move this earlier
						if (nextOutBal == outEdge
								or nextOutBal == nextOutEdge
								or outBal == outEdge
								or outBal == nextOutEdge
								or balVertex == v // redundant with test put earlier
								) 
							continue;

						std::cout << "merging: " << v << " "<< outEdge << " " << nextOutEdge << " "
											<< edges[outEdge].length << " " << edges[nextOutEdge].length << std::endl;

						//						out.PrintlnSeq(std::cout);
						//						next.PrintlnSeq(std::cout);
						ssize_t nextOutDest = edges[nextOutEdge].dest;
						graph.MergeOutEdges(v, outEdge, nextOutEdge);
					
						removedEdges.push_back(nextOutEdge);
						if (edges[outEdge].dest != edges[nextOutEdge].dest)
							removedVertices.push_back(edges[nextOutEdge].dest);
						cout << "After adding " << edges[nextOutEdge].dest << " for removal: "
								 << vertices[nextOutDest].InDegree() << " "
								 << vertices[nextOutDest].OutDegree() << endl;

						// Perform the corresponding balanced edge merge if possible.
						balVertex = edges[outBal].dest;

						ssize_t nextOutSrc = edges[nextOutBal].src;
						std::cout << "merging balanced: " << nextOutBal << " " << outBal<< std::endl;
						graph.MergeInEdges(balVertex, outBal, nextOutBal);
						
						removedEdges.push_back(nextOutBal);
						if (edges[nextOutBal].src != edges[outBal].src)
							removedVertices.push_back(edges[nextOutBal].src);

						cout << "After adding " << edges[nextOutBal].src << " for removal: "
								 << vertices[nextOutSrc].InDegree() << " "
								 << vertices[nextOutSrc].OutDegree() << endl;

//						CheckEdges(graph.vertices, graph.edges);
					}
					else {
						// Did not remove 'nextOutEdge', so try and merge the current edge with it.
						// I guess I should try all pairs of out edges...
						outEdgeIndex = nextOutEdgeIndex;
					}					
					// Advance to the next out edge.
					nextOutEdgeIndex = vertices[v].NextOut(nextOutEdgeIndex);
				}
			}
		}
		//		assert(graph.CheckGraphStructureBalance());
		cout << "removing " << removedVertices.size() << " v and " << removedEdges.size()
				 << " e." << endl;
		graph.Prune(removedVertices, removedEdges);
		ssize_t ncondensed;
		ncondensed = graph.CondenseSimplePaths();
		graph.RemoveWhirls(maxLength);
		ncondensed += graph.CondenseSimplePaths();
		assert(graph.CheckGraphStructureBalance());
		CheckEdges(graph.vertices, graph.edges);
		cout << "after everything, merged: " << ncondensed << endl;
	}
	while(removedEdges.size() > 0);

	ssize_t e;
	for (e =0 ; e < graph.edges.size(); e++ ){
		assert(graph.edges[e].length == graph.edges[e].seq.length);
	}

	ssize_t numZeroDegree = RemoveZeroDegreeVertices(graph.vertices, graph.edges);
	cout << "removed " << numZeroDegree << " zero degree vertices." << endl;
	std::string edgeOutName, gvzOutName, pathOutName, euGraphName;
	edgeOutName = outGraphName + ".edge";
	gvzOutName  = outGraphName + ".dot";
	pathOutName = outGraphName + ".path";
	euGraphName = outGraphName + ".graph";
	
	BVertex::WriteFixedDegree = 0;
	CondenseEdgeLists(graph.vertices);
	WriteIntervalGraph(outGraphName, graph, vertexSize, report);
  PrintGraph(graph.vertices, graph.edges, euGraphName, report);
  PrintEdges(graph.vertices, graph.edges, edgeOutName, report);
  GVZPrintBGraph(graph.vertices, graph.edges, gvzOutName, report);
	WriteReadPaths(pathOutName, graph.paths, graph.pathLengths, report);

	EndReport(report);
	report.close();
}

