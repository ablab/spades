/***************************************************************************
 * Title:          PrintGraphSummary.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/21/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "DeBruijnGraph.h"
#include "IntervalGraph.h"
#include "RepeatSearch.h"
#include <string>
#include "IntegralTupleStatic.h"

#include <fstream>
#include "utils.h"

void PrintSimpleRepeats(IntervalGraph &graph) {
  ssize_t e;
	ssize_t src,dest;
  for (e = 0; e < graph.edges.size(); e++) {
		src = graph.edges[e].src; 
		dest = graph.edges[e].dest;
		if (graph.vertices[src].InDegree() == 
				graph.vertices[dest].OutDegree() and
				graph.vertices[src].OutDegree() == 1 and
				graph.vertices[dest].InDegree() == 1) {
			std::cout << e << " " << graph.edges[e].length << " " 
								<< graph.edges[e].intervals->size() << " ";
			ssize_t edgeIndex, edgeIndexI;
			for (edgeIndexI = graph.vertices[src].FirstIn();
					 edgeIndexI < graph.vertices[src].EndIn();
					 edgeIndexI = graph.vertices[src].NextIn(edgeIndexI)) {
				edgeIndex = graph.vertices[src].in[edgeIndexI];
				std::cout << graph.edges[edgeIndex].length << " ";
			}
			for (edgeIndexI = graph.vertices[dest].FirstOut();
					 edgeIndexI < graph.vertices[dest].EndOut();
					 edgeIndexI = graph.vertices[dest].NextOut(edgeIndexI)) {
				edgeIndex = graph.vertices[dest].out[edgeIndexI];
				std::cout << graph.edges[edgeIndex].length << " ";
			}
			std::cout << std::endl;
		}
	}
}

void PrintCycles(IntervalGraph &graph, ssize_t maxCycleLength) {
  //UNUSED// ssize_t v;
  std::set<ssize_t> cycleEdges, cycleSearch, cycleVertices;
  ssize_t e;
  for (e = 0; e < graph.edges.size(); e++) {
    if (graph.SearchForUndirectedCycle2(graph.edges[e].src, graph.edges[e].dest, cycleEdges, 
					graph.edges[e].length,
					maxCycleLength, 0, cycleSearch, cycleVertices)) {
			// TODO: check -- does it print as a side-effect?
		}
  }
}

void PrintEarlyEnds(IntervalGraph &g, ssize_t readLength ) {
	ssize_t e, i;
	ssize_t destVertex;
	ssize_t edgeLength = 0;
	ssize_t maxOverpass;
	ssize_t readEnd, overpass;
	for (e = 0; e < g.edges.size(); e++ ) {
		destVertex = g.edges[e].dest;
		if (g.vertices[destVertex].OutDegree() == 0) {
			edgeLength = g.edges[e].length;
			maxOverpass = 0;
			for (i = 0; i < g.edges[e].intervals->size(); i++) {
				readEnd = (*g.edges[e].intervals)[i].edgePos + readLength - (*g.edges[e].intervals)[i].readPos;
				overpass = readEnd - edgeLength;
				if (overpass > maxOverpass) {
					maxOverpass = overpass;
				}
			}																					
			if (maxOverpass > 0) {
				std::cout << e << " " << g.edges[e].length 
									<< " " << maxOverpass << std::endl;
			}
		}
	}
}

void PrintEdgeMults(IntervalGraph &g) {

	ssize_t e;
	ssize_t newMult;
	for (e = 0; e < g.edges.size(); e++) {
		newMult = g.CountReadsContainedInEdge(e) +
			g.CountReadsPassingThroughEdge(e) +
			g.CountReadsExtendingIntoEdge(e, 40);
		std::cout << e << " " << g.edges[e].length << " " << g.edges[e].intervals->size()
							<< " " << newMult << std::endl;
	}
}

void PrintRepeatLengths(TVertexList &vertices, TEdgeList &edges, std::string repeatFileName) {
  ssize_t v;
  ssize_t outEdge, outEdgeIndex;
  std::ofstream repeatOut; 
  if (repeatFileName != "" ) {
    openck(repeatFileName, repeatOut, std::ios::out);
  }
  ssize_t repeatNumber = 0;
  for (v = 0; v < vertices.size(); v++) {
    if (vertices[v].InDegree() > 1 and
				vertices[v].OutDegree() == 1) {
      outEdgeIndex = vertices[v].FirstOut();
      outEdge = vertices[v].out[outEdgeIndex];
      std::cout << vertices[v].InDegree() << " " << edges[outEdge].length << std::endl;
      if ( repeatFileName != "" ) {
				std::stringstream titleStrm;
				titleStrm << "repeat_" << repeatNumber << " multiplicity " << vertices[v].InDegree() 
									<< " " << edges[outEdge].length;
				edges[outEdge].seq.PrintSeq(repeatOut, titleStrm.str());
				std::cout << std::endl;
      }
      ++repeatNumber;
    }
  }
  if (repeatFileName != "")
    repeatOut.close();
}

void PrintVertexStatistics(TVertexList &vertices, TEdgeList &edges) {
  std::vector<ssize_t> inDeg;

  ssize_t v;
  ssize_t prevSize;
  ssize_t c;
  ssize_t inDegree;
  for (v = 0; v < vertices.size(); v++ ) {
    inDegree = vertices[v].InDegree();
    if (inDegree > inDeg.size()) {
      prevSize = inDeg.size();
      inDeg.resize(inDegree+1);
      for (c = prevSize; c < inDegree; c++) {
				inDeg[c] = 0;
      }
      inDeg[inDegree]++;
    }
  }
}


void PrintUsage() {
	std::cout << "usage: printGraphSummary file.bgraph [options]" << std::endl;
	std::cout << " where options includes: " << std::endl
						<< "  -components          Print a summary of the components in the graph " << std::endl
						<< "  -vstats              Print a summary of the in/out statistics of a vertex" << std::endl
						<< "  -repeats             Print a summary of suspected repeat lengths" << std::endl
						<< "  -repeatFile file     Print a summary of suspected repeat lengths into 'file'" << std::endl
						<< "  -sources             Print the lengths of sources " << std::endl
						<< "  -sinks               Print lengths of the sinks " << std::endl
						<< "  -edges edgefile      Optionally read in an edgefile. Useful if sequences" << std::endl
						<< "                       are needed." << std::endl
						<< "  -repeatDist          Print the minimum distance between two repeats. " << std::endl
						<< "  -earlyEnds readlen   Print edges that have reads that should go past the" << std::endl
						<< "                       end of the edge." << std::endl
						<< "  -simpleRepeats       Print edges that follow canonical tangle. " << std::endl
						<< "  -printCycles maxCycleLength   Print cycles (TODO)" << std::endl
						<< "  -printEM             Print edge multiplicities" << std::endl
						<< "  -vertexSize v        Vertex size (default 20)" << std::endl
						<< std::endl;
}


int main(int argc, char* argv[]) {
  std::string edgeFile, repeatFile;
  edgeFile = "";
  repeatFile = "";
  if (argc < 2) {
    PrintUsage();
    exit(1);
  }
  std::string graphBase = argv[1];
	
	std::string graphFile = graphBase + ".bgraph";
	std::string intvFile  = graphBase + ".intv";
	std::string pathFile  = graphBase + ".path";
  int argi = 2;
  ssize_t printComponents = 0;
  ssize_t printVStats     = 0;
  ssize_t printRepeats    = 0;
  ssize_t printSources    = 0;
  ssize_t printSinks      = 0;
  ssize_t printCycles     = 0;
  ssize_t printRepeatDist = 0;
	ssize_t printEdgeMult   = 0;
	ssize_t printEarlyEnds  = 0;
	ssize_t printSimpleRepeats = 0;
  int vertexSize = 20;
  while (argi < argc) {
    if (strcmp(argv[argi], "-components")== 0)
      printComponents = 1;
    else if (strcmp(argv[argi], "-vstats") == 0)
      printVStats = 1;
    else if (strcmp(argv[argi], "-repeats") == 0)
      printRepeats = 1;
    else if (strcmp(argv[argi], "-sources") == 0)
      printSources = 1;
    else if (strcmp(argv[argi], "-sinks") == 0)
      printSinks = 1;
    else if (strcmp(argv[argi], "-edges") == 0)
      edgeFile = argv[++argi];
    else if (strcmp(argv[argi], "-repeatFile") == 0)
      repeatFile = argv[++argi];
    else if (strcmp(argv[argi], "-repeatDist") == 0)
      printRepeatDist = 1;
    else if (strcmp(argv[argi], "-printCycles") == 0)
      printCycles = atoi(argv[++argi]);
		else if (strcmp(argv[argi], "-printEM") == 0)
			printEdgeMult = 1;
    else if (strcmp(argv[argi], "-vertexSize") == 0)
      vertexSize = atoi(argv[++argi]);
		else if (strcmp(argv[argi], "-earlyEnds") == 0) {
			printEarlyEnds = atoi(argv[++argi]);
		}
		else if (strcmp(argv[argi], "-simpleRepeats") == 0) {
			printSimpleRepeats = 1;
		}
    else {
      PrintUsage();
			std::cout << "bad opeion: " << argv[argi]<< std::endl;
      exit(1);
    }
    argi++;
  }
  
  IntervalGraph graph;
  graph.vertexSize = vertexSize;
	//  graph.ReadIntervalGraph(graphFile, intvFile, pathFile, 0);
	ReadIntervalGraph(graphBase, graph, vertexSize, 1);
  Unmark(graph.vertices);
  Unmark(graph.edges);
  
  if ( edgeFile != "" ) {
    ReadSequences(edgeFile, graph.edges);
  }
  if (printComponents) {
    PrintComponents(graph.vertices, graph.edges);
  }
  if (printVStats) {
    PrintVertexStatistics(graph.vertices, graph.edges);
  }
  if (printRepeats or repeatFile != "") {
    PrintRepeatLengths(graph.vertices, graph.edges, repeatFile);
  }
  if (printCycles) {
    PrintCycles(graph, printCycles );
  }
	if (printEdgeMult) {
		PrintEdgeMults(graph);
	}
	if (printEarlyEnds) {
		PrintEarlyEnds(graph, printEarlyEnds);
	}
	if (printSimpleRepeats) {
		PrintSimpleRepeats(graph);
	}
  if (printSources ) {
    ssize_t v;
    std::cout << "sources " << std::endl;
    for (v = 0; v < graph.vertices.size(); v++ ) {
      if (graph.vertices[v].InDegree() == 0 and
					graph.vertices[v].OutDegree() == 1) {
				std::cout << v << " " << graph.edges[graph.vertices[v].out[graph.vertices[v].FirstOut()]].length << std::endl;
      }
    }
  }
  if (printSinks ) {
    ssize_t v;
    std::cout << "sinks: " << std::endl;
    for (v = 0; v < graph.vertices.size(); v++ ) {
      if (graph.vertices[v].OutDegree() == 0 and
					graph.vertices[v].InDegree() == 1) {
				std::cout << v << " " << graph.edges[graph.vertices[v].in[graph.vertices[v].FirstIn()]].length << std::endl;
      }
    }
  }

  if (printRepeatDist) {
    FindShortestDistanceBetweenTwoRepeats(graph.vertices, graph.edges);
  }


  return 0;
}
