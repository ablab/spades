/***************************************************************************
 * Title:          DeBruijn.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/28/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <vector>
#include <iostream>
#include <algorithm>

// sequence utilities
#include "SeqReader.h"
#include "SeqUtils.h"
#include "SimpleSequence.h"
#include "DNASequence.h"
#include "utils.h"

// Definitions for the de bruijn graph
#include "Vertex.h"
#include "Edge.h"
#include "DeBruijnGraph.h"

// Definitions for the branching graph.
#include "BEdge.h"
#include "BVertex.h"

#include "IntegralTupleStatic.h"

ssize_t allocedMem;


void CondensePaths(VertexList &vertices, SimpleSequenceList &sequences, int tupleSize);
void SimplePrintVertex(Vertex &v, SimpleSequenceList &sequences, int tupleSize);
void SimplePrintGraph(VertexList &graph, SimpleSequenceList &sequences, int tupleSize);

void EnumeratePositions(SimpleSequenceList &sequences, ReadPositions &readPositions, 
												int tupleSize );

void AppendDisjointGraph(BVertexList &vSetA, BEdgeList &eSetA,
												 BVertexList &vSetB, BEdgeList &eSetB);

ssize_t CountPositions(SimpleSequenceList &sequences, int tupleSize);

template<class T>
void PrintPositions(SimpleSequenceList &sequences, std::vector<T> &readPositions, int tupleSize) {
  DNASequence tmpSeq;
  ssize_t p;
  for (p = 0; p < readPositions.size(); p++ ) {
    std::cout << readPositions[p].read << " " << readPositions[p].pos << " ";
    PrintTuple(sequences, readPositions[p], tupleSize);
    std::cout << std::endl;
  }
}
void CondenseSimplePath(ssize_t vertex, ssize_t out, VertexList &vertices, 
												SimpleSequenceList &sequences, int tupleSize);


void FindCycles(VertexList &deBruijn, SimpleSequenceList &sequences, int tupleSize, ssize_t balanceGraph,
								BVertexList &cycleVertices, BEdgeList &cycleEdges);


ssize_t CountUniqueTuples(SimpleSequenceList &sequences, ReadPositions &readPositions,
											int tupleSize);


void AllocateBGraph(VertexList &deBruijn,
										BVertexList &bVertices,
										BEdgeList &bEdges,
										SimpleSequenceList &sequences, int tupleSize);

void CreateBGraph(VertexList &deBruijn, 
									BVertexList &bGraph, BEdgeList &edges,
									SimpleSequenceList &sequences, int tupleSize);


void BalanceGraph(BVertexList &graph, BEdgeList &edges,SimpleSequenceList &sequences, int tupleSize);


void CreateGraph(SimpleSequenceList &sequences, 
								 VertexList &vertices, int tupleSize);


void GVZPrintGraph(VertexList &graph, 
									 SimpleSequenceList &sequences, int tupleSize,std::string &graphOutName, std::ostream &report);

void CountDegrees(VertexList &vertices, ssize_t counts[]);


void CountSimpleSize(VertexList &vertices, ssize_t &numVertices, ssize_t &numEdges);

ssize_t AdvanceSimplePath(ssize_t vertex, ssize_t out, VertexList &vertices, 
											ssize_t &lastVertex, ssize_t &nextToLastVertex, 
											unsigned char* seq=NULL );

void AdvanceCycle(ssize_t curVertex, VertexList &deBruijn, 
									SimpleSequenceList &sequences, 
									int tupleSize,
									BVertexList &cycleVertices, 
									BEdgeList &cycleEdges);

void ClearSimplePath(ssize_t vertex, ssize_t out, VertexList &vertices);

void PrintUsage() {
  std::cout << std::endl 
						<< "Build a de Bruijn graph using a set of vertices and overlaps defined" << std::endl
						<< "by a set of sequences." << std::endl << std::endl;
  std::cout << "usage: debruijn seqFile vertexListFile graphFile [options]" 
						<< std::endl
						<< "  [-vertexSize w]       Build graph with vertices as w." << std::endl
						<< "  [-noCondense]         Do not condense simple paths (in degree = out degree = 1)." << std::endl
						<< "  [-edgeOutFile filename]  Output edges to filename." << std::endl
		//						<< "  [-intvFile filename]  Output read intervals to filename." << std::endl
		// TODO: Investigate why there is no "-intvFile" option implemented.
						<< "  [-noGraph]            Suppress printing of .graph file (default is on)." << std::endl
						<< "  [-noEdge]             Suppress printing of .edge file (default is on)." << std::endl
						<< "  [-bGraph filename]    Print bgraph to filename (default seqFile.bgraph)." << std::endl
		//						<< "  [-trimShort N]        Trim short edges (< N) from graph." << std::endl // TODO: implement
						<< "  [-graphFile filename] Print graph to filename (default seqFile.graph)." << std::endl
						<< std::endl << std::endl;
}

int main(int argc, char*argv[]) {
  allocedMem = 0;
	std::string reportFileName;
  std::string sequencesFileName;
  std::string vertexListFileName;
  std::string gvzGraphOutName;
  std::string graphFileName;
  std::string baseOut;
  std::string edgeOutputFileName;
  std::string bgraphOutName;
  
  ssize_t printEdges = 1;
  ssize_t printGraph = 1;
  int tupleSize = 4;
  int vertexSize = 3;
  ssize_t condensePaths = 1;

  int argi = 1;
  if (argc < 4) {
    PrintUsage();
    exit(1);
  }
  sequencesFileName    = argv[argi++];
  vertexListFileName   = argv[argi++];
  gvzGraphOutName      = argv[argi++];
  graphFileName        = sequencesFileName + ".graph";
  edgeOutputFileName   = sequencesFileName + ".edge";
  bgraphOutName        = sequencesFileName + ".bgraph";
	reportFileName       = FormReportName(sequencesFileName);

  while(argi < argc) {
    if (strcmp(argv[argi], "-vertexSize") == 0) {
      ++argi;
      vertexSize = atoi(argv[argi]);
      tupleSize = vertexSize + 1;
    }
    else if (strcmp(argv[argi], "-bGraph") == 0) {
      bgraphOutName = argv[argi++];
    }
    else if (strcmp(argv[argi], "-noCondense") == 0) {
      condensePaths = 0;
    }
    else if (strcmp(argv[argi], "-edgeOutFile") == 0) {
      ++argi;
      edgeOutputFileName = argv[argi];
      printEdges   = 1;
    }
    else if (strcmp(argv[argi], "-noGraph") == 0 ) {
      printGraph = 0;
    }
    else if (strcmp(argv[argi], "-noEdge") == 0) {
      printEdges = 0;
    }
    else if (strcmp(argv[argi], "-graphFile") == 0) {
      ++argi;
      graphFileName = argv[argi];
      printEdges   = 1;
    }
    else {
      std::cout << "bad argument: " << argv[argi] << std::endl;
      PrintUsage();
      exit(1);
    }
    argi++;
  }
	std::ofstream report;
	openck(reportFileName, report, std::ios::app, report);
	BeginReport(argc, argv, report);
	report << "vertex\t" << vertexSize << std::endl;
  // Read in and store every sequence (and its reverse complement)
  DNASequence seq, seqRC;
  SimpleSequenceList sequences;

  std::ifstream seqIn;
  openck(sequencesFileName, seqIn, std::ios::in, report);

  SimpleSequence simpleSeq;
  DNASequence newSeq;
  ReadSimpleSequences(sequencesFileName, sequences, report);
  AppendReverseComplements(sequences);
  _SZT_ seqMemory;
  seqMemory = GetSeqListMemoryUsage(sequences);
  std::cout << "sequences require " << seqMemory << std::endl;
  // now create a list of l-1 mers
  ReadPositions tuples;
  ReadPositions overlaps;
  VertexList vertices;

  // The list of vertices should have been previously computed
  std::ifstream vertexIn, edgeIn;
  openck(vertexListFileName, vertexIn, std::ios::in, report);
  vertexIn >> vertices;
  vertexIn.close();

  // Build the graph for every vertex
  std::cout << "Creating the graph." << std::endl;
  CreateGraph(sequences, vertices, tupleSize);
  std::cout << "The graph has " << vertices.size() << " vertices (" 
						<< sizeof(Vertex) * vertices.size() << " bytes)." << std::endl;

  //  SimplePrintGraph(vertices, sequences, tupleSize-1);

  BVertexList bVertices;
  BEdgeList bEdges;
  if (condensePaths) {
    std::cout << "Creating bgraph." << std::endl;
    CreateBGraph(vertices, bVertices, bEdges, sequences, tupleSize);
    BalanceGraph(bVertices, bEdges, sequences, tupleSize);
    
    // Cycles are non-branching, cycles in the graph that contain more 
    // than one vertex but have no degree-0 in or out vertices.
    BVertexList cycleVertices;
    BEdgeList cycleEdges;
    ssize_t balanceGraph = 1;
    FindCycles(vertices, sequences, tupleSize, balanceGraph, cycleVertices, cycleEdges );
    AppendDisjointGraph(bVertices, bEdges, cycleVertices, cycleEdges);
  }

  ssize_t nRemoved;
  nRemoved = RemoveZeroDegreeVertices(bVertices, bEdges);
  std::cout << "removed " << nRemoved << " zero degree vertices " << std::endl;
  if (! condensePaths ) {
    GVZPrintGraph(vertices, sequences, tupleSize, gvzGraphOutName, report);
  }
  else {
    GVZPrintBGraph(bVertices, bEdges, gvzGraphOutName, report);
  }
  if (printEdges) 
    PrintEdges(bVertices, bEdges, edgeOutputFileName, report);
  if (printGraph)
    PrintGraph(bVertices, bEdges, graphFileName, report);

  PrintBGraph(bVertices, bEdges, tupleSize - 1, bgraphOutName, report);

	EndReport(report);
	report.close();

  return 0;
}



ssize_t CountPositions(SimpleSequenceList &sequences, int tupleSize) {
  ssize_t s;
  ssize_t numPos = 0;
  for (s = 0; s < sequences.size(); s++ ) {
    if (sequences[s].length < tupleSize) {
      std::cout << "Error, repeat size is less than the read length " << std::endl;
      exit(0);
    }
    numPos += sequences[s].length - tupleSize + 1;
  }
  return numPos;
}


void EnumeratePositions(SimpleSequenceList &sequences, ReadPositions &readPositions, int tupleSize ) {
  ssize_t numPositions;
  numPositions = CountPositions(sequences, tupleSize);
  readPositions.resize(numPositions);
  allocedMem += sizeof (ReadPos) * numPositions;
  std::cout << "allocated " << sizeof (ReadPos) * numPositions 
						<< " (" << allocedMem <<")"<< std::endl;
  ssize_t pos = 0;
  ssize_t p, s;
  for (s = 0; s < sequences.size(); s++ ) {
    for (p = 0; p < sequences[s].length - tupleSize + 1; p++) {
      readPositions[pos].read = s;
      readPositions[pos].pos  = p;
      pos++;
    }
  }
  CompareTuples<SimpleSequenceList> comp;
  comp.sequencesPtr = &sequences;
  comp.length = tupleSize;
  //  std::cout << "before sorting " << std::endl;
  //  PrintPositions(sequences, readPositions,tupleSize);
  //  exit(0);
  std::sort(readPositions.begin(), readPositions.end(), comp);
}



ssize_t CountUniqueTuples(SimpleSequenceList &sequences, ReadPositions &readPositions,
											int tupleSize) {
  ssize_t i = 0;
  ssize_t numUnique = 0;
  //  std::cout << "counting unique tuples in " << readPositions.size() << " read positions "<< std::endl;
  while (i < readPositions.size()) {
    ++numUnique;
    i++;
    while (i < readPositions.size() and
					 CompareTupleSeq(readPositions[i], readPositions[i-1], sequences, tupleSize) == 0) {
      i++;
    }
  }
  return numUnique;
}



void CreateGraph(SimpleSequenceList &sequences, 
								 VertexList &vertices, int tupleSize) {
  //UNUSED// ssize_t numVertices;
  ssize_t v;
  ssize_t p;
  ssize_t s;
  //UNUSED// unsigned char* tupleSeq;
  ssize_t toVertex=-1, fromVertex=-1;
  unsigned char outChar, inChar;

  int vertexSize = tupleSize - 1;

  // Initialize the vertex list
  v = 0;
  p = 0;
  ssize_t lastN = -1;
  ssize_t n;
  for (s = 0; s < sequences.size(); s++ ) {
    // Make sure there isn't an 'N' in the first tuple
		if (sequences[s].length < vertexSize  + 1)
			continue;
    lastN = -1;
    for (n = vertexSize - 1; n >= 0; n-- ) {
      if (unmasked_nuc_index[sequences[s].seq[n]] >= 4) {
				lastN = n;
				fromVertex = -1;
				break;
      }
    }
    ssize_t p = 1;
    // if there are N's in the tuples, move past them
    if (lastN < 0) {
      fromVertex = LocateTuple(sequences, vertices, vertexSize,(char*) &(sequences[s].seq[p-1]));
      assert(fromVertex >= 0);
    }
    else {
			// if lastN >= 0, there is an N in the first vertex, and so it shouldn't be considered
      assert(fromVertex == -1);
    }

    for (; p < sequences[s].length - vertexSize+ 1; p++) {
      if (unmasked_nuc_index[sequences[s].seq[p + vertexSize  - 1]] >= 4) {
				lastN = p + vertexSize - 1;
				// Do not save the last vertex, since it is skipped by an 'N'
				fromVertex = -1;
      }
      // If we are able to lookup the position of this vertex (lastN < p), do so
      if (lastN < p) {
				toVertex = LocateTuple(sequences, vertices, vertexSize, (char*) &(sequences[s].seq[p])); 
				// If we have a previous vertex (there wasn't an 'N'), assign it here.
				// Otherwise, we skip this first step, and continue after assigning fromVertex
				// at the end of the if statement.
				if (fromVertex > -1) {
					assert(toVertex >= 0);
	  
					// Store the out vertices
					outChar  = sequences[s].seq[p + vertexSize - 1];
					inChar   = sequences[s].seq[p-1];

					if (unmasked_nuc_index[outChar] >= 4 or 
							unmasked_nuc_index[inChar] >= 4) {
						std::cout << "This read makes bad things happen at " << p << std::endl;
						DNASequence tmp;
						tmp = sequences[s];
						tmp.PrintSeq(std::cout);
						std::cout << std::endl;
					}
					assert(unmasked_nuc_index[outChar] < 4);
					assert(unmasked_nuc_index[inChar] < 4);
					vertices[fromVertex].out[(unsigned char) unmasked_nuc_index[outChar]] = toVertex;
					vertices[toVertex  ].in[ (unsigned char) unmasked_nuc_index[inChar]]  = fromVertex;
				}
				fromVertex = toVertex;
      }
    }
  }
}

void SimplePrintVertex(Vertex &v, SimpleSequenceList &sequences, int tupleSize) {
  ssize_t e;
  PrintTuple(sequences, v, tupleSize);
  std::cout << " pos: " << v.pos 
						<< " read: " << v.read << std::endl;
  std::cout << "   in edge ";
  for (e = 0; e < 4; e++) {
    std::cout << "\t" << e << " to " << v.in[e];
  }
  std::cout << std::endl;
  std::cout << "  out edge: ";
  for (e = 0; e < 4; e++) {
    std::cout << "\t" << e << " to " << v.out[e];
  }
  std::cout << std::endl;
}

void SimplePrintGraph(VertexList &graph, SimpleSequenceList &sequences, int tupleSize) {
  ssize_t v;
  for (v = 0; v < graph.size(); v++) {
    std::cout << "vertex: " << v << " ";
    SimplePrintVertex(graph[v], sequences, tupleSize);
  }
}
		      
		

void GVZPrintGraph(VertexList &graph, 
									 SimpleSequenceList &sequences, int tupleSize,std::string &graphOutName,
									 std::ostream &report) {
  std::ofstream graphOut;
  openck(graphOutName, graphOut, std::ios::out, report);

  ssize_t v, e;
  graphOut << "digraph G {"<<std::endl
					 << "\tsize=\"8,8\";"<<std::endl;
  for (v = 0; v < graph.size(); v++ ) {
    // If this vertex has been marked as removed, continue 
    // without printing it
    //    std::cout << "vertex " << v << " degree: " << graph[v].Degree() << std::endl;
    if (graph[v].OutDegree() > 0) {
      for (e = 0; e < 4; e++ ) {
				if (graph[v].out[e] != -1) {
					graphOut << "\t" << v << " -> " << graph[v].out[e] << std::endl;
				}
      }
      graphOut << "\t" << v << " [label=\"";
      if (tupleSize <= 6) {
				PrintTuple(sequences, graph[v], tupleSize-1, graphOut);
				graphOut << ", ";
      }
      graphOut << " " << v << "\"];" << std::endl;
    }
  }
  graphOut << "}" << std::endl;
}


void CountDegrees(VertexList &vertices, ssize_t counts[]) {
  ssize_t v;
  ssize_t d;
  counts[0] = counts[1] = counts[2] = counts[3] = 0;
  for (v = 0; v < vertices.size(); v++ ) {
    d = vertices[v].OutDegree();
    if (d > 4) {
      std::cout << "outdegree: " << d << std::endl;
      assert(d <= 4);
    }
    if (!vertices[v].Singleton()) 
      counts[d]++;
  }
}

void ClearSimplePath(ssize_t vertex, ssize_t out, VertexList &vertices) {
  ssize_t prevVertex;
  ssize_t outIndex;
	ssize_t firstVertex = vertex; // TODO: check this is the correct initialization
	ssize_t prevOutIndex;
  do {
    prevVertex = vertex;
		prevOutIndex = out;
		vertex = vertices[vertex].out[out];
    outIndex = vertices[prevVertex].FirstOutIndex();
		assert(prevVertex >= 0);
		assert(prevOutIndex >= 0 && prevOutIndex < 4);
    vertices[prevVertex].out[prevOutIndex] = -1;
		out = vertices[vertex].FirstOutIndex();
  }
  while (vertex != firstVertex and
				 vertex >= 0 and
				 vertices[vertex].InDegree() == 1 and
				 vertices[vertex].OutDegree() == 1);
}
		    

ssize_t AdvanceSimplePath(ssize_t vertex, ssize_t out, VertexList &vertices, 
											ssize_t &lastVertex, ssize_t &nextToLastVertex,
											unsigned char* seq ) {
  // move forward along a path until a branching vertex has been reached
  ssize_t firstVertex = vertex;
  nextToLastVertex = vertex;
  vertex = vertices[vertex].out[out];
  ssize_t length = 1;
  while (vertex >= 0 and
				 vertices[vertex].InDegree() == 1 and
				 vertices[vertex].OutDegree() == 1) {
    if (seq != NULL) {
      seq[length - 1] = nuc_char[out];
    }
    out = vertices[vertex].DegreeOneOutEdge();
    nextToLastVertex = vertex;
    assert(out < 4);
    vertex = vertices[vertex].out[out];
    ++length;
    if (vertex == firstVertex)
      break;
  }
  if (seq != NULL) {
    seq[length - 1] = nuc_char[out];
  }
  lastVertex = vertex;

  return length;
}


void AllocateBGraph(VertexList &deBruijn,
										BVertexList &bVertices,
										BEdgeList &bEdges,
										SimpleSequenceList &sequences, int tupleSize) {
  ssize_t numBVertices, numBEdges;
	// Count the number of edges and vertices in the graph
	// that has all simple paths condensed into one
  CountSimpleSize(deBruijn, numBVertices, numBEdges);
	std::cout << "branching graph:" << numBVertices << " vertices " 
						<< numBEdges << " edges." << std::endl;
  bVertices.resize(numBVertices);
  bEdges.resize(numBEdges);
  
  //UNUSED+// ssize_t e;
  ssize_t v ;
  ssize_t bv = 0;
  for (v = 0; v < deBruijn.size(); v++ ) {
    if (deBruijn[v].IsBranch()) {
      bVertices[bv].index = v;
      bVertices[bv].read  = deBruijn[v].read;
      bVertices[bv].pos   = deBruijn[v].pos;	
      bv++;
    }
  }
  // Make the list of vertices binary-searchable
  CompareTuples<SimpleSequenceList> compare;
  compare.sequencesPtr = &sequences;
  compare.length = tupleSize - 1;
  std::sort(bVertices.begin(), bVertices.end(), compare);
}

void AdvanceCycle(ssize_t curVertex, VertexList &deBruijn, 
									SimpleSequenceList &sequences, 
									int tupleSize,
									BVertexList &cycleVertices, 
									BEdgeList &cycleEdges) {
  ssize_t inIndex, outIndex;
  cycleVertices.push_back(BVertex());
  cycleEdges.push_back(BEdge());
  ssize_t curCycleVertex = cycleVertices.size()-1;
  ssize_t curCycleEdge   = cycleEdges.size()-1;
	  
  ssize_t readPos, readIndex;
  char *tuplePtr;

  ssize_t cycleVertex;
  ssize_t lastVertex, nextToLastVertex;

  // Find the index of the first in for cur vertex
  inIndex = deBruijn[curVertex].FirstInIndex();
  assert(inIndex < 4);

  // now add this cycle
  readIndex = deBruijn[curVertex].read;
  readPos   = deBruijn[curVertex].pos;
  tuplePtr  = (char*) &(sequences[readIndex].seq[readPos]);
  cycleVertex = LocateTuple(sequences, deBruijn, tupleSize-1, (char*) tuplePtr);
  assert(cycleVertex >= 0);
	
  outIndex = deBruijn[curVertex].FirstOutIndex();

  ssize_t pathLength;
  pathLength = AdvanceSimplePath(cycleVertex, outIndex, deBruijn, lastVertex, nextToLastVertex);

  cycleEdges[curCycleEdge].seq.seq = (unsigned char*) new char[tupleSize - 1 + pathLength];
  cycleEdges[curCycleEdge].seq.length = tupleSize - 1 + pathLength;
  cycleEdges[curCycleEdge].length = cycleEdges[curCycleEdge].seq.length;
  cycleEdges[curCycleEdge].multiplicity = 1;
  cycleEdges[curCycleEdge].index = curCycleEdge;
  strncpy((char*) cycleEdges[curCycleEdge].seq.seq, 
					(const char*) &(sequences[readIndex].seq[readPos]), tupleSize -1) ;

  AdvanceSimplePath(cycleVertex, outIndex, deBruijn, lastVertex, nextToLastVertex, 
										&(cycleEdges[curCycleEdge].seq.seq[tupleSize-1]));

  std::cout << "cycle is of length " << cycleEdges[curCycleEdge].seq.length << std::endl;
  std::cout << std::endl;
	
  // this cycle is done. clear it so it is not traversed again
  ClearSimplePath(cycleVertex, outIndex, deBruijn);

  //UNUSED// ssize_t lastCycleVertex = cycleVertices.size()-1;
  cycleVertices[curCycleVertex].out[outIndex] = curCycleEdge;
  cycleEdges[curCycleEdge].dest = curCycleVertex;
  cycleEdges[curCycleEdge].src  = curCycleVertex;

  unsigned char inChar = sequences[readIndex].seq[readPos];
  cycleVertices[curCycleVertex].in[unmasked_nuc_index[inChar]] = curCycleEdge;
}

void FindCycles(VertexList &deBruijn, SimpleSequenceList &sequences, int tupleSize, ssize_t balanceGraph,
								BVertexList &cycleVertices, BEdgeList &cycleEdges) {
  ssize_t v;
  ssize_t curVertex, nextVertex;
  char *tuplePtr;
  char* rcTuplePtr;
  //UNUSED+// ssize_t  curCycleVertex;
  ssize_t curCycleEdge;
  for (v = 0; v <deBruijn.size(); v++ ) {
    if (deBruijn[v].OutDegree() == 1 and deBruijn[v].InDegree() == 1) {
      curVertex = v;
      nextVertex = deBruijn[curVertex].FirstOut();
      while (nextVertex != curVertex and deBruijn[nextVertex].OutDegree() == 1) {
				nextVertex = deBruijn[nextVertex].FirstOut();
      }
      if (nextVertex == curVertex) {
				//				std::cout << "Found a cycle starting/ending arunt: " << curVertex << std::endl;
				AdvanceCycle(curVertex, deBruijn, sequences, tupleSize, cycleVertices, cycleEdges);
				if (balanceGraph) {
					// Hopefully this is the only code that requires a specific check for balancing other
					// than the sanity checks.
					ssize_t balancedVertex;
					tuplePtr = (char*) &(sequences[deBruijn[curVertex].read].seq[deBruijn[curVertex].pos]);
					MakeRC((char*)tuplePtr, tupleSize-1, (unsigned char*&) rcTuplePtr); // TODO: fix compiler warning
					balancedVertex = LocateTuple(sequences, deBruijn, tupleSize-1, (char*) rcTuplePtr);
					delete [] rcTuplePtr;
					assert(balancedVertex >= 0);
					AdvanceCycle(balancedVertex, deBruijn, sequences, tupleSize, cycleVertices, cycleEdges);
					curCycleEdge = cycleEdges.size()-1;
					cycleEdges[curCycleEdge].balancedEdge = curCycleEdge-1;
					cycleEdges[curCycleEdge-1].balancedEdge = curCycleEdge;
				}
      }
    }
  }
}
		

void AppendDisjointGraph(BVertexList &vSetA, BEdgeList &eSetA,
												 BVertexList &vSetB, BEdgeList &eSetB) {
  ssize_t totalVertices = vSetA.size() + vSetB.size();
  ssize_t nAVertex = vSetA.size();
  ssize_t nAEdge   = eSetA.size();
  ssize_t totalEdges = eSetA.size() + eSetB.size();

  vSetA.resize(totalVertices);
  eSetA.resize(totalEdges);
  
  ssize_t v, e;
  // Update edge/vertex indices in set b
  for (v = 0; v < vSetB.size(); v++ ) {
    for (e = 0; e < 4; e++) {
      if (vSetB[v].out[e] >= 0) {
				vSetB[v].out[e] += nAEdge;
      }
      if (vSetB[v].in[e] >= 0) {
				vSetB[v].in[e] += nAEdge;
      }
    }
    vSetA[nAVertex + v] = vSetB[v];
  }
  for (e = 0; e < eSetB.size(); e++) { 
    eSetB[e].src += nAVertex;
    eSetB[e].dest += nAVertex;
    eSetB[e].balancedEdge += nAEdge;
    eSetA[nAEdge + e] = eSetB[e];
  }
}

void CreateBGraph(VertexList &deBruijn, 
									BVertexList &bVertices, 
									BEdgeList   &bEdges,
									SimpleSequenceList &sequences, int tupleSize) {

  ssize_t v, e;
  ssize_t bv;
  ssize_t readIndex;
  ssize_t readPos;
  char* tuplePtr;
  ssize_t pathLength;
  ssize_t lastVertex;
  ssize_t lastRead, lastPos;
  ssize_t lastBVertex;
  ssize_t edgeIndex;
  edgeIndex = 0;
  // Allocate the vertices
  AllocateBGraph(deBruijn, bVertices, bEdges, sequences, tupleSize);
  ssize_t nextToLastVertex;
  for (bv = 0; bv < bVertices.size(); bv++ ) {
    // First locate this vertex in the deBruijn graph.
    readIndex = bVertices[bv].read;
    readPos   = bVertices[bv].pos;
    tuplePtr  = (char*) &(sequences[readIndex].seq[readPos]);

    v = LocateTuple(sequences, deBruijn, tupleSize-1, (char*) tuplePtr);
    // Since we created this vertex from the other graph we must find it.
    assert(v >= 0);

    // Now look to see how long each edge is leaving this vertex
    for ( e = 0; e < 4; e++ ) {
      if ( deBruijn[v].out[e] >= 0 ) {
				// Traverse simple path to compute its length
				pathLength = AdvanceSimplePath(v, e, deBruijn, lastVertex, nextToLastVertex );
				bVertices[bv].out[e]         = edgeIndex;
				bEdges[edgeIndex].seq.seq    = (unsigned char*) new char[tupleSize - 1 + pathLength];
				bEdges[edgeIndex].seq.length = tupleSize - 1 + pathLength;
				/*
				std::cout << "path from " << bv << " out: " << e << " " << v 
									<<" " << lastVertex <<" ";
				PrintTuple(sequences, deBruijn[v], tupleSize-1);				
				std::cout << "  ...  ";
				PrintTuple(sequences, deBruijn[lastVertex], tupleSize-1);
				std::cout << std::endl;
				*/
				strncpy((char*) bEdges[edgeIndex].seq.seq, 
								(const char*) &(sequences[readIndex].seq[readPos]), tupleSize -1) ;

				// Now traverse the simple path again, but copy in the sequence
				AdvanceSimplePath(v, e, deBruijn, lastVertex, nextToLastVertex, 
													&(bEdges[edgeIndex].seq.seq[tupleSize-1]));
	
				// Make sure we don't traverse this path again (when finding cycles)
				ClearSimplePath(v, e, deBruijn);
				// Now find the last vertex in the branching graph
				lastRead = deBruijn[lastVertex].read;
				lastPos  = deBruijn[lastVertex].pos;
				tuplePtr = (char*) &(sequences[lastRead].seq[lastPos]);
				lastBVertex = LocateTuple(sequences, bVertices, tupleSize-1, tuplePtr);
				if (lastBVertex < 0) { 
					std::cout << "searching for tuple of len: " << tupleSize-1 << " failed " << std::endl;
					PrintTuple(sequences, deBruijn[lastVertex], tupleSize-1);
				}
				assert(lastBVertex >= 0);
				// store the destination
				bEdges[edgeIndex].dest = lastBVertex;

				// point the destination back at this vertex
				if (nextToLastVertex == -1) {
					// There was no next to last vertex, so there was a self-loop
					nextToLastVertex = lastVertex;
				}
				unsigned char inChar;
				ssize_t nextToLastRead = deBruijn[nextToLastVertex].read;
				ssize_t nextToLastPos  = deBruijn[nextToLastVertex].pos;
				inChar = sequences[nextToLastRead].seq[nextToLastPos];
				bEdges[edgeIndex].src = bv;
				bVertices[lastBVertex].in[unmasked_nuc_index[inChar]] = edgeIndex;
				bEdges[edgeIndex].length = bEdges[edgeIndex].seq.length;

				// Initialize unused fields
				bEdges[edgeIndex].multiplicity = 0;
				bEdges[edgeIndex].index = -1;
				bEdges[edgeIndex].flagged = 0;

				edgeIndex++;
      }
    }
  }
}


void CountSimpleSize(VertexList &vertices, ssize_t &numVertices, ssize_t &numEdges) {

  ssize_t v, e;
  numEdges = 0;
  numVertices = 0;
  for (v = 0; v < vertices.size(); v++ ) {

    // this vertex creates a new edge if there 
    // is a simple edge originating from it (indegree == 0)
    // or there are multiple edges exiting it (outdegree > 1)
    if (vertices[v].OutDegree() != 1 or 
				vertices[v].InDegree() != 1 ) {
      // this vertex will be in the simple graph if it has 
      // multiple edges going into it, or
      // one edge in and no out, or one out and no in (end of a path)
      numVertices++;
      for (e = 0; e < 4; e++ ) {
				if (vertices[v].out[e] >= 0)
					numEdges++;
      }
    }
  }
}


void BalanceGraph(BVertexList &vertices, BEdgeList &edges,
									SimpleSequenceList &sequences, int tupleSize) {
  // for each edge there MUST exist a reverse complement edge.

  ssize_t v, e;

  ssize_t balancedVertex;
  char *tuplePtr;
  unsigned char* rcTuplePtr;
  //UNUSED// ssize_t vertexRead, vertexPos;
  DNASequence rc;
  ssize_t edgeIndex;
  for (v = 0; v < vertices.size(); v++ ) {
    // Step 1, get the reverse complement of this 
    tuplePtr = (char*) &(sequences[vertices[v].read].seq[vertices[v].pos]);
    MakeRC(tuplePtr, tupleSize-1, rcTuplePtr);
    balancedVertex = LocateTuple(sequences, vertices, tupleSize-1, (char*) rcTuplePtr);
    delete [] rcTuplePtr;
    // The balanced vertex MUST exist
    if (balancedVertex < 0) {
      std::cout << "did not find vertex: " << rcTuplePtr << std::endl;
    }

    assert(balancedVertex >= 0);

    // Consider all edges that leave this vertex. Any edge that leaves
    // this vertex should have a balanced edge that enters the balanced 
    // vertex, with the complementary character to that of the char leaving this
    // vertex
    // So if the out is at 'a' (1), the in is at 't' (3)
    ssize_t balancedEdge;

    for (e = vertices[v].FirstOut();
				 e < vertices[v].EndOut();
				 e = vertices[v].NextOut(e)) {
      if (vertices[v].out[e] >= 0) {
				edgeIndex = vertices[v].out[e];
				// Find what slot the edge should be going into the balanced vertex
				unsigned char inChar;
				inChar = comp_bin[e];
				balancedEdge = vertices[balancedVertex].in[inChar];
				// Find the balanced edge going into the balanced vertex
				if (vertices[balancedVertex].in[inChar] == -1) {
					std::cout << "vertex: " << v << " ";
					PrintTuple(sequences, vertices[v], tupleSize-1);
					std::cout << " balanced " << balancedVertex << " ";
					PrintTuple(sequences, vertices[balancedVertex], tupleSize - 1);
					std::cout << " char: " << (ssize_t) inChar << std::endl;
					std::cout << vertices[v].read << " " << vertices[v].pos << std::endl;
				}

				assert(vertices[balancedVertex].in[inChar] != -1);
				assert (edges[balancedEdge].dest == balancedVertex);

				// Store the indices of the balanced edges
				edges[balancedEdge].balancedEdge = edgeIndex;
				edges[edgeIndex].balancedEdge = balancedEdge;
				/*
					std::cout << "edges " << balancedEdge << " and " 
					<< edgeIndex << " are balanced " << std::endl;
				*/
      }
    }
  }
}

