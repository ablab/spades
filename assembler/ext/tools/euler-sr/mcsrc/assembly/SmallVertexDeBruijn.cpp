/***************************************************************************
 * Title:          SmallVertexDeBruijn.cpp 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  02/28/2010
 *
 * Copyright (c) 2008-2010 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>

// sequence utilities
#include "SeqReader.h"
#include "SeqUtils.h"
#include "SimpleSequence.h"
#include "DNASequence.h"
#include "utils.h"
#include "IntegralTupleStatic.h"
#include "DeBruijnGraph.h"
#include "BEdge.h"



using namespace std;

class BranchingVertex : public CountedIntegralTuple {
public:
	ssize_t in[4];
	ssize_t out[4];
	unsigned char    outCount[4];
	static ssize_t length;
	

	// TODO: Can we get rid of the flags to save some memory?
	// HAS BITFIELD
	bool flagged:1;
	//	bool marked:1; // TODO: keep this or not?  It's only for GVZPrintBGraph(), which could be run as a separate program afterwards.



	BranchingVertex() {
		out[0] = out[1] = out[2] = out[3] = -1;
		in[0] = in[1] = in[2] = in[3] = -1;
		//		index = 0;
		//index = -1;
		flagged = 0; //TESTING//
	}

  ssize_t FirstIn() {
		ssize_t i = 0;
		while(i < 4 and in[i] == -1) i++;
		return i;
  }
	ssize_t EndIn() {
		return 4;
	}
	ssize_t NextIn(ssize_t i) {
		while(i < 4 and in[i] == -1) i++;
		return i; // TODO: check
	}

	ssize_t FirstOut() {
		//UNUSED// ssize_t i = 0;
		return NextOut(-1);
	}
	ssize_t EndOut() {
		return 4;
	}
	ssize_t NextOut(ssize_t i) {
		i++;
		while (i < 4 and out[i] == -1) i++;
		return i;
	}

  ssize_t InDegree() {
    ssize_t i, d;
    d = 0;
    for (i = 0; i < EndIn(); i++ ) {
      if (in[i] != -1) d++;
    }
    return d;
  }
  ssize_t OutDegree() {
    ssize_t i, d;
    d = 0;
    for (i = 0; i < EndOut(); i++ ) {
      if (out[i] != -1) d++;
    }
    return d;
  }

	bool operator<(const BranchingVertex &b) const{ 
		return cmp(b) < 0;
		//		return tuple < b.tuple;
	}
	bool operator<(const CountedIntegralTuple &t) const{
		return cmp(t) < 0;
		//		return tuple < t.tuple;
	}

	bool operator>(const BranchingVertex &b) const{
		return cmp(b) > 0;
		//		return tuple > b.tuple;
	}
	bool operator==(const BranchingVertex &b) const{
		return cmp(b) == 0;
		//		return tuple == b.tuple;
	}
	bool operator!=(const BranchingVertex &b) const {
		return cmp(b) != 0;
		//		return tuple != b.tuple;
	}
	BranchingVertex &operator=(const BranchingVertex &b) {
		if (this != &b) {
			//			tuple = b.tuple;
			CopyTuple(b);
			memcpy(in, b.in, sizeof(ssize_t)*4);
			memcpy(out, b.out, sizeof(ssize_t)*4);
		}
		return *this;
	}
	void Write(ostream &outs) {
		outs << -1 << " "; // index is not yet defined, so use -1
		//		outs << index << " ";
		//UNUSED// ssize_t o, i;

		outs << out[0] << " ";
		outs << out[1] << " ";
		outs << out[2] << " ";
		outs << out[3] << " ";
		outs << in[0] << " ";
		outs << in[1] << " ";
		outs << in[2] << " ";
		outs << in[3];
	}

	//	void Unmark() {
	//		marked = GraphVertex::NotMarked;
	//	}
};

ssize_t BranchingVertex::length = 0;

// TODO: This is overkill.
// If trimEnd==0, the current code in this file doesn't use coverage at all.
// If trimEnd>0, it only checks coverage <2 vs. >=2.
#define COVERAGE_BITS 20

class Adjacency {
public:
	// HAS BITFIELD
	_UINT_ coverage:COVERAGE_BITS; // never more than 1M coverage makes sense
	bool destG:1;
	bool destA:1;
	bool destC:1;
	bool destT:1;
	bool srcG:1;
	bool srcA:1;
	bool srcC:1;
	bool srcT:1;

	unsigned char destCoverage[4];

	void Clear() {
		destG = destA = destC = destT = 0;
		srcG = srcA = srcC = srcT = 0;
		destCoverage[0] = destCoverage[1] = destCoverage[2] = destCoverage[3] = 0;
	}

	Adjacency() {
		Clear();
	}

	void IncrementDest(ssize_t index) {
		if (destCoverage[index] == 255) return;
		else destCoverage[index]++;
	}
			
	ssize_t DegreeOneOutIndex() {
		if (destA) return 0;
		if (destC) return 1;
		if (destG) return 2;
		if (destT) return 3;
		// Make sure this isn't called on degree 0 vertices
		assert(0);
		return -1;
	}
	
	ssize_t FirstOut() {
		if (destA) return 0;
		if (destC) return 1;
		if (destG) return 2;
		if (destT) return 3;
		return 4; // TODO: check
	}

	ssize_t EndOut() {
		return 4;
	}

	ssize_t NextOut(ssize_t i) {
		// reached the end of the list
		if (i == 3) return 4;
		if (i == 2 and destT) return 3;
		if (i == 1 and destG) return 2;
		if (i == 1 and destT) return 3;
		if (i == 0 and destC) return 1;
		if (i == 0 and destG) return 2;
		if (i == 0 and destT) return 3;
		return 4;
	}
	
	ssize_t DegreeOneInIndex() {
		if (srcA) return 0;
		if (srcC) return 1;
		if (srcG) return 2;
		if (srcT) return 3;
		// Make sure this isn't called on degree 0 vertices
		assert(0);
		return -1;
	}

	ssize_t GetSrc(ssize_t index) {
		switch (index) {
		case (0):
			return srcA;
		case (1):
			return srcC;
		case (2):
			return srcG;
		case (3):
			return srcT;
		default :
			std::cout << "error, cannot index more than 4 nucleotides" << std::endl;
			assert(0);
		}
		return 0; // Never reaches here. Quiet compiler warnings
	}

	// TODO: change types to int or unsigned char
	ssize_t SetSrc(ssize_t index, ssize_t value) {
		switch (index) {
		case (0):
			return srcA = value;
		case (1):
			return srcC = value;
		case (2):
			return srcG = value;
		case (3):
			return srcT = value;
		default :
			std::cout << "error, cannot index more than 4 nucleotides" << std::endl;
			assert(0);
		}
		return 0; // Never reaches here. Quiet compiler warnings
	}

	ssize_t GetDest(ssize_t index) {
		switch (index) {
		case (0):
			return destA;
		case (1):
			return destC;
		case (2):
			return destG;
		case (3):
			return destT;
		default :
			std::cout << "error, cannot index more than 4 nucleotides" << std::endl;
			assert(0);
		}
		return 0; // Never reaches here. Quiet compiler warnings
	}

	ssize_t SetDest(ssize_t index, ssize_t value) {
		switch (index) {
		case (0):
			return destA = value;
		case (1):
			return destC = value;
		case (2):
			return destG = value;
		case (3):
			return destT = value;
		default :
			std::cout << "error, cannot index more than 4 nucleotides" << std::endl;
			assert(0);
			return 0; // Never reaches here. Quiet compiler warnings
		}
	}

	ssize_t OutDegree() {
		return destG + destA + destC + destT;
	}

	ssize_t InDegree() {
		return srcG + srcA + srcC + srcT;
	}

	ssize_t IsBranch() {
		return !(InDegree() == 1 and OutDegree() == 1);
	}

};

typedef std::vector<BranchingVertex> BVertexList;

void AppendDisjointGraph(BVertexList &vSetA, BEdgeList &eSetA,
												 BVertexList &vSetB, BEdgeList &eSetB);

void FindCycles(CountedIntegralTuple *deBruijn, Adjacency *dbAdj, ssize_t numVertices,
								ssize_t balanceGraph,
								BVertexList &cycleVertices, BEdgeList &cycleEdges);

void AllocateBGraph(CountedIntegralTuple *deBruijn, Adjacency *deBruijnAdj, ssize_t numVertices,
										BVertexList &bVertices,	BEdgeList &bEdges);

void CreateBGraph(CountedIntegralTuple *deBruijn, Adjacency * deBruijnAdj, ssize_t numVertices,
									BVertexList &bVertices, BEdgeList   &bEdges);

void BalanceGraph(BVertexList &vertices, BEdgeList &edges);

void CreateGraph(std::string seqFileName,
								 CountedIntegralTuple *vertices, Adjacency* adj, 
								 ssize_t numVertices, int tupleSize, ssize_t allowGapped,
								 std::ostream &report
								 );

ssize_t RemoveLowCoverageEdges(CountedIntegralTuple *deBruijn,
															 Adjacency *deBruijnAdj,
															 ssize_t numVertices,
															 int minCoverage);
/*
	void GVZPrintGraph(VertexList &graph, 
	SimpleSequenceList &sequences, int tupleSize,std::string &graphOutName );
*/

void CountSimpleSize(CountedIntegralTuple *vertices, Adjacency *adj, ssize_t numVertices,
										 ssize_t &numBVertices, ssize_t &numBEdges);

ssize_t AdvanceSimplePath(ssize_t curVertexIndex, ssize_t outIndex, 
											CountedIntegralTuple *vertices, Adjacency *adj, ssize_t numVertices,
											ssize_t &lastVertexIndex, unsigned char &srcIndex,	unsigned char* seq );

void EraseSimplePath(ssize_t curVertexIndex, ssize_t firstEdge, ssize_t destVertexIndex,
										 CountedIntegralTuple *vertices, Adjacency *adj, ssize_t numVertices);

ssize_t RemoveZeroDegreeDBVertices(CountedIntegralTuple *deBruijn, Adjacency *deBruijnAdj, 
															 ssize_t numVertices);

ssize_t LookupBVertex(BVertexList &vertices, CountedIntegralTuple &vertex);
ssize_t LookupVertex(CountedIntegralTuple *vertices, ssize_t numVertices, CountedIntegralTuple &vertex);
void PrintSmallEdges(CountedIntegralTuple *deBruijn, Adjacency *deBruijnAdj, ssize_t numVertices);
void TrimSmallEdges(CountedIntegralTuple *deBruijn, Adjacency *deBruijnAdj, ssize_t numVertices, ssize_t trimEnd);

void PrintUsage() {
  std::cout << std::endl 
						<< "Build a de Bruijn graph using a set of vertices and overlaps defined" << std::endl
						<< "by a set of sequences." << std::endl << std::endl;
  std::cout << "usage: svdeBruijn seqFile vertexSize [options]" 
						<< std::endl
						<< "  [-noCondense]         Do not condense simple paths (in degree=out degree = 1)" << std::endl
						<< "  [-edgeOutFile filename] Output edges to filename." << std::endl
						<< "  [-intvFile filename]  Output read intervals to filename." << std::endl
						<< "  [-noGraph]            Suppress printing of .graph file (default is on)." << std::endl
						<< "  [-noEdge]             Suppress printing of .edge file (default is on)." << std::endl
						<< "  [-bGraph filename]    Print bgraph to filename" << std::endl
						<< "  [-trimShort N]        Trim short edges (length < N) from graph." << std::endl
						<< "  [-noAllowGapped]      Enforce that all k-1 mers from sequences have a vertex" << std::endl
						<< "                        in the graph.  Default is to allow gapped." << std::endl
						<< "  [-minCov c]           Remove edges from de Bruijn graph if they appear fewer" << std::endl
						<< "                        than 'c' times in reads.  c must be <= 255." << std::endl
						<< std::endl << std::endl;
	// TODO: Investigate why "-intvFile" and "-noAllowGapped" are not implemented,
	// and why "-graphFile" is implemented but not listed here
}


ssize_t VERTEX_SIZE;




int main(int argc, char*argv[]) {
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
	ssize_t allowGapped   = 1;
	int minCoverage = 0;
  int argi = 1;
  if (argc < 2) {
    PrintUsage();
    exit(1);
  }
  sequencesFileName    = argv[argi++];
	vertexSize = atoi(argv[argi++]);
	vertexListFileName   = sequencesFileName + ".spect";
  gvzGraphOutName      = sequencesFileName + ".dot";
  graphFileName        = sequencesFileName + ".graph";
  edgeOutputFileName   = sequencesFileName + ".edge";
  bgraphOutName        = sequencesFileName + ".bgraph";
	reportFileName       = FormReportName(sequencesFileName);
	ssize_t trimShort = 0;
	tupleSize = vertexSize + 1;
  while(argi < argc) {
    if (strcmp(argv[argi], "-bGraph") == 0) {
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
		else if (strcmp(argv[argi], "-trimShort") == 0) {
			trimShort = atoi(argv[++argi]);
		}
		else if (strcmp(argv[argi], "-minCov") == 0) {
			minCoverage = atoi(argv[++argi]);
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
	CountedIntegralTuple::SetTupleSize(vertexSize);
	BranchingVertex::length = vertexSize;

	// The following convention is used:
	//   Sequences are stored as integers so that sequence[i] = bit[i].
	//   That way, we shift left to move backwards in a sequence, 
	//   and shift right to move forwards. 

	// Set up the masks for transforming in/out vertices
	// These could be hard-coded for more effort.
	//UNUSED// ssize_t i;

  DNASequence seq, seqRC;


  std::ifstream seqIn;

  SimpleSequence simpleSeq;
  DNASequence newSeq;

  //UNUSED// ssize_t seqMemory;

  // now create a list of l-1 mers
  
  // The list of vertices should have been previously computed
  std::ifstream vertexIn, edgeIn;

  openck(vertexListFileName, vertexIn, std::ios::in | std::ios::binary, report);
	CountedIntegralTuple *dbVertices;
	ssize_t numVertices;

	// Read the list of vertices

	// TODO: tupleSize was added to the file format late; reconcile with command-line
	ssize_t tupleSize_SSZT;
	vertexIn.read((char*) &tupleSize_SSZT, sizeof(tupleSize_SSZT));
	if (vertexSize != tupleSize_SSZT) {
		std::cerr << "Warning: command-line parameter vertexSize=" << vertexSize
							<< " disagrees with binary file tupleSize=" << tupleSize_SSZT << std::endl;
	}
	CountedIntegralTuple::SetTupleSize(vertexSize);


	vertexIn.read((char*) &numVertices, sizeof(ssize_t));
	dbVertices = new CountedIntegralTuple[numVertices];
	vertexIn.read((char*) dbVertices, sizeof(CountedIntegralTuple)*numVertices);
  vertexIn.close();

	Adjacency *dbAdj = new Adjacency[numVertices];

  // Build the graph for every vertex
  std::cout << "Creating the graph." << std::endl;
	CreateGraph(sequencesFileName, 
							dbVertices, dbAdj, numVertices, IntegralTuple::tupleSize+1, allowGapped,
							report);

	if (minCoverage > 0) {
		ssize_t numRemoved;
		numRemoved = RemoveLowCoverageEdges(dbVertices, dbAdj, numVertices, minCoverage);
		cout << "removed " << numRemoved << " low coverage edges." << endl;
	}

  std::cout << "The graph has " << numVertices << " vertices (" 
						<< sizeof(CountedIntegralTuple) * numVertices << " bytes)." << std::endl;

	numVertices =	RemoveZeroDegreeDBVertices(dbVertices, dbAdj, numVertices);
	cout << "and after zero degree removal: " << numVertices << endl;
  BVertexList bVertices;
  BEdgeList bEdges;
	if (trimShort) {
		TrimSmallEdges(dbVertices, dbAdj, numVertices, trimShort);
		numVertices = RemoveZeroDegreeDBVertices(dbVertices,dbAdj, numVertices);
	}
  if (condensePaths) {
    std::cout << "Creating the branching bgraph." << std::endl;

		CreateBGraph(dbVertices, dbAdj, numVertices, 
								 bVertices, bEdges);

    BalanceGraph(bVertices, bEdges);
    
    // Cycles are non-branching, cycles in the graph that contain more
    // than one vertex but have no degree-0 in or out vertices.  These
    // are not picked up as branching vertices. They can probably be
    // removed from most assemblies, but they are included here just
    // in case.

    BVertexList cycleVertices;
    BEdgeList cycleEdges;
    ssize_t balanceGraph = 1;
    FindCycles(dbVertices, dbAdj, numVertices, 
							 balanceGraph, cycleVertices, cycleEdges );
    AppendDisjointGraph(bVertices, bEdges, cycleVertices, cycleEdges);
  }

  ssize_t nRemoved;
  nRemoved = RemoveZeroDegreeVertices(bVertices, bEdges);


  if (! condensePaths ) {
		//    GVZPrintGraph(vertices, sequences, tupleSize, gvzGraphOutName);
  }
  else {
		//    GVZPrintBGraph(bVertices, bEdges, gvzGraphOutName, report);
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



void CreateGraph(std::string seqFileName,
								 CountedIntegralTuple *vertices, Adjacency* adj, 
								 ssize_t numVertices, int tupleSize, ssize_t allowGapped,
								 std::ostream &report
								 ) {
  ssize_t v;
  ssize_t p;
  //UNUSED// ssize_t s;
  //UNUSED// unsigned char* tupleSeq;
  ssize_t toVertexIndex, fromVertexIndex;
  unsigned char outChar, inChar;

  int vertexSize = tupleSize - 1;

  // Initialize the vertex list
  v = 0;
  p = 0;
  ssize_t nextN = -1;
  ssize_t n;
	DNASequence read, readRC;
	std::ifstream readsIn;
	openck(seqFileName, readsIn, std::ios::in, report);
	DNASequence reads[2];
	CountedIntegralTuple tuple;
	CountedIntegralTuple *vertexPtr;
	ssize_t readIndex = 0;
	while(SeqReader::GetSeq(readsIn, read, SeqReader::noConvert)) {
		readIndex++;
		MakeRC(read, readRC);
		reads[0].seq = read.seq;
		reads[0].length = read.length;
		reads[1].seq = readRC.seq;
		reads[1].length = read.length;

		// Check to see if there is a masked nucleotide in this read.
		// If so, completely discard it.
		ssize_t readHasMask = 0;
		for (n = 0; n < read.length; n++ ) {
			if (unmasked_nuc_index[read.seq[n]] >= 4) {
				readHasMask = 1;
				break;
			}
		}
		if (readHasMask) {
			continue;
		}

		ssize_t strand;
		for (strand = 0; strand < 2; strand++) {
			// assume no vertices found.
			fromVertexIndex = -1;
			toVertexIndex   = -1;

			// Make sure there is enough sequence in this read to 
			// create an edge
			if (reads[strand].length < vertexSize  + 1)
				continue;
			

			// Make sure the entire first two tuples are not masked
			nextN = -1;
			for (n = 0; n < vertexSize + 1; n++ ) {
				if (unmasked_nuc_index[reads[strand].seq[n]] >= 4) {
					nextN = n;
					break;
				}
			}

			// Create edges for every pair of adjacent tuples in the read.
			ssize_t p;
			
			for (p = 0; p < reads[strand].length - vertexSize + 1; p++) {

				// the next tuple has a masked character, so it is not valid.
				if (unmasked_nuc_index[reads[strand].seq[p + vertexSize  - 1]] >= 4) {
					nextN = p + vertexSize - 1;
					// Do not save the last vertex, since it is skipped by an 'N'
					fromVertexIndex = -1;
					toVertexIndex   = -1;
				}

				// If this vertex does not contain any masked positions.
				if (nextN < p) {
					
					// The first time entering this loop, only the from
					// vertex will be set and no edges are created.
					// The next time through the loop, if the next tuple is valid
					// the to vertex is found and an edge is created.
					if (fromVertexIndex == -1) {
						tuple.StringToTuple(&(reads[strand].seq[p]));
						vertexPtr = std::lower_bound(vertices, vertices + numVertices, tuple);
						
						if (vertexPtr == vertices+ numVertices or 
								*vertexPtr != tuple) {
							if (allowGapped) {
								fromVertexIndex = -1;
								//								cout << "skipping gapped vertex. " << endl;
							}
							else {
								std::cout << "Error, a read contains a vertex that is not represented" << std::endl;
								std::cout << "in the vertex list." << std::endl;
								read.PrintlnSeq(std::cout);
								assert(0);
							}
						}		
						else {
							fromVertexIndex = vertexPtr - vertices;
						}
					}
					else {
						// A from vertex is already assigned, find the 
						// to vertex, and create the edge
						tuple.StringToTuple(&(reads[strand].seq[p]));
						vertexPtr = std::lower_bound(vertices, vertices + numVertices, tuple);
						if (vertexPtr == vertices+numVertices or 
								*vertexPtr != tuple) {
							if (allowGapped) {
								toVertexIndex = -1;
								fromVertexIndex = -1;
							}
							else {
								std::cout << "Error, a read contains a vertex that is not represented" << std::endl;
								std::cout << "in the vertex list." << std::endl;
								assert(0);
							}
						}	
						else {
							toVertexIndex   = vertexPtr - vertices;
							outChar  = reads[strand].seq[p + vertexSize - 1];
							inChar   = reads[strand].seq[p-1];
							
							if (unmasked_nuc_index[outChar] >= 4 or 
									unmasked_nuc_index[inChar] >= 4) {
								std::cout << "Part of a read is masked, but shouldn't be " << p << std::endl;
								reads[strand].PrintlnSeq(std::cout);
								assert(unmasked_nuc_index[outChar] < 4);
								assert(unmasked_nuc_index[inChar] < 4);
							}
						
							adj[fromVertexIndex].SetDest(unmasked_nuc_index[outChar], 1);
							adj[toVertexIndex].SetSrc(unmasked_nuc_index[inChar], 1);
							adj[fromVertexIndex].IncrementDest(unmasked_nuc_index[outChar]);

							adj[fromVertexIndex].coverage++; // not used for now.
							if (adj[fromVertexIndex].coverage == 0) {
								// Deal w/overflow
								// TODO: We really only use coverage = 0, 1, or >=2, so could optimize this.
								adj[fromVertexIndex].coverage--;
							}

							fromVertexIndex = toVertexIndex;
						}
					}
				}
			}
		}
		read.Reset();
		readRC.Reset();
  }
}

ssize_t LookupBVertex(BVertexList &vertices, CountedIntegralTuple &bvertex) {
	BVertexList::iterator vertexIt;
	vertexIt = std::lower_bound(vertices.begin(), vertices.end(), bvertex );
	if (vertexIt == vertices.end() || 
			bvertex.cmp(*vertexIt) != 0 // compare tuples
			//			IntegralTuple::cmp((*vertexIt).tuple, bvertex.tuple) != 0 // compare tuples
			//			(*vertexIt).tuple != bvertex.tuple
			) {
		std::cout << "ERROR, This should only look for existing branching vertices."<< std::endl;
		assert(0);
	}
	return vertexIt - vertices.begin();
}

ssize_t LookupVertex(CountedIntegralTuple *vertices, ssize_t numVertices, CountedIntegralTuple &vertex) {
	CountedIntegralTuple *vertexPtr;
	CountedIntegralTuple query;
	vertexPtr = std::lower_bound(vertices, vertices + numVertices, vertex);
	if (vertexPtr == vertices + numVertices || 
			*vertexPtr != vertex // compare tuples
			//			vertexPtr->tuple != vertex.tuple
			) {
		std::cout << "ERROR. Every vertex considered must be present in the vertex list." << std::endl;
		//		assert(0);
		return -1;
	}
	return vertexPtr - vertices;
}

ssize_t GetDestVertexIndex(CountedIntegralTuple* vertices, ssize_t numVertices,
											 ssize_t vertexIndex, char outIndex) {
	CountedIntegralTuple vertex = vertices[vertexIndex];
	CountedIntegralTuple destVertex;

	ForwardNuc(vertex, outIndex, destVertex);
	return LookupVertex(vertices, numVertices, destVertex);
}

ssize_t GetSrcVertexIndex(CountedIntegralTuple *vertices, ssize_t numVertices,
											ssize_t vertexIndex, char inIndex) {
	CountedIntegralTuple vertex = vertices[vertexIndex];
	CountedIntegralTuple srcVertex;

	BackwardsNuc(vertex, inIndex, srcVertex);
	return LookupVertex(vertices, numVertices, srcVertex);
}


	
void AllocateBGraph(CountedIntegralTuple *deBruijn, Adjacency *deBruijnAdj, ssize_t numVertices,
										BVertexList &bVertices,	BEdgeList &bEdges) {

  ssize_t numBVertices, numBEdges;

	// Count the number of branching vertices and edges.
  CountSimpleSize(deBruijn, deBruijnAdj, numVertices,
									numBVertices, numBEdges);

	// Create the more heavy-weight branching vertices and edges

	// These contain edge indices.
  bVertices.resize(numBVertices);

	// The edges here contain sequences.
  bEdges.resize(numBEdges);
  
  //UNUSED+// ssize_t e;
  ssize_t v ;
  ssize_t bv = 0;
  for (v = 0; v < numVertices; v++ ) {
    if (deBruijnAdj[v].IsBranch()) {
			//			bVertices[bv].tuple = deBruijn[v].tuple;
			bVertices[bv].CopyTuple(deBruijn[v]);
      bv++;
    }
  }

  // Make the list of vertices binary-searchable
	// But this should already be in order since the de Bruijn list is 
	// in order.
  std::sort(bVertices.begin(), bVertices.end());
}


void FindCycles(CountedIntegralTuple *deBruijn, Adjacency *deBruijnAdj, ssize_t numVertices,
								ssize_t balanceGraph,
								BVertexList &cycleVertices, BEdgeList &cycleEdges) {
  ssize_t v;
  //UNUSED// char *tuplePtr;
	//	char *rcTuplePtr;
  ssize_t curCycleVertex, curCycleEdge;
	/*
		At this point in time all simplepaths in the de Bruijn graph 
		should have been removed from the graph.  The only things that are
		left are rings in the graph.  

		Each ring should be replaced by a single vertex and edge that is a loop. 
	*/
	unsigned char srcIndex, destIndex, inNuc;
	CountedIntegralTuple curVertex, nextVertex;
	//UNUSED+// ssize_t nextVertexIndex;
	ssize_t curVertexIndex ;
	ssize_t pathLength;
	//UNUSED// unsigned char *cycleSeq, *cycleSeqRC;
  for (v = 0; v < numVertices; v++ ) {
    if (deBruijnAdj[v].OutDegree() == 1 and deBruijnAdj[v].InDegree() == 1) {
			// Found a ring.
      curVertexIndex = v;


			// Break open the ring so that the traversal ends.
			inNuc  = deBruijnAdj[curVertexIndex].DegreeOneInIndex();
			deBruijnAdj[curVertexIndex].SetSrc(inNuc, 0);

			// 
      destIndex = deBruijnAdj[curVertexIndex].DegreeOneOutIndex();
			curVertex = deBruijn[curVertexIndex];

			// Allocate the new vertices for the graph.
			cycleVertices.push_back(BranchingVertex());
			cycleEdges.push_back(BEdge());

			ssize_t lastVertexIndex;

			pathLength =  AdvanceSimplePath(curVertexIndex, destIndex,
																			deBruijn, deBruijnAdj, numVertices,
																			lastVertexIndex, srcIndex, 
																			// just count the length this time
																			(unsigned char*) NULL); 
			curCycleEdge = cycleEdges.size() - 1;
			cycleEdges[curCycleEdge].seq.seq = new unsigned char[pathLength + CountedIntegralTuple::tupleSize];
			cycleEdges[curCycleEdge].seq.length = pathLength + CountedIntegralTuple::tupleSize;

			std::string curVertexString;
			deBruijn[curVertexIndex].ToString(curVertexString);
			memcpy(cycleEdges[curCycleEdge].seq.seq, 
						 curVertexString.c_str(), curVertexString.size());
			AdvanceSimplePath(curVertexIndex, destIndex,
												deBruijn, deBruijnAdj, numVertices,
												lastVertexIndex, srcIndex, 
												// just count the length this time
									&cycleEdges[curCycleEdge].seq.seq[CountedIntegralTuple::tupleSize]);  

			// Link in the cycle.
			curCycleVertex = cycleVertices.size() - 1;
			curCycleEdge   = cycleEdges.size() - 1;
			cycleVertices[curCycleVertex].out[destIndex] = curCycleEdge;
			cycleVertices[curCycleVertex].in[srcIndex]   = curCycleEdge;
			cycleEdges[curCycleEdge].src = curCycleVertex;
			cycleEdges[curCycleEdge].dest = curCycleVertex;

			//			AdvanceCycle(curVertex, deBruijn, sequences, tupleSize, cycleVertices, cycleEdges);
			if (balanceGraph) {
				// Hopefully this is the only code that requires a specific check for balancing other
				// than the sanity checks.
				ssize_t balancedVertexIndex;
				CountedIntegralTuple balancedVertex;
				deBruijn[curVertexIndex].MakeRC(balancedVertex);
				balancedVertexIndex = LookupVertex(deBruijn, numVertices, balancedVertex);
				assert(balancedVertexIndex >= 0);

				if ((balancedVertexIndex == curVertexIndex) or
						(deBruijnAdj[balancedVertexIndex].OutDegree() == 0 and
						 deBruijnAdj[balancedVertexIndex].InDegree() == 1) ) {
					// If the vertex is it's own RC, or the RC is along the 
					// path of this vertex, that will create problems, because the
					// RC will not map to the mirror position as the forward strand.
					// The solution is to just not store the sequence.
					
					cycleEdges.pop_back();
					cycleVertices.pop_back();
				}
				else {
					cycleVertices.push_back(BranchingVertex());
					cycleEdges.push_back(BEdge());
					// Balance the graph.
					curCycleEdge = cycleEdges.size()-1;
					curCycleVertex = cycleVertices.size() - 1;
					
					cycleEdges[curCycleEdge].seq.seq = 
						new unsigned char[pathLength + CountedIntegralTuple::tupleSize];
					cycleEdges[curCycleEdge].seq.length = pathLength + CountedIntegralTuple::tupleSize;
					

					destIndex  = deBruijnAdj[balancedVertexIndex].FirstOut();
					inNuc = deBruijnAdj[balancedVertexIndex].DegreeOneInIndex();
					deBruijnAdj[balancedVertexIndex].SetSrc(inNuc, 0);

					pathLength =  AdvanceSimplePath(balancedVertexIndex, destIndex,
																					deBruijn, deBruijnAdj, numVertices,
																					lastVertexIndex, srcIndex, 
																					// just count the length this time
																					&cycleEdges[curCycleEdge].seq.seq[CountedIntegralTuple::tupleSize]); 
				
					std::string balancedSeq;
					deBruijn[balancedVertexIndex].ToString(balancedSeq);
					memcpy(cycleEdges[curCycleEdge].seq.seq, 
								 balancedSeq.c_str(), balancedSeq.size());

				
					cycleVertices[curCycleVertex].out[destIndex] = curCycleEdge;
					cycleVertices[curCycleVertex].in[srcIndex]   = curCycleEdge;
					cycleEdges[curCycleEdge].src = curCycleVertex;
					cycleEdges[curCycleEdge].dest = curCycleVertex;

					
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


void CreateBGraph(CountedIntegralTuple *deBruijn, Adjacency * deBruijnAdj, ssize_t numVertices,
									BVertexList &bVertices, BEdgeList   &bEdges) {

  ssize_t v, e;
  ssize_t bv;
  //UNUSED// char* tuplePtr;
  ssize_t pathLength;
  ssize_t edgeIndex = 0;

  // Allocate the vertices
  AllocateBGraph(deBruijn, deBruijnAdj, numVertices, bVertices, bEdges);
  //UNUSED// ssize_t nextToLastVertex;
	CountedIntegralTuple tuple;
	ssize_t lastVertexIndex;
	//UNUSED+// unsigned char destIndex;
	unsigned char srcIndex ;
  for (bv = 0; bv < bVertices.size(); bv++ ) {
    // First locate this vertex in the deBruijn graph.
		//    tuple.tuple = bVertices[bv].tuple;
		tuple.CopyTuple(bVertices[bv]);
		v = LookupVertex(deBruijn, numVertices, tuple);
		
    // Since we created this vertex from the other graph we must find it.
    assert(v >= 0);

    // Now look to see how long each edge is leaving this vertex
    for ( e = 0; e < 4; e++ ) {
      if ( deBruijnAdj[v].GetDest(e) > 0 ) {
				// Traverse simple path to count its length.
				pathLength =  AdvanceSimplePath(v,e,
																				deBruijn, deBruijnAdj, numVertices,
																				lastVertexIndex, srcIndex, 
																				// just count the length this time
																				(unsigned char*) NULL); 

				bVertices[bv].out[e]         = edgeIndex;
				bEdges[edgeIndex].seq.seq    = (unsigned char*) new char[CountedIntegralTuple::tupleSize + pathLength];
				bEdges[edgeIndex].seq.length = CountedIntegralTuple::tupleSize + pathLength;

				// Initialize the first tuple sequence of the edge
				CountedIntegralTuple bVertex;
				//				bVertex.tuple = bVertices[bv].tuple;
				bVertex.CopyTuple(bVertices[bv]);
				std::string   bVertexSeq;
				bVertex.ToString(bVertexSeq);
				strncpy((char*) bEdges[edgeIndex].seq.seq, 
								(const char*) bVertexSeq.c_str(), CountedIntegralTuple::tupleSize);

				// Now traverse the simple path again, but copy in the sequence
				pathLength =  AdvanceSimplePath(v,e,
																				deBruijn, deBruijnAdj, numVertices, lastVertexIndex, srcIndex, 
																				// just count the length this time
																				&(bEdges[edgeIndex].seq.seq[CountedIntegralTuple::tupleSize]));

				// Make sure we don't traverse this path again (when finding cycles)
				//				ClearSimplePath(v, e, deBruijn);

				// Now find the last vertex in the branching graph
				ssize_t destBranchingIndex;
				destBranchingIndex = LookupBVertex(bVertices, deBruijn[lastVertexIndex]);
				
				// Connect the edge
				bEdges[edgeIndex].dest = destBranchingIndex;
				bEdges[edgeIndex].src  = bv;

				// Connect the dest vertex to the edge
				bVertices[destBranchingIndex].in[srcIndex] = edgeIndex;
				bEdges[edgeIndex].length = bEdges[edgeIndex].seq.length;

				// Initialize unused fields
				bEdges[edgeIndex].multiplicity = 0;
				bEdges[edgeIndex].index = -1;
				bEdges[edgeIndex].flagged = 0; //TESTING//

				edgeIndex++;
      }
    }
  }
}


void CountSimpleSize(CountedIntegralTuple *vertices, Adjacency *adj, ssize_t numVertices,
										 ssize_t &numBVertices, ssize_t &numBEdges) {

  ssize_t v, e;
  numBEdges = 0;
  numBVertices = 0;
  for (v = 0; v < numVertices; v++ ) {
    // this vertex creates a new edge if there 
    // is a simple edge originating from it (indegree == 0)
    // or there are multiple edges exiting it (outdegree > 1)
		
    if (adj[v].OutDegree() != 1 or 
				adj[v].InDegree()  != 1 ) {
      // this vertex will be in the simple graph if it has 
      // multiple edges going into it, or
      // one edge in and no out, or one out and no in (end of a path)
      numBVertices++;
      for (e = 0; e < 4; e++ ) {
				if (adj[v].GetDest(e) != 0)
					numBEdges++;
      }
    }
  }
}

void BalanceGraph(BVertexList &vertices, BEdgeList &edges) {
  // for each edge there MUST exist a reverse complement edge.

  ssize_t v, e;

  ssize_t balancedVertex;
  //UNUSED// char *tuplePtr;
	//  unsigned char* rcTuplePtr;
  //UNUSED// ssize_t vertexRead, vertexPos;
  DNASequence rc;
  ssize_t edgeIndex;
	std::string tupleSeq;
	CountedIntegralTuple vertex, vertexRC;
	//UNUSED// ssize_t rcVertexIndex;
  for (v = 0; v < vertices.size(); v++ ) {
    // Step 1, get the reverse complement of this 
		//		vertex.tuple = vertices[v].tuple;
		vertex.CopyTuple(vertices[v]);
		vertex.MakeRC(vertexRC);
		balancedVertex = LookupBVertex(vertices, vertexRC);
		assert(balancedVertex >= -1);
    // The balanced vertex MUST exist
    if (balancedVertex < 0) {
			//      std::cout << "did not find vertex: " << rcTuplePtr << std::endl;
      std::cout << "did not find vertex: " << vertexRC << std::endl;
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
					std::string tupleString, balancedString;
					vertices[v].ToString(tupleString);
					vertices[balancedVertex].ToString(balancedString);
					std::cout << tupleString << " balanced " << balancedString << std::endl;
					std::cout << " char: " << (ssize_t) inChar << std::endl;
				}

				assert(vertices[balancedVertex].in[inChar] != -1);
				assert(edges[balancedEdge].dest == balancedVertex);

				// Store the indices of the balanced edges
				edges[balancedEdge].balancedEdge = edgeIndex;
				edges[edgeIndex].balancedEdge = balancedEdge;
      }
    }
  }
}

void EraseSimplePath(ssize_t curVertexIndex, ssize_t outIndex, ssize_t destVertexIndex,
										 CountedIntegralTuple *vertices, Adjacency *adj, ssize_t numVertices) {
	/*	string srcstr, deststr;
	vertices[curVertexIndex].ToString(srcstr);
	vertices[destVertexIndex].ToString(deststr);

		cout << "erasing from " << curVertexIndex << " (" << vertices[curVertexIndex].tuple  << ") to " << destVertexIndex 
			 << " ( " << vertices[destVertexIndex].tuple << ")" << endl;
	cout << srcstr << " "  << deststr << endl;
	*/
	ssize_t firstVertexIndex = curVertexIndex;
	CountedIntegralTuple curVertex;
	CountedIntegralTuple nextVertex;
	string tupstr;
	// clear the first edge
	adj[curVertexIndex].SetDest(outIndex, 0);
	
	curVertex   = vertices[curVertexIndex];
	vertices[curVertexIndex].ToString(tupstr);
	//	cout << "start: " << tupstr << " ";
	ForwardNuc(curVertex, outIndex, nextVertex);
	
	// unlink the first edge

	ssize_t nextVertexIndex;
	nextVertexIndex = LookupVertex(vertices, numVertices, nextVertex);
	
	// this is already following an edge, so the min length is 1
  ssize_t length = 1;

	assert(nextVertexIndex >= 0);
	assert(nextVertexIndex < numVertices);

	// Follow this until the next branching vertex
	CountedIntegralTuple rctup;
  while (nextVertexIndex != destVertexIndex) {
		assert(adj[nextVertexIndex].OutDegree() == 1);
		// Record the character here
		
		// Record where we have been.
		curVertex        = nextVertex;
		curVertexIndex   = nextVertexIndex;
		/*		vertices[curVertexIndex].ToString(tupstr);
		vertices[curVertexIndex].MakeRC(rctup);
		cout << tupstr << " " << curVertexIndex << ", " << vertices[curVertexIndex].tuple << " . (rc " << rctup.tuple << ") ";
		*/
		// Which nucleotide leaves this vertex
    outIndex = adj[curVertexIndex].DegreeOneOutIndex();

		// unlink this vertex.
		adj[curVertexIndex].Clear();

		// Lookup the next vertex.
    assert(outIndex < 4);
		ForwardNuc(curVertex, outIndex, nextVertex);
		nextVertexIndex = LookupVertex(vertices, numVertices, nextVertex);
		assert(nextVertexIndex != -1);
		// Increase the length of this edge
    ++length;

		// If a cycle was hit, stop.  This is necessary because
		// when the edges are unlinked after traversal, the vertex may
		// stop being a branching vertex.  
    if (nextVertexIndex == firstVertexIndex)
      break;
	}
	// Erase the last edge on this path.
	ssize_t srcIndex;
	vertices[nextVertexIndex].MakeRC(rctup);
	vertices[nextVertexIndex].ToString(tupstr);
	/*	cout << tupstr << " " << nextVertexIndex << " . " << vertices[nextVertexIndex].tuple << " rc (" << rctup.tuple << ") " << endl;
	 */
	//	srcIndex = curVertex.tuple & 3;
	srcIndex = curVertex.GetNucleotide0();
	adj[nextVertexIndex].SetSrc(srcIndex, 0);
	//	cout << endl;
}

ssize_t AdvanceSimplePath(ssize_t curVertexIndex, ssize_t outIndex, 
											CountedIntegralTuple *vertices, Adjacency *adj, ssize_t numVertices,
											ssize_t &lastVertexIndex, unsigned char &srcIndex,	unsigned char* seq ) {
  // Move forward along a path until a branching vertex has been reached
	// Start at vertex 'vertex' in the direction of 'out'.
	// 
	//  vertex = vertices[vertex].out[out];

	ssize_t firstVertexIndex = curVertexIndex;
	CountedIntegralTuple curVertex;
	CountedIntegralTuple nextVertex;
	curVertex   = vertices[curVertexIndex];
	ForwardNuc(curVertex, outIndex, nextVertex);
	ssize_t nextVertexIndex;
	nextVertexIndex = LookupVertex(vertices, numVertices, nextVertex);
	
	// this is already following an edge, so the min length is 1
  ssize_t length = 1;

	assert(nextVertexIndex >= 0);
	assert(nextVertexIndex < numVertices);

	// Follow this until the next branching vertex
  while (adj[nextVertexIndex].InDegree() == 1 and
				 adj[nextVertexIndex].OutDegree() == 1) {
		// Record the character here
    if (seq != NULL) {
      seq[length - 1] = nuc_char[outIndex];
    }

		// Record where we have been.
		curVertex        = nextVertex;
		curVertexIndex   = nextVertexIndex;
		// Which nucleotide leaves this vertex
    outIndex = adj[curVertexIndex].DegreeOneOutIndex();

		if (seq != NULL) {
			// Unlink this edge once it has been traversed so that it is no longer
			// used.
			adj[curVertexIndex].SetDest(outIndex, 0);
		}

		// Lookup the next vertex.
    assert(outIndex < 4);
		ForwardNuc(curVertex, outIndex, nextVertex);
		nextVertexIndex = LookupVertex(vertices, numVertices, nextVertex);
		assert(nextVertexIndex != -1);
		// Increase the length of this edge
    ++length;

		// If a cycle was hit, stop.  This is necessary because
		// when the edges are unlinked after traversal, the vertex may
		// stop being a branching vertex.  
    if (nextVertexIndex == firstVertexIndex)
      break;
  }

	// Append the character of the last vertex
  if (seq != NULL) {
    seq[length - 1] = nuc_char[outIndex];
  }
	
	// Store where this ended
  lastVertexIndex = nextVertexIndex;

	// Record which slot in the dest vertex this should arrive.
	//	srcIndex = curVertex.tuple & 3;
	srcIndex = curVertex.GetNucleotide0();
  return length;
}


ssize_t RemoveZeroDegreeDBVertices(CountedIntegralTuple *deBruijn, Adjacency *deBruijnAdj, 
																ssize_t numVertices) {
	ssize_t v;
	ssize_t c = 0;
	for (v = 0; v < numVertices; v++ ) {
		if (deBruijnAdj[v].OutDegree() != 0 or
				deBruijnAdj[v].InDegree() != 0) {
			deBruijnAdj[c] = deBruijnAdj[v];
			deBruijn[c] = deBruijn[v];
			c++;
		}
	}
	return c;
}
		

void TrimSmallEdges(CountedIntegralTuple *deBruijn, Adjacency *deBruijnAdj, ssize_t numVertices, ssize_t trimEnd) {
	ssize_t v;
	ssize_t numTrimmed = 0;

	for (v = 0; v < numVertices; v++) {
		if (deBruijnAdj[v].InDegree() == 0 and deBruijnAdj[v].OutDegree() == 1) {
			ssize_t lastVertexIndex;
			unsigned char srcIndex;
			unsigned char e = deBruijnAdj[v].DegreeOneOutIndex();
			ssize_t pathLength;

			pathLength =  AdvanceSimplePath(v,e,
																			deBruijn, deBruijnAdj, numVertices,
																			lastVertexIndex, srcIndex, (unsigned char*) NULL);

			ssize_t sourceRC, destRC, edgeRC;
			CountedIntegralTuple tupRC;
			deBruijn[v].MakeRC(tupRC);
			sourceRC = LookupVertex(deBruijn, numVertices, tupRC);
			deBruijn[lastVertexIndex].MakeRC(tupRC);
			destRC = LookupVertex(deBruijn, numVertices, tupRC);
			edgeRC = comp_bin[srcIndex];

			ssize_t lastRCIndex;
			unsigned char rcSrcIndex;
			pathLength =  AdvanceSimplePath(destRC, edgeRC,
																			deBruijn, deBruijnAdj, numVertices,
																			lastRCIndex, rcSrcIndex, (unsigned char*) NULL);
			assert(lastRCIndex == sourceRC);

			// The section of code below was disabled.
			// I'm re-enabling it; to suppress it, use command line option "-trimShort 0"
			// -- GT

			/* BEGIN previously disabled code */
			if (pathLength < trimEnd and deBruijnAdj[v].coverage < 2) {
				//				numTrimmed+=2;

				// This path is short, remove it from the assembly.
				//				cout << "erase for: " << endl;
				EraseSimplePath(v, e, lastVertexIndex, 
												deBruijn, deBruijnAdj, numVertices);
				numTrimmed++;
				if (deBruijn[v] != deBruijn[destRC] and // compare tuples
						deBruijn[lastVertexIndex] != deBruijn[sourceRC]) { // compare tuples
					//				if (deBruijn[v].tuple != deBruijn[destRC].tuple and
					//						deBruijn[lastVertexIndex].tuple != deBruijn[sourceRC].tuple) {
					string firstFor, firstRev;
					deBruijn[v].ToString(firstFor);
					deBruijn[sourceRC].ToString(firstRev);
					EraseSimplePath(destRC, edgeRC, sourceRC,
													deBruijn, deBruijnAdj, numVertices);
					numTrimmed++;
				}
			}
			/* END previously disabled code */
		}
	}
	cout << "trimmed " << numTrimmed << endl;

}

void PrintSmallEdges(CountedIntegralTuple *deBruijn, Adjacency *deBruijnAdj, ssize_t numVertices) {
	ssize_t v;
	ssize_t numShort = 0;
	for (v = 0; v < numVertices; v++) {
		if (deBruijnAdj[v].InDegree() == 0) {
			ssize_t lastVertexIndex;
			unsigned char srcIndex;
			unsigned char e = deBruijnAdj[v].DegreeOneOutIndex();
			ssize_t pathLength;
			unsigned char *ptr = NULL;
			pathLength =  AdvanceSimplePath(v,e,
																			deBruijn, deBruijnAdj, numVertices,
																			lastVertexIndex, srcIndex, ptr);
			if (pathLength < 10) {
				numShort++;
			}
		}
	}
	cout << "counted: " << numShort << " short paths" << endl;
	exit(0);
}


ssize_t RemoveLowCoverageEdges(CountedIntegralTuple *deBruijn, Adjacency *deBruijnAdj, ssize_t numVertices, int minCoverage){

	ssize_t v;
	ssize_t outEdge;
	ssize_t numRemoved = 0;
	for (v = 0; v < numVertices; v++ ){ 
		int e;
		for (e = 0; e < 4 ; e++ ) {
			if (deBruijnAdj[v].GetDest(e) and deBruijnAdj[v].destCoverage[e] < minCoverage ) {
				// Remove this adjacency.
				deBruijnAdj[v].destCoverage[e] = 0;
				deBruijnAdj[v].SetDest(e,0);
				ssize_t destVertexIndex;
				CountedIntegralTuple destVertexTuple, *destVertexPtr;
				ForwardNuc(deBruijn[v], e, destVertexTuple);
				destVertexPtr = std::lower_bound(deBruijn, deBruijn + numVertices, destVertexTuple);
				assert(*destVertexPtr == destVertexTuple);
				ssize_t destVertexAdj = destVertexPtr - deBruijn;
				unsigned char toNuc = deBruijn[v].GetNucleotide(0);
				deBruijnAdj[destVertexAdj].SetSrc(toNuc, 0);
				++numRemoved;
			}
		}
	}
	return numRemoved;
}
