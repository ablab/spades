/***************************************************************************
 * Title:          IntegralPrintReadIntervals.cpp
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  02/28/2010
 *
 * Copyright (c) 2007-2010 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <iostream>
#include <string>
#include <vector>
#include "compatibility.h"
#include "DeBruijnGraph.h"
#include "ReadIntervals.h"
#include "SeqReader.h"
#include "SeqUtils.h"
#include "SimpleSequence.h"
#include "IntegralTupleStatic.h"
#include "IntervalGraph.h"


using namespace std;

void PrintUsage() {
  std::cout << "usage: integralPrintReadIntervals sequences " << endl << endl
						<< "[-mapGappedReads]    Try and map reads even when not every " << endl
						<< "                     k+1-mer in the read maps." << endl << endl 
						<< "[-printCoverage out] Just print the multiplicity of " << endl
						<< "                     edges to the branching graph." << endl << endl
            << "[-altEdgeOvp ovp]    A list of k+1-mers " << endl
						<< "                     from the alternative edge list." << endl;
}


// class MappedRead {
// public:
// 	ssize_t readIndex;
// 	ssize_t readPos;
// 	MappedRead(ssize_t i, ssize_t p) {
// 		readIndex = i;
// 		readPos   = p;
// 	}
// 	MappedRead &operator=(const MappedRead &rhs) {
// 		if (this != &rhs) {
// 			readIndex = rhs.readIndex;
// 			readPos   = rhs.readPos;
// 		}
// 		return *this;
// 	}
// };

ssize_t ThreadReadAlternativeEdge(DNASequence &read, 
													    ssize_t &readPos, 
															IntervalGraph &g, 
															ssize_t altEdgeIndex, 
															ssize_t altEdgePos,
															ssize_t readIndex);

ssize_t ThreadReadExact(DNASequence &read, ssize_t &readPos, IntervalGraph &g, ssize_t edgeIndex, ssize_t edgePos,
										ssize_t &curIntv,
										std::vector<ssize_t> &readPosList, 
										std::vector<ssize_t> &edgeIndexList, 
										std::vector<ssize_t> &edgePosList,
										std::vector<ssize_t> &lengthList);


//typedef vector<MappedRead> MappedReadList;
//typedef vector<MappedReadList> MappedReadMatrix;

// TODO: Should this routine be here or in another file?
std::string HighlightString(unsigned char *seq, ssize_t seqLength,
														ssize_t start, ssize_t windowLength);

int main(int argc, char* argv[] ) {

  std::string overlapFileName;
	std::string edgeFileName;
  std::string sequenceListFileName;
  std::string readIntervalFileName;
	std::string pathFileName;
	std::string bgraphFileName;
  std::string altEdgePosTupleFileName;
  _INT_ vertexSize, overlapSize;

  if (argc < 2) {
    PrintUsage();
    exit(1);
  }

  int argi = 1;
  sequenceListFileName = argv[argi++];



	bgraphFileName = sequenceListFileName + ".bgraph";
	overlapFileName = sequenceListFileName + ".iovp";
	edgeFileName = sequenceListFileName + ".edge";
	readIntervalFileName = sequenceListFileName + ".intv";
	pathFileName  = sequenceListFileName + ".path";
	
	ssize_t mapGappedReads = 0;
	ssize_t printCoverageOnly = 0;
  ssize_t useAltEdgeList = 0;
	while (argi < argc) {
		if (strcmp(argv[argi], "-mapGappedReads") == 0){ 
			mapGappedReads = 1;
		}
		else if (strcmp(argv[argi], "-printCoverage") == 0) {
			printCoverageOnly = 1;
		}
    else if (strcmp(argv[argi], "-altEdgeOvp") == 0) {
       useAltEdgeList = 1;
       altEdgePosTupleFileName = argv[++argi];
    }
		else {
			PrintUsage();
			cout << "bad option: " << argv[argi] << endl;
			exit(0);
		}
		++argi;
	}
	
	std::string reportFileName = FormReportName(sequenceListFileName);
	std::ofstream reportFile;
	openck(reportFileName, reportFile, std::ios::app, reportFile);
	BeginReport(argc, argv, reportFile);

	std::ifstream readsIn;
	openck(sequenceListFileName, readsIn, std::ios::in, reportFile);

	//	std::ifstream edgeIn;
	//	openck(edgeFileName, edgeIn, std::ios::in, reportFile);

	
	IntervalGraph graph;
	_INT_ skipIntervals = 1; // Flag to skip reading intervals since they are stored here.

	ReadIntervalGraph(sequenceListFileName, graph, vertexSize, skipIntervals, reportFile);
  overlapSize = vertexSize + 1;

	//  ReadSequences(edgeFileName, graph.edges, reportFile); // Seems to be redundant with ReadIntervalGraph() above
	//  cout << "read edges." << endl;

  // Read in the vertex list

	// Read the overlaps.
	ssize_t numOverlaps;
  std::ifstream overlapIn;
	EdgePosIntegralTuple *edgePosTupleList;


	openck(overlapFileName, overlapIn, std::ios::in | std::ios::binary, reportFile);
	overlapIn.read((char*) &numOverlaps, sizeof(ssize_t));
	edgePosTupleList = new EdgePosIntegralTuple[numOverlaps];
	overlapIn.read((char*) edgePosTupleList, sizeof(EdgePosIntegralTuple)*numOverlaps);
	std::cout << "read: " << numOverlaps << " overlaps." << std::endl;
  //
  // Read in the overlaps for alternative (deleted) edges 
	// if they are specified.	
  //
  EdgePosIntegralTuple *altEdgePosTupleList = (EdgePosIntegralTuple *) NULL;	
	ssize_t nAltEdgePosTuples;
  if (useAltEdgeList) {
    // read in the alt edge tuples. 
		cout << "READING ALT EDGE TUPLES." << endl;
    ifstream altEdgeTupleFile;
	  openck(altEdgePosTupleFileName, altEdgeTupleFile, std::ios::in | std::ios::binary, reportFile);	
		// TODO: should we add tupleSize to start of file for consistency with .spect files?
    altEdgeTupleFile.read((char*) &nAltEdgePosTuples, sizeof(ssize_t));
		altEdgePosTupleList = new EdgePosIntegralTuple[nAltEdgePosTuples];
    altEdgeTupleFile.read((char*) altEdgePosTupleList, sizeof(EdgePosIntegralTuple)*nAltEdgePosTuples);
  } 
  PathIntervalList paths;
	PathLengthList pathLengths;
	std::vector<BareReadIntervalList> edgeReadIntervals;
	//	vector<MappedReadMatrix> edgeReadMap;

	//	ssize_t e;
	//	if (!printCoverageOnly) {
	//		edgeReadIntervals.resize(graph.edges.size());
	//		//		edgeFullReadIntervals.resize(edges.size());
	//		edgeReadMap.resize(graph.edges.size());
	//		for (e = 0; e < graph.edges.size(); e++ ){ 
	//			edgeReadMap[e].resize(graph.edges[e].length - vertexSize);
	//		}
	//	}

	if (!printCoverageOnly) {
		edgeReadIntervals.resize(graph.edges.size());
	}
	
	std::cout << "storing read intervals. " << std::endl;
	ssize_t skipGapped = 1;


	DNASequence read, readRC;
	PathIntervalList path;
	ssize_t readIndex;
	DNASequence* readPair[2], *readPtr;
	readPair[0] = &read;
	readPair[1] = &readRC;
	readIndex = 0;

	EdgePosIntegralTuple readTuple;
	//	EdgePosIntegralTuple::tupleSize = overlapSize;
	EdgePosIntegralTuple::SetTupleSize(overlapSize);
	ssize_t rn = 0;

	// How intervals from the current read map to the graph,
	// in order of increasing positions within the read:
	// interval i:
	//     nucleotides readPosList[i] through readPosList[i]+lengthList[i]-1
	// maps to edge number edgeList[i], starting at position edgePosList[i]
	//
	// rc-versions for dual read, and the intervals will go in opposite order:
	// interval i of forward read corresponds to numIntervals-1-i for rc

	std::vector<ssize_t> edgeList, lengthList, readPosList, edgePosList;
	std::vector<ssize_t> rcEdgeList, rcLengthList, rcReadPosList, rcEdgePosList;	


	while (SeqReader::GetSeq(readsIn, read, SeqReader::noConvert)) {

		if (skipGapped) {
			ssize_t seqIsGapped = 0;
			ssize_t curPos;
			for (curPos = 0; curPos < read.length; curPos++ ){ 
				if (unmasked_nuc_index[read.seq[curPos]] >= 4) {
					seqIsGapped = 1;
					break;
				}
			}
			if (seqIsGapped) {
				paths.push_back(NULL);
				pathLengths.push_back(0);
				paths.push_back(NULL);
				pathLengths.push_back(0);				
				rn++;
				continue;
			}
		}

		MakeRC(read, readRC);
		ssize_t edgeIndex = -1,  edgePos = -1;
    ssize_t altEdgeIndex = -1,  altEdgePos = -1;
		//UNUSED// ssize_t r, s;
		ssize_t overlapIndex;
		ssize_t nextN;
		ReadIntervalList readPath;
		//UNUSED// ssize_t pathPos;


		// simplify the code by an RC
		ssize_t spacing = 100000;
		PrintStatus(rn, spacing);

		//UNUSED// ssize_t forPathLength  = -1;
		
		if (read.length < overlapSize) {
			// add placeholders for forward/reverse strand
			paths.push_back(NULL);
			pathLengths.push_back(0);
			paths.push_back(NULL);
			pathLengths.push_back(0);
			rn++;
			continue;
		}
		//		for (r = 0; r < 2; r++) {
		//    std::cout << "processing sequence " << s << std::endl;
		readPath.clear();
		readPtr = readPair[0];
		nextN   = -1;
		ssize_t p;
		//UNUSED// ssize_t prevEdge = -1;
		ssize_t fullReadMapped = 1;
		//UNUSED// ssize_t gapSpansEdges = 0;

		// 
		ssize_t readPos;
		//UNUSED// ssize_t npos = 0;
		//UNUSED// ssize_t nEdges = 0;
		ssize_t listSize = readPtr->length - overlapSize + 1;
		if (readPosList.size() < listSize) {
			readPosList.resize(listSize);
			edgePosList.resize(listSize);
			edgeList.resize(listSize);
			lengthList.resize(listSize);
			rcReadPosList.resize(listSize);
			rcEdgePosList.resize(listSize);
			rcEdgeList.resize(listSize);
			rcLengthList.resize(listSize);
		}

		for (p = std::min((_SSZT_) (overlapSize-1), readPtr->length); p >= 0; p--) { // TODO: remove cast
			if (numeric_nuc_index[read.seq[p]] >= 4) {
				nextN = p;
				fullReadMapped = 0;
				break;
			}
		}

		ssize_t numIntervals, curIntv;
		readPos = 0;
		numIntervals = 0;
		curIntv = 0;
		while (readPos < listSize) {
				
			// 
			// Locate a position on the read starting at readPos that maps
			// to an edge in the graph.
			//
			edgeIndex = -1;
			for (; readPos < listSize; readPos++ ) {
				//
        //  If the read has a masked nucleotide, skip past this position in the mapping
	      //
				if (numeric_nuc_index[readPtr->seq[readPos + overlapSize - 1]] >= 4) {
					nextN = readPos + overlapSize - 1;
					fullReadMapped = 0;
				}
        //
        //  Otherwise, this position may be a valid tuple, try and locate it in the
			  //  list of edges.
        //
				if (readPos > nextN) {
					readTuple.StringToTuple(&readPtr->seq[readPos]);
					overlapIndex = LookupBinaryTuple(edgePosTupleList, numOverlaps, readTuple);
				  //
 				  // Did not map the read to an edge in the regular edge list.
          // try and map it to the list of alternative edges.
          //
					if (overlapIndex < 0) {
            if (useAltEdgeList) {
              overlapIndex = LookupBinaryTuple(altEdgePosTupleList, nAltEdgePosTuples, readTuple);
							if (overlapIndex >= 0) {
								altEdgeIndex = altEdgePosTupleList[overlapIndex].edge;
								altEdgePos   = altEdgePosTupleList[overlapIndex].pos;
                edgeIndex = -1;
                edgePos   = -1;


								ssize_t  destEdge = graph.altEdges[altEdgeIndex].destEdge;
								ssize_t  destEdgeLength = graph.edges[destEdge].length;


								// Detect if tuple from altEdge would overhang destEdge.
								// As long as initial tuple does not overhang,
								// when ThreadReadAlternativeEdge extends the initial tuple,
								// it will stop short of overhanging it.
								if (altEdgePos + graph.altEdges[altEdgeIndex].offset + IntegralTuple::tupleSize > destEdgeLength) {
									//									string ts;
									//									readTuple.ToString(ts);

									// show read, alt edge, dest edge sequences, aligned at tuple
									ssize_t destEdgePos = altEdgePos + graph.altEdges[altEdgeIndex].offset;
									ssize_t maxOffset = std::max(std::max(readPos,destEdgePos),altEdgePos);

									cout
										<< "Bad alternative edge: found tuple but it overhangs destEdge: "
										<< " readNum/2=" << rn
										<< " readPos=" << readPos
										<< " readLength=" << readPtr->length
										<< " / altEdgeIndex=" << altEdgeIndex
										<< " altEdgePos=" << altEdgePos
										<< " altEdgeLength=" << graph.altEdges[altEdgeIndex].seq.length
										<< " / destEdge=" << destEdge
										<< " destEdgePos=" << destEdgePos
										<< " destEdgeLength=" << destEdgeLength
										<< " / overlapIndex=" << overlapIndex
										<< " altEdgeOffset=" << graph.altEdges[altEdgeIndex].offset
										<< " tupleSize=" << IntegralTuple::tupleSize
										<< endl;

									//									cout << "Tuple:    " << ts << endl;
									cout
										<< "Read:     "
										<< string(maxOffset-readPos, ' ')
										<< HighlightString(read.seq,read.length,readPos,IntegralTuple::tupleSize)
										<< endl;
									cout
										<< "altEdge:  "
										<< string(maxOffset-altEdgePos, ' ')
										<< HighlightString(graph.altEdges[altEdgeIndex].seq.seq,
																			 graph.altEdges[altEdgeIndex].seq.length,
																			 altEdgePos,
																			 IntegralTuple::tupleSize)
										<< endl;
									cout
										<< "destEdge: "
										<< string(maxOffset-destEdgePos, ' ')
										<< HighlightString(graph.edges[destEdge].seq.seq,
																			 graph.edges[destEdge].seq.length,
																			 altEdgePos + graph.altEdges[altEdgeIndex].offset,
																			 
																			 IntegralTuple::tupleSize)
										<< endl;
										

									// Cancel this alternative edge
									altEdgeIndex = -1;
									altEdgePos = -1;

								} else {
									break;
								}
							}   
            }
            else {
	  					fullReadMapped = 0;
            }
					}
					else {
						edgeIndex = edgePosTupleList[overlapIndex].edge;
						edgePos   = edgePosTupleList[overlapIndex].pos;
						altEdgeIndex = -1;  altEdgePos = -1;
						break;
					}
				}
			}
			if (readPos < listSize) {
				if (edgeIndex != -1) {
					// An anchoring edge was found. 
					// Thread the rest of the read through the graph.
					//				readPosList.clear(); edgeList.clear(); edgePosList.clear(); lengthList.clear();
					//UNUSED// ssize_t startReadPos = readPos;
					numIntervals += ThreadReadExact(*readPtr, readPos, graph, edgeIndex, edgePos,
																					curIntv,
																					readPosList, edgeList, edgePosList, lengthList);
				}
				else if (altEdgeIndex != -1) {
					ssize_t altEdgeThreadLength;
					ssize_t prevReadPos = readPos;


					// Determine how many consecutive tuples from read and
					// alternative edge match (starting in readPtr at position readPos,
					// and in altEdgeIndex at altEdgePos).
					// Update readPos to point to the next non-matching nucleotide,
					// or truncate to length of read or current edge.
					altEdgeThreadLength = ThreadReadAlternativeEdge(*readPtr, readPos, graph, altEdgeIndex, altEdgePos, paths.size());


					/*
						cout << paths.size() << " threaded alternative edge of length: " << altEdgeThreadLength
							 << " to alt edge " << altEdgeIndex << " pos: " << altEdgePos 
							 << " " << prevReadPos << " " << readPos << endl;
					*/

					readPosList[curIntv] = prevReadPos;
					edgeList[curIntv]    = graph.altEdges[altEdgeIndex].destEdge;
					edgePosList[curIntv] = graph.altEdges[altEdgeIndex].offset + altEdgePos;
					lengthList[curIntv]  = readPos - prevReadPos + IntegralTuple::tupleSize;

						// Threading a read through an alternative edge adds only
								// one interval because the alternative edge sequence
								// is mapped to just one edge.
						++numIntervals;
						++curIntv;
				}
			}
			//
			// This threads to the end of a match, or the end of a read.
			// If the whole read is mapped, done anyway.  Otherwise, part
			// of the read may be missing due to an error in it.
			// Move forward with the read pos to attempt to thread the rest
			// of the read.

			++readPos;
		}
		// end trying to fit the entire read into the graph.
			
		if (numIntervals > 0) {
			// Create the edge and path intervals
			paths.push_back(new PathInterval[numIntervals]);
			pathLengths.push_back(numIntervals);
			readIndex = paths.size() - 1;
			ssize_t i;
			for (i =0 ; i < numIntervals; i++) {
				edgeReadIntervals[edgeList[i]].
					push_back(BareReadInterval(readIndex,
																		 readPosList[i],
																		 edgePosList[i],
																		 lengthList[i]));
				paths[readIndex][i].edge  = edgeList[i];
				paths[readIndex][i].index = edgeReadIntervals[edgeList[i]].size() - 1;//readPosList[i];
			}

			// 
			// Add the reverse complement path.
			paths.push_back(new PathInterval[numIntervals]);
			pathLengths.push_back(numIntervals);
			ssize_t rci = numIntervals - 1;
			for (i = 0; i < numIntervals ;i++) {
				ssize_t rcEdge =  graph.edges[edgeList[i]].balancedEdge;

				assert(readPtr->length - (readPosList[i] + lengthList[i]) >= 0);
				assert(graph.edges[edgeList[i]].length - (edgePosList[i] + lengthList[i]) >=0);
				edgeReadIntervals[rcEdge].push_back(BareReadInterval(readIndex + 1,
																														 readPtr->length - (readPosList[i] + lengthList[i]),
																														 graph.edges[edgeList[i]].length - (edgePosList[i] + lengthList[i]),
																														 lengthList[i]));
				paths[readIndex + 1][rci].edge = rcEdge;
				paths[readIndex + 1][rci].index = edgeReadIntervals[rcEdge].size() - 1;
				--rci;
			}
																															 
			/*																															 
																																			 rcReadPosList[rci] = readPtr->length - (readPosList[i] + lengthList[i]);
																																			 rcEdgePosList[rci] = graph.edges[edgeList[i]].length - (edgePosList[i] + lengthList[i]);
																																			 rcLengthList[rci]  = lengthList[i];
			*/
		}
		else {
			// The read could not be mapped back to the
			// graph.  Add a placeholder for this.
			paths.push_back(NULL);
			pathLengths.push_back(0);	
			paths.push_back(NULL);
			pathLengths.push_back(0);	
		}
			//		}
		++rn;
		assert(paths.size() == rn*2);

		read.Reset();
		readRC.Reset();
	}

	std::cout << std::endl; // Terminate status dots

	if (!printCoverageOnly) {
		std::ofstream readIntervalOut;
		openck(readIntervalFileName, readIntervalOut, std::ios::out, reportFile);

		ssize_t edgeIndex;
		ssize_t intervalIndex;
		intervalIndex = 0;
		ssize_t ri;
		for (edgeIndex = 0; edgeIndex < graph.edges.size();edgeIndex++) {
			readIntervalOut << "EDGE " << edgeIndex 
											<< " Length " << graph.edges[edgeIndex].length
											<< " Multiplicity " << edgeReadIntervals[edgeIndex].size()
											<< std::endl;
			for (ri = 0; ri < edgeReadIntervals[edgeIndex].size(); ri++) {
				readIntervalOut << "INTV " << edgeReadIntervals[edgeIndex][ri].read 
												<< " " << edgeReadIntervals[edgeIndex][ri].readPos 
												<< " " << edgeReadIntervals[edgeIndex][ri].length
												<< " " << edgeReadIntervals[edgeIndex][ri].edgePos << std::endl;
			}
		}
		readIntervalOut.close();
		// Read intervals are fixed, now determine the index of each path interval
		// into the edge interval list.
		/*
		ssize_t pathIndex, readPos;
		for (e = 0; e < edgeReadIntervals.size(); e++) {

			ssize_t allFound = 1;
			for (ri = 0; ri < edgeReadIntervals[e].size(); ri++) {
				pathIndex = edgeReadIntervals[e][ri].read;
				// The index of every path is temporarily set to 
				// equal the read position.  Since there is only
				// one read position per path and per edge, the order
				// of the read positions correponds to the order of
				// the intervals on the edge (since the read positions 
				// are assigned in the order that intervals are added).
				//
				readPos   = edgeReadIntervals[e][ri].readPos;
				ssize_t pi;
				for (pi = 0; pi < pathLengths[pathIndex]; pi++) {
					if (paths[pathIndex][pi].index == readPos and
							paths[pathIndex][pi].set == 0) {
						assert(paths[pathIndex][pi].edge == e);
						paths[pathIndex][pi].index = ri;
						paths[pathIndex][pi].set = 1;
						break;
					}
				}
				if (pi == pathLengths[pathIndex]) {
					allFound = 0;
					assert(0);
				}
			}
		}
		*/
		cout << "writing to paths: " << pathFileName << endl;
		WriteReadPaths(pathFileName, paths, pathLengths, reportFile);
	}
	EndReport(reportFile);
	return 0;
}

ssize_t ThreadReadAlternativeEdge(DNASequence &read, ssize_t &readPos, IntervalGraph &g, ssize_t altEdgeIndex, ssize_t altEdgePos, ssize_t readIndex) {

  unsigned char *altEdgeSeq   = g.altEdges[altEdgeIndex].seq.seq;
  ssize_t  altEdgeLength = g.altEdges[altEdgeIndex].seq.length;

	ssize_t  destEdge = g.altEdges[altEdgeIndex].destEdge;
	ssize_t  destEdgeLength = g.edges[destEdge].length;
	// Truncate alternative edge sequence if it would overhang destEdge
	altEdgeLength = std::min(altEdgeLength, destEdgeLength - g.altEdges[altEdgeIndex].offset);

  ssize_t  ep = altEdgePos + IntegralTuple::tupleSize;
  ssize_t  rp = readPos    + IntegralTuple::tupleSize;
	ssize_t  prevReadPos = readPos;
	while (ep < altEdgeLength and rp < read.length and read.seq[rp] == altEdgeSeq[ep]) { ep++; rp++; readPos++; }
  return rp - prevReadPos - IntegralTuple::tupleSize;
} 


ssize_t ThreadReadExact(DNASequence &read, ssize_t &readPos, IntervalGraph &g, ssize_t edgeIndex, ssize_t edgePos,
										ssize_t &curIntv,
										std::vector<ssize_t> &readPosList, 
										std::vector<ssize_t> &edgeIndexList, 
										std::vector<ssize_t> &edgePosList,
										std::vector<ssize_t> &lengthList) {
	
	// This begins assuming that read[readPos... readPos + tupleSize -1] ==
	// edges[edgeIndex][edgePos... edgePos+tuplesize]
	
	ssize_t intvLength;
	ssize_t numNewIntv = 0;
	while (readPos < read.length) {
		intvLength = 0;
		ssize_t edgeStartPos = edgePos;
		ssize_t readStartPos = readPos;
		while((readPos < read.length - IntegralTuple::tupleSize) and
					(edgePos < g.edges[edgeIndex].length - IntegralTuple::tupleSize) and
					(read.seq[readPos+IntegralTuple::tupleSize] == 
					 g.edges[edgeIndex].seq.seq[edgePos+IntegralTuple::tupleSize])) {
			++readPos;
			++edgePos;
			++intvLength;
		}
		readPosList[curIntv]   = readStartPos;
		edgeIndexList[curIntv] = edgeIndex;
		edgePosList[curIntv]   = edgeStartPos;
		lengthList[curIntv]    = intvLength + IntegralTuple::tupleSize;
		++curIntv;
		++numNewIntv;
		if ((readPos == read.length - IntegralTuple::tupleSize  ) or
				(edgePos < g.edges[edgeIndex].length - IntegralTuple::tupleSize))
			return numNewIntv;

		// This is done when the read is fully matched (first case)
		// or there is an error and the read is not fully mapped (second)

		edgeStartPos = 0;
		intvLength   = 0;

		// Determine the next edge.
		ssize_t outEdge, outEdgeIndex;
		ssize_t dest;
		dest = g.edges[edgeIndex].dest;
		ssize_t foundOut = 0;
		for (outEdgeIndex = g.vertices[dest].FirstOut();
				 outEdgeIndex != g.vertices[dest].EndOut();
				 outEdgeIndex = g.vertices[dest].NextOut(outEdgeIndex)) {
			outEdge= g.vertices[dest].out[outEdgeIndex];
			assert(g.edges[outEdge].length >= IntegralTuple::tupleSize);
			if (read.seq[readPos + IntegralTuple::tupleSize] == 
					g.edges[outEdge].seq.seq[IntegralTuple::tupleSize - 1]) {
				foundOut = 1;
				break;
			}
		}

		if (!foundOut) {
			// No possible extension worked. stop.
			return numNewIntv;
		}
		else {
			// Next interval starts at the beginning of the edge
			readPos++; 
			edgePos = 0;
			edgeIndex = outEdge;
		}
	}
	return numNewIntv;
}

// TODO: Should this routine be here or in another file?
std::string HighlightString(unsigned char *seq, ssize_t seqLength,
														ssize_t start, ssize_t windowLength) {
	std::string out((const char *) seq, seqLength);

	if (out.length() != seqLength) {
		std::cerr
			<< "HighlightString: length mismatch; seqLength=" << seqLength
			<< " but got " << out.length() << " characters" << std::endl;
	}

	// If sequence is too short, pad with dashes
	if (out.length() < start) {
		out.insert(out.length(), start - out.length(), '-');
	}

	// open region to highlight
	out.insert(start, 1, '[');

	// If sequence is too short, pad with dashes
	if (out.length() < start + windowLength + 1) {
		out.insert(out.length(), start + windowLength + 1 - out.length(), '-');
	}

	// close region to highlight
	out.insert(start + 1 + windowLength, 1, ']');
	return out;
}
