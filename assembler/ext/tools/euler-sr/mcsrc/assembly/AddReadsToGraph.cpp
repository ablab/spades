/***************************************************************************
 * Title:          AddReadsToGraph.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  10/01/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/

#include "compatibility.h"
#include "IntervalGraph.h"
#include "SeqReader.h"
#include "ReadPos.h"
#include "SortedTupleList.h"
#include "IntegralTupleStatic.h"



ssize_t CountReadPosMatches(ReadPositions &positions, ssize_t curIndex) {
	ssize_t nextIndex = curIndex;
	if (curIndex < positions.size()) { 
		nextIndex = curIndex+1;
		while (nextIndex < positions.size() and
					 positions[nextIndex] == positions[curIndex])
			nextIndex++;

	}
	return nextIndex - curIndex;
}

void PrintUsage() {
	std::cout << "addReadsToGraph- Given a list of reads, try and map them" << std::endl;
	std::cout << "  back to a graph.  Any read that does map is added to " << std::endl;
	std::cout << "  the graph." << std::endl;
	std::cout << "usage: addReadsToGraph graphBase curReads extraReads vtxSize prefixLength graphOutBase mappedReads" << std::endl;
}

int main(int argc, char* argv[]) {
	
	std::string graphBase, oldReadsName, newReadsName, 
		mappedNewReadsName, graphOutBase;
	int vertexSize, prefixLength;
	std::string edgeFile;
	int argi = 1;
	if (argc < 7) {
		PrintUsage();
		exit(0);
	}
	graphBase = argv[argi++];
	oldReadsName = argv[argi++];
	newReadsName = argv[argi++];
	vertexSize = atoi(argv[argi++]);
	prefixLength = atoi(argv[argi++]);
	graphOutBase = argv[argi++];
	mappedNewReadsName = argv[argi++];
	
	edgeFile = graphBase + ".edge";
	
	
	IntervalGraph graph;
	ReadIntervalGraph(graphBase, graph, vertexSize);
	
	std::ifstream newReads;
	openck(newReadsName, newReads, std::ios::in);
	
	// Quickly count the number of reads in the old reads file.
	std::ifstream oldReadsIn;
	openck(oldReadsName, oldReadsIn, std::ios::in);
	std::string line;
	ssize_t nOldReads = 0;

	std::ofstream mappedNewReads;
	openck(mappedNewReadsName, mappedNewReads, std::ios::out);
	char c;
	while(oldReadsIn) {
		if ((c = oldReadsIn.get()) == '>')
			nOldReads++;
		mappedNewReads << c;
	}
	oldReadsIn.close();
	

	ReadPositions positions;
	// Package the edges in a structure that may 
	// be accessed easily
	SimpleSequenceList edgeSequences;	
	
	// Make a list of query-able read positions.
	ssize_t numPos = 0;
	ssize_t e, p;
	
	// map reads according to de Bruijn edge size.
	_INT_ dbEdgeLength = vertexSize + 1; 
	edgeSequences.resize(graph.edges.size());
	for (e = 0; e < graph.edges.size(); e++) {
		if (graph.edges[e].length - dbEdgeLength + 1 > 0)
			numPos += graph.edges[e].length - dbEdgeLength + 1;
	}
	positions.resize(numPos);
	ssize_t curPos = 0;
	for (e = 0; e < graph.edges.size(); e++) {
		edgeSequences[e].seq = graph.edges[e].seq.seq;
		edgeSequences[e].length = graph.edges[e].length;
		for (p = 0; p < graph.edges[e].length - dbEdgeLength + 1; p++) {
			positions[curPos].read = e;
			positions[curPos].pos  = p;
			curPos++;
		}
	}
	positions.resize(numPos);
	ReadPos::sequences = &edgeSequences;
	ReadPos::hashLength = dbEdgeLength;
	CompareTuples<SimpleSequenceList> comp;
  comp.sequencesPtr = &edgeSequences;
  comp.length = dbEdgeLength;
  std::sort(positions.begin(), positions.end(), comp);

	DNASequence read;
	ReadPos readPos;
	ssize_t index;
	//UNUSED// ssize_t overlap = vertexSize + 1;
	//UNUSED// ssize_t numMatching;
	ssize_t curEdge, curEdgePos;
	std::vector<ssize_t> edgeList, edgePosList, edgeIntvLengthList;
	ssize_t curIntv;
	ssize_t curRead = nOldReads * 2;
	ssize_t readPosCount;
	ssize_t numAdded = 0;
	while(SeqReader::GetSeq(newReads, read, SeqReader::noConvert)) {
		// Try and map the read prefix to the graph.

		if (read.length >= dbEdgeLength) {
			index = LocateFirstTuple(edgeSequences, positions,
															 dbEdgeLength, (char*) read.seq);
			if (index != -1) {
				readPosCount = CountReadPosMatches(positions, index);
				
				if (readPosCount == 1) {
					++numAdded;
					// The first tuple maps uniquely
					edgeList.clear();
					edgePosList.clear();
					edgeIntvLengthList.clear();
					curEdge = positions[index].read;
					curEdgePos = positions[index].pos;
					edgeList.push_back(curEdge);
					edgePosList.push_back(curEdgePos);
					edgeIntvLengthList.push_back(dbEdgeLength);
					p = 1;
					curIntv = 0;
					// Try and map as much of this read as possible.

					while (readPosCount == 1 and 
								 p < prefixLength - dbEdgeLength and
								 p < read.length - dbEdgeLength+1) {
						index = LocateFirstTuple(edgeSequences, positions,
																		 dbEdgeLength, (char*) &(read.seq[p]));
						if (index != -1) 
							readPosCount = CountReadPosMatches(positions, index);
						else 
							readPosCount = 0;

						if (readPosCount == 1) {
							if (positions[index].read == curEdge and
									positions[index].pos == curEdgePos + 1) {
								edgeIntvLengthList[curIntv]++;
							}
							else {
								curEdge = positions[index].read;
								curEdgePos = positions[index].pos;
								edgeList.push_back(curEdge);
								edgePosList.push_back(curEdgePos);
								edgeIntvLengthList.push_back(dbEdgeLength);
								curIntv++;
							}						
							curEdge = positions[index].read;
							curEdgePos = positions[index].pos;
						}
						p++;
					}
					// Now add this read interval to the graph.
					ssize_t curReadPos = 0;
					ssize_t bEdge;
					ssize_t bEdgePos;
					ssize_t bReadPos;
					ssize_t edge;
					PathInterval *forPath, *revPath;
					forPath = new PathInterval[edgeList.size()];
					revPath = new PathInterval[edgeList.size()];
					graph.paths.push_back(forPath);
					graph.paths.push_back(revPath);
				
					// Add the path to the list of paths.
					graph.pathLengths.push_back(edgeList.size());
					graph.pathLengths.push_back(edgeList.size());
					//				std::cout << "new read: " << curRead << std::endl;
					for (e = 0; e < edgeList.size(); e++) {
						edge = edgeList[e];
						bEdge = graph.edges[edge].balancedEdge;

						
						graph.edges[edge].intervals->push_back(ReadInterval(edge,
																															 curRead,
																															 curReadPos,
																															 edgePosList[e],
																															 edgeIntvLengthList[e],
																															 e));

					

						forPath[e].edge = edge;
						forPath[e].index = graph.edges[edge].intervals->size() -1;
						//					std::cout << edge << " " << graph.edges[e].intervals.size()-1 << " ";

						assert(graph.edges[edge].length == graph.edges[bEdge].length);
					
						bEdgePos = graph.edges[edge].length - 
							(edgePosList[e] + edgeIntvLengthList[e]);
						bReadPos = read.length - (curReadPos + edgeIntvLengthList[e]);
					
						graph.edges[bEdge].intervals->push_back(ReadInterval(bEdge,
																																curRead+1,
																																bReadPos,
																																bEdgePos,
																																edgeIntvLengthList[e],
																																e));
					
						revPath[edgeList.size() - e - 1].edge  = bEdge;
						revPath[edgeList.size() - e - 1].index = graph.edges[bEdge].intervals->size() -1;
						curReadPos += edgeIntvLengthList[e] - vertexSize;
					}

					// Print this sequence to a new file for a record
					// of what was appended.
					read.PrintlnSeq(mappedNewReads);


					// Since this read was mapped, we've added two reads
					// to the global list.
					curRead+=2;
				}
			} // done checking 
		}
	}// end while reading reads.
	std::cout << "added: " << numAdded << " reads." << std::endl;
	WriteIntervalGraph(graphOutBase, graph, vertexSize);
	
	return 0;
}	

