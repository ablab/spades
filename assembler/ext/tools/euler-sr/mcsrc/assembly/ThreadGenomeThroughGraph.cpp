/***************************************************************************
 * Title:          ThreadGenomeThroughGraph.cpp 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  09/18/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "SeqReader.h"
#include "ReadPos.h"
#include "SeqReader.h"
#include "IntervalGraph.h"
#include "utils.h"
#include "graph/DisjointSet.h"
#include "PrintGraph.h"
#include "IntegralTupleStatic.h"

void ThreadGenome(DNASequence &genome, ssize_t genomePos, 
									IntervalGraph &g, ssize_t curEdge, ssize_t edgePos, 
									std::list<ssize_t> &maxPath, std::list<ssize_t> &curPath, 
									ssize_t &maxDepth, ssize_t maxEdgeLength,
									ssize_t depth= 0, ssize_t curLength = 0);


int main(int argc, char* argv[]) {
	ssize_t MIN_UNIQUE_LENGTH = 300;
	if (argc < 4) {
		std::cout << "usage: threadGenome graph vertexSize genome minEdgeLength" << std::endl;
		exit(0);
	}
	std::string graphName, genomeName;
	ssize_t startEdge;
	int vertexSize;
	
	graphName = argv[1];
	vertexSize = atoi(argv[2]);
	genomeName = argv[3];
	MIN_UNIQUE_LENGTH = atoi(argv[4]);
	
	/*	std::vector<int> startEdges;
	int argi = 4;
	while(argi < argc) {
		startEdges.push_back(atoi(argv[argi]));
		++argi;
	}
	*/
	IntervalGraph graph;
	ReadIntervalGraph(graphName, graph, vertexSize, 1);
	
	DNASequence genome, genomeRC;
	std::ifstream seqIn;
	openck(genomeName, seqIn, std::ios::in);
	std::cout << "reading genome." << std::endl;
	SeqReader::GetSeq(seqIn, genome, SeqReader::noConvert);
	MakeRC(genome, genomeRC);
	std::cout << "done." << std::endl;


	// 
	// Make the genome query-able.
	//
	SimpleSequenceList edges;
	edges.resize(2);
	ssize_t e, p;
	edges[0].seq = genome.seq;
	edges[0].length = genome.length;
	edges[1].seq = genomeRC.seq;
	edges[1].length = genomeRC.length;

	ReadPos::sequences  = &edges;
	ReadPos::hashLength = vertexSize;
	ReadPositions readPosList;
	ReadPositions genomeReadPosList;
	std::cout << "storing pointers into genome." << std::endl;
	for (p = 0; p < genome.length - vertexSize + 1; p++ ) {
		readPosList.push_back(ReadPos(0, p));
	}
	for (p = 0; p < genome.length - vertexSize + 1; p++ ) {
		readPosList.push_back(ReadPos(1, p));
	}
	std::cout << "making genome query-able." << std::endl;
	std::sort(readPosList.begin(), readPosList.end());
	// Now try and locate the end spot on the edge in the genome.
	std::vector<ssize_t> pathLengths;
	std::vector<std::set<ssize_t> > pathStarts;
	std::vector<std::set<ssize_t> >tangles;
	std::vector<std::list<ssize_t> >  paths;
	std::set<ssize_t> tangleEdges;
	pathLengths.resize(graph.edges.size());
	paths.resize(graph.edges.size());
	pathStarts.resize(graph.edges.size());
	for (e =0 ; e < graph.edges.size(); e++ ) {
		pathLengths[e] = 0;
	}

	for (e =0 ; e < graph.edges.size(); e++ ) {
		//		startEdge = startEdges[e];
		if (graph.edges[e].length < MIN_UNIQUE_LENGTH) {
			continue;
		}
		unsigned char *edgePtr;
		startEdge = e;
		ssize_t edgePos = graph.edges[startEdge].length - vertexSize;
		edgePtr     = &graph.edges[startEdge].seq.seq[edgePos];
		ssize_t i;
		std::cout << "looking for: ";
		for (i = 0; i < vertexSize;i++ ){
			std::cout << edgePtr[i];
		}std::cout << std::endl;
		
	
		CompareReadPos readPosCompStr;
		//		std::cout << "searching for match on edge of length: " << edge.length << std::endl;
		ReadPositions::iterator firstIt, lastIt;
		
		firstIt = std::lower_bound(readPosList.begin(), readPosList.end(),
															 (const char*) edgePtr, readPosCompStr );
		
		if (firstIt != readPosList.end()) {
			std::cout << "found last tuple on edge: " << startEdge << std::endl;
			ssize_t strand = (*firstIt).read;
			ssize_t pos    = (*firstIt).pos;
			std::list<ssize_t> path, maxPath;
			ssize_t maxDepth = 0;

			if (strand == 0) {
				ThreadGenome(genome, pos, graph, startEdge, edgePos, maxPath, path, maxDepth, MIN_UNIQUE_LENGTH);
			}
			else {
				ThreadGenome(genomeRC, pos, graph, startEdge, edgePos, maxPath, path, maxDepth, MIN_UNIQUE_LENGTH);
			}

			// Now find the distance from one long edge to another.

			std::list<ssize_t>::iterator listIt;
			ssize_t pathLength = 0;
			
			if (maxPath.size() < 3)
				continue;

			listIt = maxPath.begin();

			ssize_t firstEdgeIndex, lastEdgeIndex;
			ssize_t dest = graph.edges[*listIt].dest;
			firstEdgeIndex = *listIt;
			++listIt;
			ssize_t endFound = 0;
			for(; listIt != maxPath.end(); ++listIt) {
				if (graph.edges[*listIt].length > MIN_UNIQUE_LENGTH) {
					// Found a long edge.
					// Done looking ahead on the path.
					lastEdgeIndex = *listIt;
					endFound = 1;
					break;
				}
				dest = graph.edges[*listIt].dest;
				paths[*listIt].push_back(firstEdgeIndex);
				pathLength += graph.edges[*listIt].length - graph.vertices[dest].vertexSize;
			}
			// Now store the path.
			if (endFound) {
				pathLengths[firstEdgeIndex] = pathLength;
				listIt = maxPath.begin();
				firstEdgeIndex = *listIt;
				++listIt;
				std::cout << "set starting with: " << e 
									<< " (" << graph.edges[e].length  << ") ... " 
									<< lastEdgeIndex << " (" << graph.edges[lastEdgeIndex].length 
									<< ") " << pathLength << " ";
				
				for(; listIt != maxPath.end(); ++listIt) {
					// Store the path that starts at firstEdgeIndex.
					if (graph.edges[*listIt].length > MIN_UNIQUE_LENGTH)
						break;
					paths[firstEdgeIndex].push_back(*listIt);
					pathStarts[*listIt].insert(firstEdgeIndex);
					std::cout << *listIt << " (" << graph.edges[*listIt].length << ") ";
				}
				std::cout << std::endl;
			}
		}
		else {
			std::cout << "not found." << std::endl;
		}
	}
	std::set<ssize_t> finishedStarts;
	std::set<ssize_t> curTangleSet;
	std::set<ssize_t> queue, starts;
	std::set<ssize_t> tangleVertices;
	// For each edge, check to see if a tangle starts at that edge.
	for (e = 0; e < graph.edges.size(); e++ ){ 
		tangleEdges.clear();
		tangleVertices.clear();
		if (pathLengths[e] > 0 and finishedStarts.find(e) == finishedStarts.end()) {
			// A tangle starts at this edge.
			// Traverse the paths
			std::list<ssize_t>::iterator pathIt;
			std::set<ssize_t>::iterator pathStartIt;			
			queue.clear();
			starts.clear();
			queue.insert(e);
			ssize_t curStart;
			while (queue.size() > 0) {
				curStart = *queue.begin();
				queue.erase(curStart);
				starts.insert(curStart);
				for (pathIt = paths[curStart].begin(); 
						 pathIt != paths[curStart].end(); ++pathIt) {
					tangleEdges.insert(*pathIt);

					tangleVertices.insert(graph.edges[*pathIt].src);
					tangleVertices.insert(graph.edges[*pathIt].dest);
					// look to see if the start corresponding to the path it
					// needs to be searched or not.
					for (pathStartIt = pathStarts[*pathIt].begin();
							 pathStartIt != pathStarts[*pathIt].end();
							 ++pathStartIt) {
						if (starts.find(*pathStartIt) == starts.end()) {
							if (queue.find(*pathStartIt) == queue.end()) {
								queue.insert(*pathStartIt);
							}
						}
					}
				}
			}
			tangles.push_back(tangleVertices);
			std::cout << "tangle of size: " << starts.size() << std::endl;
			for (pathStartIt = starts.begin();
					 pathStartIt != starts.end();
					 ++pathStartIt) {
				assert(*pathStartIt >= 0 and *pathStartIt < graph.edges.size());
				std::cout << *pathStartIt << " (" << graph.edges[*pathStartIt].length << "," 
									<< pathLengths[*pathStartIt] << ") ";
				finishedStarts.insert(*pathStartIt);
			}
			std::cout << std::endl;
			std::cout << "tangle edges are: " << std::endl;
			std::set<ssize_t>::iterator setIt, setEnd;
			for (setIt = tangleEdges.begin(),
						 setEnd = tangleEdges.end();
					 setIt != setEnd;
					 ++setIt) {
				std::cout << *setIt << " (" << graph.edges[*setIt].length << ") ";
			}
			std::cout << std::endl;

		}
	}
	
	std::stringstream tangleNameStream;
	
	std::string tangleName;
	ssize_t t;
	for (t = 0; t < tangles.size(); t++) {
		if (tangles[t].size() > 2) {
			std::ofstream tangleOut;
			tangleName ="";
			tangleNameStream.str(tangleName);
			tangleNameStream << "tangle." << t << ".dot";
			openck(tangleNameStream.str(), tangleOut, std::ios::out);
			PrintGraphSubset(graph, tangles[t], MIN_UNIQUE_LENGTH, tangleOut);
			tangleOut.close();
			tangleOut.clear();
		}
	}

	return 0;
}

void ThreadGenome(DNASequence &genome, ssize_t genomePos, 
									IntervalGraph &g, ssize_t curEdge, ssize_t edgePos, 
									std::list<ssize_t> &maxPath, std::list<ssize_t> &curPath,
									ssize_t &maxDepth, ssize_t maxEdgeLength, ssize_t depth, ssize_t curLength) {

	if (depth > maxDepth) {
		maxDepth = depth;
		maxPath  = curPath;
	}
	curPath.push_back(curEdge);

	ssize_t edgeLength = g.edges[curEdge].length - edgePos;
	// adjust for edges that go past the end of the genome, possible
	// with bacteria.
	if (edgeLength + genomePos > genome.length) {
		edgeLength -= edgeLength + genomePos - genome.length;
	}
	
	ssize_t e;
	ssize_t match = 1;

	for (e = 0; e < edgeLength and match; e++ ){ 
		assert(e + genomePos < genome.length);
		if (genome.seq[e + genomePos] != g.edges[curEdge].seq.seq[edgePos + e])
			match = 0;
	}
	//UNUSED// ssize_t extensionFound = 0;
	if (match) {
		//		std::cout << curEdge << " (" << depth << "," << curLength << "," << edgeLength << ") ";
		ssize_t outEdge, outEdgeIndex;
		ssize_t dest = g.edges[curEdge].dest;
		for (outEdgeIndex = g.vertices[dest].FirstOut();
				 outEdgeIndex != g.vertices[dest].EndOut();
				 outEdgeIndex = g.vertices[dest].NextOut(outEdgeIndex)) {
			outEdge = g.vertices[dest].out[outEdgeIndex];
			//			if (g.edges[outEdge].length < maxEdgeLength)
				ThreadGenome(genome, genomePos + edgeLength - g.vertices[dest].vertexSize, g, outEdge, 0, 
										 maxPath, curPath, maxDepth, maxEdgeLength, depth+1, 
										 curLength + edgeLength - g.vertices[dest].vertexSize);
		}
	}
	curPath.pop_back();
}
