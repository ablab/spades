/***************************************************************************
 * Title:          JoinSourcesAndSinks.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  09/28/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
using namespace std;

#include "IntervalGraph.h"
#include "DeBruijnGraph.h"
#include "SimpleSequence.h"
#include "mctypes.h"
#include "ContigMap.h"
#include "MateTable.h"
#include "MateLibrary.h"
#include <map>
#include "IntegralTupleStatic.h"
#include "Scaffold.h"
#include "compatibility.h"


#define M_INCORRECT 0
#define M_CORRECT   1
#define M_AMBIGUOUS 2


ssize_t EdgesOverlap(SimpleSequence &sourceSeq,
								 SimpleSequence &sinkSeq,
								 ssize_t maxOverlap, ssize_t minOverlap);



void PrintUsage() {
	std::cout << "usage: joinsas  in_base min_overlap out_base [options]" << std::endl;
	std::cout << "  -vertexSize v  Set the vertex size of the graph. " << std::endl;
	std::cout << "  -maxEdit e     Find match with maximum edit distance e, 1 mismatch/delete" << std::endl;
	std::cout << "  -useMap mapFile Use a file to compare the joined contigs with a reference."<<std::endl
						<< "                 for benchmarking." << std::endl;
	std::cout << "  -mateFile file Read mate information from 'file'." << std::endl;
	std::cout << "  -mateRules file Read mate rules from 'file'." << std::endl;
	std::cout << "  -minMateCount c Only join edges if there are at least 'C' mates joining them." 
						<< std::endl;
	std::cout << "  -minScaffoldOverlap 'S' Only join edges if they are connected by 'c' mates" << std::endl
						<< "                      and the two sequences share at least 'S' overlapping nucleotides." << std::endl;
	// TODO: document "-mateType"
}


int main(int argc, char* argv[]) {

	if (argc < 4) {
		PrintUsage();
		exit(1);
	}

	std::string baseIn, baseOut;
	ssize_t minOverlap;
	ssize_t maxEdit;

	int argi = 1;
	baseIn = argv[argi++];
	minOverlap = atoi(argv[argi++]);
	baseOut  = argv[argi++];

	int vertexSize = 20;
	ssize_t useMap = 0;
	ssize_t useMatePairs = 0;
	ssize_t useRuleFile = 0;
	ssize_t minMateCount = 5;
	ssize_t minScaffoldOverlap = 5;
	ssize_t mateType = 0;
	maxEdit = 1;
	std::string mapFileName;
	std::string mateFileName, mateRuleFileName;
	
	while (argi < argc) {
		
		if (strcmp(argv[argi], "-vertexSize") == 0) {
			vertexSize = atoi(argv[++argi]);
		}
		else if (strcmp(argv[argi], "-maxEdit" ) == 0){ 
			maxEdit = atoi(argv[++argi]);
		}
		else if (strcmp(argv[argi], "-mateType" ) == 0){ 
			mateType = atoi(argv[++argi]);
		}
		else if (strcmp(argv[argi], "-useMap") == 0){
			mapFileName = argv[++argi];
			useMap = 1;
		}
		else if (strcmp(argv[argi], "-mateFile") == 0){
			mateFileName = argv[++argi];
			useMatePairs = 1;
		}
		else if (strcmp(argv[argi], "-mateRules") == 0){
			mateRuleFileName = argv[++argi];
			useRuleFile = 1;
		}
		else if (strcmp(argv[argi], "-minMateCount") == 0) {
			minMateCount = atoi(argv[++argi]);
		}
		else if (strcmp(argv[argi], "-minScaffoldOverlap") == 0) {
			minScaffoldOverlap = atoi(argv[++argi]);
		}
		else {
			PrintUsage();
			std::cout << "bad option: " << argv[argi] << std::endl;
			exit(1);
		}
		++argi;
	}
	
	if (useMatePairs xor useRuleFile) {
		std::cout << "You must specify both a mate and rule file to use mate pairs." 
							<< std::endl;
		exit(1);
	}
	ReadMateList mateList;
	if (useMatePairs) {
		ReadMateTable(mateFileName, mateList);
	}

	RuleList rules;
	if (useRuleFile) {
		ParseRuleFile(mateRuleFileName, rules);
	}

	IntervalGraph graph;
	ReadIntervalGraph(baseIn, graph, vertexSize);
	cout << "vertex size: " << vertexSize << endl;

	std::string bGraphOutName, intvOutName, edgeOutName, pathOutName, graphOutName, altEdgeName;
	FormGraphFileNames(baseOut, bGraphOutName, graphOutName, intvOutName, 
										 pathOutName, edgeOutName, altEdgeName);


	std::vector<ssize_t> sourceVertices, sourceEdges, 
		sinkVertices, sinkEdges;

	// create a sparse matrix of source/sink matches
	std::vector<std::vector<ssize_t> > sinkMatch, sourceMatch;
	

	std::vector<std::map<ssize_t, MateCount> > endToBegin, sourceToSinkMateCount;
	vector<ssize_t> sourcePairedEdgeCount;
	sourcePairedEdgeCount.resize(graph.edges.size());
	std::fill(sourcePairedEdgeCount.begin(), sourcePairedEdgeCount.end(), 0);

	// IF a map is specified, read it in.
	std::vector<ContigMap> contigMap;
	if (useMap) {
		std::ifstream mapFile;
		openck(mapFileName,mapFile);
		ReadUniqueContigMap(mapFile, contigMap);
	}

	// count the number of sources

	ssize_t v;
	IntMatrix scoreMat, pathMat;
	
	IntMatrix matchMat;
	CreateMatrix(matchMat, 4, 4);
	
	matchMat[0][0] = 1; matchMat[1][1] = 1; matchMat[2][2] = 1; matchMat[3][3] = 1;
	
	matchMat[0][1] = -1; matchMat[0][2] = -1; matchMat[0][3] = -1;
	matchMat[1][0] = -1; matchMat[1][2] = -1;	matchMat[1][3] = -1;
	matchMat[2][0] = -1; matchMat[2][1] = -1; matchMat[2][3] = -1;
	matchMat[3][0] = -1; matchMat[3][1] = -1;	matchMat[3][2] = -1;

	for (v=  0; v < graph.vertices.size(); v++) {
		if (graph.vertices[v].InDegree() == 0 and
				graph.vertices[v].OutDegree() == 1) {
			sourceVertices.push_back(v);
			sourceEdges.push_back(graph.vertices[v].out[graph.vertices[v].FirstOut()]);
		}
		else if (graph.vertices[v].InDegree() == 1 and
						 graph.vertices[v].OutDegree() == 0) {
			sinkVertices.push_back(v);
			sinkEdges.push_back(graph.vertices[v].in[graph.vertices[v].FirstIn()]);
		}
	}


	sinkMatch.resize(sourceVertices.size());
	sourceMatch.resize(sinkVertices.size());

	sourceToSinkMateCount.resize(graph.edges.size());
	std::vector<ssize_t> ovpLengths;
	ovpLengths.resize(sourceVertices.size());

  ssize_t sourceEdge, sinkEdge;
	ssize_t ovp;
	std::vector<ssize_t> vToRemove, eToRemove;
	ssize_t so, si;

	//UNUSED// ssize_t numCorrect = 0;
	//UNUSED// ssize_t numMisjoined = 0;
	//UNUSED// ssize_t numAmbiguous = 0;

	ssize_t perfectOvp, alignedOvp;
	ssize_t perfect[3], aligned[3], aNotP[3];
	ssize_t numPerfect, numAligned;
	numPerfect = numAligned = 0;
	perfect[0] = perfect[1] = perfect[2] = 0;
	aligned[0] = aligned[1] = aligned[2] = 0;
	aNotP[0] = aNotP[1] = aNotP[2] = 0;
	//UNUSED// ssize_t intv;
	//UNUSED// ssize_t readIndex, mateIndex, pathIndex, pathPos;
	//UNUSED// ssize_t matePathIndex;
	//UNUSED// ssize_t isRC;
	if (useMatePairs) {
		MatePairScaffoldJoinEdges(graph, mateList, mateType, minMateCount, matchMat, vToRemove);
	} // done using mate-pairs
	else {
		for (sourceEdge = 0; sourceEdge < sourceEdges.size(); sourceEdge++ ){
			for (sinkEdge = 0; sinkEdge < sinkEdges.size(); sinkEdge++) { 
				ssize_t maxScoreSourcePos, maxScoreSinkPos, maxScore;
				so = sourceEdges[sourceEdge];
				si = sinkEdges[sinkEdge];

				SimpleSequence source, sink;
				source.seq    = (unsigned char*) graph.edges[so].seq.seq;
				source.length = graph.vertices[graph.edges[so].src].vertexSize - 1;
				sink.seq      = &(graph.edges[si].seq.seq[graph.edges[si].length - 
																									graph.vertices[graph.edges[si].dest].vertexSize +1]);
				sink.length = graph.vertices[graph.edges[si].dest].vertexSize - 1;
				//UNUSED// int sinkVertexSize = graph.vertices[graph.edges[si].dest].vertexSize;

				ssize_t nMisMatch, nIndel;
				GrowMatrices(scoreMat, pathMat, source.length+1, sink.length+1);
			
				PrefixSuffixAlign(source, sink, maxScoreSourcePos, maxScoreSinkPos, maxScore, nMisMatch, nIndel,
													scoreMat, pathMat, matchMat);
			
				
				perfectOvp = alignedOvp = 0;
				// Do two tests for overlap, in the future I'll just do one,
				// but they are being compared for now.
				if ((ovp = EdgesOverlap(graph.edges[sourceEdges[sourceEdge]].seq,
																graph.edges[sinkEdges[sinkEdge]].seq,
																vertexSize, minOverlap))) {
					if (ovp >= minOverlap) {
						numPerfect++;
						perfectOvp = 1;
					}
				}

				if (maxScore  > minOverlap and nIndel + nMisMatch <= maxEdit) {
					numAligned++;
					alignedOvp = 1;
					ovp = maxScoreSourcePos;
				}
				ssize_t matchType = 0;
				if (useMap) {
					// This code is for validation purposes.  The map maps contigs back to the genome
					// so that we can make sure the correct sources and sinks have been joined. 
					ssize_t sinkEnd = contigMap[si].refEnd;
					ssize_t sourceBegin =  contigMap[so].refPos;
				
					if (sinkEnd == 0 or  sourceBegin == 0){ 
						matchType = 2;
					}
					else if (abs(sourceBegin - sinkEnd) < 100) {
						matchType = 1;
					}
					else 
						matchType = 0;
				}
			
				if (perfectOvp) {
					perfect[matchType]++;
				}
				if (alignedOvp) {
					aligned[matchType]++;
				}
			
				if (alignedOvp and !perfectOvp) {
					aNotP[matchType]++;
				}
				if (perfectOvp || alignedOvp) {
					sinkMatch[sourceEdge].push_back(sinkEdge);
					sourceMatch[sinkEdge].push_back(sourceEdge);
					ovpLengths[sourceEdge] = ovp;
				}
			}
		}
		//UNUSED// ssize_t m;
			
		//UNUSED// ssize_t source;
		//UNUSED// ssize_t middleVertex;
		ssize_t sinkLength;
		//UNUSED// ssize_t sourceEdgeIndex;
		ssize_t sinkVertex, sourceVertex;
		// merge edges that have an appropriate unique overlap
		//	for (sourceEdgeIndex = 0; sourceEdgeIndex < sourceEdges.size(); sourceEdgeIndex++ ){
		for (si = 0; si < sinkEdges.size(); si++) {
			// the sink has one match,
			sinkEdge = sourceEdge = -1;
			ssize_t overlapLength;
			if (sourceMatch[si].size() == 1 and 
					// the match to the sink just has one match
					sinkMatch[sourceMatch[si][0]].size() == 1) { 
				sinkEdge   = sinkEdges[si];
				sourceEdge = sourceEdges[sourceMatch[si][0]];
				overlapLength = ovpLengths[sourceMatch[si][0]];

				// found two edges to connect.

				sinkLength = graph.edges[sinkEdge].length;

				sourceVertex = graph.edges[sourceEdge].src;
				sinkVertex   = graph.edges[sinkEdge].dest;

				// make sure the balanced edges are also joined.

				ssize_t sinkBal, sourceBal;
				sinkBal   = graph.edges[sinkEdge].balancedEdge;
				sourceBal = graph.edges[sourceEdge].balancedEdge;

				ssize_t soe; // sink balanced edge
				// find the edge index of the balanced source
				ssize_t sinkBalIndex = SSIZE_MAX;
				//UNUSED// ssize_t sourceBalIndex;
				for (soe = 0; soe < sourceEdges.size(); soe++) {
					if (sourceEdges[soe] == sinkBal) {
						sinkBalIndex = soe;
						break; 
					}
				} 
				if (sinkBalIndex >= sourceEdges.size()) {
					std::cout << "ERROR! No balanced edge found!" << std::endl;
					exit(0);
				}

				//UNUSED// ssize_t be;
				//UNUSED// ssize_t beFound = 0;

				// Join the sink edge with the source edge.
       
				ssize_t mergedSinkEdge, mergedSourceEdge;
				ssize_t mergedSinkEdgeIndex, mergedSourceEdgeIndex;
				mergedSinkEdgeIndex   = graph.vertices[sinkVertex].FirstIn();
				mergedSinkEdge        = graph.vertices[sinkVertex].in[mergedSinkEdgeIndex];

				mergedSourceEdgeIndex = graph.vertices[sourceVertex].FirstOut();
				mergedSourceEdge      = graph.vertices[sourceVertex].out[mergedSourceEdgeIndex];
			
				// Unlink the sink edge from its vertex
				graph.vertices[sinkVertex].in[mergedSinkEdgeIndex] = -1;

				vToRemove.push_back(sinkVertex);
				// Link the sink edge to the source edge.  This will creat a simple
				// path that will be reduced to one edge later.
			
				assert(graph.vertices[sourceVertex].InDegree() == 0);
				// The sourceVertex did not have any in-edges, so the slot 0 is open.
				graph.vertices[sourceVertex].in[0] = sinkEdge;
				graph.edges[sinkEdge].dest         = sourceVertex;

				// Fix the duplicated sequence.
				graph.edges[sinkEdge].length     -= overlapLength;
				graph.edges[sinkEdge].seq.length -= overlapLength;
				std::cout << "joined edges to make new one of length: " <<
					graph.edges[sinkEdge].length + graph.edges[sourceEdge].length << std::endl;
			} // end checking to see if there is a unique overlap
		}
	}
	assert(CheckEdges(graph.vertices, graph.edges));
	std::cout << "pruning " << vToRemove.size() << " vertices and " 
						<< eToRemove.size() << " edges " << std::endl;
		
	graph.Prune(vToRemove, eToRemove);
	graph.CondenseSimplePaths();
	graph.PrintIntervalGraph(bGraphOutName, intvOutName);
	PrintGraph(graph.vertices,graph.edges, graphOutName);
	PrintEdges(graph.vertices, graph.edges, edgeOutName);
	WriteReadPaths(pathOutName, graph.paths, graph.pathLengths);

	return 0;
}


		

ssize_t EdgesOverlap(SimpleSequence &sourceSeq,
								 SimpleSequence &sinkSeq,
								 ssize_t maxOverlap, ssize_t minOverlap) {
	ssize_t ovp;
	ssize_t p;
	ssize_t overlapFound = 0;
	DNASequence tmpSrc, tmpSink;
	tmpSrc.seq = sourceSeq.seq;
	tmpSrc.length = maxOverlap;
	tmpSink.seq = &(sinkSeq.seq[sinkSeq.length - maxOverlap - 1]);
	tmpSink.length = maxOverlap;
	for (ovp = minOverlap; !overlapFound and ovp < maxOverlap ; overlapFound or ovp++) {
		overlapFound = 1;
		for (p = 0; p < ovp; p++ ) {
			if (sourceSeq.seq[p] !=
					sinkSeq.seq[sinkSeq.length - ovp + p]) {
				overlapFound = 0;
				break;
			}
		}
	}

	if (overlapFound)
		return ovp;
	else
		return 0;
}

	
						 
