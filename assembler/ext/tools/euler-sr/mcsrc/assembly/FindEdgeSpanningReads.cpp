/***************************************************************************
 * Title:          FindEdgeSpanningReads.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  10/25/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "IntervalGraph.h"
#include "parseblast/BlastParser.h"
#include "parseblast/BlastResult.h"
#include "SeqReader.h"
#include "DNASequence.h"
#include "graph/MSTAlgo.h"

ssize_t ParseEdgeIndex(std::string &edgeName) {
	std::stringstream estrm(edgeName);
	char edgeLit[4];
	if (edgeName.size() > 4) {
		estrm.get(edgeLit, 5);
	}
	ssize_t edgeIndex;
	estrm >> edgeIndex;
	return edgeIndex - 1;
}

class StoreComponentFunctor {
public:
	BEdgeList *edges;
  BVertexList *vertices;
	std::vector<ssize_t> edgeList;
	ssize_t numSources, numSinks;
	void operator()(ssize_t vertexIndex) {
		ssize_t outEdge, outEdgeIndex, inEdge, inEdgeIndex;
		//UNUSED// ssize_t sentinalEdge;
		if ((*vertices)[vertexIndex].InDegree() == 0 and 
				(*vertices)[vertexIndex].OutDegree() == 1)
			++numSinks;
		if ((*vertices)[vertexIndex].InDegree() == 1 and 
				(*vertices)[vertexIndex].OutDegree() == 0)
			++numSources;

		for (inEdgeIndex = (*vertices)[vertexIndex].FirstIn();
				 inEdgeIndex < (*vertices)[vertexIndex].EndIn();
				 inEdgeIndex = (*vertices)[vertexIndex].NextIn(inEdgeIndex)) {
			inEdge = (*vertices)[vertexIndex].in[inEdgeIndex];
			if ((*edges)[inEdge].marked == GraphEdge::NotMarked) {
				edgeList.push_back(inEdge);
      }
		}// end adding in edges
		for (outEdgeIndex = (*vertices)[vertexIndex].FirstOut();
				 outEdgeIndex < (*vertices)[vertexIndex].EndOut();
				 outEdgeIndex = (*vertices)[vertexIndex].NextIn(inEdgeIndex)) {
			outEdge = (*vertices)[vertexIndex].out[outEdgeIndex];
			if ((*edges)[outEdge].marked == GraphEdge::NotMarked) {
				edgeList.push_back(outEdge);
      }
		}// end adding in edges
	}
};
				

int main(int argc, char* argv[]) {

	std::string alignmentName, graphBase, bGraphName, intvName, edgeName, gvzGraphOut;
	ssize_t closeToEnd;
	int argi = 1;
	if (argc != 4) {
		std::cout << "usage: findESRs bGraphFile alignmentFile gvzGraphOut " << std::endl;
		std::cout << "       Find Edge Spanning Reads" << std::endl;
		std::cout << "       alignmentFile is a blast table of read alignments " << std::endl;
		exit(0);
	}
	closeToEnd = 4;
	bGraphName = argv[argi++];
	alignmentName = argv[argi++];
	gvzGraphOut = argv[argi++];
		
	BVertexList vertices;
	BEdgeList   edges;
		
	ReadBGraph(bGraphName, vertices, edges);

	ssize_t v;
	//UNUSED// ssize_t e;
	std::vector<std::set<ssize_t>* > components;
	StoreComponentFunctor storeComponentFn;
	storeComponentFn.vertices = &vertices;
	storeComponentFn.edges    = &edges;
	
	//	std::cout << "storing components " << std::endl;
	Unmark(vertices);
	Unmark(edges);
	std::vector<ssize_t> edgeComponent;
	edgeComponent.resize(edges.size());
	ssize_t componentIndex = 0;
	ssize_t c;
	std::vector<ssize_t> sources, sinks;

	for (v = 0; v < vertices.size(); v++ ) {
		storeComponentFn.edgeList.clear();
		storeComponentFn.numSources = 0;
		storeComponentFn.numSinks   = 0;
		if (vertices[v].marked == GraphVertex::NotMarked) {
			storeComponentFn.edgeList.clear();
			TraverseDFS(vertices, edges, v, storeComponentFn);

			for (c = 0; c < storeComponentFn.edgeList.size(); c++ ) {
				edgeComponent[storeComponentFn.edgeList[c]] = componentIndex;
			}
			sources.push_back(storeComponentFn.numSources);
			sinks.push_back(storeComponentFn.numSinks);
			componentIndex++;
			
		}
	}
	
	std::vector<std::map<ssize_t,ssize_t> > compGraph;
	compGraph.resize(componentIndex);

	//	std::cout << "got " << componentIndex << " components " << std::endl;
	
  BlastResult blastResult;
  
  if (!(ReadBlastTable(alignmentName, blastResult))) {
    std::cout << "Could not parse the blast file " << std::endl;
		exit(0);
  }
	//	std::cout << "parsed: " << blastResult.size() << " blast entries " << std::endl;

	ssize_t i;
	ssize_t numSpanning = 0;
	ssize_t edge1Index, edge2Index;
	ssize_t nextHit;
	ssize_t numSpanningAdjacent = 0; 
	ssize_t numSpanningComponents = 0;
	std::vector<ssize_t> spStart, spEnd;
	std::vector<std::map<ssize_t, ssize_t> > newOutEdges;
	newOutEdges.resize(vertices.size());
	
	for (i = 0; i < blastResult.size(); i++) {
		if (blastResult[i]->hsps.size() > 2) {
			edge1Index = ParseEdgeIndex(blastResult[i]->hsps[0].sbjct);
			/*
			std::cout << "checking " << blastResult[i]->queryName << " " 
								<< blastResult[i]->hsps[0].sbjct << std::endl;
			*/
			// lame code having to deal with palindromic edges
			if (edges[edge1Index].balancedEdge == edge1Index) {
				nextHit = 1;
			}
			else {
				nextHit = 2;
			}
			if (nextHit > blastResult[i]->hsps.size()) {
				std::cout << "WARNING: symmetric graph should produce symmetric results " << std::endl;
			}
			if (blastResult[i]->hsps[0].sbjct != blastResult[i]->hsps[nextHit].sbjct) {
				/*
				std::cout << "found multiple hits: " 
									<< blastResult[i]->hsps[0].sbjct << " "
									<< blastResult[i]->hsps[nextHit].sbjct << std::endl;
				*/
				ssize_t match1Start, match1End, match2Start, match2End;
				ssize_t sbj1Start, sbj1End, sbj2Start, sbj2End;

				edge2Index = ParseEdgeIndex(blastResult[i]->hsps[nextHit].sbjct);
				match1Start = blastResult[i]->hsps[0].qryPos;
				match1End   = blastResult[i]->hsps[0].qryEnd;
				match2Start = blastResult[i]->hsps[nextHit].qryPos;
				match2End   = blastResult[i]->hsps[nextHit].qryEnd;
				
				sbj1Start = blastResult[i]->hsps[0].refPos;
				sbj1End   = blastResult[i]->hsps[0].refEnd;
				sbj2Start = blastResult[i]->hsps[nextHit].refPos;
				sbj2End   = blastResult[i]->hsps[nextHit].refEnd;
				
				ssize_t tmp;
				if (sbj1Start > sbj1End ) {
					tmp = sbj1Start;
					sbj1Start = sbj1End;
					sbj1End  = tmp;
				}
				if (sbj2Start > sbj2End ) {
					tmp = sbj2Start;
					sbj2Start = sbj2End;
					sbj2End  = tmp;
				}

				ssize_t overlap = -1;
				ssize_t distance = -1;
				
				if (match2Start >= match1Start and
						match2Start <= match1End) {
					if (match2End <= match1End) {
						overlap = match2End - match2Start;
					}
					else {
						overlap = match1End - match2Start;
					}
				}
				else if (match2End > match1Start and
								 match2End < match1End)
					overlap = match2End - match1Start;
				else if (match2Start < match1Start and
								 match2End > match1End) {
					overlap = match1End - match1Start;
				}
				else {
					// the two hits do not overlap
					if (match1End < match2Start) {
						distance = match2Start - match1End;
					}
					else {
						distance = match1End - match2End;
					}
				}
				//				std::cout << "distance: " << distance << " overlap: " << overlap << std::endl;
				if ((distance != -1 or overlap != -1) and 
						((distance == -1 and overlap < 20) or
						 (overlap  == -1 and distance < 100))) {
					/*
					if ((vertices[edges[edge1Index].dest].OutDegree() == 1 and
								 vertices[edges[edge2Index].src].InDegree() == 1) or 
								(vertices[edges[edge1Index].src].InDegree() == 1 and
								 vertices[edges[edge2Index].dest].OutDegree() == 1)) {
					*/
						// Now check and see if the overlap is not in the middle of any edges
						if ((sbj1Start < closeToEnd and (edges[edge2Index].length - sbj2End < closeToEnd)) or
								(sbj2Start < closeToEnd and (edges[edge1Index].length - sbj1End < closeToEnd))) {
							
							
							if (sbj1Start < closeToEnd and (edges[edge2Index].length - sbj2End < closeToEnd)) {
							// edge goes from edge2 dest to edge1 src
								if (newOutEdges[edges[edge2Index].dest].find(edges[edge1Index].src) !=
										newOutEdges[edges[edge2Index].dest].end())
									newOutEdges[edges[edge2Index].dest][edges[edge1Index].src]++;
								else
									newOutEdges[edges[edge2Index].dest][edges[edge1Index].src] = 1;
							}
							
							
							if (sbj2Start < closeToEnd and (edges[edge1Index].length - sbj1End < closeToEnd)) {
							// edge goes from edge2 dest to edge1 src
								if (newOutEdges[edges[edge1Index].dest].find(edges[edge2Index].src) !=
										newOutEdges[edges[edge1Index].dest].end())
									newOutEdges[edges[edge1Index].dest][edges[edge2Index].src]++;
								else
									newOutEdges[edges[edge1Index].dest][edges[edge2Index].src] = 1;
								
							}
							/*
								std::cout << blastResult[i]->queryName << " " << match1Start <<" " << match1End 
								<< " " << match2Start << " " << match2End << " "
								<< blastResult[i]->hsps[0].sbjct << " " << sbj1Start << " " << sbj1End << " "
								<< blastResult[i]->hsps[1].sbjct << " " << sbj2Start << " " << sbj2End << " "
											<< edgeComponent[edge1Index] << " " << edgeComponent[edge2Index] << std::endl;
							*/
							numSpanning++;
							if (edgeComponent[edge1Index] != edgeComponent[edge2Index]) {
								numSpanningComponents++;
								if (compGraph[edge1Index].find(edge2Index) == compGraph[edge1Index].end())
									compGraph[edge1Index][edge2Index] = 1;
								else
									compGraph[edge1Index][edge2Index]++;
								
								if (compGraph[edge2Index].find(edge1Index) == compGraph[edge2Index].end())
									compGraph[edge2Index][edge1Index] = 1;
								else
									compGraph[edge2Index][edge1Index]++;
							}
							else {
								// either edge1Index can be reached from edge2, or edge 2 from edge1, check both.
								ssize_t dest;
								dest = edges[edge1Index].dest;
								if (vertices[dest].LookupOutIndex(edge2Index) >= 0) {
									numSpanningAdjacent++;
								}
								else {
									dest = edges[edge2Index].dest;
									if (vertices[dest].LookupOutIndex(edge1Index) >= 0) {
									numSpanningAdjacent++;
									}
								}
								//							}
						}
					}
				}
			}				
		}
	}
	/*
	std::cout << "found " << numSpanning << " spanning edges " << std::endl;
	std::cout << numSpanningComponents << " span different components " << std::endl;
	std::cout << numSpanningAdjacent << " span adjacent edges " << std::endl;
	*/
	std::ofstream compOut;
	openck(gvzGraphOut, compOut, std::ios::out);
	compOut << "digraph G {"  << std::endl
					<< "\tsize=\"8,8\";"<<std::endl;

	ssize_t comp;
	std::string color;
	for (comp = 0; comp < compGraph.size(); comp++ ) {
		std::map<ssize_t,ssize_t>::iterator setIt;
		for (setIt = compGraph[comp].begin(); 
				 setIt != compGraph[comp].end(); ++setIt) {
			if ((*setIt).second < 2) {
				color = "black";
			}
			else if ((*setIt).second >= 2 and (*setIt).second < 4) {
				color = "blue";
			}
			else
				color = "red";
			compOut << comp << " -> " << (*setIt).first << " [color=\""<< color << "\"];" << std::endl;
			//			std::cout << comp << " -> " << (*setIt).first << " " << (*setIt).second << std::endl;
		}
	}

	for (v = 0; v < newOutEdges.size(); v++) { 
		if (newOutEdges[v].size() > 0) {
			std::map<ssize_t, ssize_t>::iterator mapIt;
			for (mapIt = newOutEdges[v].begin();
					 mapIt != newOutEdges[v].end(); ++mapIt) {
				compOut << v << " -> " << (*mapIt).first << " [label=\"" << (*mapIt).second 
									<<  "\" color=\"purple\"];" << std::endl;
			}
		}
	}
	compOut << "};" << std::endl;
	compOut.close();
}
