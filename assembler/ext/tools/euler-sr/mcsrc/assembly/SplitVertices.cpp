/***************************************************************************
 * Title:          SplitVertices.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  01/19/2010
 *
 * Copyright (c) 2007-2010 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "IntervalGraph.h"
#include "DeBruijnGraph.h"
#include "SimpleSequence.h"
#include "graph/DisjointSet.h"
#include <iterator>
#include "IntegralTupleStatic.h"


void PrintUsage() {
	std::cout << "splitVertices - Split all internal vertices into two," 
						<< std::endl;
	std::cout << "                connect all new vertices with an edge," 
						<< std::endl;
	std::cout << "                make each vertex represent a sequence of size 1,"
						<<std::endl;
	std::cout << "                and update all intervals and edges accordingly."
						<<std::endl;
	std::cout << "usage: splitVertices graphIn  vertexSize graphOut" 
						<< std::endl;
}

ssize_t GetOverlap(TVertex &vSrc, TEdge &edge, TVertex &vDest);

typedef std::list<ssize_t> SplitPosList;

class VertexInterval {
public:
	ssize_t vertex;
	ssize_t index;
	ssize_t length;
	ssize_t pos;
	ssize_t intvVertex;
	unsigned char*seqPtr;
	VertexInterval(ssize_t v, ssize_t i, ssize_t p, ssize_t l) {
		vertex = v; index = i; pos = p; length = l;
		intvVertex = -1;
	}
};


// Define a Disjoint Vertex Interval. This is an interval on a vertex 
// which is a set of k elements where k is the vertex size.  They are 
// turned into a disjoint set so that connected components may be found
// after splitting vertices.
typedef DisjointSet::DJVertex<VertexInterval> DJVertexInterval;


// TODO: Determine whether EdgeCount vectors take a lot of space, and if
// we can pack the various fields into bitfields.

class EdgeCount {
public:
	ssize_t outDegree;
	ssize_t inDegree;
	ssize_t index;
	int vertexSize;
	// HAS BITFIELD
	bool traversed:1;
	unsigned char *seqPtr;
	EdgeCount() {
		outDegree  = inDegree = 0;
		index      = -1;
		vertexSize = 1;
		seqPtr     = NULL;
	}
};

typedef DisjointSet::DJVertex<EdgeCount> DJEdgeCount;
typedef std::vector<DJEdgeCount> DJEdgeCountVertex;
typedef std::vector<DJEdgeCountVertex> DJEdgeCountVertexList;

typedef DisjointSet::DJVertex<TVertex> DJVertex;
typedef std::vector<DJVertex> SplitDJVertex;
typedef std::vector<SplitDJVertex> SplitDJVertexList;


// Each vertex may be represented as a list of vertex intervals.
typedef std::vector<DJVertexInterval> DJVertexIntervals;

// All split vertices are stored as a list of lists of vertex intervals.
typedef std::vector<DJVertexIntervals> DJVertexIntervalsList;


void DivideReadInterval(IntervalGraph &graph,
												TVertexList &vertices,
												TEdgeList   &edges,
												TVertex &vertex, TEdge &edge,
												ssize_t intervalEdge, ssize_t intervalIndex,
												ssize_t vertexIsSrc, 
												DJVertexIntervals &vertexIntervals);

ssize_t AddSplit(SplitPosList &splits, ssize_t pos) {
	ssize_t index = 0;
	SplitPosList::iterator splitIt, prevIt;
	splitIt = splits.begin();
	if (splitIt != splits.end() and *splitIt == pos) 
		return index;
	
	for (prevIt = splitIt, ++index;  // Start on the 1st pos
			 splitIt != splits.end(); splitIt++, prevIt++) {
		if (*splitIt == pos)  
			// Signal that no splits were necessary
			return 0;
		if (*splitIt > pos) {
			break;
		}
	}
	// Pos is greater than everything in the list,
	// append it.
	splits.insert(splitIt, pos);
	// Signal that a split was added
	return 1;
}


ssize_t UnionSplitLists(SplitPosList &src, SplitPosList &dest, 
										 int vertexSize, ssize_t overlap) {
	assert(overlap < vertexSize);

	// Union the list in src to dest.
	SplitPosList::iterator destIt, srcIt;
	// Add the src splits to the dest

	ssize_t srcOriginalSize, destOriginalSize;
	srcOriginalSize = src.size();
	destOriginalSize = dest.size();

	// Add the splits that are in the source to the dest.
	
	// Skip splits that are in the src that do not overlap
	// with the dest.

	std::cout << "BEFORE dest split pos list:" << std::endl;
	std::copy(dest.begin(), dest.end(), std::ostream_iterator<ssize_t>(std::cout, "d, "));
	std::cout << std::endl;
	std::cout << "BEFORE src split pos list:" << std::endl;
	std::copy(src.begin(), src.end(), std::ostream_iterator<ssize_t>(std::cout, "s, "));
	std::cout << std::endl;

	for (srcIt = src.begin(); srcIt != src.end() and *srcIt < vertexSize - overlap; ++srcIt) ;
	std::cout << "num src: " << src.size() << std::endl;
	for (; srcIt != src.end(); ++srcIt) {
		AddSplit(dest, *srcIt - (vertexSize - overlap));
		std::cout << *srcIt - (vertexSize - overlap) << std::endl;
	}
	
	// Add the dest splits to the src;
	for(destIt = dest.begin(); destIt != dest.end() and *destIt < overlap - 1; ++destIt) {
		AddSplit(src, *destIt +  vertexSize - overlap );
		std::cout << *destIt + vertexSize - overlap << std::endl;
	}

	// Add the splits caused by the overlap of the vertices
	// Spilts happen after the index, pos - 1.

	AddSplit(src, vertexSize - overlap - 1);

	AddSplit(dest, overlap - 1);

	std::cout << "AFTER dest split pos list:" << std::endl;
	std::copy(dest.begin(), dest.end(), std::ostream_iterator<ssize_t>(std::cout, "d, "));
	std::cout << std::endl;
	std::cout << "AFTER src split pos list:" << std::endl;
	std::copy(src.begin(), src.end(), std::ostream_iterator<ssize_t>(std::cout, "s, "));
	std::cout << std::endl;

	return (srcOriginalSize < src.size() or
					destOriginalSize < dest.size());
}

void StoreVertexBalance(TVertexList &vertices, TEdgeList &edges, 
												std::vector<ssize_t> &vertexBalance) {
	//UNUSED+// ssize_t  v;
	ssize_t e, b;

	//UNUSED// ssize_t srcBal, destBal;
	vertexBalance.resize(vertices.size());
	for (e = 0; e < edges.size(); e++) {
		b = edges[e].balancedEdge;
		vertexBalance[edges[e].src] = edges[b].dest;
		vertexBalance[edges[e].dest] = edges[b].src;
	}
}

void SumVertexDegrees(DJEdgeCountVertexList &splitVertices) {
	ssize_t v, sv;
	for (v = 0; v < splitVertices.size(); v++ ){
		for (sv = 0; sv < splitVertices[v].size(); sv++) {
			DisjointSet::Find(&(splitVertices[v][sv]));
			if (splitVertices[v][sv].parent != &(splitVertices[v][sv])) {
				splitVertices[v][sv].parent->value->inDegree +=
					splitVertices[v][sv].value->inDegree;
			
				splitVertices[v][sv].parent->value->outDegree +=
					splitVertices[v][sv].value->outDegree;
			}
		}
	}
	for (v = 0; v < splitVertices.size(); v++ ){
		for (sv = 0; sv < splitVertices[v].size(); sv++) {
			splitVertices[v][sv].value->inDegree =
				splitVertices[v][sv].parent->value->inDegree;
			splitVertices[v][sv].value->outDegree =
				splitVertices[v][sv].parent->value->outDegree;
		}
	}
}

void GlueOverlappingVertices(TVertexList       &vertices,
														 TEdgeList         &edges,
														 int               vertexSize,
														 DJEdgeCountVertexList &splitVertices) {
	//UNUSED// ssize_t joinDone = 1;
	ssize_t vSrc, e, ei;
	ssize_t prevSrcI, prevDestI;
	for (vSrc = 0; vSrc < vertices.size(); vSrc++ ){
		for (ei = vertices[vSrc].FirstOut();
				 ei < vertices[vSrc].EndOut();
				 ei = vertices[vSrc].NextOut(ei)) {
			e = vertices[vSrc].out[ei];
			ssize_t vDest = edges[e].dest;
			ssize_t overlap;
			overlap = GetOverlap(vertices[vSrc], edges[e], vertices[vDest]);
			if (overlap > 0) {
				ssize_t ovpi;
				ssize_t srcStart = vertexSize - overlap;
				DJEdgeCount *srcBeginParent, *destBeginParent, *srcEndParent, *destEndParent;
				//UNUSED// DJEdgeCount *beginParent, *endParent;
				// Record the parents so that we know which vertex got merged 
				// into which other.
				srcBeginParent  = DisjointSet::Find(&splitVertices[vSrc][srcStart]);
				destBeginParent = DisjointSet::Find(&splitVertices[vDest][0]);
				srcEndParent    = DisjointSet::Find(&splitVertices[vSrc][srcStart + overlap-1]);
				destEndParent   = DisjointSet::Find(&splitVertices[vDest][overlap-1]);

				for (ovpi = 0; ovpi < overlap; ovpi++) {
					prevSrcI  = splitVertices[vSrc][ovpi + srcStart].value->index;
					prevDestI = splitVertices[vDest][ovpi].value->index;
					
					DisjointSet::Union(&(splitVertices[vSrc][ovpi + srcStart]), 
														 &(splitVertices[vDest][ovpi]));
					/*
						std::cout << "joined: " << vSrc << " " << ovpi + srcStart << " " 
										<< vDest << " " << ovpi << " "
										<< prevSrcI << ", " << prevDestI 
										<< " to " << splitVertices[vDest][ovpi].parent->value->index
										<< std::endl;
					*/
				}

				/*
					beginParent = splitVertices[vDest][0].parent;
					endParent   = splitVertices[vSrc][vertexSize - 1].parent;
				*/
				
				//				splitVertices[vDest][0].value->inDegree
				
				// Fix the edges incident to the dest vertex.
				//UNUSED// ssize_t inEdge, inEdgeIndex;
				//UNUSED// DJEdgeCount *altBeginParent, *altEndParent;

				/*				beginParent->value->inDegree = srcBeginParent->value->inDegree 
					+ destBeginParent->value->inDegree;
				*/
				//beginParent->value->inDegree++;
				//endParent->value->outDegree++;

				// Fix the edges exiting the src vertex
				/*
				endParent->value->outDegree = srcEndParent->value->outDegree +
					destEndParent->value->outDegree;
				*/
				// adjust the degrees of the vertices adjacent to the newly joined ones.
				if (srcStart > 0) {
					splitVertices[vSrc][srcStart-1].value->outDegree++;
				}
				if (overlap < vertexSize) {
					splitVertices[vDest][overlap].value->inDegree++;
				}
				Pause();
			}
		}			
	}	 
}
														 

void SplitOverlappingVertices(TVertexList &vertices,
															TEdgeList   &edges,
															int vertexSize, 
															std::vector<SplitPosList> &splitPositions) {

	ssize_t splitDone = 1;
	ssize_t v, outEdge, outEdgeIndex;
	//UNUSED// ssize_t overlapLength;
	ssize_t unionMadeSplit;
	ssize_t splitIter = 0;
	ssize_t numSplitVertices = 0;
	ssize_t overlap;
	ssize_t dest;
	ssize_t numSplits = 0;
	while (splitDone) {
		splitDone = 0;
		numSplitVertices = 0;
		for (v = 0; v < vertices.size(); v++ ) { 
			for (outEdgeIndex = vertices[v].FirstOut();
					 outEdgeIndex != vertices[v].EndOut();
					 outEdgeIndex = vertices[v].NextOut(outEdgeIndex)) {
				outEdge = vertices[v].out[outEdgeIndex];
				dest    = edges[outEdge].dest;
				overlap = GetOverlap(vertices[v], edges[outEdge], vertices[dest]);
        if (overlap > 0) {
					// The vertices for the edges overlap.
					std::cout << v << " -> " << dest << ", overlap:"  << overlap << std::endl;
					std::cout << outEdge <<" " << edges.size() << " " 
										<< edges[outEdge].dest << " " << splitPositions.size() << std::endl;
					unionMadeSplit = UnionSplitLists(splitPositions[v], 
																					 splitPositions[edges[outEdge].dest],
																					 vertexSize, overlap);
					// Look to see if anything has changed.
					splitDone |= unionMadeSplit;
					numSplitVertices += 2;
					numSplits += splitPositions[v].size() + splitPositions[edges[outEdge].dest].size();
				}
			}
		}
		std::cout << "Split overlapping vertices iter: " << splitIter 
							<< " split: " << numSplitVertices << " " 
							<< " num intvs: " << numSplits << std::endl;
		++splitIter;
	}
}


void InitVertexIntervals(TVertexList &vertices,
												 std::vector<SplitPosList> &splitPositions,
												 DJVertexIntervalsList &vertexIntervals) {
	assert(splitPositions.size() == vertices.size());
	
	vertexIntervals.resize(splitPositions.size());

	// Turn the lists of split positions into vectors of intervals.
	ssize_t v;
	ssize_t index;
	ssize_t numSplits;
	ssize_t intervalLength;
	ssize_t intervalPos;
	//UNUSED// ssize_t prevIntervalEnd;
	for (v = 0; v < vertexIntervals.size(); v++) {
		if (!splitPositions[v].empty()) {
			numSplits = splitPositions[v].size();
			if (numSplits> 0) {
				
				vertexIntervals[v].resize(numSplits + 1);
				SplitPosList::iterator splitIt, nextIt, splitEnd;
				splitEnd = splitPositions[v].end();

				nextIt = splitIt  = splitPositions[v].begin();
				++nextIt;
				index = 0;
				// All intervals begin at the beginning of the vertex
				intervalPos    = 0;
				// The intervals in split positions are addressed *after*
				// the split.
				intervalLength = *splitIt + 1;

				for (splitIt = splitPositions[v].begin(); splitIt != splitEnd; ++splitIt, ++index ) {
					// splitIt is the end of an interval, so the 
					// length is end - beginning + 1.
					intervalLength = *splitIt - intervalPos + 1;
					assert(intervalLength > 0);
					vertexIntervals[v][index].parent = &(vertexIntervals[v][index]);
					vertexIntervals[v][index].value  = new VertexInterval(v,index, 
																																	intervalPos, 
																																	intervalLength);
					intervalPos = *splitIt + 1;
				}
				// Add an interval for the last split.
				vertexIntervals[v][index].parent = &(vertexIntervals[v][index]);
				vertexIntervals[v][index].value = new VertexInterval(v, index, intervalPos, 
																														 vertices[v].vertexSize - intervalPos);
			}
		}
	}
}


ssize_t GetOverlap(TVertex &vSrc, TEdge &edge, TVertex &vDest) {
	return (vSrc.vertexSize + vDest.vertexSize - edge.length);
}

void StoreIntervalIndices(TVertexList &vertices,
													TEdgeList   &edges,
													std::vector<SplitPosList> &splitPositions,
													DJVertexIntervalsList &vertexIntervals) {

	InitVertexIntervals(vertices, splitPositions, vertexIntervals);
	std::cout << "vertex intervals: " << std::endl;
	//UNUSED+// ssize_t e;
	ssize_t v ;
	for (v = 0; v < vertexIntervals.size() ; v++) {
		std::cout << "v: " << v << " ";
		ssize_t vi;
		for (vi = 0; vi < vertexIntervals[v].size(); vi++) {
			std::cout << vertexIntervals[v][vi].value->vertex << " " << 
				vertexIntervals[v][vi].value->index << " " <<
				vertexIntervals[v][vi].value->pos << " " <<
				vertexIntervals[v][vi].value->length << ", ";
			if (vi > 0 and 
					vertexIntervals[v][vi].value->pos ==
					vertexIntervals[v][vi-1].value->pos) {
				std::cout << "error, two equal interval positions " << v << std::endl;
				assert(0);
			}
		}
		std::cout << std::endl;
	}

	ssize_t outEdgeIndex, outEdge;
	ssize_t src, dest;
	ssize_t overlap;

	//UNUSED// ssize_t curProxy = 0;
	for (v = 0; v < vertices.size(); v++) {
		for (outEdgeIndex = vertices[v].FirstOut(); 
				 outEdgeIndex < vertices[v].EndOut();
				 outEdgeIndex = vertices[v].NextOut(outEdgeIndex)) {
			outEdge = vertices[v].out[outEdgeIndex];
			dest    = edges[outEdge].dest;
			src     = v; // keep naming straight.
			overlap = GetOverlap(vertices[src], edges[outEdge], vertices[dest]);
			std::cout << src << " " << outEdge << "(" << edges[outEdge].length 
								<< ") " << dest << " " << overlap << std::endl;
			// the two vertices have overlapping sequence
			if (overlap > 0) {
				ssize_t srcLast, destFirst;
				srcLast = splitPositions[src].size() - 1;
				// Find the split position corresponding 
				// to the overlap of these two vertices
				SplitPosList::iterator srcSplitIt, srcSplitEnd;
				std::cout << "searching for " << vertices[src].vertexSize - overlap - 1 
									<< " in " << std::endl;
				std::copy(splitPositions[src].begin(), splitPositions[src].end(), 
									std::ostream_iterator<ssize_t>(std::cout, ", "));
				std::cout << std::endl;
									
				srcLast = 0;
				for (srcSplitIt = splitPositions[src].begin();
						 srcSplitIt != splitPositions[src].end(); 
						 ++srcSplitIt, ++srcLast) {
					if (*srcSplitIt == vertices[src].vertexSize - overlap - 1)
						break;
				}
				// A split position corresponding to 
				// V - overlap - 1 *must* exist in split positions
				// for src.
				assert(srcSplitIt != splitPositions[src].end());

				SplitPosList::iterator destSplitIt, destSplitEnd;
				
				/* Print some information.*/
				destFirst = 0;
				std::cout << "searching for " << overlap - 1<< " in ";
				for (destSplitIt = splitPositions[dest].begin();
						 destSplitIt != splitPositions[dest].end();
						 ++destSplitIt) {
					std::cout << *destSplitIt << " ";
				}
				std::cout << std::endl;

				
				for (destSplitIt = splitPositions[dest].begin();
						 destSplitIt != splitPositions[dest].end();
						 ++destSplitIt, ++destFirst) {
					if (overlap - 1 == *destSplitIt) 
						break;
				}

				// a split position corresponding to overlap *must*
				// exist in the positions for dest
				assert(destFirst < splitPositions[dest].size());

				// The splits that were overlapping need to be linked.
				// In this case they are represented as disjoint sets,
				// so join the sets with a union.
				ssize_t s, d;
				
				ssize_t numSrcSplit = splitPositions[src].size();
				for (s = srcLast + 1, d = 0; s < numSrcSplit; s++, d++) {
					std::cout << "before: " << dest <<  " " << d << " " 
										<< vertexIntervals[dest][d].parent 
										<< " " << src << " " << s << " " << vertexIntervals[src][s].parent
										<< " ";
 
					DisjointSet::Union(&(vertexIntervals[src][s]), &(vertexIntervals[dest][d]));
					std::cout << "after " << vertexIntervals[dest][d].parent 
										<< " " << vertexIntervals[src][s].parent 
										<< std::endl;
				}
			}
		}
	}
}


DJVertexInterval *GetParent(DJVertexInterval* vi) {
	while (vi != vi->parent) 
		vi = vi->parent;
	return vi;
}

ssize_t GetVertexIntervalIndex(DJVertexInterval &vi) {
	DJVertexInterval *parent;
	parent = GetParent(&vi);
	return parent->value->intvVertex;
}

void SetSplitVertexIndices(ssize_t curNumVertices,
													 DJVertexIntervalsList &vertexIntervals) {
	// Set the index of the vertex intervals according to their 
	// location in the full graph.
	// All overlapping split vertex intervals have been linked together
	// in the DJVertexIntervalsList structure.  Each parent of a joint-set
	// is assigned in index in graph list of vertices, so that the 
	// vertex-interval may be transformed into a vertex in the graph.

	ssize_t v, i;
	DJVertexInterval *parent;
	for (v = 0; v < vertexIntervals.size(); v++ ){ 
		for (i = 0; i < vertexIntervals[v].size(); i++) {
			parent = GetParent(&(vertexIntervals[v][i]));
			if (parent->value->intvVertex == -1) {
				parent->value->intvVertex = curNumVertices;
				curNumVertices++;
			}
			vertexIntervals[v][i].value->intvVertex = parent->value->intvVertex;
			std::cout << "interval " << v << " " << i << " references vertex " 
								<< parent->value->intvVertex << std::endl;
		}
	}
}


void CreateSplitVertexIndexMap(DJVertexIntervalsList &vertexIntervals,
															 std::map<ssize_t, ssize_t>    &splitIndexMap) {
	ssize_t v;
	for (v = 0; v < vertexIntervals.size(); v++ ) {
		// Add this to the map.
		if (vertexIntervals[v].size() > 0) {
			splitIndexMap[v] = vertexIntervals[v][0].value->intvVertex;
		}
	}
}


ssize_t CountUniqueVertexIntervals(DJEdgeCountVertexList &splitVertices) {
	ssize_t v, sv;

	for (v = 0; v < splitVertices.size(); v++) {
		for (sv = 0; sv < splitVertices[v].size(); sv++ ) {
			//parents.insert(splitVertices[v][sv].parent);
			splitVertices[v][sv].parent->value->index = -1;
		}
	}
	ssize_t numParents = 0;
	for (v = 0; v < splitVertices.size(); v++) {
		for (sv = 0; sv < splitVertices[v].size(); sv++ ) {
			//parents.insert(splitVertices[v][sv].parent);
			if (splitVertices[v][sv].parent->value->index == -1) {
				splitVertices[v][sv].parent->value->index = numParents;
				++numParents;
			}
		}
	}
	return numParents;
}

void MergeVertexIntervals(DJEdgeCountVertexList &splitVertices) {
	// This may require lots of resizing, but that's ok.
	ssize_t v, sv;
	ssize_t vertexLength;
	
	if (splitVertices.size() == 0) {
		std::cout << "ERROR, this should be ran on a graph with at least 1 vertx." <<std::endl;
		exit(0);
	}
	// Determine the sizes of all split vertices by merging simple paths.
	//	unsigned char*sequence = new unsigned char[vertices[0].vertexSize];
	for (v = 0; v < splitVertices.size(); v++) {
		vertexLength = 1;

		ssize_t cur = 0;
		if (splitVertices[v].size() == 0) {
			std::cout << "ERROR! the vertex size should be greater than 0!" << std::endl;
			exit(1);
		}
		// Special case here the first vertex is branching.
		if (splitVertices[v][0].value->outDegree > 1) {
			splitVertices[v][0].value->vertexSize = 1;
			cur = 0;
		}
		vertexLength = 1;
		for (sv = 0; sv < splitVertices[v].size() - 1; sv++) {
			// If the current vertex (sv) out degree is 0, no special
			// branches are required other than the link to the next internal
			// vertex, so the next vertex is simply removed.
			// 
			
			if (splitVertices[v][sv].value->outDegree <= 0 and
					splitVertices[v][sv+1].value->inDegree == 0 ) {
				// This vertex and the next should be merged, it is a simple edge.
				//				sequence[vertexLength] = seqPtr[sv+1];

				vertexLength++;
				// mark the next vertex for removal since it will be merged with the previous
				splitVertices[v][sv+1].value->index = -1;
				//				std::cout << "marked " << v << " " << sv + 1 << " for remvoal." << std::endl;
			}
			if (splitVertices[v][sv+1].value->inDegree > 0) {
				// The next vertex should not be merged with the current merged.
				// End the current merged vertex, and start a new one for 
				// the next vertex.
				splitVertices[v][cur].value->vertexSize = vertexLength;

				// To create the super vertex, the in degree is the same as that 
				// of the first vertex, nad the out degree is whatever the 
				// last vertex (at sv) is.
				/*				std::cout << "compare outs: " << v << " " << cur << " " << sv << " "
									<< splitVertices[v][cur].value->outDegree
									<< " " 
									<< splitVertices[v][sv].value->outDegree << std::endl;
				*/
				splitVertices[v][cur].value->outDegree = 
					splitVertices[v][sv].value->outDegree;

				vertexLength = 1;
				//				sequence[0] = seqPtr[sv+1];
				cur = sv + 1;
			}
		}
		if (cur < splitVertices[v].size()) {
			splitVertices[v][cur].value->vertexSize = vertexLength;
			splitVertices[v][cur].value->outDegree = splitVertices[v][splitVertices[v].size()-1].value->outDegree;
		}
	}


	// Fix the connectivity of the parents.

	// Step 1, reset the parents 
	for (v = 0; v < splitVertices.size(); v++ ) {
		DJEdgeCountVertex::iterator it, end;
		for (it = splitVertices[v].begin(), end = splitVertices[v].end();
				 it != end; ++it) {
			(*it).value->traversed = GraphVertex::NotMarked;
		}
	}

	// Remove the vertices marked for removal.
	for (v= 0; v < splitVertices.size(); v++ ) {
		DJEdgeCountVertex::iterator it;
		it = splitVertices[v].begin();
		while (it != splitVertices[v].end()) {
			if ( (*it).value->index == -1)
				splitVertices[v].erase(it);
			else
				++it;
		}
	}
}





ssize_t CountUniqueVertexIntervals(DJVertexIntervalsList &vertexIntervals) {
	// Count the number of disjoint (not glued) vertices
	// by counting the number of unique parents.
	ssize_t v, p;
	DJVertexInterval *parent;
	ssize_t curDisjointVertex = 0;
	std::set<DJVertexInterval*> parents;
	for (v = 0; v < vertexIntervals.size(); v++) {
		for (p = 0; p < vertexIntervals[v].size(); p++) {
			parent = GetParent(&(vertexIntervals[v][p]));
			parents.insert(parent);
		}
		if (vertexIntervals[v].size() > 0)
			++curDisjointVertex;
	}
	return parents.size();
}

ssize_t ConnectVertices(TVertexList &vertices,
										TEdgeList &edges,
										ssize_t src, ssize_t dest, 
										ssize_t &curNumEdges) {
	// Look to see if src is already connected with dest
	ssize_t outEdge, outEdgeIndex;
	ssize_t edgeExists = 0;
	for (outEdgeIndex = vertices[src].FirstOut();
			 outEdgeIndex < vertices[src].EndOut();
			 outEdgeIndex = vertices[src].NextOut(outEdgeIndex)) {
		outEdge = vertices[src].out[outEdgeIndex];
		if (edges[outEdge].dest == dest) {
			edgeExists = 1;
			break;
		}
	}
	if (!edgeExists) {
		edges[curNumEdges].src  = src;
		edges[curNumEdges].dest = dest;
		// find the first available spot to add the cur edge to src
		ssize_t i;
		for (i = 0; i < vertices[src].out.size(); i++ ) {
			if (vertices[src].out[i] == -1) {
				vertices[src].out[i] = curNumEdges;
				break;
			}
		}
		if (i == vertices[src].out.size()) {
			vertices[src].out.push_back(curNumEdges);
		}
		// Find the first available spot to add the cur edge to dest
		for (i = 0; i < vertices[dest].in.size(); i++) {
			if (vertices[dest].in[i] == -1) {
				vertices[dest].in[i] = curNumEdges;
				break;
			}
		}
		if (i == vertices[dest].in.size()) {
			vertices[dest].in.push_back(curNumEdges);
		}
		++curNumEdges;
	}
	return curNumEdges-1;
}

ssize_t VerticesAreConnected(TVertexList &vertices, TEdgeList &edges,
												 ssize_t src, ssize_t dest) {
	
	// Simply return 1 if src is connected to dest by an edge.
	// 0 otherwise.
	ssize_t outEdgeIndex, outEdge;
	for (outEdgeIndex = vertices[src].FirstOut();
			 outEdgeIndex < vertices[src].EndOut();
			 outEdgeIndex = vertices[src].NextOut(outEdgeIndex)) {
		outEdge = vertices[src].out[outEdgeIndex];
		if (edges[outEdge].dest == dest) {
			return 1;
		}
	}
	return 0;
}


void GetEdgeSequence(TVertexList &vertices, TEdgeList &edges,
										 VertexInterval &src, VertexInterval &dest, 
										 unsigned char *&edgeSeq, _SSZT_ &edgeSeqLength) {
	edgeSeqLength = src.length + dest.length;

	// Find an edge that corresponds to src.
	ssize_t outEdgeIndex, outEdge;
	ssize_t inEdgeIndex, inEdge;
	unsigned char *srcSeqPtr, *destSeqPtr;

	// This may only be called when the src vertex has an 
	// out edge.
	ssize_t srcVertex  = src.vertex;
	ssize_t destVertex = dest.vertex;

	assert(vertices[srcVertex].OutDegree() > 0 ||
				 vertices[srcVertex].InDegree() > 0);
	
	if (vertices[srcVertex].OutDegree() > 0) {
		outEdgeIndex = vertices[srcVertex].FirstOut();
		outEdge      = vertices[srcVertex].out[outEdgeIndex];
		srcSeqPtr    = edges[outEdge].seq.seq;
	}
	else {
		inEdgeIndex = vertices[srcVertex].FirstIn();
		inEdge = vertices[srcVertex].in[inEdgeIndex];
		srcSeqPtr = &(edges[inEdge].seq.seq[edges[inEdge].length - vertices[srcVertex].vertexSize]);
	}
									
	// This also may only be called when dest has an in-edge
	assert(vertices[destVertex].InDegree()> 0 || 
				 vertices[destVertex].OutDegree()> 0 );
	if (vertices[destVertex].OutDegree() > 0) {
		outEdgeIndex = vertices[destVertex].FirstOut();
		outEdge      = vertices[destVertex].out[outEdgeIndex];
		destSeqPtr   = edges[outEdge].seq.seq;
	}
	else {
		inEdgeIndex = vertices[destVertex].FirstIn();
		inEdge = vertices[destVertex].in[inEdgeIndex];
		destSeqPtr = &(edges[inEdge].seq.seq[edges[inEdge].length - vertices[destVertex].vertexSize]);
	}
	
	ssize_t s, d;
	assert(edgeSeqLength + 1 > 0);
	edgeSeq = new unsigned char[edgeSeqLength+1];
	for (s = 0; s < src.length; s++) {
		edgeSeq[s] = srcSeqPtr[s + src.pos];
	}

	for (d = 0; d < dest.length; d++) {
		edgeSeq[src.length + d] = destSeqPtr[d + dest.pos];
	}
	edgeSeq[src.length + d] = '\0';
	std::cout << src.length << " " << dest.length << " got sequence for edge: " << edgeSeq << std::endl;
}


void SetEdgeSequence(TVertexList &vertices, TEdgeList &edges, 
										 VertexInterval &src, VertexInterval &dest,
										 TEdge &newEdge) {

	GetEdgeSequence(vertices, edges, src, dest, newEdge.seq.seq, newEdge.seq.length);
	newEdge.length = newEdge.seq.length;
}

ssize_t GetVertexIndex(std::map<ssize_t, ssize_t> &intvIndexMap,
									 ssize_t vertex) {

	// Either a vertex is in the map, so it's split
	// and it should be referenced elsewhere, or 
	// it's not split, and it's left alone.
	if (intvIndexMap.find(vertex) == intvIndexMap.end())
		return vertex;
	else 
		return intvIndexMap[vertex];
}


void DivideReadIntervals(IntervalGraph &graph,
												 TVertexList &vertices, TEdgeList &edges,
												 ssize_t vertexIndex, ssize_t vertexIsSrc, // May be src or dest.
												 ssize_t edgeIndex,   // May or may not pass through vertex.
												 DJVertexIntervals &vertexIntervals) {
	ssize_t intv;
	TEdge &edge = edges[edgeIndex];
	//UNUSED// int vertexSize = vertices[vertexIndex].vertexSize;
	ssize_t lastVI  = vertexIntervals.size()-1;
	ssize_t firstVI = 0;
	ssize_t pathIndex, pathPos;
	ssize_t prevEdge, prevIntv;
	ssize_t nextEdge, nextIntv;

	//UNUSED// ssize_t prevIntvDest;
	ssize_t doDivideReadInterval = 1;
	ssize_t numIntervals = edge.intervals->size();
	for (intv = 0; intv < numIntervals; intv++) {

		pathIndex = (*edge.intervals)[intv].read;
		pathPos   = (*edge.intervals)[intv].pathPos;
		std::cout << "dividing read interval: " << edgeIndex << " " << intv << " " 
							<< pathIndex << " " << pathPos << std::endl;
		if (vertexIsSrc) { 
			// It's possible this path has already been split, check that here.
			if (pathPos > 0) {
				prevEdge = graph.paths[pathIndex][pathPos-1].edge;
				prevIntv = graph.paths[pathIndex][pathPos-1].index;
				if (edges[prevEdge].dest == vertexIntervals[lastVI].value->intvVertex) {
					// The last interval ends in the last interval in vertexIntervals.
					// That means that the path has already been split.
					doDivideReadInterval = 0;

					// Although there is no need to divide the end of this
					// read interval, its length should be updated so that
					// the divided part is cut off.
					ssize_t offset = vertices[edges[edgeIndex].src].vertexSize - 
						vertices[vertexIntervals[lastVI].value->intvVertex].vertexSize;

					(*edges[edgeIndex].intervals)[intv].readPos += offset;
					(*edges[edgeIndex].intervals)[intv].length  -= offset;
					(*edges[edgeIndex].intervals)[intv].edgePos = 0;
				}
			}
		}
		else {
			// The vertex is a dest vertex, check to make sure that 
			// the end of this interval hasn't already been split.
			if (pathPos < graph.pathLengths[pathIndex]-1) {
				// It's possible that this path has already been split
				nextEdge = graph.paths[pathIndex][pathPos+1].edge;
				nextIntv = graph.paths[pathIndex][pathPos+1].index;
				if (edges[nextEdge].src == vertexIntervals[firstVI].value->intvVertex) {
					doDivideReadInterval = 0;

					//  the end of this read interval shouldn't be divided 
					// on the vertexInterval list, but it should be 
					// truncated
					ssize_t truncLength = vertices[edges[edgeIndex].dest].vertexSize - 
						vertices[vertexIntervals[firstVI].value->intvVertex].vertexSize;

					(*edges[nextEdge].intervals)[intv].length -= truncLength;
					// the read pos shouldn't be affected here
				}
			}
		}

		if (doDivideReadInterval) {
			DivideReadInterval(graph, vertices, edges,
												 vertices[vertexIndex], edges[edgeIndex],
												 edgeIndex, intv, vertexIsSrc, vertexIntervals);
		}
	}
}


ssize_t GetStartVertexInterval(ssize_t edgePosInVertex, DJVertexIntervals &vertexIntervals) {
	ssize_t vi;
	for (vi = 0; vi < vertexIntervals.size(); vi++) {
		ssize_t vertexIntvStart, vertexIntvEnd;
		vertexIntvStart = vertexIntervals[vi].value->pos;
		vertexIntvEnd   = vertexIntvStart + vertexIntervals[vi].value->length;
		if (edgePosInVertex >= vertexIntvStart and
				edgePosInVertex < vertexIntvEnd) {
			return vi;
		}
	}
	std::cout << "ERROR! should have found part of the vertex that overlaps an edge begin" << std::endl;
	assert(0);
	return -1; // Never reach here.  Quiet the compiler warnings.
}


ssize_t GetEndVertexInterval(ssize_t edgePosInVertex, DJVertexIntervals &vertexIntervals) {
	ssize_t vi;
	for (vi = 0; vi < vertexIntervals.size(); vi++) {
		ssize_t vertexIntvStart, vertexIntvEnd;
		vertexIntvStart = vertexIntervals[vi].value->pos;
		vertexIntvEnd   = vertexIntvStart + vertexIntervals[vi].value->length;
		if (edgePosInVertex > vertexIntvStart and
				edgePosInVertex <= vertexIntvEnd) {
			return vi + 1;
		}
	}
	std::cout << "ERROR! should have found part of the vertex that overlaps an edge end" << std::endl;
	assert(0);
	return -1; // Never reach here.  Quiet the compiler warnings.
}

ssize_t LookupEdgeToDest(TVertex &vertex,
										 TEdgeList &edges,
										 ssize_t dest) {
	ssize_t outEdge, outEdgeIndex;
	for (outEdgeIndex = vertex.FirstOut();
			 outEdgeIndex < vertex.EndOut();
			 outEdgeIndex = vertex.NextOut(outEdgeIndex)) {
		outEdge = vertex.out[outEdgeIndex];
		if (edges[outEdge].dest == dest) 
			return outEdge;
	}
	return -1;
}

void StoreVertexIntervalEdges(TVertexList &vertices, 
															TEdgeList &edges,
															DJVertexIntervals &vertexIntervals,
															ssize_t startIntervalIndex,
															ssize_t endIntervalIndex,
															std::vector<ssize_t> &edgeList) {

	ssize_t vi;
	ssize_t src, dest, edge;
	for (vi = startIntervalIndex; vi < endIntervalIndex-1; vi++) {
		// The vertices corresponding to the vertex intervals 
		src  = vertexIntervals[vi].value->intvVertex;
		dest = vertexIntervals[vi+1].value->intvVertex;
		edge = LookupEdgeToDest(vertices[src], edges, dest);
		// These vertices *must* be connected.
		assert(edge != -1);
		edgeList.push_back(edge);
	}
}


void DivideReadInterval(IntervalGraph &graph,
												TVertexList &vertices,
												TEdgeList   &edges,
												TVertex &vertex, TEdge &edge,
												ssize_t intervalEdge, ssize_t intervalIndex,
												ssize_t vertexIsSrc, 
												DJVertexIntervals &vertexIntervals) {
	// Divide a read interval.
	// We have an edge:
	// ---------------------------------------------->
	// (********)                           (********)	
	// =================  (case 1, RI at beg of edge).
	//   ================== (case 2 RI overlaps src)
	// (case 3 RI overlaps dest)   ===================
	//

	ssize_t vertexInterval, lastVertexInterval;
	ssize_t vertexBeginInEdge, vertexEndInEdge;

	// Configure the boundaries of the read interval
	// overlap with this vertex

	// Use 0-based half-open boundaries: (0..10], for example
	// for an interval of length 10.
	ssize_t readIntvPosInVertex, readIntvEndInVertex;
		
	// Congfigure the boundaries of where the vertex overlaps
	// the edge.  This simply depends if the vertex is the 
	// source or dest of the edge.
	if (vertexIsSrc == 1) {
		// Vertex starts at the beginning of the edge.  Configure the 
		// length of the read interval split accordingly.
		vertexBeginInEdge = 0; 
		vertexEndInEdge   = vertex.vertexSize;
	}
	else {
		vertexBeginInEdge = edge.length - vertex.vertexSize;
		vertexEndInEdge   = edge.length;
	}

	// Get quick access to the interval that will be split.
	ReadInterval &intv  = (*edge.intervals)[intervalIndex];

	// Do some checking.  Make sure that part of the read interval 
	// overlaps with the vertex.
	if (intv.edgePos >= vertexEndInEdge or
			intv.edgePos + intv.length - 1 < vertexBeginInEdge) {
		// There is no overlap, just finish
		return;
	}

	std::cout << vertexBeginInEdge<< " " << vertexEndInEdge 
						<< " splitting path " << intervalEdge << std::endl;
	graph.PrintPath((*edge.intervals)[intervalIndex].read, std::cout);
	std::cout << "is src: " << vertexIsSrc << std::endl;
	std::cout << "along positions: ";
	ssize_t vp;
	for (vp = 0; vp < vertexIntervals.size(); vp++) {
		std::cout << vertexIntervals[vp].value->pos << " " << 
			vertexIntervals[vp].value->length << ", ";
	}
	std::cout << std::endl;


	// The read interval overlaps with the vertex, find out where in
	// the vertex the overlap begins.
	if (intv.edgePos <= vertexBeginInEdge and 
			intv.edgePos + intv.length > vertexBeginInEdge) {
		// The read interval is before or at the beginning
		// of the vertex.  Start the interval at the beginning
		// of the vertex.
		readIntvPosInVertex = 0;
	}
	else {
		// The read interval starts inside the vertex.
		readIntvPosInVertex = intv.edgePos - vertexBeginInEdge;
	}

	// Find where this interval ends in the vertex
	if (intv.edgePos + intv.length >= vertexEndInEdge)
		readIntvEndInVertex = vertex.vertexSize;
	else
		readIntvEndInVertex = intv.edgePos + intv.length - vertexBeginInEdge;

	// Insert the new edges into the path.
	
	// Find the indices of the vertex intervals.
	vertexInterval     = GetStartVertexInterval(readIntvPosInVertex, vertexIntervals);
	lastVertexInterval = GetEndVertexInterval(readIntvEndInVertex, vertexIntervals);
	std::cout << "first interval: " << vertexInterval << " last " << lastVertexInterval << std::endl;
	std::vector<ssize_t> vertexIntvEdges;

	// Each pair of adjacent vertex intervals is connected by 
	// an edge.  Store those edges.
	StoreVertexIntervalEdges(vertices, edges, vertexIntervals,
													 vertexInterval, lastVertexInterval,
													 vertexIntvEdges);

	// Now replace the path with the split intervals.
	ssize_t pathStartPos;
	ssize_t i;
	ssize_t numIntervals = lastVertexInterval - vertexInterval;
	ssize_t pathIndex, pathPos;

	pathIndex = intv.read;
	pathPos   = intv.pathPos;

	ssize_t intvIndex, intvEdge;
	ssize_t readPos;
	//UNUSED// ssize_t edgePos;
	//UNUSED+// ssize_t  length;
	ssize_t remIntvLength=0;
	if (vertexIsSrc) {
		// Fix the interval after the source vertex

		// There is an interval after the interval that was split.
		// 
		// Since there is more than one interval in the path,
		// vertex must be an internal vertex on the path.

		ssize_t sucPathEdge = graph.paths[intv.read][intv.pathPos].edge;
		ssize_t sucPathIntv = graph.paths[intv.read][intv.pathPos].index;

		ssize_t numVertexIntervals = vertexIntervals.size();
		ssize_t lastIntervalVertex = vertexIntervals[numVertexIntervals-1].value->intvVertex;


		
		ssize_t intvOvpLength;
		// The amount that the interval overlaps the source
		// vertex.
		intvOvpLength =  vertex.vertexSize - intv.edgePos;
		
		// Truncate the portion of the read interval that
		// overlaps

		// Store the original read pos so that we can make sure
		// that once the path has been split, the read pos
		// corresponding to the first interval is the same.

		ssize_t origReadPos = (*edges[sucPathEdge].intervals)[sucPathIntv].readPos;

		
		(*edges[sucPathEdge].intervals)[sucPathIntv].readPos += intvOvpLength;
		(*edges[sucPathEdge].intervals)[sucPathIntv].length  -= intvOvpLength;

		// Start the read pos at the end of the vertex interval list.
		readPos = (*edges[sucPathEdge].intervals)[sucPathIntv].readPos;
		graph.ReplacePathRange(intv.read, intervalEdge, intervalIndex,
													 intv.pathPos, intv.pathPos, vertexIntvEdges);

		//		std::cout << "after splicing path. " << vertexIntvEdges.size() << std::endl;
		//		graph.PrintPath(intv.read, std::cout);
		ssize_t srcVertexLength;
		ssize_t curVertexLength;
		curVertexLength = vertices[lastIntervalVertex].vertexSize;

		// Compute the overlap length of the original read interval
		// with the 1 .. n-1'th read intervals.
		if (intvOvpLength >= curVertexLength) {
			readPos -= curVertexLength;
			intvOvpLength -= curVertexLength;
		}
		else {
			readPos -= intvOvpLength;
			intvOvpLength = 0;
		}

		// Fix the read positions of the new intervals
		ssize_t prevVertexOverlapLength;
		for (i = numIntervals - 2; i >= 0 ; i--) {

			intvEdge = graph.paths[pathIndex][pathPos + i].edge;
			intvIndex = graph.paths[pathIndex][pathPos + i].index;
			srcVertexLength = vertices[edges[intvEdge].src].vertexSize;
			// Some sanity checks.
			assert(remIntvLength >= 0);
			assert(edges[intvEdge].length > 0);

			if (intvOvpLength < srcVertexLength) {
				readPos -= intvOvpLength;
				prevVertexOverlapLength = intvOvpLength;
			}
			else {
				readPos -= srcVertexLength;
				prevVertexOverlapLength = srcVertexLength;
			}

			
			(*edges[intvEdge].intervals)[intvIndex].readPos = readPos;

			(*edges[intvEdge].intervals)[intvIndex].length  = prevVertexOverlapLength + curVertexLength;

			curVertexLength = srcVertexLength;

			// account for the remainder of the 
			intvOvpLength -= prevVertexOverlapLength;
			(*edges[intvEdge].intervals)[intvIndex].edgePos = edges[intvEdge].length - 
				(*edges[intvEdge].intervals)[intvIndex].length;
		}
		assert(readPos == origReadPos);
		std::cout << "after fixing (src) boundaries " << std::endl;
		graph.PrintPath(intv.read, std::cout);

		
	}
	else {
		// It is necessary to truncate the interval before the 
		// replaced subpath.

		ssize_t precPathEdge = graph.paths[intv.read][intv.pathPos].edge;
		ssize_t precPathIntv = graph.paths[intv.read][intv.pathPos].index;
		ssize_t precDestVertex   = edges[precPathEdge].dest;
		ssize_t precOvpLen   = (*edges[precPathEdge].intervals)[precPathIntv].edgePos + 
			(*edges[precPathEdge].intervals)[precPathIntv].length - 
			(edges[precPathEdge].length - vertices[precDestVertex].vertexSize);


		assert(precOvpLen > 0);

		// The overlap length must be at least the length of the first interval
		// since the overlaps are used to define the vertex intervals.
		
		ssize_t firstVertex = vertexIntervals[0].value->intvVertex;
		
		// Now trim the length of the preceeding edge interval to end after the
		// first vertex interval.
		
		(*edges[precPathEdge].intervals)[precPathIntv].length -=
			precOvpLen + vertices[firstVertex].vertexSize;
		
		// the read pos starts after the end of this interval.
		readPos = (*edges[precPathEdge].intervals)[precPathIntv].readPos + 
			(*edges[precPathEdge].intervals)[precPathIntv].length;
		
		pathStartPos = intv.pathPos + 1;
		std::cout << "inserting " << vertexIntvEdges.size() << " edges."  << std::endl;
		graph.ReplacePathRange(intv.read, intervalEdge, intervalIndex,
													 pathStartPos, pathStartPos, vertexIntvEdges);

		std::cout << "after splicing path. " << std::endl;
		graph.PrintPath(intv.read, std::cout);
		remIntvLength = precOvpLen;
		pathPos = pathStartPos;
		pathIndex = intv.read;
		//UNUSED+// ssize_t  srcVertexIndex;
		ssize_t destVertexIndex;
		destVertexIndex = edges[graph.paths[pathIndex][pathPos].edge].dest;
		for (i = 0; i < numIntervals-1; i++) {
			
			intvEdge = graph.paths[pathIndex][pathPos+ i].edge;
			intvIndex = graph.paths[pathIndex][pathPos+ i].index;
			
			// Some sanity checks.
			assert(remIntvLength >= 0);
			assert(edges[intvEdge].length > 0);

			(*edges[intvEdge].intervals)[intvIndex].readPos = readPos;
			if (remIntvLength >= edges[intvEdge].length) 
				(*edges[intvEdge].intervals)[intvIndex].length  = edges[intvEdge].length;
			else
				(*edges[intvEdge].intervals)[intvIndex].length = remIntvLength;
			
			// account for the remainder of the 
			remIntvLength -= vertices[edges[intvEdge].src].vertexSize;
			(*edges[intvEdge].intervals)[intvIndex].edgePos = 0;
			readPos += vertices[vertexIntervals[vertexInterval + i].value->intvVertex].vertexSize;
		}

		std::cout << "after fixing (dest) boundaries " << std::endl;
		graph.PrintPath(intv.read, std::cout);
	}

	std::cout << "adding " << vertexInterval <<"  " << lastVertexInterval
						<<"  " << vertexIntvEdges.size() << " edges." << std::endl;

	// Fix the read positions of the replaced path.
	// Remember what you thought about on the run...
}
												 
void TrimEdgeSequence(TEdge &edge,
											ssize_t start, ssize_t end) {

	// Fix the edge sequence.
	assert(end - start > 0);
	unsigned char *newSeq = new unsigned char[end - start];
	assert(end-  start <= edge.length);
	memcpy(newSeq, &(edge.seq.seq[start]), end - start);
	
	// Update the edge length and sequence.
	edge.length = end - start;
	edge.seq.length = edge.length;
	delete[] edge.seq.seq;
	edge.seq.seq = newSeq;
}


void GlueOverlappingVertices(IntervalGraph         &graph,
														 TVertexList           &vertices,
														 TEdgeList             &edges,
														 DJVertexIntervalsList &vertexIntervals) {
	
	// The division of each vertex into vertex intervals is stored
	// in vertexIntervals.  This replaces a set of vertex intervals
	// that all have a common parent with one vertex (gluing them
	// together).  

	// 

	// First count how many components there are of 'glued' vertices.
	ssize_t numUniqueIntervals;
	numUniqueIntervals = CountUniqueVertexIntervals(vertexIntervals);
	
	// Set the split intervals to be appended to the end of vertices.
	SetSplitVertexIndices(vertices.size(), vertexIntervals);

	// Create a map that given a vertex that has been split,
	// gives the index of that vertex in the appended list.
	// Vertices that are unchanged should not be found in the map.

	std::map<ssize_t, ssize_t> intvIndexMap;
	CreateSplitVertexIndexMap(vertexIntervals, intvIndexMap);

	ssize_t origNumVertices = vertices.size();
	ssize_t origNumEdges    = edges.size();
	ssize_t v, e;
	// Make room for the new vertices.
	vertices.resize(vertices.size() + numUniqueIntervals);
	edges.resize(edges.size()       + numUniqueIntervals);


	// Initialize interval storage for the new edges.
	for(e = origNumEdges; e < edges.size(); e++ )
		edges[e].Init();
	
	// Now, connect every pair of overlapping vertex intervals, unless
	// they have already been connected.
	
	ssize_t srcIndex, destIndex;
	ssize_t intv;
	ssize_t edgeIndex = origNumEdges;
	//UNUSED// ssize_t curEdge;
	ssize_t origVertex;
	// Use each original split vertex to link glued vertices
	// 
	//UNUSED// ssize_t isSrc  = 1;
	//UNUSED// ssize_t isDest = 0;
	for (v = 0; v < origNumVertices; v++ ){
		origVertex = v;
		if (vertexIntervals[v].size() == 0)
			continue;
		
		for (intv = 0; intv < vertexIntervals[v].size() - 1; intv++) {
			
			// Find the indices of the glued vertices
			// in the vertex list.
			srcIndex   = vertexIntervals[v][intv].value->intvVertex;
			destIndex  = vertexIntervals[v][intv+1].value->intvVertex;

			// Configure the vertex lengths
			vertices[srcIndex].vertexSize = 
				vertexIntervals[v][intv].value->length;

			vertices[destIndex].vertexSize = 
				vertexIntervals[v][intv+1].value->length;

			// All vertices generated from vertex intervals should
			// be joined by edges.  Create the edge if it doesn't
			// already exist.
			if (!VerticesAreConnected(vertices, edges, srcIndex, destIndex)) {
				ssize_t newEdge;
				newEdge = ConnectVertices(vertices, edges,
																	srcIndex, destIndex, edgeIndex);
				
				std::cout << "connecting vertices " << srcIndex << " " << destIndex 
									<< " with edge: " << edgeIndex <<" "  << vertices[srcIndex].vertexSize + 
					vertices[destIndex].vertexSize << std::endl;

				// Once the connectivity of the edge is done, 
				// create the sequence.
				SetEdgeSequence(vertices, edges, 
												*vertexIntervals[v][intv].value,
												*vertexIntervals[v][intv + 1].value,
												edges[newEdge]);
			}
		} // Done joining intervals from this vertex

		ssize_t outEdgeIndex, outEdge, inEdgeIndex, inEdge;
		ssize_t overlap;
		// Get the index of the last interval for the split vertex.
		// This interval will connect to vertices as the source.
		ssize_t lastIntvVertexIndex =
			vertexIntervals[v][vertexIntervals[v].size()-1].value->intvVertex;
		ssize_t firstIntvVertexIndex =
			vertexIntervals[v][0].value->intvVertex;

		for (inEdgeIndex = vertices[origVertex].FirstIn();
				 inEdgeIndex < vertices[origVertex].EndIn();
				 inEdgeIndex = vertices[origVertex].NextIn(inEdgeIndex)) {
			inEdge  = vertices[origVertex].in[inEdgeIndex];
			// get rid of the previous in index
			vertices[origVertex].in[inEdgeIndex] = -1;
			overlap = GetOverlap(vertices[edges[inEdge].src], edges[inEdge], vertices[origVertex]);
			if (overlap <= 0) {
				// The vertices V, origVertex, connected by inEdge =(V,origVertex)
				// do not overlap.
				
				// Look to see if the source vertex is one that 
				/*				std::cout << "dividing read intervals for edge: " << inEdge 
									<< " into " << origVertex << std::endl;
										DivideReadIntervals(graph, vertices, edges,
														origVertex, isDest, 
														inEdge, vertexIntervals[v]);
				*/
				// Route the in edge to the first split vertex 
				edges[inEdge].dest = firstIntvVertexIndex;
				ssize_t edgeEndPos = edges[inEdge].length - 
					vertices[origVertex].vertexSize + 
					vertices[firstIntvVertexIndex].vertexSize;
				
				TrimEdgeSequence(edges[inEdge], 0, edgeEndPos);
			}
		}

		// Reconnect the edges (V, src) to (V, srcIndex)


		// Add out-edges from the src that are not
		// overlapping

		for (outEdgeIndex = vertices[origVertex].FirstOut();
				 outEdgeIndex < vertices[origVertex].EndOut();
				 outEdgeIndex = vertices[origVertex].NextOut(outEdgeIndex)) {
			outEdge = vertices[origVertex].out[outEdgeIndex];
			vertices[origVertex].out[outEdgeIndex] = -1;
			overlap = GetOverlap(vertices[origVertex], edges[outEdge], vertices[edges[outEdge].dest]);
			
			if (overlap <= 0) {
 				// The two vertices do not overlap, so the edge will
				// not be created when connecting vertex intervals.  Instead,
				// the edge should be connected to the new split vertex

				// Re-link the out edge
				std::cout << "Dividing read intervals for out: " << outEdge
									<< " from " << origVertex << std::endl;
				/*				DivideReadIntervals(graph, vertices, edges, 
														origVertex, isSrc,
														outEdge, vertexIntervals[v]);
				*/			
				edges[outEdge].src = lastIntvVertexIndex;

				ssize_t altDestIndex;
				altDestIndex = GetVertexIndex(intvIndexMap, edges[outEdge].dest);
				/*
				if (altDestIndex == edges[outEdge].dest) {

					// Unlink this from the original vertex (although it will 
					// be deleted anyway.
					
				}
				else {
					// The dest vertex is split and has not been updated.  The edge should be re-routed to 
					// an alternative vertex.
					ssize_t destFirstIntvIndex;
					assert(vertexIntervals[altDestIndex].size() > 0);
					destFirstIntvIndex = vertexIntervals[altDestIndex][0].value->intvVertex;

					ssize_t oldEdgeIndex = vertices[edges[outEdge].dest].LookupInIndex(outEdge);

					assert(oldEdgeIndex >= 0);
					assert(edges[outEdge].dest >= 0);
					vertices[edges[outEdge].dest].in[oldEdgeIndex] = -1;
					
					// Link in the new edge.
					edges[outEdge].dest = destFirstIntvIndex;
					std::cout << "moved src of edge: " << outEdge 
										<< " to " << srcIndex << std::endl;
				}
				*/
			}
		}		
	}
}

ssize_t FindEdgeIndex(TVertex &src, TEdgeList &edges, ssize_t dest) {
	ssize_t e;
	for (e = 0; e < src.out.size(); e++ ) {
		if (src.out[e] >= 0 and
				edges[src.out[e]].dest == dest) {
			return src.out[e];
		}
	}
	return -1;
}


int main(int argc, char* argv[]) {
	if (argc < 3) {
		PrintUsage();
		exit(1);
	}
	std::string graphInName, graphOutName;
	int vertexSize;
	
	graphInName = argv[1];
	vertexSize  = atoi(argv[2]);
	graphOutName = argv[3];

	IntervalGraph graph;
	std::cout << "reading graph " << std::endl;
	ReadIntervalGraph(graphInName, graph, vertexSize);
	std::cout << "done. " << std::endl;
	// Get a quick reference to the vertices and edges
	TVertexList &vertices = graph.vertices;
	TEdgeList   &edges    = graph.edges;

	std::vector<ssize_t> vertexBalance;
	std::vector<ssize_t> originalVertices;
	std::vector<unsigned char*> vertexLabels;
	StoreVertexBalance(vertices, edges, vertexBalance);

	// Count the number of internal vertices;
	ssize_t v, e;
	ssize_t numInternal = 0;
	for (v = 0; v < vertices.size(); v++) { 
		if (vertices[v].InDegree() != 0 and
				vertices[v].OutDegree() != 0) {
			numInternal++;
		}
	}
	
	// Store the vertex labels for constructing the
	// edge sequences later on.
	vertexLabels.resize(vertices.size());
	for (v = 0; v < vertices.size(); v++ ) {
		vertexLabels[v] = (unsigned char*) new unsigned char[vertexSize];
		if (vertices[v].OutDegree() > 0) {
			memcpy((char*)vertexLabels[v],
						 edges[vertices[v].out[vertices[v].FirstOut()]].seq.seq,
						 vertexSize);
		}
		else {
			ssize_t edgeIndex;
			edgeIndex = vertices[v].in[vertices[v].FirstIn()];
			memcpy((char*)vertexLabels[v],
						 &(edges[edgeIndex].seq.seq[edges[edgeIndex].seq.length -
																				vertexSize]),
						 vertexSize);
		}
	}
		
	//	SplitDJVertexList splitVertices;
	DJEdgeCountVertexList splitVertices;
	splitVertices.resize(vertices.size());
	for (v = 0; v < splitVertices.size(); v++) {
		splitVertices[v].resize(vertexSize);
	}

	
	// copy over the graph structure
	ssize_t first = 0;
	ssize_t last  = vertexSize - 1;
	ssize_t sv;
	ssize_t edgeIndex = edges.size();


	
	// Make room for the new edge list.


	// Copy the previous connectivity.
	edgeIndex = edges.size();

	std::vector<ssize_t> origSrcList, origDestList, origEdgeLengthList;
	origSrcList.resize(edges.size());
	origDestList.resize(edges.size());
	origEdgeLengthList.resize(edges.size());
	for (e = 0; e < edges.size(); e++ ) {
		origSrcList[e] = edges[e].src;
		origDestList[e] = edges[e].dest;
		origEdgeLengthList[e] = edges[e].length;
	}
	
	// Initialize the split vertices.
	for (v = 0; v < vertices.size(); v++) {
		for (sv = 0; sv < vertexSize ; sv++ ) {
			splitVertices[v][sv].value = new EdgeCount;
			splitVertices[v][sv].parent = &(splitVertices[v][sv]);
		}
	}
	
	// Initialize connectivity.
	for (v = 0; v < vertices.size(); v++) {
		(*splitVertices[v][first].value).inDegree = vertices[v].InDegree();
		(*splitVertices[v][last].value).outDegree = vertices[v].OutDegree();
	}
	
	// Initialize seqeunce pointers.
	for (v = 0; v < vertices.size(); v++) {
		unsigned char *vertexPtr;
		if (vertices[v].OutDegree() != 0) {
			vertexPtr = edges[vertices[v].out[vertices[v].FirstOut()]].seq.seq;
		}
		else {
			assert(vertices[v].InDegree() > 0);
			ssize_t edgeIndex = vertices[v].in[vertices[v].FirstIn()];
			
			vertexPtr = &(edges[edgeIndex].seq.seq[edges[edgeIndex].seq.length - 
																						 vertices[v].vertexSize]);
		}
		for (sv = 0; sv < vertexSize; sv++ ){
			(*splitVertices[v][sv].value).seqPtr = vertexLabels[v] + sv;
		}
	}

	ssize_t vertexIndex = 0;
	for (v = 0; v < vertices.size(); v++) {
		for (sv = 0; sv < vertexSize; sv++ ) {
			splitVertices[v][sv].value->index = vertexIndex++;
		}
	}

	// Glue overlapping vertices here.
	GlueOverlappingVertices(vertices, edges, vertexSize, splitVertices);
	SumVertexDegrees(splitVertices);
	// Set the parent pointers for all disjoint vertices so that there is just 
	// 1 link to the parent.

	for (v = 0; v < vertices.size(); v++) {
		std::cout << v << " ";
		for (sv = 0; sv < splitVertices[v].size(); sv++ ){
			std::cout << splitVertices[v][sv].value->inDegree << " " << splitVertices[v][sv].value->outDegree << ",";
		}
		std::cout << std::endl;
	}
	MergeVertexIntervals(splitVertices);
	std::cout << " after merging: " << std::endl;
	for (v = 0; v < vertices.size(); v++) {
		std::cout << v << " ";
		for (sv = 0; sv < splitVertices[v].size(); sv++ ){
			std::cout << splitVertices[v][sv].value->inDegree << " " << splitVertices[v][sv].value->outDegree << ",";
		}
		std::cout << std::endl;
	}

	ssize_t numUniqueVertices;
	numUniqueVertices = CountUniqueVertexIntervals(splitVertices);

	
	std::cout << "there are: " << numUniqueVertices << " unique vertices." << std::endl;
	// make space that actually holds connectivity.
	TVertexList tSplitVertices;
	tSplitVertices.resize(numUniqueVertices);

	// Restore the edge connectivity of the original edges.
	//UNUSED// ssize_t nOrigEdges = edges.size();
	
	for (e = 0; e< edges.size(); e++) {
		ssize_t origSrc, origDest;
		origSrc  = edges[e].src;
		origDest = edges[e].dest;
		
		ssize_t newSrc, newDest;
		ssize_t overlap;
		overlap = GetOverlap(vertices[origSrc], edges[e], vertices[origDest]);
		ssize_t srcEnd = splitVertices[origSrc].size() - 1;
		newSrc = splitVertices[origSrc][srcEnd].parent->value->index;
		ssize_t srcLen, destLen;
		//UNUSED// ssize_t s;

		if (overlap <= 0) {
			// Simply connect the end of src to the beginning of dest.
			newDest = splitVertices[origDest][0].parent->value->index;
			srcLen  = splitVertices[origSrc][srcEnd].parent->value->vertexSize;
			destLen = splitVertices[origDest][0].parent->value->vertexSize;
			ssize_t offset = vertexSize - srcLen;
			ssize_t newLength = edges[e].length - (vertexSize - srcLen) - (vertexSize - destLen);
			// shift the sequence back 
			memmove(edges[e].seq.seq,
							(char*) &(edges[e].seq.seq[offset]),
							newLength);
			edges[e].length = newLength;
			edges[e].seq.length = newLength;
		}
		else {
			// The two eges overlap, so it's not so easy as linking the
			// end of the src to the beginning of the dest.

			// Find where the new src liks in.
			// This is just a sanity check for now.
			ssize_t destPos   = 0;
			ssize_t destStart = overlap;
			for (sv = 0; sv < splitVertices[origDest].size() and
						 destPos != destStart; sv++ ) {
				destPos += splitVertices[origDest][sv].parent->value->vertexSize;
			}
			// Make sure the source was found.
			assert(destPos < vertices[origDest].vertexSize);
			assert(sv < splitVertices[origDest].size());
			newDest = splitVertices[origDest][sv].parent->value->index;
			
			memcpy(edges[e].seq.seq,
						 splitVertices[origSrc][srcEnd].parent->value->seqPtr,
						 splitVertices[origSrc][srcEnd].parent->value->vertexSize);

			memcpy(edges[e].seq.seq + splitVertices[origSrc][srcEnd].parent->value->vertexSize,
						 splitVertices[origDest][0].parent->value->seqPtr,
						 splitVertices[origDest][0].parent->value->vertexSize);

			edges[e].length = splitVertices[origSrc][srcEnd].parent->value->vertexSize + 
				splitVertices[origDest][0].parent->value->vertexSize;
			
			edges[e].seq.length = edges[e].length;
		}

		edges[e].src  = newSrc;
		edges[e].dest = newDest;
		tSplitVertices[newSrc].AddOutEdge(e);
		tSplitVertices[newDest].AddInEdge(e);
	}

	// Configure the sizes of the split vertices.
	for (v = 0; v < splitVertices.size(); v++ )  {
		for (sv = 0; sv < splitVertices[v].size(); sv++ ) {
			tSplitVertices[splitVertices[v][sv].parent->value->index].vertexSize =
				splitVertices[v][sv].parent->value->vertexSize;
		}
	}

	//UNUSED// ssize_t nNewEdges = 0;
	// Count the number of split vertices that have not been linked by an edge.
	ssize_t outEdge, outEdgeIndex;
	std::vector<TEdge> newEdges;
	for (v = 0; v < splitVertices.size(); v++) {
		for (sv = 0; sv < splitVertices[v].size() - 1; sv++) {
			ssize_t src, dest;
			src = splitVertices[v][sv].parent->value->index;
			dest = splitVertices[v][sv+1].parent->value->index;
			ssize_t edgeFound = 0;
			for (outEdgeIndex = tSplitVertices[src].FirstOut();
					 outEdgeIndex != tSplitVertices[src].EndOut();
					 outEdgeIndex = tSplitVertices[src].NextOut(outEdgeIndex)) { 
				outEdge = tSplitVertices[src].out[outEdgeIndex];
				if ((outEdge < edges.size() and edges[outEdge].dest == dest) or
						(newEdges[outEdge - edges.size()].dest == dest)) {
					edgeFound = 1;
					break;
				}
			}
			if (edgeFound == 0) {
				// This pair of vertices requires a new edge.
				ssize_t lastEdge = newEdges.size();
				newEdges.push_back(TEdge());
				newEdges[lastEdge].Init();

				newEdges[lastEdge].src = src;
				newEdges[lastEdge].dest = dest;
				std::cout << "new edge " << src << " " << dest << " " << v << " " << sv << std::endl;

				tSplitVertices[src].AddOutEdge(lastEdge + edges.size());
				tSplitVertices[dest].AddInEdge(lastEdge + edges.size());

				newEdges[lastEdge].seq.length = 
					newEdges[lastEdge].length = tSplitVertices[src].vertexSize + 
					tSplitVertices[dest].vertexSize;

				newEdges[lastEdge].seq.seq = new unsigned char[newEdges[lastEdge].length];

				// The edge sequence is simply the combination of the two 
				// vertex sequences.
				memcpy((char*) newEdges[lastEdge].seq.seq,
							 splitVertices[v][sv].parent->value->seqPtr,
							 splitVertices[v][sv].parent->value->vertexSize);

				memcpy((char*) (newEdges[lastEdge].seq.seq +
												splitVertices[v][sv].parent->value->vertexSize),
							 splitVertices[v][sv+1].parent->value->seqPtr,
							 splitVertices[v][sv+1].parent->value->vertexSize);
			}
		}
	}
					
	std::cout << "created: " << newEdges.size() << " new edges. " << std::endl;
	edges.insert(edges.end(),newEdges.begin(), newEdges.end());
	
	// Assign the index of each split vertex to that of the parent for quick reference.


	// Determine the balance of hte new edges.

	ssize_t nAssignedBalance = 0;
 
	for (v = 0; v < vertices.size(); v++ ) {
		ssize_t bv = vertexBalance[v];
		assert(splitVertices[v].size() == splitVertices[bv].size());
		ssize_t svl = splitVertices[v].size();
#if _DEBUG_
		for (sv = 0; sv < svl; sv++ ) {
			assert(splitVertices[v][sv].value->inDegree ==
						 splitVertices[bv][svl - sv - 1].value->outDegree);
			assert(splitVertices[v][sv].value->outDegree ==
						 splitVertices[bv][svl - sv - 1].value->inDegree);
		}
#endif
		for (sv = 0; sv < svl - 1; sv++ ){
			// Find the edge connecting this vertex and the next
			ssize_t src, dest;
			ssize_t bSrc, bDest;
			src = splitVertices[v][sv].parent->value->index;
			dest = splitVertices[v][sv+1].parent->value->index;
			
			bSrc = splitVertices[bv][svl - sv - 2].parent->value->index;
			bDest = splitVertices[bv][svl - sv - 1].parent->value->index;
			
			//UNUSED// ssize_t edgeIndex, bEdgeIndex;
			ssize_t edge, bEdge;
			edge  = -1;
			bEdge = -1;
			edge  = FindEdgeIndex(tSplitVertices[src], edges, dest);
			bEdge = FindEdgeIndex(tSplitVertices[bSrc], edges, bDest);
			
			// We must have found both edges.
			assert(edge >= 0);
			assert(bEdge >= 0);

			edges[edge].balancedEdge = bEdge;
			edges[bEdge].balancedEdge = edge;
			nAssignedBalance++;
		}
	}
	std::cout << "assigned balance of: " << nAssignedBalance << std::endl;
	
	std::cout << "The new graph is:" << std::endl;
	WriteVertexList(std::cout, tSplitVertices);
	WriteEdgeList(std::cout, edges);

	ssize_t p, pi;
	PathLengthList &pathLengths = graph.pathLengths;
	PathIntervalList &paths = graph.paths;
	

	for (p = 0; p < paths.size(); p++ ) {

		// Make sure there is a path to process
		if (pathLengths[p] <= 0) 
			continue;

		// Indices for referencing the path edge, and path edge index.
		ssize_t pe, pei;
		ssize_t sv;
		pe  = paths[p][0].edge;
		pei = paths[p][0].index;

		// The position on the read for the split path.
		ssize_t splitReadPos = (*edges[pe].intervals)[pei].readPos;
		ssize_t splitEdgePos = (*edges[pe].intervals)[pei].edgePos;
		ssize_t splitPathPos = 0;
		ssize_t src, dest;
		ssize_t origSrc, origDest;
		ReadInterval intv;
		
		src = edges[pe].src;
		dest = edges[pe].dest;
		origSrc = origSrcList[pe];
		origDest = origDestList[pe];

		if (splitEdgePos < vertexSize) {
			// The path begins in the vertex sequence.  It may
			// be necessary to split it.
			ssize_t edgeSrcOffset;
			edgeSrcOffset = vertexSize - tSplitVertices[src].vertexSize;

			// Determine the connectivity of the edge to the split vertices
			src  = edges[pe].src;
			dest = edges[pe].dest;
			// Find the original vertices so that we know how much of it is split.
			origSrc  = origSrcList[pe];
			origDest = origDestList[pe];
			
			// Find the splitvertex containing this one.
			ssize_t pos = 0;
			ssize_t pathStart = -1;
			for (sv = 0; sv < splitVertices[origSrc].size(); sv++) {
				if (pos + splitVertices[origSrc][sv].parent->value->vertexSize > splitEdgePos) {
					pathStart = sv;
					break;
				}
				pos += splitVertices[origSrc][sv].parent->value->vertexSize;
			}
			assert(pathStart >= 0);
			if (pathStart < splitVertices[origSrc].size() - 1) {
				// This path starts in a vertex that has multiple splits.
				// It needs to be split itself.
				
				splitEdgePos -= pos;

				for (sv = pathStart; sv < splitVertices[origSrc].size() - 1; sv++ ) {
					ssize_t splitSrc  = splitVertices[origSrc][sv].parent->value->index;
					ssize_t splitDest = splitVertices[origSrc][sv+1].parent->value->index;
					//UNUSED// ssize_t foundDest = 0;
					ssize_t edgeIndex;
					// Find the out edge that connects src to dest
					edgeIndex = FindEdgeIndex(tSplitVertices[splitSrc], edges, splitDest);
					assert(edgeIndex != -1);

					// Create the new split interval.
					intv.edgePos = splitEdgePos;
					// This edge connects two internal vertices, so the length of 
					// the interval is guaranteed
					// to be the length of the two internal vertices combined.
					intv.length  = edges[edgeIndex].length;
					intv.readPos = splitReadPos;
					intv.pathPos = splitPathPos;
					intv.read    = p;
					
					std::cout << "path: " << p << " " << intv.length << " " << intv.readPos << " " << intv.pathPos << std::endl;
					// Add this to the split edge (later).

					// Add it to the path (later).
					
					// Prepare for the next interval.

					splitEdgePos = 0;
					splitReadPos += splitVertices[origSrc][sv].parent->value->vertexSize;
					splitPathPos++;
				}
			}
		}
	
		// Now loop over the intervals joining edges.
		
	
		// Step 1.  Process the first vertex on this path, if the path
		// starts out in a split vertex.
	
		for (pi = 0; pi < pathLengths[p] - 1; pi++) {
			// cache the edge and index of this path
			pe  = paths[p][pi].edge;
			pei = paths[p][pi].index;

			// This connects two split vertices
			src  = edges[pe].src;
			dest = edges[pe].dest;
			origSrc  = origSrcList[pe];
			origDest = origDestList[pe];

			assert(pe < origEdgeLengthList.size());			
			
			// This is an external edge (it connects the ends of two verties that were connected before
			// splitting.  Adjust the length of it.
			//			(*edges[pe].intervals)[pei].edgePos -= (vertexSize - tSplitVertices[src].vertexSize);
			
			// If this is the last path interval, don't try and automatically add
			// all of the intervals at the dest vertex, because the path interval may not
			// even reach it.
			//if (pi == pathLengths[p] - 1) 
			//				break;

			// Adjust the path interval length to represent the length of the split dest vertex.
			(*edges[pe].intervals)[pei].length -= 
				((vertexSize - tSplitVertices[src].vertexSize) +
				 (vertexSize - tSplitVertices[dest].vertexSize));

			std::cout << "connecting two." << std::endl;
			std::cout << (*edges[pe].intervals)[pei].edgePos << " " << (vertexSize - tSplitVertices[src].vertexSize)
								<< " " << (*edges[pe].intervals)[pei].length << " " 
								<< (vertexSize - tSplitVertices[src].vertexSize) +
				(vertexSize - tSplitVertices[dest].vertexSize) << std::endl;


			splitReadPos += (*edges[pe].intervals)[pei].length;

			// Add edges corresponding to the split vertices
			std::cout << "splitting internal vertex: " << origDest << " " << splitVertices[origDest].size()
								<< std::endl;
			for (sv = 0; sv < splitVertices[origDest].size()-1; sv++ ) {
				ssize_t intSrc, intDest;
				intSrc = splitVertices[origDest][sv].parent->value->index;
				intDest = splitVertices[origDest][sv+1].parent->value->index;

									
				edgeIndex = FindEdgeIndex(tSplitVertices[intSrc], edges, intDest);
				assert(edgeIndex >= 0);

				// Create the new split interval.
				intv.edgePos = 0;
				// This edge connects two internal vertices, so the length of 
				// the interval is guaranteed
				// to be the length of the two internal vertices combined.
					
				intv.length  = edges[edgeIndex].length;
				intv.readPos = splitReadPos;
				intv.pathPos = splitPathPos;
				intv.read    = p;
				
				std::cout << "internal: " << p << " " << intv.length << " " << intv.readPos << " " << intv.pathPos 
									<< std::endl;
				splitReadPos += tSplitVertices[intSrc].vertexSize;
				splitPathPos++;

				// Add this new interval to the path (later).
				// .
			}
		}

		//
		// Now adjust the very last interval.
		//

		// cache the edge and index of this path
		assert(pi == (pathLengths[p]-1));
		pe  = paths[p][pi].edge;
		pei = paths[p][pi].index;

		src  = edges[pe].src;
		dest = edges[pe].dest;
		origSrc  = origSrcList[pe];
		origDest = origDestList[pe];

		ssize_t lengthInDestVertex;
		lengthInDestVertex = (*edges[pe].intervals)[pei].edgePos + (*edges[pe].intervals)[pei].length
			- (edges[pe].length - tSplitVertices[dest].vertexSize);
		
		// This interval doees not overlap with the dest, done modifying it.
		if (lengthInDestVertex <= 0 or lengthInDestVertex <= tSplitVertices[dest].vertexSize) {
			splitReadPos += (*edges[pe].intervals)[pei].length;
			continue;
		}

		// This interval extends a bit into the dest.  

		// If the last interval extends past the end of the split dest, 
		// trim it.
		if (lengthInDestVertex > tSplitVertices[dest].vertexSize) {
			(*edges[pe].intervals)[pei].length -= lengthInDestVertex;
			(*edges[pe].intervals)[pei].length += tSplitVertices[dest].vertexSize;
		}

		// Find out where the vertex interval ends in the next vertex
		ssize_t destPos = 0;
		ssize_t destEnd = -1;
		for (sv = 0; sv < splitVertices[origDest].size(); sv++ ) {
			if (destPos + splitVertices[origDest][sv].parent->value->vertexSize >= lengthInDestVertex) {
				destEnd = sv;
				break;
			}
			destPos += splitVertices[origDest][sv].parent->value->vertexSize;
		}
		
		assert(destEnd != -1);
		// Now add intervals corresponding to the remainder of the path in the last vertex
		for (sv = 0; sv < destEnd; sv++ ) {
			assert(lengthInDestVertex >= 0);
			ssize_t intSrc, intDest;
			intSrc = splitVertices[origDest][sv].parent->value->index;
			intDest = splitVertices[origDest][sv+1].parent->value->index;
			ssize_t edgeIndex;
			edgeIndex = FindEdgeIndex(tSplitVertices[intSrc], edges, intDest);
			assert(edgeIndex >= 0);
			
			intv.readPos = splitReadPos;
			intv.edgePos = 0;
			intv.pathPos = splitPathPos;
			if (lengthInDestVertex < edges[edgeIndex].length)
				intv.length = lengthInDestVertex;
			else
				intv.length = edges[edgeIndex].length;
			
			lengthInDestVertex -= edges[edgeIndex].length;
		}
	}

	// Now fix the intervals that cross past vertices.

	exit(0);

	std::vector<SplitPosList> splitPositions;
	splitPositions.resize(vertices.size());
	SplitOverlappingVertices(vertices, edges, vertexSize, splitPositions);

	DJVertexIntervalsList vertexIntervals;
	StoreIntervalIndices(vertices, edges, splitPositions, vertexIntervals);
	ssize_t numIntervals = CountUniqueVertexIntervals(vertexIntervals);
	
	std::cout << "got: " << numIntervals << " overlapping intervals." << std::endl;

	GlueOverlappingVertices(graph, vertices, edges, vertexIntervals);
	
	//	exit(0);
	for (v = 0; v < splitPositions.size(); v++ ){
		if (vertices[v].OutDegree() > 0) {
			SplitPosList::iterator splitIt;
			std::cout << v << " : ";
			for (splitIt = splitPositions[v].begin();
					 splitIt != splitPositions[v].end();
					 ++splitIt) {
				std::cout << *splitIt << " ";
			}
			std::cout << std::endl;
			ssize_t firstOut;
			ssize_t balOut;
			ssize_t balVertex;
			firstOut = vertices[v].FirstOut();
			balOut = edges[vertices[v].out[firstOut]].balancedEdge;
			balVertex = edges[balOut].dest;
			std::cout << balVertex << " : ";
			for (splitIt = splitPositions[balVertex].begin();
					 splitIt != splitPositions[balVertex].end();
					 ++splitIt) {
				std::cout << *splitIt << " ";
			}
			std::cout << std::endl;
		}
	}
	
	// Each vertex splitting operation gives rise to one 
	// new edge and one new vertex.  Make room for them.
	
	std::cout << "growing graph by " << numInternal 
						<< " vertices and edges." << std::endl;
	

	// Assign edge balance
	//UNUSED// ssize_t firstOutIndex, balOutEdge;
	//UNUSED// ssize_t splitEdge, splitEdgeIndex;
	//UNUSED// ssize_t balSplitEdge, balSplitIndex;
	
	// vertex interval index, balanced vertex interval index
	ssize_t vi, bvi;
	// current vertex index, current balanced vertex index.
	ssize_t cvi, cbvi;
	// next vertex index, next balanced vertex index
	ssize_t nvi, nbvi;

	//UNUSED// ssize_t splitVertexIndex;
	//UNUSED// ssize_t balancedSplitVertexIndex, balancedIntervalIndex;
	//UNUSED+// ssize_t  intvEdgeIndex;
	ssize_t intvEdge, balIntvEdge;
	for (vi = 0; vi < vertexIntervals.size(); vi++) { 
		// don't do anything to vertices that have not been split
		if (vertexIntervals[vi].size() == 0) 
			continue;

		// Join the balance of the edges that connect the new internal vertices
		assert(vertexIntervals[vi].size() > 0);
		
		// Look up the balanced vertex index.
		bvi = vertexBalance[vi];


		assert(vertexIntervals[bvi].size() > 0);

		
		assert(vertexIntervals[vi].size() == vertexIntervals[bvi].size());
		numIntervals = vertexIntervals[vi].size();
		ssize_t i, ei;
		for (i = 0; i < vertexIntervals[vi].size() - 1; i++) {
			cvi = vertexIntervals[vi][i].value->intvVertex;
			nvi = vertexIntervals[vi][i+1].value->intvVertex;
			intvEdge = -1;
			
			for (ei = vertices[cvi].FirstOut();
					 ei != vertices[cvi].EndOut();
					 ei = vertices[cvi].NextOut(ei)) {
				e = vertices[cvi].out[ei];
				if (edges[e].dest == nvi) {
					// Found the out edge that links the two split intervals.
					intvEdge = e;
					break;
				}
			}
			balIntvEdge = -1;
			cbvi = vertexIntervals[bvi][numIntervals - i - 1].value->intvVertex;
			nbvi = vertexIntervals[bvi][numIntervals - i - 2].value->intvVertex;
			for (ei = vertices[cbvi].FirstIn(); 
					 ei != vertices[cbvi].EndIn();
					 ei = vertices[cbvi].NextIn(ei)) {
				e = vertices[cbvi].in[ei];
				if (edges[e].src == nbvi) {
					balIntvEdge = e;
				}
			}
			assert(intvEdge != -1);
			assert(balIntvEdge != -1);
			edges[intvEdge].balancedEdge = balIntvEdge;
			edges[balIntvEdge].balancedEdge = balIntvEdge;
		}
	}
	std::vector<ssize_t> verticesToRemove, edgesToRemove;
	for (v = 0; v < vertexIntervals.size(); v++ ) {
		if (vertexIntervals[v].size() > 0) {
			verticesToRemove.push_back(v);
		}
	}
	graph.Prune(verticesToRemove, edgesToRemove);

	WriteIntervalGraph(graphOutName, graph, 1);
}


