/***************************************************************************
 * Title:          ReadIntervals.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef READ_INTERVALS_H_
#define READ_INTERVALS_H_

#include "compatibility.h"
#include "DeBruijnGraph.h"
#include "ReadPaths.h"
#include "BEdge.h"
#include <algorithm>
#include <string>
#include <list>

using namespace std;
class EdgeInterval {
 public:
  ssize_t edge;
  ssize_t edgePos;
  ssize_t length;
  ssize_t readPos;
  EdgeInterval() {
    edge = -1;
    edgePos  = -1;
    readPos  = -1;
    length = 0;
    
  }
  EdgeInterval& operator=(const EdgeInterval &ei) {
		if (this != &ei) {
			edge    = ei.edge;
			edgePos = ei.edgePos;
			readPos = ei.readPos;
			length  = ei.length;
		}
		return *this;
  }
  void Print(std::ostream &out) {
    std::cout << "edge: " << edge << " edgePos: " << edgePos 
							<< " readPos: " << readPos << " length: " << length;
  }
};

// Read intervals are for edges to map which reads pass through them.

class BareReadInterval {
 public:
  ssize_t read;
  ssize_t readPos;
	ssize_t edgePos;
  ssize_t length;
	BareReadInterval(ssize_t readP, ssize_t readPosP, ssize_t edgePosP, ssize_t lengthP) {
		read    = readP;
		readPos = readPosP;
		edgePos = edgePosP;
		length  = lengthP;
	}
	BareReadInterval() {
		read = readPos = edgePos = length = -1;
	}
	BareReadInterval & operator=(const BareReadInterval &ri) {
		if (this != &ri) {
			read = ri.read;
			readPos = ri.readPos;
			edgePos = ri.edgePos;
			length  = ri.length;
		}
		return *this;
	}
};
typedef std::vector<BareReadInterval> BareReadIntervalList;
typedef std::pair<BareReadInterval*, ssize_t> BareReadIntervalIndex;
typedef std::list<BareReadIntervalIndex> BareReadIntervalIndexList;

class ReadInterval : public BareReadInterval {
 public:
  /* This could probably be 10-12. 
		 The path length is bounded by the genome length (assuming an
		 interval at each position), but is typically much shorter.
  */
	// HAS BITFIELD
  //	int pathPos:27;
  _INT_ pathPos:(_INT_BITS_-5);

	//  
	// There is comfortably enough space for 
	// 12 flags here.
	//
	bool markedForDeletion:1;
	bool traversed:1;
	bool detached:1;
											 

	ReadInterval(ssize_t edgeP, ssize_t readP, ssize_t readPosP, 
							 ssize_t edgePosP, ssize_t lengthP, ssize_t pathPosP) : 
	BareReadInterval(readP, readPosP, edgePosP, lengthP) {
		//		edge = edgeP;
		pathPos = pathPosP;
	}

  ReadInterval() {
		//    edge    = -1;
    read    = -1;
    readPos = -1;
    edgePos = -1;
    length  = -1;
		pathPos = -1;
		traversed = 0;
		markedForDeletion = 0;
		detached = 0;
  }
  void Print(std::ostream &out) {
		//    out << " edge: " << edge << " edgePos " << edgePos << " read " << read 
    out << "edgePos " << edgePos << " read " << read 
				<< " read pos " << readPos << " len: " << length;
  }

  ReadInterval &operator=(const EdgeInterval &ei) {
		//    edge    = ei.edge;
    readPos = ei.readPos;
    edgePos = ei.edgePos;
    length  = ei.length;
    return *this;
  }    

  ReadInterval &operator=(const ReadInterval &ri) {
		if (this != &ri) {
			//    edge    = ri.edge;
			read    = ri.read;
			readPos = ri.readPos;
			edgePos = ri.edgePos;
			length  = ri.length;
			pathPos = ri.pathPos;
			markedForDeletion = ri.markedForDeletion;
			traversed = ri.traversed;
			detached  = ri.detached;
		}
    return *this;
  }
	
	ssize_t IsMarkedForDeletion() {
		return markedForDeletion == 1;
	}
	ssize_t MarkForDeletion() {
		markedForDeletion = 1;
		return markedForDeletion;
	}
};

class CompareBareReadIntervalIndicesByRead {
 public:
	BareReadIntervalList* readIntervalPtr;
	// Return true if interval a < interval b
	ssize_t operator()(const ssize_t &a, const ssize_t &b) {
		if ((*readIntervalPtr)[a].read != (*readIntervalPtr)[b].read) 
			return (*readIntervalPtr)[a].read < (*readIntervalPtr)[b].read;
		if ((*readIntervalPtr)[a].readPos != (*readIntervalPtr)[b].readPos) 
			return (*readIntervalPtr)[a].readPos < (*readIntervalPtr)[b].readPos;
		// interval a = interval b, so < is false
		return 0;
	}
};


class CompareReadIntervalsByReadPos {
 public:
  ssize_t operator()(const ReadInterval &a, const ReadInterval &b) {
    // first sort on edges, if either edge is marked as don't care
		// skip the comparison
		//    if (a.edge != b.edge and a.edge != -1 and b.edge !=-1) return a.edge < b.edge;
    
    // next sort on reads
    if (a.read != b.read) return a.read < b.read;

    // the intervals are from the same read, sort according to the
    // order of the reads
    
    return a.readPos < b.readPos;

    // the edges are equal
    if (a.edgePos != b.edgePos) return a.edgePos < b.edgePos;
  }
};
	

class CompareBareReadIntervalsByEdgePos {
 public:
	ssize_t operator()(const BareReadInterval &a, const BareReadInterval &b) const {
    // next sort on reads
    if (a.read != b.read) return a.read < b.read;

    // the intervals are from the same read, sort according to the
    // order of the reads
		if (a.readPos != b.readPos) return a.readPos < b.readPos;

    // the edges are equal
    return a.edgePos < b.edgePos;
	}
};

class CompareBareReadIntervalsByReadPos {
 public:
	ssize_t operator()(const BareReadInterval &a, const BareReadInterval &b) const {
		if (a.read != b.read) return a.read < b.read;
		
		if (a.edgePos != b.edgePos) return a.edgePos < b.edgePos;
		
		return a.readPos < b.readPos;
	}
};

class CompareReadIntervals {
 public:
  ssize_t operator()(const ReadInterval &a, const ReadInterval &b) {
    // first sort on edges, if either edge is marked as don't care
		// skip the comparison
		//    if (a.edge != b.edge and a.edge != -1 and b.edge !=-1) return a.edge < b.edge;
    
		// next sort on reads
    if (a.read != b.read) return a.read < b.read;

    // the edges are equal
    if (a.edgePos != b.edgePos) return a.edgePos < b.edgePos;

    // the intervals are from the same read, sort according to the
    // order of the reads
    
    return a.readPos < b.readPos;
  }
};

class EdgeIntervalList {
 public:
  static int tupleSize;
  std::vector<EdgeInterval> edgeIntervals;
  EdgeIntervalList() {
    edgeIntervals.resize(0);
  }
  
  ssize_t IncrementEdgeInterval(ssize_t edge, ssize_t edgePos, ssize_t readPos, ssize_t &edgeMultiplicity) {
    ssize_t cur = 0;
    cur  = edgeIntervals.size() - 1;
    if (edgeIntervals.size() > 0 and 
				edgeIntervals[cur].edge == edge and 
				edgeIntervals[cur].edgePos + edgeIntervals[cur].length - tupleSize == edgePos - 1) {
      edgeIntervals[cur].length++;
    }
    else {
      ++cur;
      edgeIntervals.resize(cur + 1);
      edgeIntervals[cur].edge    = edge;
      edgeIntervals[cur].edgePos = edgePos;
      edgeIntervals[cur].length  = tupleSize;
      edgeIntervals[cur].readPos = readPos;
      ++edgeMultiplicity;
    }
		return cur;
  }
};


typedef std::vector<ReadInterval> ReadIntervalList;
typedef std::vector<BareReadInterval> BareReadIntervalList;
typedef std::vector<BareReadIntervalList> BareReadIntervalListList;
typedef std::vector<EdgeIntervalList> EdgeIntervalListList;

ssize_t IncrementReadIntervalList(ssize_t edge, ssize_t edgePos, ssize_t read, ssize_t prevPos, ssize_t readPos,
															ReadIntervalList &readIntervals, int tupleSize,
															ssize_t prevEdge);

class IntervalEdge : public BEdge {
 public:
  ReadIntervalList *intervals;
  IntervalEdge() : BEdge() {}
	void Init() {
		seq.seq = NULL;
		seq.length = 0;
		intervals = new ReadIntervalList;
	}
  void Nullify() {
		intervals = new ReadIntervalList;
    BEdge::Nullify();
  }
  IntervalEdge &operator=(const IntervalEdge &edge) {
		if (this != &edge) {
			// copy parent type information
			BEdge::operator=(edge);
			//    intervals->clear();
			// copy over an entire list of intervals
			intervals = edge.intervals;
		}
		return *this;
  }
	void Clear() {
		if (intervals != NULL) {
			intervals->clear();
			std::vector<ReadInterval>().swap(*intervals);
			delete intervals;
			intervals = NULL;
		}
		((BEdge*)this)->Nullify();
	}
};

void PathReadPosToIntervalIndex(std::vector<BareReadIntervalList> &edgeReadIntervals,
																std::vector<std::vector<ssize_t> > &edgeIntervalIndices,
																PathIntervalList &paths,
																std::vector<ssize_t> &pathLengths);

typedef std::vector<IntervalEdge> IntervalEdgeList;


// sure, this is a pretty strict templatization, but I'm not
// sure how else to do it when I'm using subclasses
template<typename E>
void PrintReadIntervals(std::vector<E> &edges,
												std::string &readIntervalFileName,
												std::ostream &report = std::cout
												) {

  std::ofstream readIntervalOut;
  openck(readIntervalFileName, readIntervalOut, std::ios::out, report);

  //UNUSED// ssize_t v, e;
  ssize_t edgeIndex;
  ssize_t intervalIndex;
  intervalIndex = 0;
  for (edgeIndex = 0; edgeIndex < edges.size(); edgeIndex++ ) {
    if (edges[edgeIndex].IsNullified() == 1)
      continue;
    readIntervalOut << "EDGE " << edgeIndex 
										<< " Length " << edges[edgeIndex].length
										<< " Multiplicity " << edges[edgeIndex].intervals->size()
										<< std::endl;
    for (intervalIndex = 0; intervalIndex < edges[edgeIndex].intervals->size();
				 intervalIndex++) {
      readIntervalOut << "INTV" 
											<< " " << (*edges[edgeIndex].intervals)[intervalIndex].read 
											<< " " << (*edges[edgeIndex].intervals)[intervalIndex].readPos 
											<< " " << (*edges[edgeIndex].intervals)[intervalIndex].length
											<< " " << (*edges[edgeIndex].intervals)[intervalIndex].edgePos << std::endl;
    }
  }
  readIntervalOut.close();
}


template<typename V, typename E>
ssize_t StoreIntervalPathForward(E &edge, ssize_t intervalIndex, 
														 std::vector<V> &vertices, std::vector<E> &edges,
														 ssize_t vertexLength, std::vector<ssize_t> &path) {
  /* 
     Given an edge and a read interval on that edge, search
     back for the entire length of the path.  Add each edge that gets traced back 
     to 'path'.
  */

  ssize_t edgeOut;
  ssize_t destVertex;
  ssize_t nextEdge;
  ssize_t curEdge = edge.index;
  ssize_t nextIntervalIndex;

  /*
    std::cout << "edge: " << edge.index << " interval: " << intervalIndex 
    << " re:" << edge.intervals[intervalIndex].read
    << " rp:" << edge.intervals[intervalIndex].readPos
    << " le:" << edge.intervals[intervalIndex].length
    << " ep:" << edge.intervals[intervalIndex].edgePos << std::endl;
  */
  do {
    destVertex = edges[curEdge].dest;
    // Look through each in edge for the continuation of this read path
    nextIntervalIndex = -1;
    for (edgeOut = vertices[destVertex].FirstOut();
				 edgeOut < vertices[destVertex].EndOut();
				 edgeOut = vertices[destVertex].NextOut(edgeOut)) {
      nextEdge = vertices[destVertex].out[edgeOut];
      //  std::cout << "  search succ edge " << nextEdge << std::endl;
      nextIntervalIndex = GetSuccEdgeReadIntvIndex((*edges[curEdge].intervals)[intervalIndex],
																									 edges[nextEdge], vertexLength);
      if (nextIntervalIndex != -1) {
				path.push_back(nextEdge);
				break;
      }
    }
    curEdge = nextEdge;
    intervalIndex = nextIntervalIndex;
  }
  while (intervalIndex != -1 and edgeOut < 4);
  return path.size();
}

template<typename E>
ssize_t GetSuccEdgeReadIntvIndex(ReadInterval &intv, E &edge, ssize_t vertexLength) {
  ssize_t succIntvIndex;
  ssize_t intervalEnd;
	succIntvIndex = GetReadFirstIndex(intv.read, edge);
  if (succIntvIndex >= 0) {
    while ((*edge.intervals)[succIntvIndex].read == intv.read) {
      intervalEnd = intv.readPos + intv.length;
      if (intervalEnd - vertexLength == (*edge.intervals)[succIntvIndex].readPos )
				return succIntvIndex;
      succIntvIndex++;
    }
  }
  // no succeeding read index found
  return -1;
}



template<typename V, typename E>
	ssize_t StoreIntervalPathReverse(E &edge, ssize_t intervalIndex, 
															 std::vector<V> &vertices, std::vector<E> &edges,
															 ssize_t vertexLength, std::vector<ssize_t> &path) {
  /* 
     Given an edge and a read interval on that edge, search
     back for the entire length of the path.  Add each edge that gets traced back 
     to 'path'.
  */
  ssize_t edgeIn;
  ssize_t srcVertex;
  ssize_t prevEdge, curEdge;
  ssize_t prevIntervalIndex;
  curEdge = edge.index;
  if ((*edge.intervals)[intervalIndex].readPos == 0)
    return 0;

  do {
    srcVertex = edges[curEdge].src;
    prevIntervalIndex = -1;
    // Look through each in edge for the continuation of this read path
    for (edgeIn = vertices[srcVertex].FirstIn();
				 edgeIn < vertices[srcVertex].EndIn();
				 edgeIn = vertices[srcVertex].NextIn(edgeIn)) {
      prevEdge = vertices[srcVertex].in[edgeIn];
      prevIntervalIndex = GetPrecEdgeReadIntvIndex((*edges[curEdge].intervals)[intervalIndex],
																									 edges[prevEdge], vertexLength);
      if (prevIntervalIndex != -1) {
				path.push_back(prevEdge);
				break;
      }
    }
    curEdge = prevEdge;
    intervalIndex = prevIntervalIndex;
  }
  while (prevIntervalIndex != -1 and
				 edgeIn < 4 and
				 (*edges[curEdge].intervals)[intervalIndex].readPos != 0);
  return path.size();
}

template<typename E>
ssize_t GetPrecEdgeReadIntvIndex(ReadInterval &intv, E &edge, ssize_t vertexLength) {
  ssize_t precIntvIndex;
  ssize_t intervalEnd;
  //stopfn(intv.edgePos);
  precIntvIndex = GetReadFirstIndex(intv.read, edge);

  if (precIntvIndex >= 0) {
    while ((*edge.intervals)[precIntvIndex].read == intv.read) {
      intervalEnd = (*edge.intervals)[precIntvIndex].readPos + 
				(*edge.intervals)[precIntvIndex].length;
      if (intervalEnd - vertexLength == intv.readPos )
				return precIntvIndex;
      precIntvIndex++;
    }
  }
  // no preceeding read index found
  return -1;
}



ssize_t CountEdges(std::string &readIntervalFileName);

template<typename E>
ssize_t BinaryPrintReadIntervals(std::string &intervalFileName,
														 std::vector<E> &edges,
														 std::ostream &report = std::cout
														 ) {
	std::ifstream intvIn;
	openck(intervalFileName, intvIn, std::ios::in | std::ios::binary, report);
	ssize_t numEdges;
	intvIn.read((char*) &numEdges, sizeof(ssize_t));
	ssize_t e;
	for (e = 0; e < numEdges; e++) {
		intvIn.read((char*) &edges[e].length, sizeof(ssize_t));
		intvIn.read((char*) &edges[e].multiplicity, sizeof(ssize_t));
		(*edges[e].intervals).resize(edges[e].multiplicity);
		intvIn.read((char*) &((*edges[e].intervals)[0]), sizeof(ReadInterval) * edges[e].multiplicity);
	}
	return 1;
}


template<typename E>
ssize_t ReadReadIntervals(std::string &readIntervalFileName,
											std::vector<E> &edges,
											std::ostream &report = std::cout
											) {
  std::ifstream readIntervalIn;
  openck(readIntervalFileName, readIntervalIn, std::ios::in, report);
  std::string word;
  ssize_t value;
  ssize_t curEdge = -1;
  ssize_t curIntv = -1;
	ssize_t maxReadIndex = -1;
	//UNUSED// char *linePtr;
	//UNUSED// ssize_t read, readPos, length, edgePos;
	string line;
	//UNUSED// ssize_t edgeLength, multiplicity;
  while(readIntervalIn) {
    if (! (readIntervalIn >> word)) break;
    if (word == "EDGE") {
			if (curEdge % 10000 ==9999) {
				std::cout << "edge: " << curEdge + 1 << std::endl;
			}
      ++curEdge;
      if (! (readIntervalIn >> value)) { 
				std::cout << "Should have edge number " << std::endl;
				exit(1);
      }
      if (! (readIntervalIn >> word)) {
				std::cout << "Should have 'Length' token " << std::endl;
				exit(1);
      }
      if (word != "Length" ){ 
				std::cout << "Badly formatted interval file " << word 
									<< " should be 'Length'"<<std::endl;
				exit(1);
      }
      if (!(readIntervalIn >> edges[curEdge].length)) {
				std::cout << "Should specify length " << std::endl;
				exit(1);
      }
      if (!(readIntervalIn >> word)) {
				std::cout << "Should have 'Multiplicity' token" << std::endl;
				exit(1);
      }
      if (word != "Multiplicity") {
				std::cout << "Badly formatted interval file " << word
									<< " should be 'Multiplicity'" << std::endl;
				exit(1);
      }
      readIntervalIn >> edges[curEdge].multiplicity;
      (*edges[curEdge].intervals).resize(edges[curEdge].multiplicity);
      curIntv = 0;
    }
    else if (word == "INTV") {
      assert(curIntv < edges[curEdge].multiplicity);
			ssize_t read, readPos, length, edgePos;
			readIntervalIn >> read;
			readIntervalIn >> readPos;
			readIntervalIn >> length;
			readIntervalIn >> edgePos;
      (*edges[curEdge].intervals)[curIntv].read = read;
			(*edges[curEdge].intervals)[curIntv].readPos = readPos;
			(*edges[curEdge].intervals)[curIntv].length = length;
			(*edges[curEdge].intervals)[curIntv].edgePos = edgePos;
			/*			
      readIntervalIn >> (*edges[curEdge].intervals)[curIntv].read
										 >> (*edges[curEdge].intervals)[curIntv].readPos 
										 >> (*edges[curEdge].intervals)[curIntv].length
										 >> (*edges[curEdge].intervals)[curIntv].edgePos;
			*/
			assert((*edges[curEdge].intervals)[curIntv].length >= 0);
			if ((*edges[curEdge].intervals)[curIntv].read > maxReadIndex)
				maxReadIndex = (*edges[curEdge].intervals)[curIntv].read;
      curIntv++;
    }
    else {
      std::cout << "Badly formatted interval file: " << word << std::endl;
      exit(1);
    }
  }
	if (maxReadIndex < 0) {
		std::cout << "Error, did not read in any paths, halting assembly." << std::endl;
		exit(0);
	}
	return maxReadIndex;
}

void SetReadPathIntervals(ReadIntervalList &readIntervals,
													PathIntervalList &paths,
													PathLengthList &pathLengths);

void PrintReadIntervals(ReadIntervalList &readIntervals,
												SimpleSequenceList &edges,
												std::vector<ssize_t> &mult,
												std::string &readIntervalFileName);

void StoreAllReadIntervals(ReadPositions &dbEdges, _INT_ dbEdgeSize, 
													 SimpleSequenceList &edges, 
													 std::string &sequenceListFileName, 
													 std::vector<BareReadIntervalList> &edgeReadIntervals,
													 PathIntervalList &paths,
													 std::vector<ssize_t> &pathLengths,
													 _INT_ skipGapped);

void StoreReadIntervals(EdgeIntervalListList &edgeIntervals,
												ReadIntervalList &readIntervals);

void StoreEdgeInterval(ReadPositions &overlaps, _INT_ overlapSize,
											 SimpleSequenceList &edges,
											 SimpleSequence &sequence,
											 EdgeIntervalList &edgeInterval,
											 std::vector<ssize_t> &mult);


void StoreEdgeIntervals(ReadPositions &overlaps, _INT_ overlapSize,
												SimpleSequenceList &edges,
												SimpleSequenceList &sequences,
												EdgeIntervalListList &edgeIntervals,
												PathIntervalList &paths,
												PathLengthList &pathLengths,
												std::vector<ssize_t> &mult);

void EdgeToReadIntervals(EdgeIntervalListList &edgeIntervals,
												 ReadIntervalList &readIntervals);


void SortReadIntervalsByReadPos(ReadIntervalList &list);
void SortReadIntervals(ReadIntervalList &list);
void SortBareReadIntervalsByReadPos(BareReadIntervalList &list);
void SortBareReadIntervalsByEdgePos(BareReadIntervalList &list);
void SortReadIntervalIndicesByRead(BareReadIntervalList &list, 
																	 std::vector<ssize_t> &listIndices);


class CompareReadIntervalRead {
 public:
  ssize_t operator()(const ReadInterval &read1, const ssize_t readIndex) {
    return read1.read < readIndex;
  }
};


template<typename E>
ssize_t TraceReadIntervalReverse(E &curEdge, ssize_t readInterval, E &prevEdge, int vertexSize) {
  ssize_t prevIntvIndex;
  prevIntvIndex = GetReadFirstIndex((*(curEdge.intervals))[readInterval].read, prevEdge);
  if (prevIntvIndex == -1)
    return -1;

  // Now try and match the end of this segment of the read interval with 
  // the previous read interval.

  while (prevIntvIndex < (*prevEdge.intervals).size() and
				 (*prevEdge.intervals)[prevIntvIndex].read == (*curEdge.intervals)[readInterval].read) {
    if ((*prevEdge.intervals)[prevIntvIndex].readPos +
				(*prevEdge.intervals)[prevIntvIndex].length - vertexSize ==
				(*curEdge.intervals)[readInterval].readPos)
      return prevIntvIndex;
    ++prevIntvIndex;
  }
  return -1;
}

template<typename E>
ssize_t TraceReadIntervalForward(E &curEdge, ssize_t readInterval, E &nextEdge, int vertexSize) {
  ssize_t nextIntvIndex;
  /*
    std::cout << "tracing next " << readInterval << " " << curEdge.index  
    << " " << nextEdge.index << std::endl;
  */
	// Locate the first time this read occurs in this edge.  The read intervals
	// should be sorted according to read index on the edge, so if this read maps to the edge 
	// several times the intervals will all be contiguous.

  nextIntvIndex = GetReadFirstIndex((*curEdge.intervals)[readInterval].read, nextEdge);
  //  std::cout << " next: " << nextIntvIndex << std::endl;
  if (nextIntvIndex == -1)
    return -1;

  // Now try and match the end of this segment of the read interval with 
  // the nextious read interval.
  ssize_t curEnd = (*curEdge.intervals)[readInterval].readPos + 
    (*curEdge.intervals)[readInterval].length;
  while (nextIntvIndex < (*nextEdge.intervals).size() and
				 (*nextEdge.intervals)[nextIntvIndex].read == (*curEdge.intervals)[readInterval].read) {
    if (curEnd == (*nextEdge.intervals)[nextIntvIndex].readPos + vertexSize)
      return nextIntvIndex;
    ++nextIntvIndex;
  }
  return -1;
}

template<typename E>
ssize_t GetReadFirstIndex(ssize_t read, E &edge) {
  //UNUSED// ssize_t first, last, cur;
  ReadIntervalList::iterator readPos;
  CompareReadIntervalRead comp;
  ssize_t readIndex;
  readPos = std::lower_bound(edge.intervals->begin(), edge.intervals->end(), read, comp);
  if (readPos == edge.intervals->end())
    return -1;
  else {
    readIndex = readPos - edge.intervals->begin();
    if ((*edge.intervals)[readIndex].read == read)
      return readIndex;
    else
      return -1;
  }
}


ssize_t  CountEdgeIntervals(EdgeIntervalListList &edgeIntervals);

class EdgeMappedReadPos : public ReadPos {
 public:
  ssize_t edge;
  ssize_t edgePos;
  EdgeMappedReadPos() : ReadPos() {
    edge = -1;
    edgePos  = -1;
  }
  friend std::ostream &operator<<(std::ostream &out, const EdgeMappedReadPos &rp) {
    out << (ReadPos&) rp << " " << rp.edge << " " << rp.edgePos;
    return out;
  }
  friend std::istream &operator>>(std::istream &in, EdgeMappedReadPos &rp) {
    in >> (ReadPos&) rp >> rp.edge >> rp.edgePos;
    return in;
  }
};

typedef std::vector<EdgeMappedReadPos> EdgeMappedReadPositions;



#endif
