/***************************************************************************
 * Title:          BVertex.h 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  01/19/2010
 *
 * Copyright (c) 2007-2010 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef BVERTEX_H_
#define BVERTEX_H_

#include "ReadPos.h"
#include "graph/GraphAlgo.h"
#include "SeqUtils.h"
#include <iostream>
using namespace std;

#define ALPHABET_SIZE 4

// TODO: deal with sparc compiler warning
//  "warning: `class BVertex' has virtual functions but non-virtual destructor"
class BVertex : public ReadPos, public GraphVertex {
 public:
  ssize_t index;
	static ssize_t WriteFixedDegree;
  std::vector<ssize_t> out;
  std::vector<ssize_t> in;

	// TODO: The only uses of this "length" variable are setting it in the
	// two constructors below, and this line in lib/graph/GraphAlgo.h:
	//       vertices[curVertex].length = curTime;
	// None of those are actual lengths, and the value is never really used.
	// Disabling it.
	//UNUSED?//  int length;

	// HAS BITFIELD
  bool flagged: 1;

  BVertex() {
    out.resize(ALPHABET_SIZE);
    std::fill(out.begin(), out.end(), -1);
    in.resize(ALPHABET_SIZE);
    std::fill(in.begin(), in.end(), -1);
    index = -1;
		//UNUSED?//    length = 0;
  }
  BVertex& operator=(const BVertex &b) {
		if (this != &b) {
			index = b.index;
			out   = b.out;
			in    = b.in;
			//UNUSED?//			length = b.length;
		}
    return *this;
  }

  ssize_t FirstOut() {
    return NextOut(-1);
  }
  ssize_t FirstIn() {
    return NextIn(-1);
  }

  ssize_t EndOut() {
    return out.size();
  }

  ssize_t EndIn() {
    return in.size();
  }
  ssize_t LookupOutIndex(ssize_t outEdgeIndex) {
    ssize_t outEdge = EndOut()-1;
		ssize_t firstOut = FirstOut();
    while (outEdge >= firstOut and out[outEdge] != outEdgeIndex)
      outEdge = PrevOut(outEdge);
		if (outEdge >= firstOut)
      return outEdge;
    else
      return -1;
  }
  void EraseOutEdge(ssize_t outEdgeIndex ) {
    ssize_t outEdge = LookupOutIndex(outEdgeIndex);
    assert(outEdge < EndOut());
    assert(outEdge >= 0);
    out[outEdge] = -1;
  }

  ssize_t LookupInIndex(ssize_t inEdgeIndex) {
    ssize_t inEdge = EndIn()-1;
		ssize_t firstIn = FirstIn();
    while (inEdge >= firstIn and in[inEdge] != inEdgeIndex)
      inEdge = PrevIn(inEdge);
    if (inEdge >= firstIn) 
      return inEdge;
    else 
      return -1;
  }
	
  void EraseInEdge(ssize_t inEdgeIndex ) {
    ssize_t inEdge = LookupInIndex(inEdgeIndex);
    assert(inEdge < EndIn());
    assert(inEdge >= 0);
    in[inEdge] = -1;
  }
    
	ssize_t PrevOut(ssize_t curEdge) {
		if (curEdge < 0)
			return curEdge;
		while (--curEdge >= 0 and out[curEdge] == -1) ;
		return curEdge;
	}

	ssize_t PrevIn(ssize_t curEdge) {
		if (curEdge < 0)
			return curEdge;
		while (--curEdge >= 0 and in[curEdge] == -1) ;
		return curEdge;
	}

  ssize_t NextOut(ssize_t curEdge) {
    ssize_t size = out.size();
    if (curEdge >= size)
      return size;
    while (++curEdge < size and out[curEdge] == -1 ) ; 
    return curEdge;
  }

  ssize_t NextIn(ssize_t curEdge) {
    ssize_t size = in.size();
    if (curEdge >= size) {
      return size;
    }
    while (++curEdge < size and in[curEdge] == -1) ;
    return curEdge;
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

	template<typename E>
		void GetSortedOutEdgeList(std::vector<E> &edges, 
															std::vector<ssize_t> &edgeList) {
		//UNUSED// ssize_t curLength = 0;
		std::vector<std::pair<ssize_t, ssize_t> > lengthIndexList;
		ssize_t outEdge, outEdgeIndex;
		for (outEdgeIndex = FirstOut(); outEdgeIndex < EndOut(); outEdgeIndex = NextOut(outEdgeIndex)) {
			outEdge = out[outEdgeIndex];
			lengthIndexList.push_back(std::pair<ssize_t, ssize_t>(edges[outEdge].length, outEdge));
		}
		std::sort(lengthIndexList.begin(), lengthIndexList.end());
		ssize_t i;
		for (i = 0; i < lengthIndexList.size(); i++) {
			edgeList.push_back(lengthIndexList[i].second);
		}
	}

	template<typename E>
		void GetSortedInEdgeList(std::vector<E> &edges,
														 std::vector<ssize_t> &edgeList) {
		//UNUSED// ssize_t curLength = 0;
		std::vector<std::pair<ssize_t, ssize_t> > lengthIndexList;
		ssize_t inEdge, inEdgeIndex;
		for (inEdgeIndex = FirstIn(); 
				 inEdgeIndex < EndIn(); 
				 inEdgeIndex = NextIn(inEdgeIndex)) {
			inEdge = in[inEdgeIndex];
			lengthIndexList.push_back(std::pair<ssize_t, ssize_t>(edges[inEdge].length, inEdge));
		}
		std::sort(lengthIndexList.begin(), lengthIndexList.end());
		ssize_t i;
		for (i = 0; i < lengthIndexList.size(); i++) {
			edgeList.push_back(lengthIndexList[i].second);
		}
	}
	
	ssize_t AddOutEdge(ssize_t edgeIndex) {
		out.push_back(edgeIndex);
		return out.size()-1;
	}
	
	ssize_t AddInEdge(ssize_t edgeIndex) {
		in.push_back(edgeIndex);
		return in.size()-1;
	}
	ssize_t AddOutEdge(ssize_t edgeIndex,
								 char *edgeSeq,
								 ssize_t edgeLength,
								 int vertexSize) { 
		assert(edgeLength > vertexSize + 1);
		ssize_t outIndex = unmasked_nuc_index[(ssize_t) edgeSeq[vertexSize + 1]];
		// Make sure we are not overwriting anythign
		assert(out[outIndex] == -1);
		out[outIndex] = edgeIndex;
		return outIndex;
	}

	ssize_t AddInEdge(ssize_t edgeIndex,
								char *edgeSeq,
								ssize_t edgeLength,
								int vertexSize) {
		assert(edgeLength > vertexSize + 1);
		ssize_t inIndex = unmasked_nuc_index[(ssize_t) edgeSeq[edgeLength - vertexSize - 1]];
		assert(in[inIndex] == -1);
		in[inIndex] = edgeIndex;
		return inIndex;
	}
	
  friend std::ostream &operator<<(std::ostream &out, const BVertex &vertex);
  friend std::istream &operator>>(std::istream &in, BVertex &vertex);

	ssize_t CheckUniqueEdges(std::vector<ssize_t> &edgeList) {
		std::vector<ssize_t> tmpOut;
		size_t i;
		for (i = 0; i < edgeList.size(); i++ ){ 
			if (edgeList[i] != -1) {
				tmpOut.push_back(edgeList[i]);
			}
		}
		std::sort(tmpOut.begin(), tmpOut.end());
		size_t endOut = tmpOut.size();
		if (tmpOut.size() == 0)
			return 1;
		endOut--;
		for (i = 0; i < endOut; i++) {
			if (tmpOut[i] == tmpOut[i+1]) {
				return 0;
			}
		}
		return 1;
	}
	ssize_t CheckUniqueOutEdges() {
		return CheckUniqueEdges(out);
	}
	ssize_t CheckUniqueInEdges() {
		return CheckUniqueEdges(in);
	}

	void CondenseEdgeLists() {
		ssize_t curPos = 0;
		size_t o, i;
		for (o = 0; o < out.size(); o++ ){
			if (out[o] != -1) {
				out[curPos] = out[o];
				curPos++;
			}
		}
		out.resize(curPos);
		curPos = 0;
		for (i = 0; i < in.size(); i++) {
			if (in[i] != -1) {
				in[curPos] = in[i];
				curPos++;
			}
		}
		in.resize(curPos);
	}

	void Write(ostream &outs) {
		if (WriteFixedDegree) {
			outs << index << " ";
			outs << out[0] << " ";
			outs << out[1] << " ";
			outs << out[2] << " ";
			outs << out[3] << " ";
			outs << in[0] << " ";
			outs << in[1] << " ";
			outs << in[2] << " ";
			outs << in[3];
		}
		else {
			outs << index << " ";
			outs << "out " << out.size() << " ";
			size_t o, i;
			for (o = 0; o < out.size(); o++) {
				outs << out[o] << " ";
			}
			outs << "in " << in.size() << " ";
			for (i = 0; i < in.size(); i++ ){
				outs << in[i] << " ";
			}
		}
	}
};

typedef std::vector<BVertex> BVertexList;

#endif
