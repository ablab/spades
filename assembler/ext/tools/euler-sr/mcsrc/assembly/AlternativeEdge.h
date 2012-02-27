/***************************************************************************
 * Title:          AlternativeEdge.h
 * Author:         Mark Chaisson
 * Created:        2009
 * Last modified:  11/26/2009
 *
 * Copyright (c) 2009 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef ALTERNATIVE_EDGE_H_
#define ALTERNATIVE_EDGE_H_

#include "SimpleSequence.h"
#include "DNASequence.h"
#include "utils.h"
using namespace std;
class AlternativeEdge {
 public:
	ssize_t offset;
	ssize_t destEdge;
	SimpleSequence seq;
	AlternativeEdge() {
		offset = -1;
		seq.seq = NULL;
		seq.length = 0;
	}
	AlternativeEdge(SimpleSequence &pseq, ssize_t poffset) { 
		//		seq.seq = new unsigned char[pseq.length];
		//		memcpy(seq.seq, pseq.seq, pseq.length);
		seq.seq = pseq.seq;
		seq.length = pseq.length;
		offset = poffset;
	}
	AlternativeEdge(SimpleSequence &pseq, ssize_t pdestedge, ssize_t poffset) { 
		AlternativeEdge(pseq, poffset);
	  destEdge = pdestedge;
	}
	AlternativeEdge &operator=(const AlternativeEdge &rhs) {
		if (this != &rhs) {
			this->offset = rhs.offset;
			this->destEdge = rhs.destEdge;
			this->seq.seq = new unsigned char[rhs.seq.length];
			memcpy(this->seq.seq, rhs.seq.seq, rhs.seq.length);
			this->seq.length = rhs.seq.length;
		}
		return *this;
	}
};

typedef std::vector<AlternativeEdge> AlternativeEdgeVect;

void AppendAlternativeEdges(AlternativeEdgeVect &src, 
														ssize_t offset,
														AlternativeEdgeVect &dest);
	
void AppendAlternativeEdge(AlternativeEdgeVect &vect,
													 SimpleSequence &seq, 
													 ssize_t offset);
void WriteAltEdges(AlternativeEdgeVect &vect,
													 ssize_t edgeIndex, ofstream &seqOut);
ssize_t ReadAltEdge(ifstream &seqIn,
                                                                DNASequence &seq,
                                                                ssize_t &destEdge,
                                                                ssize_t &offset);

#endif 
