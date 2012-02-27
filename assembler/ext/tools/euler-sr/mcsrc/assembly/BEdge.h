/***************************************************************************
 * Title:          BEdge.h 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  01/19/2010
 *
 * Copyright (c) 2007-2010 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef BEDGE_H_
#define BEDGE_H_

#include "graph/GraphAlgo.h"

class BEdge : public GraphEdge {
 public:
  SimpleSequence seq;
  //  ReadIntervalList readIntervals;
  ssize_t src;
  ssize_t dest;
  ssize_t multiplicity;
  ssize_t index;
  ssize_t balancedEdge;
  ssize_t length;
	// HAS BITFIELD
  bool flagged      : 1;

	// TODO: check if multiplicity or length can be shortened by one bit
	// to accomodate the one bit field "flagged".
	// Both are in read statements like this one, which could be a problem.
	//    "intvIn.read((char*) &edges[e].length, sizeof(int));"
	// Also would need to check other uses and speed.

  friend std::ostream &operator<<(std::ostream &out, const BEdge &edge) {
    out << edge.src << " " << edge.dest << " " << edge.multiplicity 
				<< " " << edge.index << " " << edge.balancedEdge << " " << edge.length;
    return out;
  }
  friend std::istream &operator>>(std::istream &in, BEdge &edge) {
    in >> edge.src >> edge.dest 
       >> edge.multiplicity >> edge.index >> edge.balancedEdge >> edge.length;
    return in;
  }
	void Init() {}
  BEdge() {
    seq.seq = NULL;
    seq.length = 0;
    src  = -1;
    dest = -1;
    multiplicity = 0;
		index = -1;
    balancedEdge = -1;
		length = 0;
    flagged = GraphVertex::NotMarked;
  }
  void Copy(const BEdge &e) {
		// TODO: Check if other fields (seq, flagged) should be copied too
    src  = e.src;
    dest = e.dest;
    multiplicity = e.multiplicity;
    index = e.index;
    balancedEdge = e.balancedEdge;
    length  = e.length;
  }

  BEdge& ShallowCopy(const BEdge &e) {
    seq.length = e.seq.length;
    seq.seq = e.seq.seq;
    Copy(e);
    return *this;
  }
  BEdge& operator=(const BEdge &e) {

		if (this == &e)
			return *this;

		// TODO: Check why this has been commented out and is incomplete

    /*
      seq.length = e.seq.length;
      seq.seq = new unsigned char[seq.length];
      memcpy(seq.seq, e.seq.seq, seq.length);
      // copy the rest
      Copy(e);
    */
    // for now, just do a shallow copy, hopefully that will speed
    // things up
    return ShallowCopy(e);
  }
  ssize_t Cost() {
		assert(multiplicity > 0);
		return -multiplicity;
  }
  void Cost(ssize_t cost) {
    multiplicity = -cost;
  }

  void Nullify(){ 
    dest = -1;
    src = -1;
    if (seq.seq != NULL) {
      delete [] seq.seq;
      seq.seq = NULL;
      seq.length = 0;
    }
    multiplicity = 0;
    length = 0;
  }
  ssize_t IsNullified() {
    if (dest == -1 and src == -1 and multiplicity == 0 and length == 0)
      return 1;
    else
      return 0;
  }
};

typedef std::vector<BEdge> BEdgeList;


#endif
