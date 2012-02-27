/***************************************************************************
 * Title:          Vertex.h 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  12/15/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef VERTEX_H_
#define VERTEX_H_

#define ALPHABET_SIZE 4
#include "ReadPos.h"

// TODO:
// In Vertex.h, FirstOut(), FirstIn(), EndOut(), EndIn() return the
// actual edge numbers, while FirstOutIndex(), etc., return the index
// of the edge on a given vertex.
// But in BVertex.h and SmallVertexDeBruijn.cpp, FirstOut(),...
// return the index of the edge on the vertex, not the actual edge number.
// This potentially could lead to problems.
// class Vertex is only used in Vertex.h and DeBruijn.cpp.

class Vertex : public ReadPos {
 public:
  ssize_t out[ALPHABET_SIZE];
  ssize_t in[ALPHABET_SIZE];

  ssize_t FirstOut() {
    ssize_t e;
    for (e = 0; e < ALPHABET_SIZE; e++ ) {
      if (out[e] >= 0) 
				return out[e];
    }
		return -1; // TODO: check
  }

  ssize_t FirstOutIndex() {
    ssize_t e;
    for (e = 0; e < ALPHABET_SIZE; e++) 
      if (out[e] >= 0)
				return e;
    return e;
  }
  ssize_t FirstInIndex() {
    ssize_t e;
    for (e = 0; e < ALPHABET_SIZE; e++) 
      if (in[e] >= 0)
				return e;
    return e;
  }

  ssize_t FirstIn() {
    ssize_t e;
    for (e = 0; e < ALPHABET_SIZE; e++) 
      if (in[e] >= 0)
				return in[e];
  }
    
  ssize_t InDegree() {
    ssize_t e, d;
    d = 0;
    for (e = 0; e < ALPHABET_SIZE; e++ )
      (in[e] >= 0) and d++;
    return d;
  }

  ssize_t OutDegree() {
    ssize_t e, d;
    d = 0;
    for (e = 0; e < ALPHABET_SIZE; e++ )
      (out[e] >= 0) and d++;
    return d;
  }

  ssize_t Singleton() {
		ssize_t totalDegree = 0;
		ssize_t i,o;
		for (o = 0; o < ALPHABET_SIZE; o++ ) {
			(out[o] != -1) ? totalDegree++ : totalDegree;
		}
		for (i = 0; i < ALPHABET_SIZE; i++) {
			(in[i] != -1) ? totalDegree++ : totalDegree;
		}
		return totalDegree == 0;
  }

  ssize_t IsBranch() {
    return (InDegree() != 1 or OutDegree() != 1);
  }

  ssize_t DegreeOneOutEdge() {
    // if there is a singe out edge, return the vertex that it goes to
    ssize_t e, e1, d;
    e1 = -1;
    d  = 0;
    for (e = 0; e < ALPHABET_SIZE ; e++ ) {
      if (out[e] >= 0) {
				e1 = e;
				d++;
      }
    }
    if (d == 1) 
      return e1;
    else
      return -1;
  }
    
  Vertex &operator=(const Vertex& val) {
		if (this != &val) {
			out[0] = val.out[0]; 
			out[1] = val.out[1];
			out[2] = val.out[2];
			out[3] = val.out[3];
			in[0]  = val.in[0];
			in[1]  = val.in[1];
			in[2]  = val.in[2];
			in[3]  = val.in[3];
			pos = val.pos;
			read = val.read;
		}
    return *this;
  }
  Vertex() : ReadPos(-1,-1) {
    out[0] = out[1] = out[2] = out[3] = -1;
    in[0] = in[1] = in[2] = in[3] = -1;
  }
};


typedef std::vector<Vertex> VertexList;

#endif
