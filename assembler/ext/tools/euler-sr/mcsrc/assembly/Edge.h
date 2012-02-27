/***************************************************************************
 * Title:          Edge.h 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  11/26/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef EDGE_H_
#define EDGE_H_

class Edge {
 public:
  ssize_t dest;
  Edge() {
    dest = -1;
  }
  Edge& operator=(const Edge &e) {
		if (this != &e) {
			dest = e.dest;
		}
		return *this;
  }
};

#endif
