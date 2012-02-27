/***************************************************************************
 * Title:          Trace.h 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  11/26/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef PATH_TRACE_H_
#define PATH_TRACE_H_

#include <vector>
#include <stdlib.h>
#include <assert.h>

class PathTrace {
 public:
	// Make resizing the matrix fast.
	std::vector<ssize_t> *edges;
	ssize_t parent;
	ssize_t count;
	ssize_t size() {
		assert(edges != NULL);
		return edges->size();
	}
	PathTrace() {
		count = 0;
		edges = new std::vector<ssize_t>;
		parent = -1;
	}
	PathTrace& operator=(const PathTrace &pt) {
		if (this != &pt) {
			count = pt.count;
			edges = pt.edges;
			parent = pt.parent;
		}
		return *this;
	}
	~PathTrace() {
		// just don't delete the edge list.
	}
};


typedef std::vector<PathTrace> PathTraceList;
class TraceMap {
 public:
	ssize_t trace;
	ssize_t pos;
	TraceMap(ssize_t t, ssize_t p) {
		trace    = t; pos = p;
	}
};

class TraceMapList {
 public:
	std::vector<TraceMap> traces;
	ssize_t numStart;
	ssize_t numInternal;
	ssize_t numEnd;
	ssize_t outResolved, inResolved;
	ssize_t startEdge, endEdge;
	ssize_t inAdjacent, outAdjacent;

	ssize_t GetFirstStartTraceIndex() {
		size_t t;
		for (t = 0; t < traces.size(); t++) { 
			if (traces[t].pos == 0)
				return t;
		}
		assert(t < traces.size());
		return traces.size(); // Never reaches here.  Quiet compiler.
	}

	ssize_t GetFirstEndTraceIndex(PathTraceList &pathTraces) {
		size_t t;
		for (t = 0; t < traces.size(); t++) {
			if ((ssize_t) pathTraces[traces[t].trace].edges->size() 
					== 
					traces[t].pos + 1) {
				return t;
			}
		}
		assert(t < traces.size());
		return traces.size(); // Never reaches here.  Quiet compiler.
	}

	
	TraceMapList() {
		numStart = 0;
		numInternal = 0;
		numEnd = 0;
		outResolved = 0;
		inResolved = 0;
		inAdjacent = 0;  outAdjacent = 0;
		startEdge = endEdge = -1;
	}
};
	
typedef std::vector<TraceMapList> TraceMapMatrix;

#endif
