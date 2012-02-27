/***************************************************************************
 * Title:          ReadPaths.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  01/19/2010
 *
 * Copyright (c) 2007-2010 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef READ_PATHS_H_
#define READ_PATHS_H_

#include <istream>
#include <fstream>
#include <assert.h>
#include "utils.h"
#include "compatibility.h"

class PathInterval {
 public:
	// HAS BITFIELD
	_SSZT_ edge;
	//	int index:31;
	_SSZT_ index:(SIZE_BITS-1);
	_SSZT_ set:1;
	PathInterval() { edge = -1; index = -1; set = 0;}
	PathInterval(ssize_t e, ssize_t i) {edge= e; index = i; }
};


class Path {
 public:
	PathInterval *interval;
	ssize_t length;
	ssize_t index;
	Path() {
		length = 0;
		interval = NULL;
		index = -1;
	}
	Path& operator=(const Path &rhs) {
		if (this != &rhs) {
			interval = rhs.interval;
			length   = rhs.length;
		}
		return *this;
	}
};

typedef std::vector<PathInterval*> PathIntervalList;
typedef std::vector<Path*> PathList;
typedef std::vector<ssize_t> PathLengthList;
typedef std::vector<std::vector<ssize_t> > PathEdgeList;


class ComparePaths {
 public:

	ssize_t operator()(const Path *pathA , const Path *pathB) {
		assert(pathA);
		assert(pathB);
		ssize_t length = pathA->length;
		if (pathA->length > pathB->length)
			length = pathB->length;

		ssize_t i;
		for (i = 0; i < length; i++ ) {
			if (pathA->interval[i].edge < pathB->interval[i].edge)
				return 1;
		}
		return 0;
	}
};

ssize_t GetLastEdgeIndex(PathIntervalList &paths, PathLengthList &pathLengths, ssize_t pathIndex,
										 ssize_t &lastEdge, ssize_t &lastIntv);

ssize_t GetFirstEdgeIndex(PathIntervalList &paths, PathLengthList &pathLengths, ssize_t pathIndex,
										 ssize_t &lastEdge, ssize_t &lastIntv);

void ReadReadPaths(std::string &pathInFile, PathIntervalList &paths, PathLengthList &pathLengths, std::ostream &report=std::cout);
void WriteReadPaths(std::string &pathOutFile, PathIntervalList &paths, PathLengthList &pathLengths, std::ostream &report=std::cout);
void SortPaths(PathList &paths);
void PathIntervalListToPathList(PathIntervalList &pathIntervals,
																PathLengthList &pathLengths,
																PathList &paths);

void PathIntervalListToPathEdgeList(PathIntervalList &pathIntervals,
																		PathLengthList &pathLengths,
																		PathEdgeList &paths);

ssize_t ArePathsEqual(PathInterval *pathA, ssize_t lengthA, PathInterval *pathB, ssize_t lengthB);

void SplicePath(PathIntervalList &paths, PathLengthList &pathLengths, 
								ssize_t path, ssize_t start, ssize_t end);

void FreePaths(PathIntervalList &paths, PathLengthList &pathLengths);

#endif
