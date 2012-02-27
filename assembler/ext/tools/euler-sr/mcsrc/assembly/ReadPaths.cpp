/***************************************************************************
 * Title:          ReadPaths.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  10/31/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "ReadPaths.h"
#include <ostream>
#include <algorithm>

void ReadReadPaths(std::string &pathInFile, PathIntervalList &paths,
									 PathLengthList &pathLengths,
									 std::ostream &report) {
	ssize_t numPaths;
	std::ifstream pathIn;
	openck(pathInFile, pathIn, std::ios::in, report);

	pathIn >> numPaths;
	paths.resize(numPaths);
	pathLengths.resize(numPaths);
	PathInterval *pi;
	ssize_t p;
	ssize_t pathLength;
	ssize_t i;
	for (p = 0; p < numPaths; p++ ) {
		pathIn >> pathLength;
		pathLengths[p] = pathLength;
		ssize_t pathIndex;
		if (pathLength > 0) {
			pi = new PathInterval[pathLength];
			paths[p] = pi;
			for (i = 0; i < pathLength; i++) { 
				pathIn >> pi[i].edge >> pathIndex;
				pi[i].index = pathIndex;
			}
		}
		else 
			paths[p] = NULL;
	}
}

void WriteReadPaths(std::string &pathOutFile, PathIntervalList &paths, PathLengthList &pathLengths,
										std::ostream &report) {
	std::ofstream pathOut;
	openck(pathOutFile, pathOut, std::ios::out, report);
	ssize_t p, i;
	pathOut << paths.size() << std::endl;
	for (p = 0; p < paths.size(); p++ ) {
		pathOut << pathLengths[p]  << " ";
		for (i = 0; i < pathLengths[p]; i++)
			pathOut << " " << paths[p][i].edge <<" " << paths[p][i].index;
		pathOut << std::endl;
	}
}

void SortPaths(PathList &paths) {
	ComparePaths comp;
	std::sort(paths.begin(), paths.begin() + paths.size(), comp);
}

void PathIntervalListToPathList(PathIntervalList &pathIntervals,
																PathLengthList &pathLengths,
																PathList &paths) {
	paths.resize(pathIntervals.size());
	ssize_t i;
	for (i = 0; i < paths.size(); i++) {
		paths[i] = new Path;
		paths[i]->interval = pathIntervals[i];
		paths[i]->index    = i;
		paths[i]->length   = pathLengths[i];
	}
}


ssize_t ArePathsEqual(PathInterval *pathA, ssize_t lengthA, PathInterval *pathB, ssize_t lengthB) {
	if (lengthA != lengthB)
		return 0;

	ssize_t i;
	for (i = 0; i < lengthA; i++) {
		if (pathA[i].edge != pathB[i].edge)
			return 0;
	}
	return 1;
}

void PathIntervalListToPathEdgeList(PathIntervalList &pathIntervals,
																		PathLengthList &pathLengths,
																		PathEdgeList &paths) {
	paths.resize(pathIntervals.size());
	ssize_t i, j;
	for (i = 0; i < paths.size(); i++ ){ 
		paths[i].resize(pathLengths[i]);
		for (j = 0; j < pathLengths[i]; j++ ) 
			paths[i][j] = pathIntervals[i][j].edge;
	}
}

ssize_t GetLastEdgeIndex(PathIntervalList &paths, PathLengthList &pathLengths, ssize_t pathIndex,
										 ssize_t &lastEdge, ssize_t &lastIntv) {
	ssize_t lastIntvIndex = pathLengths[pathIndex] - 1;
	if (lastIntvIndex >= 0) {
		lastEdge = paths[pathIndex][lastIntvIndex].edge;
		lastIntv = paths[pathIndex][lastIntvIndex].index;
		return 1;
	}
	else {
		return 0;
	}
}

ssize_t GetFirstEdgeIndex(PathIntervalList &paths, PathLengthList &pathLengths, ssize_t pathIndex,
											ssize_t &firstEdge, ssize_t &firstIntv) {
	if (pathLengths[pathIndex] > 0) {
		firstEdge = paths[pathIndex][0].edge;
		firstIntv = paths[pathIndex][0].index;
		return 1;
	}
	else {
		return 0;
	}
}


void SplicePath(PathIntervalList &paths, PathLengthList &pathLengths, 
								ssize_t path, ssize_t start, ssize_t end) {
	// Splice path using zero-based half-open coordinates.
	ssize_t s, e;
	for (s = start, e = end; e < pathLengths[path]; e++, s++) {
		paths[path][s] = paths[path][e];
	}
	pathLengths[path] -= end - start;
}

void FreePaths(PathIntervalList &paths, PathLengthList &pathLengths) {
	ssize_t p;
	for (p = 0; p < paths.size();p++) {
		if (pathLengths[p] > 0) delete[] paths[p];
	}
}


