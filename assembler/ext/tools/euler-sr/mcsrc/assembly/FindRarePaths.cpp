/***************************************************************************
 * Title:          FindRarePaths.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  12/18/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "ReadPaths.h"
#include "utils.h"
#include "IntegralTupleStatic.h"


int main(int argc, char* argv[]) {
	std::string pathsIn, pathsOut;
	ssize_t limit;

	int argi = 1;
	if (argc < 3) {
		std::cout << "usage: findRarePaths pathsIn limit pathsOut" << std::endl;
		exit(0);
	}
	pathsIn = argv[argi++];
	limit   = atoi(argv[argi++]);
	//	pathsOut = argv[argi++];

	PathIntervalList pathIntervals;
	PathLengthList pathLengths;
	PathList paths;
	std::cout << "reading paths."<< std::endl;
	ReadReadPaths(pathsIn, pathIntervals, pathLengths);
	std::cout << "converting to list."<<std::endl;
	PathEdgeList pathEdges;
	PathIntervalListToPathEdgeList(pathIntervals, pathLengths, pathEdges);

	std::sort(pathEdges.begin(), pathEdges.end());
	std::cout << "finding rare. " << pathEdges.size() << std::endl;

	ssize_t i;
	ssize_t count = 1;
	ssize_t numRare = 0;
	for (i = 1; i < pathEdges.size(); i++) {
		if (pathEdges[i] == pathEdges[i-1]) {
			count++;
		}
		else {
			if (count < limit) {
				ssize_t p;
				std::cout << i-1;
				for (p = 0; p < pathEdges[i-1].size(); p++ ) {
					std::cout << " " << pathEdges[i-1][p];
				}
				std::cout << std::endl;
			}
			count = 1;
			++numRare;
		}
	}
	std::cout << numRare << " rare paths out of " << pathEdges.size() << std::endl;
	return 0;
}

