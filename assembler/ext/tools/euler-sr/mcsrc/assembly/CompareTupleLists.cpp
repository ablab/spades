/***************************************************************************
 * Title:          CompareTupleLists.cpp
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  10/24/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/

#include "ListSpectrum.h"
#include "StringMultTuple.h"
#include "IntegralTupleStatic.h"

#include <string>

void PrintUsage() {
		std::cout << "usage: compareTupleLists readList refList" << std::endl;
		std::cout << " -min m (1)" << std::endl;
		std::cout << " -max m (200)" << std::endl;
}

int main(int argc, char* argv[]) {

	std::string readSpectName, refSpectName;
	if (argc < 3) {
		PrintUsage();
		exit(1);
	}
	readSpectName = argv[1];
	refSpectName  = argv[2];
	int argi = 3;
	ssize_t min = 1;
	ssize_t max = 200;
	while (argi < argc) {
		if (strcmp(argv[argi], "-min") == 0) {
			min = atoi(argv[++argi]);
		}
		else if (strcmp(argv[argi], "-max") == 0) {
			max =atoi(argv[++argi]);
		}
		else {
			PrintUsage();
			std::cout << "bad option: " << argv[argi]<< std::endl;
			exit(1);
		}
		++argi;
	}

	ListSpectrum<StringMultTuple> readSpect, refSpect;
	readSpect.Read(readSpectName, 0);
	refSpect.Read(refSpectName, 0);
	
	std::cout << "threshold, retained , removed, shared but low, readNotRef, refNotRead" << std::endl;

	ssize_t m;
	for (m = min; m < max; m++) {
		ssize_t a, b;
		a = 0;
		b = 0;
		ssize_t readNotRef, refNotRead, removed, retained;
		removed = readNotRef = refNotRead = retained = 0;
		ssize_t removedLow = 0;
		while (a < readSpect.size() and b < refSpect.size()) {
			if (readSpect[a] == refSpect[b]) {
				// Count overlaps 
				if (readSpect[a].GetMult() >= m) {
					retained++;
				}
				else {
					refNotRead++;
					// Keep track of how many were removed due 
					// to too low of coverage in the read set.
					removedLow++;
				}
				a++; 
				b++;
			}
			else if (readSpect[a] < refSpect[b]) {
				if (readSpect[a].GetMult() >= m)
					readNotRef++;
				else
					removed++;
				a++;
			}
			else if (refSpect[b] < readSpect[a]) {
				refNotRead++;
				b++;
			}
			else {
				std::cout << "huh? " << std::endl;
				exit(0);
			}
		}

		// Account for the remainder of a list.
		if (a < readSpect.size()) {
			while (a < readSpect.size()) {
				if (readSpect[a].GetMult() >= m)
					readNotRef++;
				else {
					refNotRead++;
					removed++;
				}
				a++;
			}
		}
		else {
			refNotRead+= (refSpect.size() - b);
		}
		std::cout << m << ", " << retained 
							<< ", " << removed << ", " 
							<< removedLow << ", " 
							<< readNotRef << ", " << refNotRead << std::endl;
	}
	return 0;
}
			
				
	
