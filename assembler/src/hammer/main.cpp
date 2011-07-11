/*
 * main.cpp
 *
 *  Created on: 08.07.2011
 *      Author: snikolenko
 */
 
#include<omp.h>
#include<cmath>
#include<string>
#include<iostream>
#include<sstream>
#include<fstream>
#include<cstdlib>
#include<vector>
#include<map>
#include<list>
#include<queue>
#include<cstdarg>
#include<algorithm>
#include<cassert>

#include "hammer_tools.hpp"

using namespace std;

int main(int argc, char * argv[]) {

	if (argc < 6) {
		cout << "Usage: ./main tau qvoffset readsFilename dirprefix nthreads\n";
		return 0;
	}

	int tau = atoi(argv[1]);
	int qvoffset = atoi(argv[2]);
	
	string readsFilename = argv[3];
	string dirprefix = argv[4];
	int nthreads = atoi(argv[5]);

	cout << "Starting work on " << readsFilename << " with " << nthreads << " threads.\n"; flush(cout);

	vector<KMerStatMap> vv;
	DoPreprocessing(tau, qvoffset, readsFilename, nthreads, vv);
	ReadStatMapContainer rmsc(vv);
	cout << vv[0].size() << "\n";
	cout << "Got RMSC of size " << rmsc.size() << "\n";
	
	vector<StringKMerVector> vs(tau+1);
	vector<KMerCount> kmers;
	DoSplitAndSort(tau, nthreads, rmsc, &vs, &kmers);
	cout << "Got " << kmers.size() << " kmers.\n";
	
	
	DoClustering(tau, nthreads, dirprefix, vs, kmers);
	
	return 0;
}


