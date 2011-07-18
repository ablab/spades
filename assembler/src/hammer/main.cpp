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
#include <unordered_set>

#include "read/ireadstream.hpp"
#include "hammer_config.hpp"
#include "hammer_tools.hpp"
#include "kmer_cluster.hpp"

using namespace std;

int main(int argc, char * argv[]) {
	if (argc < 6 || argc > 7) {
		cout << "Usage: ./main tau qvoffset readsFilename dirprefix nthreads [iterno]\n";
		return 0;
	}

	int tau = atoi(argv[1]);
	int qvoffset = atoi(argv[2]);
	
	string readsFilename = argv[3];
	string dirprefix = argv[4];
	int nthreads = atoi(argv[5]);
	
	int iterno = 1; if (argc > 6) iterno = atoi(argv[6]);

	cout << "Starting work on " << readsFilename << " with " << nthreads << " threads.\n"; flush(cout);

	vector<KMerStatMap> vv;
	vector<ReadStat> * rv = ireadstream::readAllNoValidation(readsFilename);
	cout << "All reads read to memory." << endl;
	/*for (uint64_t i = 0; i < rv->size(); ++i) {
		cout << "@" << rv->at(i).read.getName() << endl << rv->at(i).read.getSequenceString().data() << endl;
	}*/

	
	for (int iter_count = 0; iter_count < iterno; ++iter_count) {
	
	DoPreprocessing(tau, qvoffset, readsFilename, nthreads, &vv, rv);
	ReadStatMapContainer rmsc(vv);
	cout << vv[0].size() << "\n";
	cout << "Got RMSC of size " << rmsc.size() << "\n";
	
	vector<StringKMerVector> vs(tau+1);
	vector<KMerCount> kmers;
	DoSplitAndSort(tau, nthreads, rmsc, &vs, &kmers, rv);
	cout << "Got " << kmers.size() << " kmers.\n";
	// free up memory
	for (uint32_t i=0; i < vv.size(); ++i) vv[i].clear(); vv.clear();
	
	
	// maps to be prepared for reconstruction
	map<KMer, KMer, KMer::less2> changes;
	unordered_set<KMer, KMer::hash> good;

	KMerClustering kmc(kmers, nthreads, tau);
	// prepare the maps
	kmc.process(dirprefix, vs, &changes, &good);
	// free up memory
	kmc.clear();
	cout << "Finished clustering." << endl;
	
	// Now for the reconstruction step; we still have the reads in rv, correcting them in place.
	#pragma omp parallel for shared(changes, good, rv) num_threads(nthreads)
	for (int i = 0; i < rv->size(); ++i) {
		CorrectRead(changes, good, &(rv->at(i).read), false);
	}

	cout << "Correction done. Printing out reads." << endl;
	
	ofstream outf; outf.open(dirprefix + "/" + (char)((int)'0' + iter_count) + ".reads.corrected");
	for (uint64_t i = 0; i < rv->size(); ++i) {
		outf << "@" << rv->at(i).read.getName() << endl << rv->at(i).read.getSequenceString().data() << endl << "+" << rv->at(i).read.getName() << endl << rv->at(i).read.getPhredQualityString(qvoffset) << endl;
	}
	outf.close();
	
	cout << "Iteration " << iter_count << " done." << endl;

	} // iterations
	
	return 0;
}


