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

	cout << "Starting work on " << readsFilename << " with " << nthreads << " threads, K=" << K << endl;

	vector<ReadStat> * rv = ireadstream::readAllNoValidation(readsFilename);
	cout << "All reads read to memory." << endl;
	/*for (uint64_t i = 0; i < rv->size(); ++i) {
		cout << "@" << rv->at(i).read.getName() << endl << rv->at(i).read.getSequenceString().data() << endl;
	}*/

	
	for (int iter_count = 0; iter_count < iterno; ++iter_count) {
	
	vector<KMerStatMap> vv;
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
	vector<ofstream *> outfv; vector<bool> changed;
	for (int i=0; i<nthreads; ++i) {
		outfv.push_back(new ofstream(dirprefix + "/" + (char)((int)'0' + iter_count) + ".reconstruct." + (char)((int)'0' + i)));
		changed.push_back(false);
	}

	#pragma omp parallel for shared(changes, good, rv, changed, outfv) num_threads(nthreads)
	for (int i = 0; i < rv->size(); ++i) {
		bool res = CorrectRead(changes, good, kmers, &(rv->at(i)), outfv[omp_get_thread_num()]);
		changed[omp_get_thread_num()] = changed[omp_get_thread_num()] || res;
	}
	bool res = false;
	for (int i=0; i<nthreads; ++i) {
		outfv[i]->close(); delete outfv[i];
		res = res || changed[i];
	}

	cout << "Correction done. Printing out reads." << endl;
	
	ofstream outf; outf.open(dirprefix + "/" + (char)((int)'0' + iter_count) + ".reads.corrected");
	for (uint64_t i = 0; i < rv->size(); ++i) {
		outf << "@" << rv->at(i).read.getName() << endl << rv->at(i).read.getSequenceString().data() << endl << "+" << rv->at(i).read.getName() << endl << rv->at(i).read.getPhredQualityString(qvoffset) << endl;
		// prepare the reads for next iteration
		rv->at(i).kmers.clear(); rv->at(i).kmers_rev.clear();
	}
	outf.close();
	
	cout << "Iteration " << iter_count << " done." << endl;
	if (!res) {
		cout << "Nothing has changed in this iteration. Exiting." << endl;
		break;
	}

	} // iterations
	delete rv;
	return 0;
}


