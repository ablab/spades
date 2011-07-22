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
#include "position_kmer.hpp"


using namespace std;

std::vector<ReadStat> * PositionKMer::rv = NULL;
uint64_t PositionKMer::revNo = 0;

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
	
	PositionKMer::rv = rv;
	PositionKMer::revNo = rv->size();
	for (uint64_t i = 0; i < PositionKMer::revNo; ++i) {
		ReadStat rs; rs.read = !(rv->at(i).read);
		rv->push_back( rs );
	}
	
	for (int iter_count = 0; iter_count < iterno; ++iter_count) {
		cout << "\n     === ITERATION " << iter_count << " ===" << endl;
	
		vector<KMerStatMap> vv;
		DoPreprocessing(tau, qvoffset, readsFilename, nthreads, &vv, rv);
		ReadStatMapContainer rmsc(vv);
		cout << "Got RMSC of size " << rmsc.size() << "\n";
		
		vector< vector<uint64_t> > vs(tau+1);
		vector<KMerCount> kmers;
		DoSplitAndSort(tau, nthreads, rmsc, &vs, &kmers, rv);
		cout << "Got " << kmers.size() << " kmers.\n";
		// free up memory
		for (uint32_t i=0; i < vv.size(); ++i) vv[i].clear(); vv.clear();
		
		KMerClustering kmc(kmers, nthreads, tau);
		// prepare the maps
		kmc.process(dirprefix, vs);
		// free up memory
		kmc.clear();
		cout << "Finished clustering." << endl;

		// Now for the reconstruction step; we still have the reads in rv, correcting them in place.
		vector<ofstream *> outfv; vector<bool> changed;
		for (int i=0; i<nthreads; ++i) {
			outfv.push_back(new ofstream(dirprefix + "/" + (char)((int)'0' + iter_count) + ".reconstruct." + (char)((int)'0' + i)));
			changed.push_back(false);
		}

		#pragma omp parallel for shared(rv, changed, outfv) num_threads(nthreads)
		for (int i = 0; i < PositionKMer::revNo; ++i) {
			bool res = CorrectRead(kmers, &(rv->at(i)), &(rv->at(i + PositionKMer::revNo)), outfv[omp_get_thread_num()]);
			changed[omp_get_thread_num()] = changed[omp_get_thread_num()] || res;
		}
		bool res = false;
		for (int i=0; i<nthreads; ++i) {
			outfv[i]->close(); delete outfv[i];
			res = res || changed[i];
		}		

		cout << "Correction done. Printing out reads." << endl;
	
		ofstream outf; outf.open(dirprefix + "/" + (char)((int)'0' + iter_count) + ".reads.corrected");	
		for (uint64_t i = 0; i < PositionKMer::revNo; ++i) {
			outf << "@" << rv->at(i).read.getName() << endl << rv->at(i).read.getSequenceString().data() << endl << "+" << rv->at(i).read.getName() << endl << rv->at(i).read.getPhredQualityString(qvoffset) << endl;		
		}
		outf.close();

		// prepare the reads for next iteration
		// delete consensuses, clear kmer data, and restore correct revcomps
		rv->resize( PositionKMer::revNo );
		cout << rv->size() << "." << endl;
		for (uint64_t i = 0; i < PositionKMer::revNo; ++i) {
			rv->at(i).kmers.clear();
			ReadStat rs; rs.read = !(rv->at(i).read);
			rv->push_back( rs );
		}
		cout << "Reads restored." << endl;
	
		if (!res) {
			cout << "Nothing has changed in this iteration. Exiting." << endl;
			break;
		}

	} // iterations

	delete rv;
	return 0;
}


