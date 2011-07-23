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
char * PositionKMer::blob = NULL;

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

	uint64_t totalReadSize;
	PositionKMer::rv = ireadstream::readAllNoValidation(readsFilename, &totalReadSize);
	cout << "All reads read to memory." << endl;

	PositionKMer::blob = new char[ (uint64_t)(totalReadSize * ( 2 + CONSENSUS_BLOB_MARGIN)) ];
	cout << "Allocated blob of size " << (uint64_t)(totalReadSize * ( 2 + CONSENSUS_BLOB_MARGIN)) << endl;
	
	PositionKMer::revNo = PositionKMer::rv->size();
	for (uint64_t i = 0; i < PositionKMer::revNo; ++i) {
		ReadStat rs; rs.read = !(PositionKMer::rv->at(i).read);
		PositionKMer::rv->push_back( rs );
	}
	cout << "Reverse complementary reads added.\n";

	uint64_t curpos = 0;
	for (uint64_t i = 0; i < PositionKMer::rv->size(); ++i) {
		PositionKMer::rv->at(i).blobpos = curpos;
		for (uint32_t j=0; j < PositionKMer::rv->at(i).read.size(); ++j) 
			PositionKMer::blob[ curpos + j ] = PositionKMer::rv->at(i).read.getSequenceString()[j];
		curpos += PositionKMer::rv->at(i).read.size();
	}
	cout << "Filled up blob." << endl;
	
	for (int iter_count = 0; iter_count < iterno; ++iter_count) {
		cout << "\n     === ITERATION " << iter_count << " ===" << endl;
	
		vector<KMerStatMap> vv;
		DoPreprocessing(tau, qvoffset, readsFilename, nthreads, &vv);
		ReadStatMapContainer rmsc(vv);
		cout << "Got RMSC of size " << rmsc.size() << "\n";
		
		vector< vector<uint64_t> > vs(tau+1);
		vector<KMerCount> kmers;
		DoSplitAndSort(tau, nthreads, rmsc, &vs, &kmers);
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

		#pragma omp parallel for shared(changed, outfv) num_threads(nthreads)
		for (int i = 0; i < PositionKMer::revNo; ++i) {
			bool res = CorrectRead(kmers, &(PositionKMer::rv->at(i)), &(PositionKMer::rv->at(i + PositionKMer::revNo)), outfv[omp_get_thread_num()]);
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
			outf << "@" << PositionKMer::rv->at(i).read.getName() << endl << PositionKMer::rv->at(i).read.getSequenceString().data() << endl << "+" << PositionKMer::rv->at(i).read.getName() << endl << PositionKMer::rv->at(i).read.getPhredQualityString(qvoffset) << endl;		
		}
		outf.close();

		// prepare the reads for next iteration
		// delete consensuses, clear kmer data, and restore correct revcomps
		PositionKMer::rv->resize( PositionKMer::revNo );
		cout << PositionKMer::rv->size() << "." << endl;
		for (uint64_t i = 0; i < PositionKMer::revNo; ++i) {
			PositionKMer::rv->at(i).kmers.clear();
			ReadStat rs; rs.read = !(PositionKMer::rv->at(i).read);
			PositionKMer::rv->push_back( rs );
		}
		cout << "Reads restored." << endl;
	
		if (!res) {
			cout << "Nothing has changed in this iteration. Exiting." << endl;
			break;
		}

	} // iterations

	//delete [] PositionKMer::rv;
	//delete [] PositionKMer::blob;
	return 0;
}


