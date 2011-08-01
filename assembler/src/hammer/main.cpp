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
#include <boost/bind.hpp>

#include "read/ireadstream.hpp"
#include "hammer_tools.hpp"
#include "kmer_cluster.hpp"
#include "position_kmer.hpp"


using namespace std;

std::vector<Read> * PositionKMer::rv = NULL;
std::vector<PositionRead> * PositionKMer::pr = NULL;
hint_t PositionKMer::revNo = 0;
hint_t PositionKMer::blob_size = 0;
hint_t PositionKMer::blob_max_size = 0;
char * PositionKMer::blob = NULL;
hint_t * PositionKMer::blobkmers = NULL;
std::vector<uint32_t> * PositionKMer::subKMerPositions = NULL;

struct KMerStatCount {
	PositionKMer km;

	KMerStatCount (uint32_t cnt, hint_t cng) : count(cnt), changeto(cng), km(cng) { }
	uint32_t count;
	hint_t changeto;

	bool isGood() const { return changeto == KMERSTAT_GOOD; }
	bool change() const { return changeto < KMERSTAT_CHANGE; }
};


int main(int argc, char * argv[]) {
	if (argc < 6 || argc > 7) {
		cout << "Usage: ./main tau qvoffset readsFilename dirprefix nthreads [iterno]\n";
		return 0;
	}

	cout << "sizeof( hint_t ) = " << sizeof(hint_t) << endl;
	cout << "sizeof( KMerStat ) = " << sizeof(KMerStat) << endl;
	cout << "sizeof( PositionRead ) = " << sizeof(PositionRead) << endl;
	cout << "sizeof( PositionKMer ) = " << sizeof(PositionKMer) << endl;
	cout << "sizeof( KMerCount ) = " << sizeof(KMerCount) << endl;
	cout << "sizeof( new KMerCount ) = " << sizeof(KMerStatCount) << endl;

	int tau = atoi(argv[1]);
	int qvoffset = atoi(argv[2]);
	
	string readsFilename = argv[3];
	string dirprefix = argv[4];
	int nthreads = atoi(argv[5]);
	
	int iterno = 1; if (argc > 6) iterno = atoi(argv[6]);

	// initialize subkmer positions
	PositionKMer::subKMerPositions = new std::vector<uint32_t>(tau + 2);
	for (uint32_t i=0; i < tau+1; ++i) PositionKMer::subKMerPositions->at(i) = (uint32_t)(i * K / (tau+1) );
	PositionKMer::subKMerPositions->at(tau+1) = K;
	cout << "SubKMer positions: "; for (uint32_t i=0; i < tau+2; ++i) cout << PositionKMer::subKMerPositions->at(i) << " "; cout << endl;

	cout << "Starting work on " << readsFilename << " with " << nthreads << " threads, K=" << K << endl;

	hint_t totalReadSize;
	PositionKMer::rv = ireadstream::readAllNoValidation(readsFilename, &totalReadSize);
	cout << "All reads read to memory." << endl;

	PositionKMer::blob = new char[ (hint_t)(totalReadSize * ( 2 + CONSENSUS_BLOB_MARGIN)) ];
	PositionKMer::blobkmers = new hint_t[ (hint_t)(totalReadSize * ( 2 + CONSENSUS_BLOB_MARGIN)) ];	
	cout << "Allocated blob of size " << (hint_t)(totalReadSize * ( 2 + CONSENSUS_BLOB_MARGIN)) << endl;
	PositionKMer::blob_size = totalReadSize;
	PositionKMer::blob_max_size = (hint_t)(totalReadSize * ( 2 + CONSENSUS_BLOB_MARGIN));
	std::fill( PositionKMer::blobkmers, PositionKMer::blobkmers + PositionKMer::blob_max_size, -1 );
	
	PositionKMer::revNo = PositionKMer::rv->size();
	for (hint_t i = 0; i < PositionKMer::revNo; ++i) {
		Read revcomp = !(PositionKMer::rv->at(i));
		PositionKMer::rv->push_back( revcomp );
	}
	cout << "Reverse complementary reads added.\n";
	
	for (int iter_count = 0; iter_count < iterno; ++iter_count) {
		cout << "\n     === ITERATION " << iter_count << " ===" << endl;

		PositionKMer::pr = new vector<PositionRead>();
		hint_t curpos = 0;

		for (hint_t i = 0; i < PositionKMer::rv->size(); ++i) {
			PositionRead pread(curpos, PositionKMer::rv->at(i).size(), i);
			PositionKMer::pr->push_back(pread);
			for (uint32_t j=0; j < PositionKMer::rv->at(i).size(); ++j) 
				PositionKMer::blob[ curpos + j ] = PositionKMer::rv->at(i).getSequenceString()[j];
			curpos += PositionKMer::rv->at(i).size();
		}
		cout << "Filled up blob. Real size " << curpos << "." << endl;
		PositionKMer::blob_size = curpos;
	
		vector<KMerNo> vv;
		DoPreprocessing(tau, qvoffset, readsFilename, nthreads, &vv);
		cout << "Got " << vv.size() << " kmer positions. Starting parallel sort." << endl;
		
		vector<KMerCount> kmers;
		ParallelSortKMerNos( &vv, &kmers, nthreads );

		// sort ( vv.begin(), vv.end(), KMerNo::less );
		cout << "KMer positions sorted." << endl;

		vector< vector<hint_t> > vs(tau+1);
		vector<SubKMerPQ> vskpq;
		DoSplitAndSort(tau, nthreads, vv, &vs, &kmers, &vskpq);
		vv.clear();

		KMerClustering kmc(kmers, nthreads, tau);
		// prepare the maps
		kmc.process(dirprefix, &vskpq);
		cout << "Finished clustering." << endl;

		// Now for the reconstruction step; we still have the reads in rv, correcting them in place.
		vector<ofstream *> outfv; vector<bool> changed;
		for (int i=0; i<nthreads; ++i) {
			outfv.push_back(new ofstream(dirprefix + "/" + (char)((int)'0' + iter_count) + ".reconstruct." + (char)((int)'0' + i)));
			changed.push_back(false);
		}

		#pragma omp parallel for shared(changed, outfv) num_threads(nthreads)
		for (int i = 0; i < PositionKMer::revNo; ++i) {
			bool res = CorrectRead(kmers, i, outfv[omp_get_thread_num()]);
			changed[omp_get_thread_num()] = changed[omp_get_thread_num()] || res;
		}
		bool res = false;
		for (int i=0; i<nthreads; ++i) {
			outfv[i]->close(); delete outfv[i];
			res = res || changed[i];
		}

		cout << "Correction done." << endl;
		cout << "Printing out reads." << endl;
	
		ofstream outf; outf.open(dirprefix + "/" + (char)((int)'0' + iter_count) + ".reads.corrected");	
		for (hint_t i = 0; i < PositionKMer::revNo; ++i) {
			PositionKMer::pr->at(i).print(outf, qvoffset);
		}
		outf.close();

		// prepare the reads for next iteration
		// delete consensuses, clear kmer data, and restore correct revcomps
		kmers.clear();
		delete PositionKMer::pr;
		std::fill( PositionKMer::blobkmers, PositionKMer::blobkmers + PositionKMer::blob_max_size, -1 );

		PositionKMer::rv->resize( PositionKMer::revNo );
		cout << PositionKMer::rv->size() << ".  " << endl;
		for (hint_t i = 0; i < PositionKMer::revNo; ++i) {
			PositionKMer::rv->push_back( !(PositionKMer::rv->at(i)) );
		}
		cout << "Reads restored." << endl;
	
		if (!res) {
			cout << "Nothing has changed in this iteration. Exiting." << endl;
			break;
		}

	} // iterations

	PositionKMer::subKMerPositions->clear();
	delete PositionKMer::subKMerPositions;
	PositionKMer::rv->clear();
	delete PositionKMer::rv;
	delete [] PositionKMer::blobkmers;
	delete [] PositionKMer::blob;
	return 0;
}


