/*
 * main.cpp
 *
 *  Created on: 08.07.2011
 *      Author: snikolenko
 */
 
#include <omp.h>
#include <cmath>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <vector>
#include <map>
#include <list>
#include <queue>
#include <cstdarg>
#include <algorithm>
#include <cassert>
#include <unordered_set>
#include <boost/bind.hpp>

#include "read/ireadstream.hpp"
#include "hammer_tools.hpp"
#include "kmer_cluster.hpp"
#include "position_kmer.hpp"
#include "subkmers.hpp"

using namespace std;

std::vector<Read> * PositionKMer::rv = NULL;
std::vector<bool> * PositionKMer::rv_bad = NULL;
std::vector<PositionRead> * PositionKMer::pr = NULL;
hint_t PositionKMer::revNo = 0;
hint_t PositionKMer::blob_size = 0;
hint_t PositionKMer::blob_max_size = 0;
char * PositionKMer::blob = NULL;
char * PositionKMer::blobquality = NULL;
hint_t * PositionKMer::blobkmers = NULL;
std::vector<uint32_t> * PositionKMer::subKMerPositions = NULL;

struct KMerStatCount {
	PositionKMer km;
	uint32_t count;
	hint_t changeto;

	KMerStatCount (uint32_t cnt, hint_t cng) :  km(cng), count(cnt), changeto(cng) { }

	bool isGood() const { return changeto == KMERSTAT_GOOD; }
	bool change() const { return changeto < KMERSTAT_CHANGE; }
};

int main(int argc, char * argv[]) {
	if (argc < 6 || argc > 7) {
		cout << "Usage: ./main tau qvoffset readsFilename dirprefix nthreads [iterno]\n";
		return 0;
	}

/*	cout << "sizeof( hint_t ) = " << sizeof(hint_t) << endl;
	cout << "sizeof( KMerStat ) = " << sizeof(KMerStat) << endl;
	cout << "sizeof( PositionRead ) = " << sizeof(PositionRead) << endl;
	cout << "sizeof( PositionKMer ) = " << sizeof(PositionKMer) << endl;
	cout << "sizeof( KMerCount ) = " << sizeof(KMerCount) << endl;
	cout << "sizeof( new KMerCount ) = " << sizeof(KMerStatCount) << endl;
	cout << "sizeof( long double ) = " << sizeof( long double ) << endl;*/

	int tau = atoi(argv[1]);
	int qvoffset = atoi(argv[2]);
	
	string readsFilename = argv[3];
	string dirprefix = argv[4];
	int nthreads = atoi(argv[5]);
	ostringstream tmp;
	
	int iterno = 1; if (argc > 6) iterno = atoi(argv[6]);

	// initialize subkmer positions
	PositionKMer::subKMerPositions = new std::vector<uint32_t>(tau + 2);
	for (uint32_t i=0; i < (uint32_t)(tau+1); ++i) PositionKMer::subKMerPositions->at(i) = (i * K / (tau+1) );
	PositionKMer::subKMerPositions->at(tau+1) = K;

	TIMEDLN("Starting work on " << readsFilename << " with " << nthreads << " threads, K=" << K);

	hint_t totalReadSize;
	PositionKMer::rv = ireadstream::readAllNoValidation(readsFilename, &totalReadSize);
	PositionKMer::rv_bad = new std::vector<bool>(PositionKMer::rv->size(), false);

	PositionKMer::blob_size = totalReadSize;
	PositionKMer::blob_max_size = (hint_t)(totalReadSize * ( 2 + CONSENSUS_BLOB_MARGIN));

	PositionKMer::blob = new char[ PositionKMer::blob_max_size ];
	PositionKMer::blobquality = new char[ PositionKMer::blob_max_size ];
	PositionKMer::blobkmers = new hint_t[ PositionKMer::blob_max_size ];
	TIMEDLN("Max blob size as allocated is " << PositionKMer::blob_max_size);

	std::fill( PositionKMer::blobkmers, PositionKMer::blobkmers + PositionKMer::blob_max_size, -1 );
	std::fill( PositionKMer::blobquality, PositionKMer::blobquality + PositionKMer::blob_max_size, (char)(qvoffset + 2) );
	
	PositionKMer::revNo = PositionKMer::rv->size();
	for (hint_t i = 0; i < PositionKMer::revNo; ++i) {
		Read revcomp = !(PositionKMer::rv->at(i));
		PositionKMer::rv->push_back( revcomp );
		PositionKMer::rv_bad->push_back(false);
	}
	TIMEDLN("All reads read to memory. Reverse complementary reads added.");
	
	for (int iter_count = 0; iter_count < iterno; ++iter_count) {
		cout << "\n     === ITERATION " << iter_count << " begins ===" << endl;

		PositionKMer::pr = new vector<PositionRead>();
		hint_t curpos = 0;

		for (hint_t i = 0; i < PositionKMer::rv->size(); ++i) {
			PositionRead pread(curpos, PositionKMer::rv->at(i).size(), i, PositionKMer::rv_bad->at(i));
			PositionKMer::pr->push_back(pread);
			for (uint32_t j=0; j < PositionKMer::rv->at(i).size(); ++j) {
				PositionKMer::blob[ curpos + j ] = PositionKMer::rv->at(i).getSequenceString()[j];
				PositionKMer::blobquality[ curpos + j ] = PositionKMer::rv->at(i).getQualityString()[j];
			}
			curpos += PositionKMer::rv->at(i).size();
		}
		TIMEDLN("Filled up blob. Real size " << curpos << ". " << PositionKMer::pr->size() << " reads.");
		PositionKMer::blob_size = curpos;
	
		vector<KMerNo> vv;
		DoPreprocessing(tau, qvoffset, readsFilename, nthreads, &vv);
		TIMEDLN("Preprocessing done. Got " << vv.size() << " kmer positions. Starting parallel sort.");
		
		vector<KMerCount> kmers;
		ParallelSortKMerNos( &vv, &kmers, qvoffset, nthreads );
		TIMEDLN("KMer positions sorted. In total, we have " << kmers.size() << " kmers.");
		vv.clear();

		SubKMerSorter * skmsorter = new SubKMerSorter( kmers.size(), &kmers, nthreads, tau, SubKMerSorter::SorterTypeStraight );
		skmsorter->runSort();
		TIMEDLN("Auxiliary subvectors sorted. Starting split kmer processing in " << min(nthreads, tau+1) << " effective threads.");

		tmp.str(""); tmp << dirprefix.data() << "/" << std::setfill('0') << std::setw(2) << iter_count << ".kmers.solid";
		ofstream ofkmers( tmp.str().data() );
		tmp.str(""); tmp << dirprefix.data() << "/" << std::setfill('0') << std::setw(2) << iter_count << ".kmers.bad";
		ofstream ofkmers_bad( tmp.str().data() );
		KMerClustering kmc(kmers, nthreads, tau);
		// prepare the maps
		kmc.process(dirprefix, skmsorter, &ofkmers, &ofkmers_bad);
		ofkmers.close();
		ofkmers_bad.close();
		delete skmsorter;
		TIMEDLN("Finished clustering. Starting reconstruction.");

		// Now for the reconstruction step; we still have the reads in rv, correcting them in place.
		vector<ofstream *> outfv; vector<bool> changed;
		for (int i=0; i<nthreads; ++i) {
			//tmp.str(""); tmp << dirprefix.data() << "/" << std::setfill('0') << std::setw(2) << iter_count << ".reconstruct." << i;
			// outfv.push_back(new ofstream( tmp.str().data() ));
			outfv.push_back(NULL);
			changed.push_back(false);
		}

		#pragma omp parallel for shared(changed, outfv) num_threads(nthreads)
		for (size_t i = 0; i < PositionKMer::revNo; ++i) {
			bool res = CorrectRead(kmers, i, outfv[omp_get_thread_num()]);
			changed[omp_get_thread_num()] = changed[omp_get_thread_num()] || res;
		}
		bool res = false;
		for (int i=0; i<nthreads; ++i) {
			if (outfv[i] != NULL) { outfv[i]->close(); delete outfv[i]; }
			res = res || changed[i];
		}

		TIMEDLN("Correction done. Printing out reads.");
	
		tmp.str(""); tmp << dirprefix << "/" << std::setfill('0') << std::setw(2) << iter_count << ".reads.corrected";
		ofstream outf( tmp.str().data() );
		tmp.str(""); tmp << dirprefix << "/" << std::setfill('0') << std::setw(2) << iter_count << ".reads.bad";
		ofstream outf_bad( tmp.str().data() );
		for (hint_t i = 0; i < PositionKMer::revNo; ++i) {
			if (PositionKMer::rv_bad->at(i)) {
				PositionKMer::pr->at(i).print(outf_bad, qvoffset);
			} else {
				PositionKMer::pr->at(i).print(outf, qvoffset);
			}
		}
		outf.close(); outf_bad.close();

		// prepare the reads for next iteration
		// delete consensuses, clear kmer data, and restore correct revcomps
		kmers.clear();
		delete PositionKMer::pr;
		std::fill( PositionKMer::blobkmers, PositionKMer::blobkmers + PositionKMer::blob_max_size, -1 );

		PositionKMer::rv->resize( PositionKMer::revNo );
		for (hint_t i = 0; i < PositionKMer::revNo; ++i) {
			PositionKMer::rv->push_back( !(PositionKMer::rv->at(i)) );
		}
		TIMEDLN("Reads restored.");

		if (!res) {
			TIMEDLN("Nothing has changed in this iteration. Exiting.");
			break;
		}

	} // iterations

	PositionKMer::subKMerPositions->clear();
	delete PositionKMer::subKMerPositions;
	PositionKMer::rv->clear();
	delete PositionKMer::rv;
	delete PositionKMer::rv_bad;
	delete [] PositionKMer::blobkmers;
	delete [] PositionKMer::blob;
	delete [] PositionKMer::blobquality;
	return 0;
}


