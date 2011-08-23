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

#include "config_struct_hammer.hpp"
#include "read/ireadstream.hpp"
#include "hammer_tools.hpp"
#include "kmer_cluster.hpp"
#include "position_kmer.hpp"
#include "subkmers.hpp"
#include "globals.hpp"

using namespace std;

std::vector<Read> * PositionKMer::rv = NULL;
std::vector<bool> * PositionKMer::rv_bad = NULL;
std::vector<PositionRead> * PositionKMer::pr = NULL;
hint_t PositionKMer::revNo = 0;
hint_t PositionKMer::lastLeftNo = 0;
hint_t PositionKMer::blob_size = 0;
hint_t PositionKMer::blob_max_size = 0;
char * PositionKMer::blob = NULL;
char * PositionKMer::blobquality = NULL;
hint_t * PositionKMer::blobkmers = NULL;
std::vector<uint32_t> * PositionKMer::subKMerPositions = NULL;

double Globals::error_rate = 0.01;
int Globals::blocksize_quadratic_threshold = 100;
double Globals::good_cluster_threshold = 0.95;
double Globals::blob_margin = 0.25;
int Globals::qvoffset = 64;
bool Globals::paired_reads = false;
int Globals::trim_quality = -1;

struct KMerStatCount {
	PositionKMer km;
	uint32_t count;
	hint_t changeto;

	KMerStatCount (uint32_t cnt, hint_t cng) :  km(cng), count(cnt), changeto(cng) { }

	bool isGood() const { return changeto == KMERSTAT_GOOD; }
	bool change() const { return changeto < KMERSTAT_CHANGE; }
};

string getFilename( const string & dirprefix, const string & suffix ) {
	ostringstream tmp;
	tmp.str(""); tmp << dirprefix.data() << "/" << suffix.data();
	return tmp.str();
}

string getFilename( const string & dirprefix, int iter_count, const string & suffix ) {
	ostringstream tmp;
	tmp.str(""); tmp << dirprefix.data() << "/" << std::setfill('0') << std::setw(2) << iter_count << "." << suffix.data();
	return tmp.str();
}

int main(int argc, char * argv[]) {
	
	string config_file = CONFIG_FILENAME;
	if (argc > 1) config_file = argv[1];
	TIMEDLN("Loading config from " << config_file.c_str());

	cfg::create_instance(config_file);
	string dirprefix = cfg::get().working_dir;
	string readsFilename = cfg::get().reads;
	int tau = cfg::get().tau;
	Globals::qvoffset = cfg::get().quality_offset;
	int nthreads = cfg::get().num_threads;
	int iterno = cfg::get().num_iterations;
	string blobFilename, kmersFilename;
	bool readBlobAndKmers = cfg::get().read_blob_and_kmers;
	bool writeBlobAndKmers = cfg::get().write_blob_and_kmers;
	bool exitAfterWritingBlobAndKMers = false;
	if (readBlobAndKmers || writeBlobAndKmers) {
		blobFilename = cfg::get().blob;
		kmersFilename = cfg::get().kmers;
		exitAfterWritingBlobAndKMers = cfg::get().exit_after_writing_blob_and_kmers;
	}
	Globals::error_rate = cfg::get().error_rate;
	Globals::blocksize_quadratic_threshold = cfg::get().blocksize_quadratic_threshold;
	Globals::good_cluster_threshold = cfg::get().good_cluster_threshold;
	Globals::blob_margin = cfg::get().blob_margin;
	Globals::trim_quality = cfg::get().trim_quality;

	Globals::paired_reads = cfg::get().paired_reads;
	string readsFilenameLeft, readsFilenameRight;
	if (Globals::paired_reads) {
		readsFilenameLeft = cfg::get().reads_left;
		readsFilenameRight = cfg::get().reads_right;
		TIMEDLN("Starting work on " << readsFilenameLeft << " and " << readsFilenameRight << " with " << nthreads << " threads, K=" << K);
	} else {
		TIMEDLN("Starting work on " << readsFilename << " with " << nthreads << " threads, K=" << K);
	}

	// initialize subkmer positions
	PositionKMer::subKMerPositions = new std::vector<uint32_t>(tau + 2);
	for (uint32_t i=0; i < (uint32_t)(tau+1); ++i) {
		PositionKMer::subKMerPositions->at(i) = (i * K / (tau+1) );
	}
	PositionKMer::subKMerPositions->at(tau+1) = K;

	hint_t totalReadSize = 0;

	if (!Globals::paired_reads) {
		PositionKMer::rv = new std::vector<Read>();
		ireadstream::readAllNoValidation(PositionKMer::rv, readsFilename, &totalReadSize, Globals::qvoffset, Globals::trim_quality);
		PositionKMer::lastLeftNo = PositionKMer::rv->size();
	} else {
		PositionKMer::rv = new std::vector<Read>();
		ireadstream::readAllNoValidation(PositionKMer::rv, readsFilenameLeft, &totalReadSize, Globals::qvoffset, Globals::trim_quality);
		PositionKMer::lastLeftNo = PositionKMer::rv->size();
		hint_t rightSize = 0;
		ireadstream::readAllNoValidation(PositionKMer::rv, readsFilenameRight, &rightSize, Globals::qvoffset, Globals::trim_quality);
		totalReadSize += rightSize;
	}

	PositionKMer::blob_size = totalReadSize + 1;
	PositionKMer::blob_max_size = (hint_t)(PositionKMer::blob_size * ( 2 + Globals::blob_margin));

	PositionKMer::blob = new char[ PositionKMer::blob_max_size ];
	PositionKMer::blobquality = new char[ PositionKMer::blob_max_size ];
	PositionKMer::blobkmers = new hint_t[ PositionKMer::blob_max_size ];
	TIMEDLN("Max blob size as allocated is " << PositionKMer::blob_max_size);

	std::fill( PositionKMer::blobkmers, PositionKMer::blobkmers + PositionKMer::blob_max_size, -1 );

	PositionKMer::revNo = PositionKMer::rv->size();
	for (hint_t i = 0; i < PositionKMer::revNo; ++i) {
		string seq = PositionKMer::rv->at(i).getSequenceString();
		Read revcomp = !(PositionKMer::rv->at(i));
		PositionKMer::rv->push_back( revcomp );
	}
	PositionKMer::rv_bad = new std::vector<bool>(PositionKMer::rv->size(), false);

	TIMEDLN("All reads read to memory. Reverse complementary reads added.");

	if (readBlobAndKmers) {
		PositionKMer::readBlob( getFilename(dirprefix, blobFilename.c_str() ).c_str() );
	}

	for (int iter_count = 0; iter_count < iterno; ++iter_count) {
		cout << "\n     === ITERATION " << iter_count << " begins ===" << endl;

		PositionKMer::pr = new vector<PositionRead>();
		hint_t curpos = 0;
		for (hint_t i = 0; i < PositionKMer::rv->size(); ++i) {
			PositionRead pread(curpos, PositionKMer::rv->at(i).size(), i, PositionKMer::rv_bad->at(i));
			PositionKMer::pr->push_back(pread);
			if (!readBlobAndKmers || iter_count > 1) {
				for (uint32_t j=0; j < PositionKMer::rv->at(i).size(); ++j) {
					PositionKMer::blob[ curpos + j ] = PositionKMer::rv->at(i).getSequenceString()[j];
					PositionKMer::blobquality[ curpos + j ] = (char)(Globals::qvoffset + PositionKMer::rv->at(i).getQualityString()[j]);
				}
			}
			curpos += PositionKMer::rv->at(i).size();
		}
		TIMEDLN("Blob done, filled up PositionReads. Real size " << curpos << ". " << PositionKMer::pr->size() << " reads.");

		vector<KMerCount> kmers;
		if (!readBlobAndKmers || iter_count > 0) {
			TIMEDLN("Doing honest preprocessing.");
			PositionKMer::blob_size = curpos;	
			vector<KMerNo> vv;
			DoPreprocessing(tau, readsFilename, nthreads, &vv);
			TIMEDLN("Preprocessing done. Got " << vv.size() << " kmer positions. Starting parallel sort.");
		
			
			ParallelSortKMerNos( &vv, &kmers, nthreads );
			TIMEDLN("KMer positions sorted. In total, we have " << kmers.size() << " kmers.");
			vv.clear();
		} else {
			TIMEDLN("Reading kmers from " << kmersFilename.c_str() );
			PositionKMer::readKMerCounts( getFilename( dirprefix, kmersFilename.c_str() ).c_str(), &kmers );
			TIMEDLN("Kmers read from " << kmersFilename.c_str());
		}


		if ( writeBlobAndKmers && iter_count == 0 ) { // doesn't make sense to overwrite the first blob
			PositionKMer::writeBlob( getFilename(dirprefix, blobFilename.c_str() ).data() );
			PositionKMer::writeKMerCounts( getFilename(dirprefix, kmersFilename.c_str() ).data(), kmers );
			TIMEDLN("Blob and kmers written.");
			if ( exitAfterWritingBlobAndKMers ) break;
		}

		SubKMerSorter * skmsorter = new SubKMerSorter( kmers.size(), &kmers, nthreads, tau, SubKMerSorter::SorterTypeStraight );
		skmsorter->runSort();
		TIMEDLN("Auxiliary subvectors sorted. Starting split kmer processing in " << min(nthreads, tau+1) << " effective threads.");

		ofstream ofkmers( getFilename(dirprefix, iter_count, "kmers.solid").data() );
		ofstream ofkmers_bad( getFilename(dirprefix, iter_count, "kmers.bad").data() );
		KMerClustering kmc(kmers, nthreads, tau);
		// prepare the maps
		kmc.process(dirprefix, skmsorter, &ofkmers, &ofkmers_bad);
		ofkmers.close();
		ofkmers_bad.close();
		delete skmsorter;
		TIMEDLN("Finished clustering. Starting reconstruction.");

		// Now for the reconstruction step; we still have the reads in rv, correcting them in place.
		vector<ofstream *> outfv; vector<hint_t> changedReads; vector<hint_t> changedNucleotides;
		for (int i=0; i<nthreads; ++i) {
			//tmp.str(""); tmp << dirprefix.data() << "/" << std::setfill('0') << std::setw(2) << iter_count << ".reconstruct." << i;
			// outfv.push_back(new ofstream( tmp.str().data() ));
			outfv.push_back(NULL);
			changedReads.push_back(0);
			changedNucleotides.push_back(0);
		}

		#pragma omp parallel for shared(changedReads, changedNucleotides, outfv) num_threads(nthreads)
		for (size_t i = 0; i < PositionKMer::revNo; ++i) {
			bool res = CorrectRead(kmers, i, outfv[omp_get_thread_num()]);
			changedNucleotides[omp_get_thread_num()] += res;
			if (res) ++changedReads[omp_get_thread_num()];
		}
		hint_t totalReads = 0; hint_t totalNucleotides = 0;
		for (int i=0; i<nthreads; ++i) {
			if (outfv[i] != NULL) { outfv[i]->close(); delete outfv[i]; }
			totalReads += changedReads[i];
			totalNucleotides += changedNucleotides[i];
		}

		TIMEDLN("Correction done. Changed " << totalNucleotides << " bases in " << totalReads << " reads. Printing out reads.");

		if (!Globals::paired_reads) {
			outputReads( false, getFilename(dirprefix, iter_count, "reads.corrected").c_str(),
					    getFilename(dirprefix, iter_count, "reads.bad").c_str() );
		} else {
			outputReads( true,  getFilename(dirprefix, iter_count, "reads.left.corrected").c_str(),
					    getFilename(dirprefix, iter_count, "reads.left.bad").c_str(),
					    getFilename(dirprefix, iter_count, "reads.right.corrected").c_str(),
					    getFilename(dirprefix, iter_count, "reads.right.bad").c_str(),
					    getFilename(dirprefix, iter_count, "reads.left.unpaired").c_str(),
					    getFilename(dirprefix, iter_count, "reads.right.unpaired").c_str() );
		}

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

		if (totalReads < 10) {
			TIMEDLN("Too few reads have changed in this iteration. Exiting.");
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


