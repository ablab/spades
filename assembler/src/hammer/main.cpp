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

#include "google/sparse_hash_map"

#include "config_struct_hammer.hpp"
#include "read/ireadstream.hpp"
#include "hammer_tools.hpp"
#include "kmer_cluster.hpp"
#include "position_kmer.hpp"
#include "subkmers.hpp"
#include "globals.hpp"

#include "memory_limit.hpp"

using std::string;
using std::vector;
using std::map;

std::vector<Read> * Globals::rv = NULL;
std::vector<bool> * Globals::rv_bad = NULL;
std::vector<PositionRead> * Globals::pr = NULL;
hint_t Globals::revNo = 0;
hint_t Globals::lastLeftNo = 0;
hint_t Globals::blob_size = 0;
hint_t Globals::blob_max_size = 0;
char * Globals::blob = NULL;
char * Globals::blobquality = NULL;
KMerNoHashMap Globals::hm = KMerNoHashMap();
std::vector<uint32_t> * Globals::subKMerPositions = NULL;

double Globals::error_rate = 0.01;
int Globals::blocksize_quadratic_threshold = 100;
double Globals::good_cluster_threshold = 0.95;
double Globals::blob_margin = 0.25;
int Globals::qvoffset = 64;
bool Globals::paired_reads = false;
int Globals::trim_quality = -1;
bool Globals::trim_left_right = false;
bool Globals::use_iterative_reconstruction = false;
bool Globals::reconstruction_in_full_iterations = false;
double Globals::iterative_reconstruction_threshold = 0.995;
int Globals::max_reconstruction_iterations = 1;

bool Globals::read_kmers_after_clustering = false;
bool Globals::write_kmers_after_clustering = false;
string Globals::kmers_after_clustering = "";

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

string getFilename( const string & dirprefix, int iter_count, const string & suffix, int suffix_num ) {
	ostringstream tmp;
	tmp.str(""); tmp << dirprefix.data() << "/" << std::setfill('0') << std::setw(2) << iter_count << "." << suffix.data() << "." << suffix_num;
	return tmp.str();
}

int main(int argc, char * argv[]) {

    const size_t GB = 1 << 30;
    limit_memory(120 * GB);
	
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
	Globals::trim_left_right = cfg::get().trim_left_right;
	Globals::use_iterative_reconstruction = cfg::get().use_iterative_reconstruction;
	Globals::iterative_reconstruction_threshold = cfg::get().iterative_reconstruction_threshold;
	Globals::max_reconstruction_iterations = cfg::get().max_reconstruction_iterations;
	Globals::reconstruction_in_full_iterations = cfg::get().reconstruction_in_full_iterations;
	Globals::read_kmers_after_clustering = cfg::get().read_kmers_after_clustering;
	Globals::write_kmers_after_clustering = cfg::get().write_kmers_after_clustering;
	Globals::kmers_after_clustering = cfg::get().kmers_after_clustering;

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
	Globals::subKMerPositions = new std::vector<uint32_t>(tau + 2);
	for (uint32_t i=0; i < (uint32_t)(tau+1); ++i) {
		Globals::subKMerPositions->at(i) = (i * K / (tau+1) );
	}
	Globals::subKMerPositions->at(tau+1) = K;

	hint_t totalReadSize = 0;

	if (!Globals::paired_reads) {
		Globals::rv = new std::vector<Read>();
		ireadstream::readAllNoValidation(Globals::rv, readsFilename, &totalReadSize, Globals::qvoffset, Globals::trim_quality);
		Globals::lastLeftNo = Globals::rv->size();
	} else {
		Globals::rv = new std::vector<Read>();
		ireadstream::readAllNoValidation(Globals::rv, readsFilenameLeft, &totalReadSize, Globals::qvoffset, Globals::trim_quality);
		Globals::lastLeftNo = Globals::rv->size();
		hint_t rightSize = 0;
		ireadstream::readAllNoValidation(Globals::rv, readsFilenameRight, &rightSize, Globals::qvoffset, Globals::trim_quality);
		totalReadSize += rightSize;
	}

	Globals::blob_size = totalReadSize + 1;
	Globals::blob_max_size = (hint_t)(Globals::blob_size * ( 2 + Globals::blob_margin));

	Globals::blob = new char[ Globals::blob_max_size ];
	Globals::blobquality = new char[ Globals::blob_max_size ];
	TIMEDLN("Max blob size as allocated is " << Globals::blob_max_size);

	Globals::revNo = Globals::rv->size();
	for (hint_t i = 0; i < Globals::revNo; ++i) {
		string seq = Globals::rv->at(i).getSequenceString();
		Read revcomp = !(Globals::rv->at(i));
		Globals::rv->push_back( revcomp );
	}
	Globals::rv_bad = new std::vector<bool>(Globals::rv->size(), false);

	TIMEDLN("All reads read to memory. Reverse complementary reads added.");

	if (readBlobAndKmers) {
		Globals::readBlob( getFilename(dirprefix, blobFilename.c_str() ).c_str() );
	}

	for (int iter_count = 0; iter_count < iterno; ++iter_count) {
		cout << "\n     === ITERATION " << iter_count << " begins ===" << endl;

		Globals::pr = new vector<PositionRead>();
		hint_t curpos = 0;
		for (hint_t i = 0; i < Globals::rv->size(); ++i) {
			PositionRead pread(curpos, Globals::rv->at(i).size(), i, Globals::rv_bad->at(i));
			Globals::pr->push_back(pread);
			if (!readBlobAndKmers || iter_count > 1) {
				for (uint32_t j=0; j < Globals::rv->at(i).size(); ++j) {
					Globals::blob[ curpos + j ] = Globals::rv->at(i).getSequenceString()[j];
					Globals::blobquality[ curpos + j ] = (char)(Globals::qvoffset + Globals::rv->at(i).getQualityString()[j]);
				}
			}
			curpos += Globals::rv->at(i).size();
		}
		Globals::blob_size = curpos;
		TIMEDLN("Blob done, filled up PositionReads. Real size " << Globals::blob_size << ". " << Globals::pr->size() << " reads.");

		vector<KMerCount*> kmers;
		Globals::hm.clear();
		#ifdef BOOST_UNORDERED_MAP
			Globals::hm.rehash(Globals::blob_size);
		#endif
		#if defined GOOGLE_SPARSE_MAP
			Globals::hm.resize(Globals::blob_size);
		#endif

		if (!readBlobAndKmers || iter_count > 0) {
			TIMEDLN("Doing honest preprocessing.");
			DoPreprocessing(tau, readsFilename, nthreads, &kmers, &Globals::hm);
			TIMEDLN("Preprocessing done. Got " << Globals::hm.size() << " kmers.");
		} else {
			TIMEDLN("Reading kmers from " << kmersFilename.c_str() );
			Globals::readKMerCounts( getFilename( dirprefix, kmersFilename.c_str() ).c_str(), &kmers );
			TIMEDLN("Kmers read from " << kmersFilename.c_str());
		}

		if ( !Globals::read_kmers_after_clustering && writeBlobAndKmers && iter_count == 0 ) { // doesn't make sense to overwrite the first blob
			Globals::writeBlob( getFilename(dirprefix, blobFilename.c_str() ).data() );
			Globals::writeKMerCounts( getFilename(dirprefix, kmersFilename.c_str() ).data(), kmers );
			TIMEDLN("Blob and kmers written.");
			if ( exitAfterWritingBlobAndKMers ) break;
		}

		if ( Globals::read_kmers_after_clustering ) {
			TIMEDLN("Reading clustering results");
			Globals::readKMerHashMap( Globals::kmers_after_clustering.c_str(), &Globals::hm, &kmers );
			TIMEDLN("Clustering results read.");
		} else {
			TIMEDLN("Starting subvectors sort.");
			SubKMerSorter * skmsorter = new SubKMerSorter( kmers.size(), &kmers, nthreads, tau, SubKMerSorter::SorterTypeStraight );
			skmsorter->runSort();
			TIMEDLN("Auxiliary subvectors sorted. Starting split kmer processing in " << min(nthreads, tau+1) << " effective threads.");

			KMerClustering kmc(kmers, nthreads, tau);
			// prepare the maps
			ofstream ofkmersnum( getFilename(dirprefix, iter_count, "kmers.num").data() );
			ofkmersnum << kmers.size() << endl;
			ofkmersnum.close();
			ofstream ofkmers( getFilename(dirprefix, iter_count, "kmers.solid").data() );
			ofstream ofkmers_bad( getFilename(dirprefix, iter_count, "kmers.bad").data() );
			kmc.process(dirprefix, skmsorter, &ofkmers, &ofkmers_bad);
			ofkmers.close();
			ofkmers_bad.close();
			delete skmsorter;
			TIMEDLN("Finished clustering.");
		}

		if ( Globals::write_kmers_after_clustering && iter_count == 0 ) {
			TIMEDLN("Writing k-mers hash after clustering.");
			Globals::writeKMerHashMap( getFilename(dirprefix, iter_count, "kmers.hash").data(), Globals::hm);
			TIMEDLN("K-mers hash written.");
		}

		if ( Globals::use_iterative_reconstruction ) {
			for ( int iter_no = 0; iter_no < Globals::max_reconstruction_iterations; ++iter_no ) {
				//ofstream ofiter( getFilename(dirprefix, iter_count, "reconstruct", iter_no).data() );
				//size_t res = IterativeReconstructionStep(nthreads, &ofiter);
				//ofiter.close();
				size_t res = IterativeReconstructionStep(nthreads, kmers);
				TIMEDLN("Solid k-mers iteration " << iter_no << " produced " << res << " new k-mers.");
				if ( res < 10 ) break;
			}
			TIMEDLN("Solid k-mers finalized.");
		}

		// Now for the reconstruction step; we still have the reads in rv, correcting them in place.
		vector<ofstream *> outfv; vector<hint_t> changedReads; vector<hint_t> changedNucleotides;
		for (int i=0; i<nthreads; ++i) {
			//outfv.push_back(new ofstream( getFilename(dirprefix, iter_count, "reconstruct", i ).data() ));
			outfv.push_back(NULL);
			changedReads.push_back(0);
			changedNucleotides.push_back(0);
		}

		#pragma omp parallel for shared(changedReads, changedNucleotides, outfv) num_threads(nthreads)
		for (size_t i = 0; i < Globals::revNo; ++i) {
			bool res = CorrectRead(Globals::hm, kmers, i, outfv[omp_get_thread_num()]);
			changedNucleotides[omp_get_thread_num()] += res;
			if (res) ++changedReads[omp_get_thread_num()];
			if (res && outfv[omp_get_thread_num()] != NULL) {
				#pragma omp critical
				{
				*(outfv[omp_get_thread_num()]) << "Final result again:  size=" << Globals::rv->at(i).size() << "\n" << Globals::rv->at(i).getSequenceString().c_str() << "\n" << Globals::rv->at(i).getPhredQualityString(Globals::qvoffset).c_str() << endl;
				}
			}
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
		for (size_t i=0; i < kmers.size(); ++i) delete kmers[i];
		kmers.clear();
		Globals::hm.clear();
		delete Globals::pr;
		//std::fill( Globals::blobhash, Globals::blobhash + Globals::blob_max_size, -1 );

		Globals::rv->resize( Globals::revNo );
		for (hint_t i = 0; i < Globals::revNo; ++i) {
			Globals::rv->push_back( !(Globals::rv->at(i)) );
		}
		TIMEDLN("Reads restored.");

		if (totalReads < 1) {
			TIMEDLN("Too few reads have changed in this iteration. Exiting.");
			break;
		}

	} // iterations

	Globals::subKMerPositions->clear();
	delete Globals::subKMerPositions;
	Globals::rv->clear();
	delete Globals::rv;
	delete Globals::rv_bad;
	//delete [] Globals::blobhash;
	delete [] Globals::blob;
	delete [] Globals::blobquality;
	return 0;
}


