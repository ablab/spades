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

// forking
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <sys/stat.h>
#include <fcntl.h>
// file size
#include <sys/stat.h>


using std::string;
using std::vector;
using std::map;

std::vector<std::string> Globals::input_filenames = std::vector<std::string>();
std::vector<hint_t> Globals::input_file_blob_positions = std::vector<hint_t>();
std::vector<uint32_t> * Globals::subKMerPositions = NULL;
std::vector<hint_t> * Globals::kmernos = NULL;
std::vector<KMerCount *> * Globals::kmers = NULL;
int Globals::iteration_no = 0;
hint_t Globals::revNo = 0;
hint_t Globals::revPos = 0;

hint_t Globals::blob_size = 0;
hint_t Globals::blob_max_size = 0;
char * Globals::blob = NULL;
char * Globals::blobquality = NULL;
bool * Globals::blobrc = NULL;

std::vector<PositionRead> * Globals::pr = NULL;

int main(int argc, char * argv[]) {

	string config_file = CONFIG_FILENAME;
	if (argc > 1) config_file = argv[1];
	TIMEDLN("Loading config from " << config_file.c_str());
	cfg::create_instance(config_file);

	// hard memory limit
    const size_t GB = 1 << 30;
    limit_memory(cfg::get().general_hard_memory_limit * GB);

    // input files with reads
    Globals::input_filenames.clear();
    if (cfg::get().input_numfiles > 0) Globals::input_filenames.push_back(cfg::get().input_file_0);
    if (cfg::get().input_numfiles > 1) Globals::input_filenames.push_back(cfg::get().input_file_1);
    if (cfg::get().input_numfiles > 2) Globals::input_filenames.push_back(cfg::get().input_file_2);
    if (cfg::get().input_numfiles > 3) Globals::input_filenames.push_back(cfg::get().input_file_3);
    if (cfg::get().input_numfiles > 4) Globals::input_filenames.push_back(cfg::get().input_file_4);

    // if we need to change single Ns to As, this is the time
    if (cfg::get().general_change_n_to_a) {
    	if (cfg::get().count_do) TIMEDLN("Changing single Ns to As in input read files.");
    	HammerTools::ChangeNtoAinReadFiles();
    	if (cfg::get().count_do) TIMEDLN("Single Ns changed, " << Globals::input_filenames.size() << " read files written.");
    }

    // estimate total read size
    hint_t totalReadSize = HammerTools::EstimateTotalReadSize();
	TIMEDLN("Estimated total size of all reads is " << totalReadSize);

	// allocate the blob
	Globals::blob_size = totalReadSize + 1;
	Globals::blob_max_size = (hint_t)(Globals::blob_size * ( 2 + cfg::get().general_blob_margin));
	Globals::blob = new char[ Globals::blob_max_size ];
	Globals::blobquality = new char[ Globals::blob_max_size ];
	Globals::blobrc = new bool[ Globals::blob_max_size ];
	TIMEDLN("Max blob size as allocated is " << Globals::blob_max_size);


	// initialize subkmer positions
	HammerTools::InitializeSubKMerPositions();

	// now we can begin the iterations
	for (Globals::iteration_no = 0; Globals::iteration_no < cfg::get().general_max_iterations; ++Globals::iteration_no) {
		cout << "\n     === ITERATION " << Globals::iteration_no << " begins ===" << endl;
		bool do_everything = cfg::get().general_do_everything_after_first_iteration && (Globals::iteration_no > 0);

		// read input reads into blob
		Globals::pr = new vector<PositionRead>();
		HammerTools::ReadAllFilesIntoBlob();

		// count k-mers
		if ( cfg::get().count_do || do_everything ) {
			HammerTools::CountKMersBySplitAndMerge();
		} else {
			HammerTools::CountAndSplitKMers(false);
		}

		// sort k-mers
		pid_t pIDsortKmerTotalsFile = 0;
		if ( cfg::get().sort_do || do_everything) {
			pIDsortKmerTotalsFile = vfork();
			if (pIDsortKmerTotalsFile == 0) {
				TIMEDLN("  [" << getpid() << "] Child process for sorting the kmers.total file starting.");
				if (cfg::get().general_gzip) {
					string systemcall = string("gunzip -c ") +
							HammerTools::getFilename(cfg::get().input_working_dir, Globals::iteration_no, "kmers.total").c_str() +
							string(" | sort -k2 -T ") + cfg::get().input_working_dir.c_str() +
							string(" | gzip -1 > ") +
							HammerTools::getFilename(cfg::get().input_working_dir, Globals::iteration_no, "kmers.total.sorted").c_str();
					TIMEDLN("  [" << getpid() << "] " << systemcall);
					if (std::system(systemcall.c_str())) {
						TIMEDLN("  [" << getpid() << "] System error with sorting kmers.total!");
					}
				} else {
					execlp("sort", "sort", "-k2",
							"-o", HammerTools::getFilename(cfg::get().input_working_dir, Globals::iteration_no, "kmers.total.sorted").c_str(),
							"-T", cfg::get().input_working_dir.c_str(),
							HammerTools::getFilename(cfg::get().input_working_dir, Globals::iteration_no, "kmers.total").data(), (char *) 0 );
				}
				_exit(0);
			}
		}

		// initialize k-mer structures
		if (Globals::kmernos) Globals::kmernos->clear(); else Globals::kmernos = new std::vector<hint_t>;
		if (Globals::kmers) Globals::kmers->clear(); else Globals::kmers = new std::vector<KMerCount *>;

//		cout << "   outputting blob of size " << Globals::blob_size << " max size " << Globals::blob_max_size << endl;
//		cout << string(Globals::blob, Globals::revPos) << endl;
//		cout << string(Globals::blob + Globals::revPos, Globals::blob_size - Globals::revPos) << endl;

		// fill in sorted unclustered k-mers
		if ( do_everything || cfg::get().subvectors_do || cfg::get().hamming_do || cfg::get().bayes_do ) {
			if (pIDsortKmerTotalsFile != 0) {
				int childExitStatus;
				waitpid(pIDsortKmerTotalsFile, &childExitStatus, 0);
			}
			TIMEDLN("K-mer sorting done, reading k-mer info from the sorted file.");
			HammerTools::ReadKmerNosFromFile( HammerTools::getFilename(cfg::get().input_working_dir, Globals::iteration_no, "kmers.total.sorted"), Globals::kmernos );
			TIMEDLN("K-mer info read.");
		}

		// fill in already prepared k-mers
		if ( !do_everything && cfg::get().input_read_solid_kmers ) {
			TIMEDLN("Loading k-mers from " << cfg::get().input_solid_kmers );
			HammerTools::ReadKmersWithChangeToFromFile( cfg::get().input_solid_kmers, Globals::kmers, Globals::kmernos );
			TIMEDLN("K-mers loaded.");
		}

		// extract and sort subvectors
		SubKMerSorter * skmsorter = NULL;
		if ( cfg::get().subvectors_do || do_everything ) {
			skmsorter = new SubKMerSorter(Globals::kmernos->size(), Globals::kmernos, 1, cfg::get().general_tau, SubKMerSorter::SorterTypeFileBasedStraight);
			skmsorter->runSort(HammerTools::getFilename(cfg::get().input_working_dir, Globals::iteration_no, "kmers.total.sorted"));
			TIMEDLN("All subvectors sorted.");
		}

		// cluster and subcluster the Hamming graph
		// TODO: refactor into two different procedures
		if ( cfg::get().hamming_do || cfg::get().bayes_do || do_everything ) {
			int clustering_nthreads = min( cfg::get().general_max_nthreads, cfg::get().bayes_nthreads);
			KMerClustering kmc(Globals::kmers, Globals::kmernos, clustering_nthreads, cfg::get().general_tau);
			boost::shared_ptr<FOStream> ofkmers = cfg::get().hamming_write_solid_kmers ?
				FOStream::init(HammerTools::getFilename(cfg::get().input_working_dir, Globals::iteration_no, "kmers.solid").c_str()) : boost::shared_ptr<FOStream>();
			boost::shared_ptr<FOStream> ofkmers_bad = cfg::get().hamming_write_bad_kmers ?
				FOStream::init(HammerTools::getFilename(cfg::get().input_working_dir, Globals::iteration_no, "kmers.bad").c_str()) : boost::shared_ptr<FOStream>();
			kmc.process((cfg::get().hamming_do || do_everything),
					cfg::get().input_working_dir, skmsorter, ofkmers, ofkmers_bad);
			TIMEDLN("Finished clustering.");
		}

		// expand the set of solid k-mers
		if ( cfg::get().expand_do || do_everything ) {
			int expand_nthreads = min( cfg::get().general_max_nthreads, cfg::get().expand_nthreads);
			TIMEDLN("Starting solid k-mers expansion in " << expand_nthreads << " threads.");
			for ( int expand_iter_no = 0; expand_iter_no < cfg::get().expand_max_iterations; ++expand_iter_no ) {
				size_t res = HammerTools::IterativeExpansionStep(expand_iter_no, expand_nthreads, *Globals::kmers);
				TIMEDLN("Solid k-mers iteration " << expand_iter_no << " produced " << res << " new k-mers.");
				if ( res < 10 ) break;
			}
			TIMEDLN("Solid k-mers finalized.");
		}

		// write the final set of k-mers
		if ((cfg::get().expand_write_kmers_result && !cfg::get().input_read_solid_kmers) || do_everything ) {
			ofstream ofkmerstotal(HammerTools::getFilename(cfg::get().input_working_dir, Globals::iteration_no, "kmers.result").data());
			boost::iostreams::filtering_ostream kmer_res;
			if (cfg::get().general_gzip)
				kmer_res.push(boost::iostreams::gzip_compressor(boost::iostreams::gzip_params(boost::iostreams::gzip::best_speed, boost::iostreams::gzip::deflated, 15, 9)));
			kmer_res.push(ofkmerstotal);
			HammerTools::PrintKMerResult( kmer_res, *Globals::kmers );
			TIMEDLN("KMers with changetos written.");
		}

		hint_t totalReads = 0;
		// reconstruct and output the reads
		if ( cfg::get().correct_do || do_everything ) {
			totalReads = HammerTools::CorrectAllReads();
		}

		// prepare the reads for next iteration
		// delete consensuses, clear kmer data, and restore correct revcomps
		for (size_t i=0; i < Globals::kmers->size(); ++i) delete Globals::kmers->at(i);
		Globals::kmers->clear();
		delete Globals::pr;

		if (totalReads < 1) {
			TIMEDLN("Too few reads have changed in this iteration. Exiting.");
			break;
		}
	}

	// clean up
	Globals::subKMerPositions->clear();
	delete Globals::subKMerPositions;
	delete [] Globals::blob;
	delete [] Globals::blobquality;

	TIMEDLN("All done. Exiting.");
	return 0;
}


