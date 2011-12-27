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

hint_t Globals::blob_size = 0;
hint_t Globals::blob_max_size = 0;
char * Globals::blob = NULL;
char * Globals::blobquality = NULL;

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
    if (cfg::get().general_change_n_to_a && cfg::get().count_do) {
    	TIMEDLN("Changing single Ns to As in input read files.");
    	HammerTools::ChangeNtoAinReadFiles();
    	TIMEDLN("Single Ns changed, " << Globals::input_filenames.size() << " read files written.");
    }

    // estimate total read size
    hint_t totalReadSize = HammerTools::EstimateTotalReadSize();
	TIMEDLN("Estimated total size of all reads is " << totalReadSize);

	// allocate the blob
	Globals::blob_size = totalReadSize + 1;
	Globals::blob_max_size = (hint_t)(Globals::blob_size * ( 2 + cfg::get().general_blob_margin));
	Globals::blob = new char[ Globals::blob_max_size ];
	Globals::blobquality = new char[ Globals::blob_max_size ];
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
		}

		// sort k-mers
		pid_t pIDsortKmerTotalsFile = 0;
		if ( cfg::get().sort_do || do_everything) {
			pIDsortKmerTotalsFile = vfork();
			if (pIDsortKmerTotalsFile == 0) {
				TIMEDLN("  [" << getpid() << "] Child process for sorting the kmers.total file starting.");
				execlp("sort", "sort", "-k2",
						"-o", HammerTools::getFilename(cfg::get().input_working_dir, Globals::iteration_no, "kmers.total.sorted").c_str(),
						"-T", cfg::get().input_working_dir.c_str(),
						HammerTools::getFilename(cfg::get().input_working_dir, Globals::iteration_no, "kmers.total").data(), (char *) 0 );
				_exit(0);
			}
		}

		// initialize k-mer structures
		if (Globals::kmernos) Globals::kmernos->clear(); else Globals::kmernos = new std::vector<hint_t>;
		if (Globals::kmers) Globals::kmers->clear(); else Globals::kmers = new std::vector<KMerCount *>;

		// fill in sorted unclustered k-mers
		if ( do_everything || cfg::get().subvectors_do || cfg::get().hamming_do || cfg::get().bayes_do ) {
			if (pIDsortKmerTotalsFile != 0) {
				int childExitStatus;
				waitpid(pIDsortKmerTotalsFile, &childExitStatus, 0);
			}
<<<<<<< HEAD
		}

		if ( Globals::write_kmers_after_clustering && iter_count == 0 ) {
			TIMEDLN("Writing k-mers hash after clustering.");
			Globals::writeKMerHashMap( getFilename(Globals::working_dir, iter_count, "kmers.hash").data(), Globals::hm);
			TIMEDLN("K-mers hash written.");
		}

		if ( Globals::use_iterative_reconstruction && (Globals::skip_iterative < 0) ) {
			if (Globals::conserve_memory) Globals::kmernos = &kmernos;
			for ( int iter_no = 0; iter_no < Globals::max_reconstruction_iterations; ++iter_no ) {
				// ofstream ofs( getFilename(Globals::working_dir, iter_count, "kmers.iterative", iter_no) );
				size_t res = IterativeReconstructionStep(nthreads, kmers, NULL);
				//ofs.close();
				TIMEDLN("Solid k-mers iteration " << iter_no << " produced " << res << " new k-mers.");

				if ( Globals::write_each_iteration_kmers ) {
					ofstream oftmp( getFilename(Globals::working_dir, iter_count, "goodkmers", iter_no ).data() );
					for ( hint_t n = 0; n < kmers.size(); ++n ) {
						if ( kmers[n]->second.isGoodForIterative() ) {
							oftmp << kmers[n]->first.str() << "\n>" << kmers[n]->first.start()
							      << "  cnt=" << kmers[n]->second.count << "  tql=" << (1-kmers[n]->second.totalQual) << "\n";
						}
					}
				}
=======
			TIMEDLN("K-mer sorting done, reading k-mer info from the sorted file.");
			HammerTools::ReadKmerNosFromFile( HammerTools::getFilename(cfg::get().input_working_dir, Globals::iteration_no, "kmers.total.sorted"), Globals::kmernos );
			TIMEDLN("K-mer info read.");
		}

		// fill in already prepared k-mers
		if ( cfg::get().input_read_solid_kmers ) {
			TIMEDLN("Loading k-mers from " << cfg::get().input_solid_kmers );
			HammerTools::ReadKmersAndNosFromFile( cfg::get().input_solid_kmers, Globals::kmers, Globals::kmernos );
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
			ofstream * ofkmers = cfg::get().hamming_write_solid_kmers ?
				new ofstream(HammerTools::getFilename(cfg::get().input_working_dir, Globals::iteration_no, "kmers.solid").c_str()) : NULL;
			ofstream * ofkmers_bad = cfg::get().hamming_write_bad_kmers ?
				new ofstream(HammerTools::getFilename(cfg::get().input_working_dir, Globals::iteration_no, "kmers.bad").c_str()) : NULL;
			kmc.process((cfg::get().hamming_do || do_everything),
					cfg::get().input_working_dir, skmsorter, ofkmers, ofkmers_bad);
			if (ofkmers) ofkmers->close();
			if (ofkmers_bad) ofkmers_bad->close();
			TIMEDLN("Finished clustering.");
		}

		// expand the set of solid k-mers
		if ( cfg::get().expand_do || do_everything ) {
			int expand_nthreads = min( cfg::get().general_max_nthreads, cfg::get().expand_nthreads);
			for ( int expand_iter_no = 0; expand_iter_no < cfg::get().expand_max_iterations; ++expand_iter_no ) {
				size_t res = HammerTools::IterativeExpansionStep(expand_iter_no, expand_nthreads, *Globals::kmers);
				TIMEDLN("Solid k-mers iteration " << expand_iter_no << " produced " << res << " new k-mers.");
>>>>>>> -- first stage of a large refactor: removed extra code, cleaned up main and hammer_tools, changed config file structure
				if ( res < 10 ) break;
			}
			TIMEDLN("Solid k-mers finalized.");
		}

		// write the final set of k-mers
		if ( cfg::get().expand_write_kmers_result || do_everything ) {
			ofstream ofkmerstotal(HammerTools::getFilename(cfg::get().input_working_dir, Globals::iteration_no, "kmers.result").data());
			HammerTools::PrintKMerResult( &ofkmerstotal, *Globals::kmers );
			ofkmerstotal.close();
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


