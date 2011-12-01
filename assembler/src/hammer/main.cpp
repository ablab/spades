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
std::vector<hint_t> * Globals::kmernos = NULL;

std::string Globals::working_dir = "";
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
double Globals::special_nonsingleton_threshold = 0;
int Globals::max_reconstruction_iterations = 1;

bool Globals::conserve_memory = false;
int Globals::num_of_tmp_files = 500;
int Globals::iteration_no = 0;
bool Globals::skip_to_clustering = false;
bool Globals::skip_to_subvectors = false;
bool Globals::skip_sorting_subvectors = false;
bool Globals::skip_cluster_merging = false;
int Globals::hamming_class_buffer = 10000;
int Globals::skip_iterative = -1;
bool Globals::unload_blob_before_merge = false;

bool Globals::likelihood_e_step = false;
bool Globals::subtract_simplex_volume = false;
bool Globals::change_n_to_random = false;

bool Globals::debug_output_clustering = false;
bool Globals::debug_output_likelihood = false;

bool Globals::read_kmers_after_clustering = false;
bool Globals::write_kmers_after_clustering = false;
bool Globals::write_each_iteration_kmers = false;
bool Globals::regular_threshold_for_correction = false;
bool Globals::discard_only_singletons = false;
bool Globals::use_true_likelihood = false;
string Globals::kmers_after_clustering = "";


void readFileIntoBlob(const string & readsFilename, hint_t & curpos, hint_t & cur_read, bool reverse_complement) {
	TIMEDLN("Reading single reads file " << readsFilename);
	ireadstream irs(readsFilename, Globals::qvoffset);
	VERIFY(irs.is_open());
	Read r;
	while (irs.is_open() && !irs.eof()) {
		irs >> r;
		size_t read_size = r.trimNsAndBadQuality(Globals::trim_quality);
		if (read_size < K) continue;
		if ( reverse_complement ) r = !r;
		PositionRead pread(curpos, read_size, cur_read, false);
		Globals::pr->push_back(pread);
		for (uint32_t j = 0; j < read_size; ++j) {
			Globals::blob[curpos + j] = r.getSequenceString()[j];
			Globals::blobquality[curpos + j] = (char) (Globals::qvoffset + r.getQualityString()[j]);
		}
		curpos += read_size;
		++cur_read;
	}
	irs.close();
}

bool correctAndUpdateOneRead( vector<ofstream *> outfv, const vector<KMerCount*> & kmers, vector<hint_t> & changedReads,
		vector<hint_t> & changedNucleotides, hint_t readno, Read & r, size_t i ) {
	bool isGood = false;
	bool res = CorrectRead(Globals::hm, kmers, readno, r, isGood, outfv[i]);
	if (!Globals::conserve_memory) Globals::rv_bad->at(readno) = isGood;
	changedNucleotides[i] += res;
	if (res) ++changedReads[i];
	if (res && outfv[i] != NULL) {
			*(outfv[i]) << "Final result again:  size=" << r.size() << "\n"
					<< r.getSequenceString().c_str() << "\n"
					<< r.getPhredQualityString(Globals::qvoffset).c_str() << endl;
	}
	return isGood;
}

void correctAndUpdateReadFile( const string & readsFilename, vector<ofstream *> outfv, const vector<KMerCount*> & kmers, vector<hint_t> & changedReads,
		vector<hint_t> & changedNucleotides, ofstream *outf_good, ofstream *outf_bad ) {
	ireadstream irs(readsFilename, Globals::qvoffset);
	VERIFY(irs.is_open());
	hint_t readno = 0;
	while (irs.is_open() && !irs.eof()) {
		Read r;
		irs >> r;
		size_t read_size = r.trimNsAndBadQuality(Globals::trim_quality);
		if (read_size < K) {
			if (outfv[0] != NULL) { *(outfv[0]) << r.getName() << "  is very bad: not a single " << K << "-mer" << endl; }
			continue;
		}
		if ( correctAndUpdateOneRead(outfv, kmers, changedReads, changedNucleotides, readno, r, 0) ) {
			if (outfv[0] != NULL) {
				*(outfv[0]) << " good, so again: " << r.getSequenceString().c_str() << "\n";
			}
			r.print(*outf_good, Globals::qvoffset);
		} else {
			r.print(*outf_bad, Globals::qvoffset);
		}
		++readno;
	}
	irs.close();
}

void correctAndUpdatePairedReadFiles( const string & readsFilenameLeft, const string & readsFilenameRight, vector<ofstream *> outfv,
		const vector<KMerCount*> & kmers, vector<hint_t> & changedReads, vector<hint_t> & changedNucleotides,
		ofstream * ofbadl, ofstream * ofcorl, ofstream * ofunpl, ofstream * ofbadr, ofstream * ofcorr, ofstream * ofunpr ) {
	ireadstream irsl(readsFilenameLeft, Globals::qvoffset);
	ireadstream irsr(readsFilenameRight, Globals::qvoffset);
	VERIFY(irsl.is_open());
	VERIFY(irsr.is_open());
	hint_t readno_left = 0;
	hint_t readno_right = Globals::lastLeftNo;
	cout << "  readno_right init=" << Globals::lastLeftNo << " size=" << Globals::pr->size() << endl;
	while (irsl.is_open() && !irsl.eof()) {
		Read l; Read r;
		irsl >> l; irsr >> r;
		if (outfv[0] != NULL)
			(*outfv[0]) << "\n " << r.getName() << "\n" << r.getSequenceString().data() << "\n";
	
		size_t read_size_left = l.trimNsAndBadQuality(Globals::trim_quality);
		size_t read_size_right = r.trimNsAndBadQuality(Globals::trim_quality);
		if (read_size_left < K && read_size_right < K) {
			if (outfv[0] != NULL) (*outfv[0]) << "very bad, size < " << K << "\n";
			continue;
		}
		bool left_res = false;
		bool right_res = false;
		if (read_size_left >= K) {
			left_res = correctAndUpdateOneRead(outfv, kmers, changedReads, changedNucleotides, readno_left, l, 0);
			++readno_left;
		}
		if (read_size_right >= K) {
			right_res = correctAndUpdateOneRead(outfv, kmers, changedReads, changedNucleotides, readno_right, r, 0);
			++readno_right;
		}
		if (  !left_res ) l.print(*ofbadl, Globals::qvoffset);
		if ( !right_res ) r.print(*ofbadr, Globals::qvoffset);
		if (  left_res && !right_res ) l.print(*ofunpl, Globals::qvoffset);
		if ( !left_res &&  right_res ) r.print(*ofunpr, Globals::qvoffset);
		if (  left_res &&  right_res ) {
			l.print(*ofcorl, Globals::qvoffset);
			r.print(*ofcorr, Globals::qvoffset);
		}
	}
	irsl.close(); irsr.close();
}

int main(int argc, char * argv[]) {

    const size_t GB = 1 << 30;
    limit_memory(300 * GB);
	
	string config_file = CONFIG_FILENAME;
	if (argc > 1) config_file = argv[1];
	getGlobalConfigParameters(config_file);

	int tau = cfg::get().tau;
	string readsFilename = cfg::get().reads;
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

	TIMEDLN("sizeof: PositionKMer=" << sizeof(PositionKMer) << " QualBitSet=" << sizeof(QualBitSet) << " KMerStat=" << sizeof(KMerStat) << " KMerCount=" << sizeof(KMerCount));
	TIMEDLN("");

	string readsFilenameLeft, readsFilenameRight;
	if (Globals::paired_reads) {
		readsFilenameLeft = cfg::get().reads_left;
		readsFilenameRight = cfg::get().reads_right;
		TIMEDLN("Starting work on " << readsFilenameLeft << " and " << readsFilenameRight << " with " << nthreads << " threads, K=" << K);
	} else {
		TIMEDLN("Starting work on " << readsFilename << " with " << nthreads << " threads, K=" << K);
	}

	if (Globals::change_n_to_random && !Globals::skip_to_clustering && !Globals::skip_to_subvectors && (Globals::skip_iterative < 0)) {
		TIMEDLN("Preprocessing: change single Ns to As with quality 2");
		if (Globals::paired_reads) {
			pid_t pIDsubstN1 = vfork();
			if (pIDsubstN1 == 0) {
				TIMEDLN("  [" << getpid() << "] Child process for substituting Ns in " << readsFilenameLeft << " starting.");
				int exitcode = system((string("sed \'n;s/\\([ACGT]\\)N\\([ACGT]\\)/\\1A\\2/g;n;n\' ") + readsFilenameLeft + " > " + getFilename(Globals::working_dir, "reads.left.input")).c_str() );
				if (exitcode != 0) {
					TIMEDLN("  [" << getpid() << "] ERROR: finished with non-zero exit code " << exitcode);
				}
				_exit(0);
			}
			pid_t pIDsubstN2 = vfork();
			if (pIDsubstN2 == 0) {
				TIMEDLN("  [" << getpid() << "] Child process for substituting Ns in " << readsFilenameRight << " starting.");
				int exitcode = system((string("sed \'n;s/\\([ACGT]\\)N\\([ACGT]\\)/\\1A\\2/g;n;n\' ") + readsFilenameRight + " > " + getFilename(Globals::working_dir, "reads.right.input")).c_str() );
				if (exitcode != 0) {
					TIMEDLN("  [" << getpid() << "] ERROR: finished with non-zero exit code " << exitcode);
				}
				_exit(0);
			}
			int childExitStatus;
			waitpid(pIDsubstN1, &childExitStatus, 0);
			waitpid(pIDsubstN2, &childExitStatus, 0);

			readsFilenameLeft = getFilename(Globals::working_dir, "reads.left.input");
			readsFilenameRight = getFilename(Globals::working_dir, "reads.right.input");
			TIMEDLN("Single Ns changed, reads written to " << readsFilenameLeft << " and " << readsFilenameRight);
		} else {
			pid_t pIDsubstN = vfork();
			if (pIDsubstN == 0) {
				TIMEDLN("  [" << getpid() << "] Child process for substituting Ns in " << readsFilename << " starting.");
				string newReads = getFilename(Globals::working_dir, "reads.input");
				int exitcode = system((string("sed \'n;s/\\([ACGT]\\)N\\([ACGT]\\)/\\1A\\2/g;n;n\' ") + readsFilename + " > " + newReads).c_str() );
				if (exitcode != 0) {
					TIMEDLN("  [" << getpid() << "] ERROR: finished with non-zero exit code " << exitcode);
				}
				_exit(0);
			}
			int childExitStatus;
			waitpid(pIDsubstN, &childExitStatus, 0);
			readsFilename = getFilename(Globals::working_dir, "reads.input").data();
			TIMEDLN("Single Ns changed, reads written to " << readsFilename);
		}
	}

	if (Globals::change_n_to_random) {
		if (Globals::paired_reads) {
			readsFilenameLeft = getFilename(Globals::working_dir, "reads.left.input");
			readsFilenameRight = getFilename(Globals::working_dir, "reads.right.input");
		} else {
			readsFilename = getFilename(Globals::working_dir, "reads.input").data();
		}
	}


	// initialize subkmer positions
	Globals::subKMerPositions = new std::vector<uint32_t>(tau + 2);
	for (uint32_t i=0; i < (uint32_t)(tau+1); ++i) {
		Globals::subKMerPositions->at(i) = (i * K / (tau+1) );
	}
	Globals::subKMerPositions->at(tau+1) = K;

	hint_t totalReadSize = 0;

	// in memory conservation mode, we don't keep Globals::rv and Globals::rv_bad at all
	if (Globals::conserve_memory) {
		// we estimate the needed memory by combining file sizes and dividing them in two (bases+quality)
		struct stat st;
		if (!Globals::paired_reads) {
			stat(readsFilename.c_str(), &st);
			totalReadSize += st.st_size;
		} else {
			stat(readsFilenameLeft.c_str(), &st);
			totalReadSize += st.st_size;
			stat(readsFilenameRight.c_str(), &st);
			totalReadSize += st.st_size;
		}
		totalReadSize = totalReadSize / (2.5);
		TIMEDLN("Estimated total size of all reads is " << totalReadSize);
	} else {
		if (!Globals::paired_reads) {
			Globals::rv = new std::vector<Read>();
			ireadstream::readAllNoValidation(Globals::rv, readsFilename,
					&totalReadSize, Globals::qvoffset, Globals::trim_quality);
			Globals::lastLeftNo = Globals::rv->size();
		} else {
			Globals::rv = new std::vector<Read>();
			ireadstream::readAllNoValidation(Globals::rv, readsFilenameLeft,
					&totalReadSize, Globals::qvoffset, Globals::trim_quality);
			Globals::lastLeftNo = Globals::rv->size();
			hint_t rightSize = 0;
			ireadstream::readAllNoValidation(Globals::rv, readsFilenameRight,
					&rightSize, Globals::qvoffset, Globals::trim_quality);
			totalReadSize += rightSize;
		}
		TIMEDLN("Total size of all reads is " << totalReadSize);
		Globals::revNo = Globals::rv->size();
		for (hint_t i = 0; i < Globals::revNo; ++i) {
			string seq = Globals::rv->at(i).getSequenceString();
			Read revcomp = !(Globals::rv->at(i));
			Globals::rv->push_back( revcomp );
		}
		Globals::rv_bad = new std::vector<bool>(Globals::rv->size(), false);
		TIMEDLN("All reads read to memory. Reverse complementary reads added.");
	}


	Globals::blob_size = totalReadSize + 1;
	Globals::blob_max_size = (hint_t)(Globals::blob_size * ( 2 + Globals::blob_margin));

	Globals::blob = new char[ Globals::blob_max_size ];
	Globals::blobquality = new char[ Globals::blob_max_size ];
	TIMEDLN("Max blob size as allocated is " << Globals::blob_max_size);

	if (readBlobAndKmers) {
		Globals::readBlob( getFilename(Globals::working_dir, blobFilename.c_str() ).c_str() );
	}

	for (int iter_count = 0; iter_count < iterno; ++iter_count) {
		cout << "\n     === ITERATION " << iter_count << " begins ===" << endl;
		Globals::iteration_no = iter_count;

		Globals::pr = new vector<PositionRead>();
		hint_t curpos = 0;
		hint_t cur_read = 0;

		// again, if we're conserving memory, we don't keep Globals::rv, so here we actually read the files
		if (Globals::conserve_memory) {
			if (!Globals::paired_reads) {
				readFileIntoBlob(readsFilename, curpos, cur_read, false);
				TIMEDLN("  readsTotal=" << cur_read << "\tblob=" << curpos);
				Globals::lastLeftNo = cur_read;
				Globals::revNo = cur_read;
				readFileIntoBlob(readsFilename, curpos, cur_read, true);
				TIMEDLN("  readsTotalWithRevComp=" << cur_read << "\tblob=" << curpos);
			} else {
				readFileIntoBlob(readsFilenameLeft,  curpos, cur_read, false);
				Globals::lastLeftNo = cur_read;
				readFileIntoBlob(readsFilenameRight, curpos, cur_read, false);
				Globals::revNo = cur_read;
				readFileIntoBlob(readsFilenameLeft,  curpos, cur_read,  true);
				readFileIntoBlob(readsFilenameRight, curpos, cur_read,  true);
			}
		} else {
			for (hint_t i = 0; i < Globals::rv->size(); ++i) {
				PositionRead pread(curpos, Globals::rv->at(i).size(), i, Globals::rv_bad->at(i));
				Globals::pr->push_back(pread);
				if (!readBlobAndKmers || iter_count > 1) {
					for (uint32_t j = 0; j < Globals::rv->at(i).size(); ++j) {
						Globals::blob[curpos + j] = Globals::rv->at(i).getSequenceString()[j];
						Globals::blobquality[curpos + j] = (char) (Globals::qvoffset + Globals::rv->at(i).getQualityString()[j]);
					}
				}
				curpos += Globals::rv->at(i).size();
				cur_read = Globals::rv->size();
			}
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

		pid_t pIDsortKmerTotalsFile = 0;

		if (!readBlobAndKmers || iter_count > 0) {
			if (Globals::conserve_memory) {
				if (!Globals::skip_to_clustering && !Globals::skip_to_subvectors && (Globals::skip_iterative < 0)) {
					TIMEDLN("Splitting kmer instances into files.");
					SplitToFiles(Globals::working_dir, iter_count);
					TIMEDLN("Kmer instances split. Starting merge in " << nthreads << " threads.");
					ofstream kmerno_file( getFilename(Globals::working_dir, iter_count, "kmers.total") );
					hint_t kmer_num = 0;

					for ( int iFile=0; iFile < Globals::num_of_tmp_files;  ) {

						std::vector<KMerNoHashMap> khashmaps(nthreads);

						#pragma omp parallel for shared(kmerno_file, kmer_num) num_threads(nthreads)
						for ( int j = 0; j< nthreads; ++j) {
							if ( j + iFile > Globals::num_of_tmp_files) continue;
							ifstream inStream( getFilename( Globals::working_dir, iter_count, "tmp.kmers", iFile+j ) );
							ProcessKmerHashFile( &inStream, khashmaps[j] );
						}

						for ( int j = 0; j< nthreads; ++j) {
							if ( j + iFile > Globals::num_of_tmp_files) continue;
							PrintProcessedKmerHashFile( &kmerno_file, kmer_num, khashmaps[j] );
						}

						iFile += nthreads;

					}

					kmerno_file.close();
					pIDsortKmerTotalsFile = vfork();
					if (pIDsortKmerTotalsFile == 0) {
						TIMEDLN("  [" << getpid() << "] Child process for sorting the kmers.total file starting.");
						execlp("sort", "sort", "-k2", "-o", getFilename(Globals::working_dir, iter_count, "kmers.total.sorted").data(),
								"-T", Globals::working_dir.c_str(), getFilename(Globals::working_dir, iter_count, "kmers.total").data(), (char *) 0 );
						_exit(0);
					}
					string cmd = "rm -rf " + getFilename(Globals::working_dir, iter_count, "tmp.kmers.*");
					if ( system(cmd.c_str()) != 0 ) { TIMEDLN("Some error with removing temporary files. Proceeding nevertheless."); }
					TIMEDLN("Merge done. There are " << kmer_num << " kmers in total.");
				} else if (Globals::skip_to_subvectors) {
					TIMEDLN("Skipping directly to subvectors, reading sorted kmers from " << getFilename(Globals::working_dir, iter_count, "kmers.total.sorted"));
				}
			} else {
				TIMEDLN("Doing honest preprocessing.");
				DoPreprocessing(tau, readsFilename, nthreads, &kmers, &Globals::hm);
				TIMEDLN("Preprocessing done. There are " << kmers.size() << " kmers in total.");
			}
		} else {
			TIMEDLN("Reading kmers from " << kmersFilename.c_str() );
			Globals::readKMerCounts( getFilename( Globals::working_dir, kmersFilename.c_str() ).c_str(), &kmers );
			TIMEDLN("Kmers read from " << kmersFilename.c_str());
		}

		if ( !Globals::read_kmers_after_clustering && writeBlobAndKmers && iter_count == 0 ) { // doesn't make sense to overwrite the first blob
			Globals::writeBlob( getFilename(Globals::working_dir, blobFilename.c_str() ).data() );
			Globals::writeKMerCounts( getFilename(Globals::working_dir, kmersFilename.c_str() ).data(), kmers );
			TIMEDLN("Blob and kmers written.");
			if ( exitAfterWritingBlobAndKMers ) break;
		}

		std::vector<hint_t> kmernos;

		if ( Globals::read_kmers_after_clustering ) {
			TIMEDLN("Reading clustering results");
			Globals::readKMerHashMap( Globals::kmers_after_clustering.c_str(), &Globals::hm, &kmers );
			TIMEDLN("Clustering results read.");
		} else {
			if (Globals::conserve_memory && Globals::skip_iterative < 0) {
				SubKMerSorter * skmsorter = NULL;

				if (!Globals::skip_to_clustering) {
					int childExitStatus;
					waitpid(pIDsortKmerTotalsFile, &childExitStatus, 0);
				} else {
					TIMEDLN("Skipping straight to clustering.");
				}

				fillInKmersFromFile( getFilename(Globals::working_dir, iter_count, "kmers.total.sorted"), &kmernos );

				if (!Globals::skip_cluster_merging) {
					TIMEDLN("KMer indices filled, starting subvector sorting.");
					skmsorter = new SubKMerSorter(kmernos.size(), &kmernos, nthreads, tau, SubKMerSorter::SorterTypeFileBasedStraight);
					skmsorter->runSort(getFilename(Globals::working_dir, iter_count, "kmers.total.sorted"));
					TIMEDLN("All subvector sorting done, starting clustering.");
				} else {
					TIMEDLN("KMer indices filled, skipping straight to merging.");
				}

				if (Globals::skip_iterative < 0) {
					KMerClustering kmc(&kmers, &kmernos, nthreads, tau);
					// prepare the maps
					ofstream ofkmersnum(getFilename(Globals::working_dir, iter_count, "kmers.num").data());
					ofkmersnum << kmers.size() << endl;
					ofkmersnum.close();
					ofstream ofkmers(getFilename(Globals::working_dir, iter_count,	"kmers.solid").data());
					ofstream ofkmers_bad(getFilename(Globals::working_dir, iter_count,	"kmers.bad").data());
					kmc.process(Globals::working_dir, skmsorter, &ofkmers, &ofkmers_bad);
					ofkmers.close();
					ofkmers_bad.close();
					TIMEDLN("Finished clustering.");
				}

			} else if (Globals::conserve_memory && Globals::skip_iterative >= 0) {
				/*TIMEDLN("Reading kmers from " << getFilename(Globals::working_dir, iter_count, "kmers.total.sorted"));
				fillInKmersAndNosFromFile( getFilename(Globals::working_dir, iter_count, "kmers.total.sorted"), &kmers, &kmernos );
				Globals::kmernos = &kmernos;
				TIMEDLN("KMer indices filled, reading solid kmers.");*/

			} else if (!Globals::conserve_memory) {
				TIMEDLN("Starting subvectors sort.");
				SubKMerSorter * skmsorter = new SubKMerSorter(kmers.size(),
						&kmers, nthreads, tau,
						SubKMerSorter::SorterTypeStraight);
				skmsorter->runSort();
				TIMEDLN("Auxiliary subvectors sorted. Starting split kmer processing in " << min(nthreads, tau+1) << " effective threads.");

				KMerClustering kmc(&kmers, nthreads, tau);
				// prepare the maps
				ofstream ofkmersnum(getFilename(Globals::working_dir, iter_count, "kmers.num").data());
				ofkmersnum << kmers.size() << endl;
				ofkmersnum.close();
				ofstream ofkmers(getFilename(Globals::working_dir, iter_count,	"kmers.solid").data());
				ofstream ofkmers_bad(getFilename(Globals::working_dir, iter_count,	"kmers.bad").data());
				kmc.process(Globals::working_dir, skmsorter, &ofkmers, &ofkmers_bad);
				TIMEDLN("Finished clustering. Closing and deleting");
				ofkmers.close();
				ofkmers_bad.close();
				TIMEDLN("Finished clustering. Closed");
				delete skmsorter;
				TIMEDLN("Finished clustering. Deleted");
				TIMEDLN("Finished clustering.");
			}
		}

		if ( Globals::write_kmers_after_clustering && iter_count == 0 ) {
			TIMEDLN("Writing k-mers hash after clustering.");
			Globals::writeKMerHashMap( getFilename(Globals::working_dir, iter_count, "kmers.hash").data(), Globals::hm);
			TIMEDLN("K-mers hash written.");
		}

		if ( Globals::use_iterative_reconstruction && (Globals::skip_iterative < 0) ) {
			if (Globals::conserve_memory) Globals::kmernos = &kmernos;
			for ( int iter_no = 0; iter_no < Globals::max_reconstruction_iterations; ++iter_no ) {
				ofstream ofs( getFilename(Globals::working_dir, iter_count, "kmers.iterative", iter_no) );
				size_t res = IterativeReconstructionStep(nthreads, kmers, &ofs);
				ofs.close();
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
				if ( res < 10 ) break;
			}
			TIMEDLN("Solid k-mers finalized.");

			// Writing final set of k-mers
			ofstream ofkmerstotal(getFilename(Globals::working_dir, iter_count, "kmers.result").data());
			PrintKMerFileWithChangeTo( &ofkmerstotal, kmers );
			ofkmerstotal.close();
			TIMEDLN("KMers with changetos written.");

		} else if (Globals::skip_iterative >= 0) {
			TIMEDLN("Loading iterative results from " << getFilename(Globals::working_dir, iter_count, "goodkmers", Globals::skip_iterative));
			// fillInKmersAnNosFromFile( getFilename(Globals::working_dir, iter_count, "kmers.bad"), &kmers, &kmernos );
			TIMEDLN("Reading kmers");
			fillInKmersWithChangeToFromFile( getFilename(Globals::working_dir, iter_count, "kmers.result"), &kmers, &kmernos );
			Globals::kmernos = &kmernos;
			TIMEDLN("Iterative results loaded.");
		}

		// Now for the reconstruction step; we still have the reads in rv, correcting them in place.
		vector<ofstream *> outfv; vector<hint_t> changedReads; vector<hint_t> changedNucleotides;
		for (int i=0; i < (Globals::conserve_memory ? 1 : nthreads); ++i) {
			outfv.push_back(new ofstream( getFilename(Globals::working_dir, iter_count, "reconstruct", i ).data() ));
			//outfv.push_back(NULL);
			changedReads.push_back(0);
			changedNucleotides.push_back(0);
		}

		// Read reconstruction and output
		if (Globals::conserve_memory) {
			TIMEDLN("Starting read reconstruction");
			if (!Globals::paired_reads) {
				ofstream ofgood(getFilename(Globals::working_dir, iter_count, "reads.corrected").c_str());
				ofstream ofbad( getFilename(Globals::working_dir, iter_count, "reads.bad").c_str());
				correctAndUpdateReadFile(readsFilename, outfv, kmers, changedReads, changedNucleotides, &ofgood, &ofbad );
				ofgood.close(); ofbad.close();
			} else {
				ofstream ofbadl(getFilename(Globals::working_dir, iter_count, "reads.left.bad").c_str());
				ofstream ofcorr(getFilename(Globals::working_dir, iter_count, "reads.right.corrected").c_str());
				ofstream ofbadr(getFilename(Globals::working_dir, iter_count, "reads.right.bad").c_str());
				ofstream ofunpl(getFilename(Globals::working_dir, iter_count, "reads.left.unpaired").c_str());
				ofstream ofunpr(getFilename(Globals::working_dir, iter_count, "reads.right.unpaired").c_str());
				ofstream ofcorl(getFilename(Globals::working_dir, iter_count, "reads.left.corrected").c_str());
				correctAndUpdatePairedReadFiles(readsFilenameLeft, readsFilenameRight, outfv, kmers, changedReads, changedNucleotides, &ofbadl, &ofcorl, &ofunpl, &ofbadr, &ofcorr, &ofunpr );
				ofbadl.close(); ofbadr.close();
				ofcorl.close(); ofcorr.close();
				ofunpl.close(); ofunpr.close();
			}
		} else {
			#pragma omp parallel for shared(changedReads, changedNucleotides, outfv) num_threads(nthreads)
			for (size_t i = 0; i < Globals::revNo; ++i) {
				correctAndUpdateOneRead(outfv, kmers, changedReads, changedNucleotides, i, Globals::rv->at(i), omp_get_thread_num());
			}
		}
		hint_t totalReads = 0; hint_t totalNucleotides = 0;
		for (int i=0; i < (Globals::conserve_memory ? 1 : nthreads); ++i) {
			if (outfv[i] != NULL) { outfv[i]->close(); delete outfv[i]; }
			totalReads += changedReads[i];
			totalNucleotides += changedNucleotides[i];
		}

		TIMEDLN("Correction done. Changed " << totalNucleotides << " bases in " << totalReads << " reads.");

		if ( !Globals::conserve_memory ) {
			if (!Globals::paired_reads) {
				outputReads(false, getFilename(Globals::working_dir, iter_count, "reads.corrected").c_str(), getFilename(
						Globals::working_dir, iter_count, "reads.bad").c_str());
			} else {
				outputReads(true, getFilename(Globals::working_dir, iter_count,
						"reads.left.corrected").c_str(),
						getFilename(Globals::working_dir, iter_count, "reads.left.bad").c_str(),
						getFilename(Globals::working_dir, iter_count, "reads.right.corrected").c_str(),
						getFilename(Globals::working_dir, iter_count, "reads.right.bad").c_str(),
						getFilename(Globals::working_dir, iter_count, "reads.left.unpaired").c_str(),
						getFilename(Globals::working_dir, iter_count, "reads.right.unpaired").c_str());
			}
		}

		// prepare the reads for next iteration
		// delete consensuses, clear kmer data, and restore correct revcomps
		for (size_t i=0; i < kmers.size(); ++i) delete kmers[i];
		kmers.clear();
		Globals::hm.clear();
		delete Globals::pr;

		if (Globals::conserve_memory) {
			if (!Globals::paired_reads) {
				readsFilename = getFilename(Globals::working_dir, iter_count, "reads.corrected");
			} else {
				readsFilenameLeft = getFilename(Globals::working_dir, iter_count, "reads.left.corrected");
				readsFilenameRight = getFilename(Globals::working_dir, iter_count, "reads.right.corrected");
			}
			// cannot skip anywhere on subsequent iterations
			Globals::skip_to_clustering = false;
		} else {
			Globals::rv->resize( Globals::revNo );
			for (hint_t i = 0; i < Globals::revNo; ++i) {
				Globals::rv->push_back( !(Globals::rv->at(i)) );
			}
			TIMEDLN("Reads restored.");
		}

		if (totalReads < 1) {
			TIMEDLN("Too few reads have changed in this iteration. Exiting.");
			break;
		}

	} // iterations

	Globals::subKMerPositions->clear();
	delete Globals::subKMerPositions;
	if (!Globals::conserve_memory) {
		Globals::rv->clear();
		delete Globals::rv;
		delete Globals::rv_bad;
	}
	delete [] Globals::blob;
	delete [] Globals::blobquality;
	return 0;
}


