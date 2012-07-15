//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/*
 * hammer_tools.cpp
 *
 *  Created on: 08.07.2011
 *      Author: snikolenko
 */

#include "standard.hpp"
#include <iostream>
#include <fstream>
#include <queue>
#include <unordered_map>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/format.hpp>

#include <time.h>
#include <sys/resource.h>
#include <iomanip>
#include <algorithm>
#ifdef USE_GLIBCXX_PARALLEL
#include <parallel/algorithm>
#endif

#include "read/ireadstream.hpp"
#include "defs.hpp"
#include "mathfunctions.hpp"
#include "valid_kmer_generator.hpp"
#include "position_kmer.hpp"
#include "globals.hpp"
#include "config_struct_hammer.hpp"
#include "subkmers.hpp"
#include "hammer_tools.hpp"
#include "mmapped_writer.hpp"

// forking
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <sys/stat.h>
#include <fcntl.h>

using namespace std;

string encode3toabyte (const string & s)  {
	string retval;
	char c = 48;
	int weight = 16;
	size_t i;
	for (i = 0; i < s.length(); i += 1) {
		if (i % 3 == 0) {
			c= 48;
			weight = 16;
		}
		c += weight * nt2num(s[i]);
		weight /= 4;
		if (i % 3 == 2) retval += c;
	}
	if (i % 3 != 0) retval += c;
	return retval;
}

void HammerTools::ChangeNtoAinReadFiles() {
	std::vector<pid_t> pids;
	for (size_t iFile=0; iFile < Globals::input_filenames.size(); ++iFile) {
		string cur_filename = HammerTools::getFilename(cfg::get().input_working_dir, Globals::input_filename_bases[iFile]);
		pid_t cur_pid = vfork();
		if (cur_pid == 0) {
			TIMEDLN("  [" << getpid() << "] Child process for substituting Ns in " << Globals::input_filenames[iFile] << " starting.");
			string cmd = string("sed \'n;s/\\([ACGT]\\)N\\([ACGT]\\)/\\1A\\2/g;n;n\' ") + Globals::input_filenames[iFile].c_str() + " > " + cur_filename.c_str();
			int exitcode = system( cmd.c_str() );
			if (exitcode != 0) {
				TIMEDLN("  [" << getpid() << "] ERROR: finished with non-zero exit code " << exitcode);
			}
			_exit(0);
		}
		pids.push_back(cur_pid);
		Globals::input_filenames[iFile] = cur_filename;
	}
	int childExitStatus;
	for (size_t iFile=0; iFile < Globals::input_filenames.size(); ++iFile) {
		waitpid(pids[iFile], &childExitStatus, 0);
	}
}

void HammerTools::DecompressIfNeeded() {
	struct stat st;
	char f_cont[2];
	vector<pid_t> pIDs(Globals::input_filenames.size(), 0);
	for (size_t iFile=0; iFile < Globals::input_filenames.size(); ++iFile) {
		if (boost::filesystem::extension(boost::filesystem::path(Globals::input_filenames[iFile])) != ".gz")
			continue;

		stat(Globals::input_filenames[iFile].c_str(), &st);
		std::ifstream ifs(Globals::input_filenames[iFile], std::ios::in | std::ios::binary );
		ifs.read(f_cont, 2);
		ifs.close();
		if ( (f_cont[0] == (char)0x1f) && (f_cont[1] == (char)0x8b) ) {
			string newFilename = HammerTools::getFilename(cfg::get().input_working_dir, Globals::input_filename_bases[iFile]);
			string oldFilename = Globals::input_filenames[iFile];
			Globals::input_filenames[iFile] = newFilename;
			Globals::input_filename_bases[iFile] = boost::filesystem::basename(boost::filesystem::path(Globals::input_filename_bases[iFile]));
			pIDs[iFile] = vfork();
			if (pIDs[iFile] == 0) {
				string systemcall = string("gunzip -c ") + oldFilename + string(" > ") + newFilename;
				TIMEDLN("  [" << getpid() << "] " << systemcall);
				if (system(systemcall.c_str())) {
					TIMEDLN("  [" << getpid() << "] System error with unzipping input files!");
				}
			_exit(0);
			}
		}
	}
	for (size_t iFile=0; iFile < Globals::input_filenames.size(); ++iFile) {
		if (pIDs[iFile] != 0) {
			int childExitStatus;
			waitpid(pIDs[iFile], &childExitStatus, 0);
		}
	}
}

hint_t HammerTools::EstimateTotalReadSize() {
	struct stat st;
	hint_t totalReadSize = 0;
	for (size_t iFile=0; iFile < Globals::input_filenames.size(); ++iFile) {
		stat(Globals::input_filenames[iFile].c_str(), &st);
		totalReadSize += st.st_size;
	}
	totalReadSize = totalReadSize / (2.5);
	return totalReadSize;
}


void HammerTools::InitializeSubKMerPositions() {
	ostringstream log_sstream;
	log_sstream.str("");
	Globals::subKMerPositions = new std::vector<uint32_t>(cfg::get().general_tau + 2);
	for (uint32_t i=0; i < (uint32_t)(cfg::get().general_tau + 1); ++i) {
		Globals::subKMerPositions->at(i) = (i * K / (cfg::get().general_tau + 1) );
		log_sstream << Globals::subKMerPositions->at(i) << " ";
	}
	Globals::subKMerPositions->at(cfg::get().general_tau + 1) = K;
	TIMEDLN("Hamming graph threshold tau=" << cfg::get().general_tau << ", k=" << K << ", subkmer positions = [ " << log_sstream.str() << "]" );
}

size_t HammerTools::ReadFileIntoBlob(const string & readsFilename, hint_t & curpos, hint_t & cur_read, bool reverse_complement) {
  TIMEDLN("Reading input file " << readsFilename);
  int trim_quality = cfg::get().input_trim_quality;
  ireadstream irs(readsFilename, cfg::get().input_qvoffset);
  Read r;
  size_t reads = 0;
  while (irs.is_open() && !irs.eof()) {
    irs >> r;
    if (reverse_complement) r = !r;
    size_t read_size = r.trimNsAndBadQuality(trim_quality);

    PositionRead pread(curpos, read_size, cur_read, false);
    if (read_size >= K) {
      pread.set_ltrim(r.ltrim());
    }
    Globals::pr->push_back(pread);

    const std::string &s = r.getSequenceString();
    memcpy(Globals::blob + curpos, s.data(), read_size);

    if (!Globals::use_common_quality) {
      const std::string &q = r.getQualityString();
      memcpy(Globals::blobquality + curpos, q.data(), read_size);
    }

    curpos += read_size;
    reads += 1;
    cur_read += 1;
  }
  irs.close();

  return reads;
}

void HammerTools::ReadAllFilesIntoBlob() {
	if (Globals::pr) Globals::pr->clear(); else Globals::pr = new vector<PositionRead>();
	hint_t curpos = 0;
	hint_t cur_read = 0;
	Globals::input_file_blob_positions.clear();
  Globals::input_file_sizes.clear();
	Globals::input_file_blob_positions.push_back(0);
	for (size_t iFile=0; iFile < Globals::input_filenames.size(); ++iFile) {
		size_t reads = ReadFileIntoBlob(Globals::input_filenames[iFile], curpos, cur_read, false);
		Globals::input_file_blob_positions.push_back(cur_read);
    Globals::input_file_sizes.push_back(reads);
	}
	Globals::revNo = cur_read;
	for (size_t iFile=0; iFile < Globals::input_filenames.size(); ++iFile) {
		ReadFileIntoBlob(Globals::input_filenames[iFile], curpos, cur_read, true);
	}
  TIMEDLN("All files were read. Used " << curpos << " bytes out of " << Globals::blob_max_size << " allocated.");
}

void HammerTools::findMinimizers( vector< pair<hint_t, pair< double, size_t > > > & v, int num_minimizers, vector< hint_t > & mmers, int which_first ) {
	if (which_first == 0) {
		sort(v.begin(), v.end(), PositionKMer::compareSubKMersGreaterSimple);
	} else if (which_first == 1) {
		sort(v.begin(), v.end(), PositionKMer::compareSubKMersLessSimple);
	} else if (which_first == 2) {
		sort(v.begin(), v.end(), PositionKMer::compareSubKMersGFirst);
	} else {
		sort(v.begin(), v.end(), PositionKMer::compareSubKMersCFirst);
	}
	vector< pair<hint_t, pair< double, size_t > > >::iterator it = v.begin();
	vector< int > kmers_in_mmer(mmers.size(), 0); // count kmers in mmers

	for ( ; it != v.end(); ) {
/*		char c = Globals::blob[ it->first ];
		if ( (which_first == 0) && ( (c == 'A') || (c == 'C') ) ) break;
		if ( (which_first == 1) && ( (c == 'G') || (c == 'T') ) ) break;
		if ( (which_first == 2) && ( (c == 'C') || (c == 'T') ) ) break;
		if ( (which_first == 3) && ( (c == 'A') || (c == 'G') ) ) break;*/
		bool erase = true;
		for (size_t j=0; j < mmers.size(); ++j ) {
			if ( (it->first >= mmers[j]) && (it->first <= mmers[j] + M - K) && (kmers_in_mmer[j] < num_minimizers) ) {
				erase = false;
				++kmers_in_mmer[j];
			}
		}
		if ( erase ) it = v.erase(it); else ++it;
	}
}

size_t my_hash(hint_t pos) {
	size_t hash = 877;
	for (size_t i = 0; i < K; i++) {
		hash = ((hash << 5) - hash) + ((int)Globals::blob[pos + i]) * 13;
	}
	return hash;
}

void HammerTools::SplitKMers() {
	int numfiles = cfg::get().count_numfiles;
	int count_num_threads = min(cfg::get().count_merge_nthreads, cfg::get().general_max_nthreads);
	int num_minimizers = 0;
	if (cfg::get().general_num_minimizers) {
		num_minimizers = *cfg::get().general_num_minimizers;
	}
	bool use_minimizers = HammerTools::doingMinimizers();
	int which_first = Globals::iteration_no % 4;

	TIMEDLN("Splitting kmer instances into files in " << count_num_threads << " threads. This takes a while");

  MMappedWriter* ostreams = new MMappedWriter[numfiles];
  for (unsigned i = 0; i < (unsigned int)numfiles; ++i) {
    std::string filename = getFilename(cfg::get().input_working_dir, Globals::iteration_no, "tmp.kmers", i);
    ostreams[i].open(filename);
	}
	Seq<K>::hash hash_function;
	size_t readbuffer = cfg::get().count_split_buffer;

	vector< vector< vector< KMerNo > > > tmp_entries(count_num_threads);
	for (int i=0; i < count_num_threads; ++i) {
		tmp_entries[i].resize(numfiles);
		for (int j=0; j < numfiles; ++j) {
			tmp_entries[i][j].reserve((int)( 1.25 * readbuffer / count_num_threads));
		}
	}

	size_t cur_i = 0, cur_limit = 0;
	int cur_fileindex = 0;

	while (cur_i < Globals::pr->size()) {
		cur_limit = min(cur_limit + readbuffer, Globals::pr->size());

		#pragma omp parallel for shared(tmp_entries) num_threads(count_num_threads)
		for (size_t i = cur_i; i < cur_limit; ++i) {
      const PositionRead &pr = Globals::pr->at(i);
      // Skip opaque reads
      if (!pr.valid())
        continue;

      size_t cpos = pr.start(), csize = pr.size();
      string s(Globals::blob + cpos, csize);
			string q;
			if (Globals::use_common_quality) {
				q = string(csize, (char)Globals::common_quality);
			} else {
				q = string(Globals::blobquality + cpos, csize);
			}
			ValidKMerGenerator<K> gen(s, q);
			if (use_minimizers ) {

				// M is the larger k-mer size -- every m-mer should have several minimizer k-mers inside
				ValidKMerGenerator<M> gen_m(s, q);
				vector< hint_t > mmers;
				//cout << s << endl;
				while (gen_m.HasMore()) {
					//for (size_t j=0; j<gen_m.pos(); ++j) cout << " ";
					//cout << string(Globals::blob + cpos + gen_m.pos() - 1, M) << endl;
					mmers.push_back( cpos + gen_m.pos() - 1 );
					gen_m.Next();
				}

				vector< pair<hint_t, pair< double, size_t > > > kmers;
				while (gen.HasMore()) {
					kmers.push_back( make_pair(cpos + gen.pos() - 1,
                                     make_pair(1 - gen.correct_probability(), my_hash(cpos + gen.pos() - 1) ) ));
					gen.Next();
				}
				HammerTools::findMinimizers( kmers, num_minimizers, mmers, which_first );

				/*for (size_t i=0; i<kmers.size(); ++i) {
					for (size_t j=0; j<kmers[i].first - cpos; ++j) cout << " ";
					cout << string(Globals::blob + kmers[i].first, K) << endl;
				}
				cout << endl;*/
				for (vector< pair<hint_t, pair< double, size_t > > >::const_iterator it = kmers.begin(); it != kmers.end(); ++it ) {
					tmp_entries[omp_get_thread_num()][it->second.second % numfiles].push_back(KMerNo(it->first, it->second.first));
				}
			} else {
				//cout << s << endl;
				while (gen.HasMore()) {
					//cout << gen.kmer().str() << endl;
					tmp_entries[omp_get_thread_num()][hash_function(gen.kmer()) % numfiles].push_back(KMerNo(cpos + gen.pos() - 1,
                                                                                                   1 - gen.correct_probability()));
					gen.Next();
				}
			}
		}
		cur_i = cur_limit;

		++cur_fileindex;

		#pragma omp parallel for shared(tmp_entries) num_threads(count_num_threads)
		for (int k=0; k < numfiles; ++k) {
      size_t sz = 0;
      for (size_t i = 0; i < (size_t)count_num_threads; ++i)
        sz += tmp_entries[i][k].size() * sizeof(tmp_entries[i][k][0]);

      ostreams[k].reserve(sz);
			for (size_t i = 0; i < (size_t)count_num_threads; ++i) {
        ostreams[k].write(&tmp_entries[i][k][0], tmp_entries[i][k].size() * sizeof(tmp_entries[i][k][0]));
			}
		}

		for (int i=0; i < count_num_threads; ++i) {
			tmp_entries[i].clear();
			tmp_entries[i].resize(numfiles);
			for (int j=0; j < numfiles; ++j) {
				tmp_entries[i][j].reserve((int)( 1.25 * readbuffer / count_num_threads));
			}
		}
	}
  delete[] ostreams;
}

void HammerTools::FillMapWithMinimizers( KMerMap & m ) {
	int num_minimizers = 0;
	if (cfg::get().general_num_minimizers) {
		num_minimizers = *cfg::get().general_num_minimizers;
	}
	int which_first = Globals::iteration_no % 4;
	for (int i = 0; i < (int)Globals::pr->size(); ++i ) {
    const PositionRead &pr = Globals::pr->at(i);
    // Skip opaque reads
    if (!pr.valid())
      continue;

    string s(Globals::blob        + pr.start() , pr.size());
		string q;

		if (Globals::use_common_quality) {
			q = string(pr.size(), (char)Globals::common_quality);
		} else {
			q = string(Globals::blobquality + pr.start(), pr.size());
		}
		ValidKMerGenerator<K> gen(s, q);
		vector< pair<hint_t, pair< double, size_t > > > kmers;
		unordered_map<hint_t, Seq<K> > seqs;
		while (gen.HasMore()) {
			hint_t cur_pos = pr.start() + gen.pos() - 1;
			kmers.push_back( make_pair(cur_pos, make_pair( 1 - gen.correct_probability(), my_hash(cur_pos) ) ));
			seqs[ cur_pos ] = gen.kmer();
			gen.Next();
		}
		ValidKMerGenerator<M> gen_m(s, q);
		vector< hint_t > mmers;
		while (gen_m.HasMore()) {
			mmers.push_back( Globals::pr->at(i).start() + gen_m.pos() - 1 );
			gen_m.Next();
		}

		HammerTools::findMinimizers( kmers, num_minimizers, mmers, which_first );
		for ( vector< pair<hint_t, pair< double, size_t > > >::const_iterator it = kmers.begin(); it != kmers.end(); ++it ) {
			Seq<K> km = seqs[it->first];
			KMerMap::iterator mit = m.find(km);
			if (mit == m.end()) {
				m[km] = make_pair(PositionKMer(it->first), KMerStat(Globals::use_common_quality, 1, KMERSTAT_GOODITER, it->second.first));
			} else {
				if (mit->second.second.count == 1) {
					QualBitSet qbs(K);
					for (uint32_t j=0; j<K; ++j) {
						qbs.set(j, Globals::blobquality[mit->second.first.start() + j]);
					}
					mit->second.second.qual = qbs;
				}
				mit->second.second.count++;
				mit->second.second.totalQual *= it->second.first;
				for (uint32_t j=0; j<K; ++j) {
					mit->second.second.qual.set(j,
							min( MAX_SHORT, mit->second.second.qual[j] + Globals::blobquality[it->first + j]) );
				}
			}
		}
	}
}

void HammerTools::CountKMersBySplitAndMerge() {
	vector<KMerCount> vec;
	KMerMap m;

	// split
	/*if ( HammerTools::doingMinimizers() ) {
		TIMEDLN("Filling map");
		FillMapWithMinimizers( m );
		Globals::number_of_kmers = m.size();
		TIMEDLN("Map filled, " << Globals::number_of_kmers << " kmers in total");
	} else {*/
    if (cfg::get().count_do) {
      HammerTools::SplitKMers();
    }

    int count_num_threads = min( cfg::get().count_merge_nthreads, cfg::get().general_max_nthreads );
    int numfiles = cfg::get().count_numfiles;

    TIMEDLN("Kmer instances split. Starting merge in " << count_num_threads << " threads.");
//   int merge_nthreads = min( cfg::get().general_max_nthreads, cfg::get().count_merge_nthreads);

    std::vector< std::vector<KMerCount> > kmcvec(numfiles);

    for (int iFile=0; iFile < numfiles; ++iFile) {
      std::string fname = getFilename(cfg::get().input_working_dir, Globals::iteration_no, "tmp.kmers", iFile);
      ProcessKmerHashFile(fname, kmcvec[iFile]);
      TIMEDLN("Processed file " << iFile);
    }

    TIMEDLN("Concat vectors");
    size_t totalsize = 0; for ( int iFile=0; iFile < numfiles; ++iFile ) totalsize += kmcvec[iFile].size();
    vec.reserve(totalsize);
    vector< size_t > enditers;
    for (int iFile=0; iFile < numfiles; ++iFile ) {
      vec.insert(vec.end(), kmcvec[iFile].begin(), kmcvec[iFile].end());
      enditers.push_back(vec.size());
      kmcvec[iFile].clear();
    }

    TIMEDLN("Merge in place");
    KMerNo::is_less_kmercount fun_ilk;
    for (int iFile=1; iFile < numfiles; ++iFile ) {
      std::inplace_merge(vec.begin(), vec.begin() + enditers[iFile-1], vec.begin() + enditers[iFile], fun_ilk );
    }
	//}

	TIMEDLN("Extracting kmernos");
	Globals::kmernos->clear();
	Globals::kmernos->reserve(vec.size());
	/*if (HammerTools::doingMinimizers()) {
    for (KMerMap::const_iterator mit = m.begin(); mit != m.end(); ++mit ) {
      Globals::kmernos->push_back(mit->second.first.start());
    }
	} else {*/
    for (size_t i=0; i < vec.size(); ++i) {
      Globals::kmernos->push_back(vec[i].first.start());
    }
	//}

	TIMEDLN("Writing subvectors");
	int tau_ = cfg::get().general_tau;
	vector< boost::shared_ptr<FOStream> > ostreams(tau_+1);
	ostreams[0] =
      FOStream::init_buf(getFilename(cfg::get().input_working_dir, Globals::iteration_no, "subkmers.sorted", 0),
                         1 << cfg::get().general_file_buffer_exp);
	for (int j = 1; j < tau_+1; ++j) {
		ostreams[j] =
        FOStream::init_buf(getFilename(cfg::get().input_working_dir, Globals::iteration_no, "subkmers", j),
                           1 << cfg::get().general_file_buffer_exp);
	}
	/*if (HammerTools::doingMinimizers()) {
		size_t num_kmer = 0;
		for ( KMerMap::const_iterator mit = m.begin(); mit != m.end(); ++mit ) {
			const hint_t & pos = mit->second.first.start();
			for (int j=0; j < tau_+1; ++j) {
				(ostreams[j]->fs) << string(Globals::blob + pos + Globals::subKMerPositions->at(j), Globals::subKMerPositions->at(j+1) - Globals::subKMerPositions->at(j))
						<< "\t" << num_kmer << "\n";
			}
			++num_kmer;
		}
	} else {*/
    for (size_t i=0; i < vec.size(); ++i) {
      const hint_t & pos = vec[i].first.start();
      for (int j=0; j < tau_+1; ++j) {
        (ostreams[j]->fs) << string(Globals::blob + pos + Globals::subKMerPositions->at(j), Globals::subKMerPositions->at(j+1) - Globals::subKMerPositions->at(j))
                          << "\t" << i << "\n";
      }
    }
	//}
	ostreams.clear();

	int sort_nthreads = cfg::get().general_max_nthreads;
	TIMEDLN("Starting child processes for sorting subvector files.");
	vector< pid_t > pids(tau_+1);
	for (int j=1; j < tau_+1; ++j) {
		pids[j] = vfork();
		if ( pids[j] == 0 ) {
			TIMEDLN("  [" << getpid() << "] Child process " << j << " for sorting subkmers starting.");
			ostringstream arg4; arg4 << "--parallel=" << sort_nthreads;
			string arg2 = string("-T") + cfg::get().input_working_dir;
			string arg1 = string("-o" + HammerTools::getFilename(cfg::get().input_working_dir, Globals::iteration_no, "subkmers.sorted", j));
			string arg3 = HammerTools::getFilename(cfg::get().input_working_dir, Globals::iteration_no, "subkmers", j);
			char *const arg[5] = {const_cast<char*>("-k1"),
					const_cast<char*>(arg1.c_str()),
					const_cast<char*>(arg2.c_str()),
					const_cast<char*>(arg3.c_str()),
					//const_cast<char*>(arg4.str().c_str()),
					NULL };
			execvp("sort", arg);
		} else if ( pids[j] < 0 ) {
			TIMEDLN("Failed to fork. Exiting.");
			exit(1);
		}
	}

  TIMEDLN("Serializing sorted kmers.");
	/*if (HammerTools::doingMinimizers()) {
    ofstream os(HammerTools::getFilename(cfg::get().input_working_dir, Globals::iteration_no, "kmers.total.ser").c_str(), ios::binary);
    for (KMerMap::const_iterator mit = m.begin(); mit != m.end(); ++mit)
      os.write((char*)&mit->second, sizeof(mit->second));
    m.clear();
    os.close();
	} else {*/
    ofstream os(HammerTools::getFilename(cfg::get().input_working_dir, Globals::iteration_no, "kmers.total.ser").c_str(), ios::binary);
    size_t sz = vec.size();
    os.write((char*)&sz, sizeof(sz));
    for (auto I = vec.begin(), E = vec.end(); I != E; ++I)
      binary_write(os, *I);
    os.close();
	//}
	if (!cfg::get().general_remove_temp_files) {
		TIMEDLN("Serializing kmernos.");
		ofstream os(HammerTools::getFilename(cfg::get().input_working_dir.c_str(), Globals::iteration_no, "kmers.numbers.ser").c_str(), ios::binary);
		size_t sz = Globals::kmernos->size();
		os.write((char*)&sz, sizeof(sz));
		os.write((char*)&(*Globals::kmernos)[0], sz*sizeof((*Globals::kmernos)[0]));
		os.close();
	}

	TIMEDLN("Waiting for subvectors to sort.");
	for (int j = 1; j < tau_+1; ++j) {
	    int status;
	    while (-1 == waitpid(pids[j], &status, 0));
	    if (!WIFEXITED(status) || WEXITSTATUS(status) != 0) {
	        TIMEDLN("Process " << j << " (pid " << pids[j] << ") failed. Exiting.");
	        exit(1);
	    }
	    HammerTools::RemoveFile(HammerTools::getFilename(cfg::get().input_working_dir, Globals::iteration_no, "subkmers", j));
	}
	TIMEDLN("Merge done. There are " << vec.size() << " kmers in total.");
}

static void fillQVec(vector<int> & qvec, hint_t index) {
	if (!Globals::use_common_quality) {
		for (uint32_t j=0; j<K; ++j) {
			qvec[j] = Globals::blobquality[index + j];
		}
	}
}

void HammerTools::KmerHashUnique(const std::vector<KMerNo> & vec, std::vector<KMerCount> & vkmc) {
	KMerCount kmc(PositionKMer(vec[0].getIndex()), KMerStat(true, 1, KMERSTAT_GOODITER, vec[0].getQual()));

	vector<int> qvec(K);
	fillQVec(qvec, vec[0].getIndex());

	bool first_occ = true;
	vkmc.push_back(kmc);
	for (size_t i=1; i < vec.size(); ++i) {
		KMerCount & curkmc = vkmc[vkmc.size() - 1];
		if (vec[i].equal(curkmc)) {
			first_occ = false;
			curkmc.second.count++;
			curkmc.second.totalQual *= vec[i].getQual();
			fillQVec(qvec, vec[i].getIndex());
		} else {
			if (!first_occ && !Globals::use_common_quality) {
				curkmc.second.qual = QualBitSet();
				for (uint32_t j=0; j<K; ++j) {
					curkmc.second.qual.set(j, min(MAX_SHORT, qvec[j]));
				}
			}
			KMerCount kmc(PositionKMer(vec[i].getIndex()), KMerStat(true, 1, KMERSTAT_GOODITER, vec[i].getQual()));
			fillQVec(qvec, vec[i].getIndex());
			vkmc.push_back(kmc);
			first_occ = true;
		}
	}
	if (!first_occ && !Globals::use_common_quality) {
		vkmc[vkmc.size() - 1].second.qual = QualBitSet();
		for (uint32_t j=0; j<K; ++j) {
			vkmc[vkmc.size() - 1].second.qual.set(j, min(MAX_SHORT, qvec[j]));
		}
	}
}

void HammerTools::ProcessKmerHashFile(const std::string &fname, std::vector<KMerCount> & vkmc) {
  std::vector<KMerNo> vec;

  // Make sure memory mapping is released as soon as possible
  {
    MMappedRecordReader<KMerNo> ins(fname, /* unlink */ true);;

    vec.resize(ins.size());
    ins.read(&vec[0], vec.size());
    VERIFY(!ins.good());
  }

#ifdef USE_GLIBCXX_PARALLEL
  // Explicitly force a call to parallel sort routine.
  __gnu_parallel::sort(vec.begin(), vec.end(), KMerNo::is_less());
#else
  std::sort(vec.begin(), vec.end(), KMerNo::is_less());
#endif
  TIMEDLN("Sorting done, starting unification.");
  if (!vec.size()) return;

  // FIXME: This needs to be parallel as well.
  HammerTools::KmerHashUnique(vec, vkmc);
}

void HammerTools::ProcessKmerHashVector( const vector< pair<hint_t, double> > & sv, KMerNoHashMap & km, std::vector<KMerCount> & vkmc ) {
	std::vector<KMerNo> vec;
	for (size_t i=0; i < sv.size(); ++i) {
		vec.push_back(KMerNo(sv[i].first, sv[i].second));
	}
	KMerNo::is_less kmerno_cmp;
	std::sort(vec.begin(), vec.end(), kmerno_cmp );
	if (!vec.size()) return;

	HammerTools::KmerHashUnique(vec, vkmc);
}

void HammerTools::PrintProcessedKmerHashFile( boost::iostreams::filtering_ostream & outf, hint_t & kmer_num, KMerNoHashMap & km ) {
	for (KMerNoHashMap::iterator it = km.begin(); it != km.end(); ++it) {
		outf << it->second->first.start() << "\t"
				<< string(Globals::blob + it->second->first.start(), K) << "\t"
				<< it->second->second.count << "\t"
				<< setw(8) << it->second->second.totalQual << "\t";
		for (size_t i=0; i < K; ++i) outf << it->second->second.qual[i] << " ";
		outf << "\n";
		delete it->second;
		++kmer_num;
	}
	km.clear();
}

void HammerTools::PrintProcessedKmerHashFile( boost::iostreams::filtering_ostream & outf, hint_t & kmer_num, std::vector<KMerCount> & km ) {
	for (std::vector<KMerCount>::const_iterator it = km.begin(); it != km.end(); ++it) {
		outf << it->first.start() << "\t"
				<< string(Globals::blob + it->first.start(), K) << "\t"
				<< it->second.count << "\t"
				<< setw(8) << it->second.totalQual << "\t";
		for (size_t i=0; i < K; ++i) outf << it->second.qual[i] << " ";
		outf << "\n";
		++kmer_num;
	}
}

void HammerTools::ReadKmersWithChangeToFromFile( const string & fname, vector<KMerCount> *kmers, vector<hint_t> *kmernos ) {
	kmernos->clear();
	kmers->clear();
	boost::shared_ptr<FIStream> fis = FIStream::init_buf(fname, 1 << cfg::get().general_file_buffer_exp);
	string buf;
	while (fis->fs.good()) {
		std::getline(fis->fs, buf);
		hint_t pos; int cnt; double qual; hint_t chg;
		sscanf(buf.c_str(), "%lu\t%*s\t%u\t%lu\t%lf", &pos, &cnt, &chg, &qual);
		kmernos->push_back(pos);
		kmers->push_back(KMerCount(PositionKMer(pos), KMerStat(Globals::use_common_quality, cnt, chg, qual)));
	}
}

string HammerTools::getFilename( const string & dirprefix, const string & suffix ) {
	ostringstream tmp;
	tmp.str(""); tmp << dirprefix.data() << "/" << suffix.data();
	return tmp.str();
}

string HammerTools::getFilename( const string & dirprefix, int iter_count, const string & suffix ) {
	ostringstream tmp;
	tmp.str(""); tmp << dirprefix.data() << "/" << std::setfill('0') << std::setw(2) << iter_count << "." << suffix.data();
	return tmp.str();
}

string HammerTools::getReadsFilename( const string & dirprefix, int read_file_no, int iter_no, const string & suffix ) {
	ostringstream tmp;
	tmp.str("");
	tmp << dirprefix.data() << "/" << Globals::input_filename_bases[read_file_no] << '.' << std::setfill('0') << std::setw(2) << iter_no << "." << suffix.data() << ".fastq";
	return tmp.str();
}

string HammerTools::getFilename( const string & dirprefix, const string & suffix, int suffix_num ) {
	ostringstream tmp;
	tmp.str(""); tmp << dirprefix.data() << "/" << suffix.data() << "." << suffix_num;
	return tmp.str();
}

string HammerTools::getFilename( const string & dirprefix, int iter_count, const string & suffix, int suffix_num ) {
	ostringstream tmp;
	tmp.str(""); tmp << dirprefix.data() << "/" << std::setfill('0') << std::setw(2) << iter_count << "." << suffix.data() << "." << suffix_num;
	return tmp.str();
}

string HammerTools::getFilename( const string & dirprefix, int iter_count, const string & suffix, int suffix_num, const string & suffix2 ) {
	ostringstream tmp;
	tmp.str(""); tmp << dirprefix.data() << "/" << std::setfill('0') << std::setw(2) << iter_count << "." << suffix.data() << "." << suffix_num << "." << suffix2.data();
	return tmp.str();
}

bool HammerTools::internalCorrectReadProcedure( const Read & r, const hint_t readno, const string & seq, const vector<KMerCount> & km,
		const PositionKMer & kmer, const uint32_t pos, const KMerStat & stat, vector< vector<int> > & v,
		int & left, int & right, bool & isGood, ofstream * ofs, bool revcomp, bool correct_threshold, bool discard_singletons ) {
	bool res = false;
	if (  stat.isGoodForIterative() || ( correct_threshold && stat.isGood() ) ) {
		isGood = true;
		if (ofs != NULL) *ofs << "\t\t\tsolid";
		//cout << "\t\t\tsolid";
		for (size_t j = 0; j < K; ++j) {
			if (!revcomp)
				v[dignucl(kmer[j])][pos + j]++;
			else
				v[complement(dignucl(kmer[j]))][K-1-pos-j]++;
		}
		if ((int) pos < left)
			left = pos;
		if ((int) pos > right)
			right = pos;
	} else {
		// if discard_only_singletons = true, we always use centers of clusters that do not coincide with the current center
		if (stat.change() && (discard_singletons
				|| km[stat.changeto].second.isGoodForIterative()
				|| ( correct_threshold && stat.isGood() ))) {
			if (ofs != NULL) *ofs << "\tchange to\n";
			//cout << "  kmer " << kmer.start() << " " << kmer.str() << " wants to change to " << stat.changeto << " " << km[stat.changeto].first.str() << endl;
			isGood = true;
			if ((int) pos < left)
				left = pos;
			if ((int) pos > right)
				right = pos;
			const PositionKMer & newkmer = km[stat.changeto].first;

			for (size_t j = 0; j < K; ++j) {
				v[dignucl(newkmer[j])][pos + j]++;
			}
			// pretty print the k-mer
			res = true;
			if (ofs != NULL) {
				for (size_t j = 0; j < pos; ++j)
					*ofs << " ";
				*ofs << newkmer.str().data();
			}
			//for (size_t j = 0; j < pos; ++j) cout << " "; cout << newkmer.str().data();
		}
	}
	return res;
}

hint_t HammerTools::IterativeExpansionStep(int expand_iter_no, int nthreads, vector<KMerCount> & kmers) {
	hint_t res = 0;

	// cycle over the reads, looking for reads completely covered by solid k-mers
	// and adding new solid k-mers on the fly
	#pragma omp parallel for shared(res) num_threads(nthreads)
	for (hint_t readno = 0; readno < Globals::revNo; ++readno) {
		PositionRead &pr = Globals::pr->at(readno);

    // skip opaque reads w/o kmers
    if (!pr.valid()) continue;
    // maybe this read has already been covered by solid k-mers
		if (pr.isDone()) continue;

		const uint32_t read_size = pr.size();
		vector<unsigned> covered_by_solid(read_size, false);
		vector<hint_t> kmer_indices(read_size, -1);

		pair<int, hint_t> it = make_pair( -1, BLOBKMER_UNDEFINED );
		while ( (it = pr.nextKMerNo(it.first)).first > -1 ) {
			kmer_indices[it.first] = it.second;
			if ( kmers[it.second].second.isGoodForIterative() ) {
				for ( size_t j = it.first; j < it.first + K; ++j )
					covered_by_solid[j] = true;
			}
		}

		bool isGood = true;
		for (size_t j = 0; j < read_size; ++j) {
			if (!covered_by_solid[j] ) { isGood = false; break; }
		}
		if (!isGood) continue;

		// ok, now we're sure that everything is covered
		// first, set this read as already done
		pr.set_done();

		// second, mark all k-mers as solid
		for (size_t j = 0; j < read_size; ++j) {
			if ( kmer_indices[j] == (hint_t)-1 ) continue;
			if ( !kmers[kmer_indices[j]].second.isGoodForIterative() &&
           !kmers[kmer_indices[j]].second.isMarkedGoodForIterative() ) {
				#pragma omp critical
				{
				++res;
				kmers[kmer_indices[j]].second.makeGoodForIterative();
				}
			}
		}
	}

	if ( cfg::get().expand_write_each_iteration ) {
		ofstream oftmp( getFilename(cfg::get().input_working_dir, Globals::iteration_no, "goodkmers", expand_iter_no ).data() );
		for ( hint_t n = 0; n < kmers.size(); ++n ) {
			if ( kmers[n].second.isGoodForIterative() ) {
				oftmp << kmers[n].first.str() << "\n>" << kmers[n].first.start()
				      << "  cnt=" << kmers[n].second.count << "  tql=" << (1-kmers[n].second.totalQual) << "\n";
			}
		}
	}

	return res;
}

void HammerTools::PrintKMerResult( boost::iostreams::filtering_ostream & outf, const vector<KMerCount> & kmers ) {
	for (vector<KMerCount>::const_iterator it = kmers.begin(); it != kmers.end(); ++it) {
		outf << it->first.start() << "\t"
			 << string(Globals::blob + it->first.start(), K) << "\t"
			 << it->second.count << "\t"
			 << it->second.changeto << "\t"
			 << setw(8) << it->second.totalQual << "\t";
		for (size_t i=0; i < K; ++i) outf << it->second.qual[i] << " ";
		outf << "\n";
	}
}

bool HammerTools::CorrectOneRead(const vector<KMerCount> & kmers,
                                 size_t & changedReads, size_t & changedNucleotides,
                                 hint_t readno, Read & r, size_t i, bool correct_threshold, bool discard_singletons, bool discard_bad ) {
	bool isGood = false;
	PositionRead & pr = Globals::pr->at(readno);

  string seq (Globals::blob + pr.start(), pr.size());
	const uint32_t read_size = pr.size();
  VERIFY(read_size >= K);

	// create auxiliary structures for consensus
	vector<int> vA(read_size, 0), vC(read_size, 0), vG(read_size, 0), vT(read_size, 0);
	vector< vector<int> > v;  // A=0, C=1, G=2, T=3
	v.push_back(vA); v.push_back(vC); v.push_back(vG); v.push_back(vT);
	isGood = false;

	// getting the leftmost and rightmost positions of a solid kmer
	int left = read_size; int right = -1;
	bool changedRead = false;
	pair<int, hint_t> it = make_pair( -1, BLOBKMER_UNDEFINED );
	while ( (it = pr.nextKMerNo(it.first)).first > -1 ) {
		const PositionKMer & kmer = kmers[it.second].first;
		const uint32_t pos = it.first;
		const KMerStat & stat = kmers[it.second].second;
		changedRead = changedRead || internalCorrectReadProcedure( r, readno, seq, kmers, kmer, pos, stat, v,
				left, right, isGood, NULL, false, correct_threshold, discard_singletons );
	}

	int left_rev = 0; int right_rev = read_size-(int)K;

	if ( left <= right && left_rev <= right_rev ) {
		left = std::min(left, (int)read_size - left_rev - (int)K);
		right = std::max(right, (int)read_size - right_rev - (int)K);
	} else if ( left > right && left_rev <= right_rev ) {
		left = (int)read_size - left_rev - (int)K;
		right = (int)read_size - right_rev - (int)K;
	}

	// at this point the array v contains votes for consensus
	/*cout << "\n" << seq << "\n";
	for (size_t k=0; k<4; ++k) {
		for (size_t j=0; j<read_size; ++j) {
			cout << v[k][j];
		}
		cout << "\n";
	}*/

	size_t res = 0; // how many nucleotides have really changed?
	// find max consensus element
	for (size_t j=0; j<read_size; ++j) {
		char cmax = seq[j]; int nummax = 0;
		for (size_t k=0; k<4; ++k) {
			if (v[k][j] > nummax) {
				cmax = nucl(k); nummax = v[k][j];
			}
		}
		if (seq[j] != cmax) ++res;
		seq[j] = cmax;
	}

	r.setSequence(seq.data(), /* preserve_trimming */ true);
	//cout << seq << "\n";

	// if discard_bad=false, we retain original sequences when needed
	if (discard_bad) {
		r.trimLeftRight(left, right+K-1);
		if ( left > 0 || right + K -1 < read_size ) changedRead = true;
	}

	changedNucleotides += res;
	if (res > 0) ++changedReads;
	return isGood;
}
  
bool HammerTools::doingMinimizers() {
  return ( cfg::get().general_minimizers && (Globals::iteration_no < 8) );
}
 
 
void HammerTools::CorrectReadsBatch(std::vector<bool> &res, std::vector<Read> &reads,
                                    size_t &changedReads, size_t &changedNucleotides,
                                    const std::vector<hint_t> &readno, const std::vector<size_t> &read_sizes,
                                    const vector<KMerCount> & kmers) {
  unsigned correct_nthreads = min(cfg::get().correct_nthreads, cfg::get().general_max_nthreads);
  bool discard_singletons = cfg::get().bayes_discard_only_singletons;
  bool correct_threshold = cfg::get().correct_use_threshold;
  bool discard_bad = cfg::get().correct_discard_bad && !cfg::get().correct_notrim;
  if (HammerTools::doingMinimizers()) discard_bad = false;
 
  vector<size_t> changedReadBuf(correct_nthreads, 0);
  vector<size_t> changedNuclBuf(correct_nthreads, 0);
 
  size_t buf_size = reads.size();
 # pragma omp parallel for shared(reads, res, readno, kmers) num_threads(correct_nthreads)
  for (size_t i = 0; i < buf_size; ++i) {
    if (read_sizes[i] >= K) {
      res[i] =
          HammerTools::CorrectOneRead(kmers,
                                      changedReadBuf[omp_get_thread_num()], changedNuclBuf[omp_get_thread_num()],
                                      readno[i], reads[i], 0, correct_threshold, discard_singletons, discard_bad);
    } else
      res[i] = false;
  }
 
  for (int i = 0; i < correct_nthreads; ++i ) {
    changedReads += changedReadBuf[i];
    changedNucleotides += changedNuclBuf[i];
  }
}

static void ConstructRead(Read &r, const PositionRead &pr) {
  int qvoffset = cfg::get().input_qvoffset;
  size_t cpos = pr.start(), csize = pr.size();
  string s(Globals::blob + cpos, csize);
  string q;
  if (Globals::use_common_quality) {
    q = string(csize, (char)Globals::common_quality);
  } else {
    q = string(Globals::blobquality + cpos, csize);
  }
  
  r.setSequence(s.c_str());
  r.setQuality(q.c_str(), 0);
  r.set_ltrim(pr.ltrim());
}
 
void HammerTools::CorrectReadFile(const vector<KMerCount> & kmers,
                                  size_t & changedReads, size_t & changedNucleotides,
                                  hint_t readno_start, size_t len,
                                  ofstream *outf_good, ofstream *outf_bad ) {
  int qvoffset = cfg::get().input_qvoffset;
 
  unsigned correct_nthreads = min(cfg::get().correct_nthreads, cfg::get().general_max_nthreads);
  size_t read_buffer_size = correct_nthreads * cfg::get().correct_readbuffer;
  vector<Read> reads(read_buffer_size);
  vector<bool> res(read_buffer_size, false);
  vector<size_t> read_size(read_buffer_size, 0);
  vector<hint_t> readno(read_buffer_size+1, 0);
  readno[0] = readno_start;
 
  int buffer_no = 0;
  size_t read = 0;
  while (read < len) {
    size_t buf_size = 0;
    for (; buf_size < read_buffer_size && read < len; ++buf_size, ++read) {
      const PositionRead &pr = Globals::pr->at(readno[buf_size]);
 
      ConstructRead(reads[buf_size], pr);
      read_size[buf_size] = reads[buf_size].trimNsAndBadQuality(cfg::get().input_trim_quality);
 
      readno[buf_size + 1] = readno[buf_size] + 1;
    }
    TIMEDLN("Prepared batch " << buffer_no << " of " << buf_size << " reads.");
 
    HammerTools::CorrectReadsBatch(res, reads,
                                   changedReads, changedNucleotides,
                                   readno, read_size, kmers);
 
    readno[0] = readno[buf_size];
 
    TIMEDLN("Processed batch " << buffer_no);
    for (size_t i = 0; i < buf_size; ++i) {
      reads[i].print(*(res[i] ? outf_good : outf_bad), qvoffset);
    }
    TIMEDLN("Written batch " << buffer_no);
    ++buffer_no;
  }
}

void HammerTools::CorrectPairedReadFiles(const vector<KMerCount> & kmers,
                                         size_t &changedReads, size_t &changedNucleotides,
                                         hint_t readno_left_start, hint_t readno_right_start, size_t len,
                                         ofstream * ofbadl, ofstream * ofcorl, ofstream * ofbadr, ofstream * ofcorr, ofstream * ofunp) {
  int qvoffset = cfg::get().input_qvoffset;

  unsigned correct_nthreads = min(cfg::get().correct_nthreads, cfg::get().general_max_nthreads);
  size_t read_buffer_size = correct_nthreads * cfg::get().correct_readbuffer;
  vector<Read> l(read_buffer_size);
  vector<Read> r(read_buffer_size);
  vector<bool> left_res(read_buffer_size, false);
  vector<bool> right_res(read_buffer_size, false);
  vector<size_t> read_size_left(read_buffer_size, 0);
  vector<size_t> read_size_right(read_buffer_size, 0);
  vector<hint_t> readno_left(read_buffer_size+1, 0);
  vector<hint_t> readno_right(read_buffer_size+1, 0);
  readno_left[0] = readno_left_start;
  readno_right[0] = readno_right_start;
 
  int buffer_no = 0;

  size_t readno = 0;
  while (readno < len) {
    size_t buf_size = 0;
    for (; buf_size < read_buffer_size && readno < len; ++buf_size, ++readno) {
      const PositionRead &prl = Globals::pr->at(readno_left[buf_size]);
      const PositionRead &prr = Globals::pr->at(readno_right[buf_size]);
 
      ConstructRead(l[buf_size], prl);
      ConstructRead(r[buf_size], prr);
 
      read_size_left[buf_size] = l[buf_size].trimNsAndBadQuality(cfg::get().input_trim_quality);
      read_size_right[buf_size] = r[buf_size].trimNsAndBadQuality(cfg::get().input_trim_quality);
 
      readno_left[buf_size + 1] = readno_left[buf_size] + 1;
      readno_right[buf_size + 1] = readno_right[buf_size] + 1;
    }
    TIMEDLN("Prepared batch " << buffer_no << " of " << buf_size << " reads.");
  
    HammerTools::CorrectReadsBatch(left_res, l,
                                   changedReads, changedNucleotides,
                                   readno_left, read_size_left, kmers);
    HammerTools::CorrectReadsBatch(right_res, r,
                                   changedReads, changedNucleotides,
                                   readno_right, read_size_right, kmers);
 
    readno_left[0] = readno_left[buf_size];
    readno_right[0] = readno_right[buf_size];
 
    TIMEDLN("Processed batch " << buffer_no);
    for (size_t i = 0; i < buf_size; ++i) {
      if (left_res[i] && right_res[i]) {
        l[i].print(*ofcorl, qvoffset);
        r[i].print(*ofcorr, qvoffset);
      } else {
        l[i].print(*(left_res[i] ? ofunp : ofbadl), qvoffset);
        r[i].print(*(right_res[i] ? ofunp : ofbadr), qvoffset);
      }
    }
    TIMEDLN("Written batch " << buffer_no);
    ++buffer_no;
  }
}

string getLargestPrefix(const string & str1, const string & str2) {
	string substr = "";
	for (size_t i = 0; i != str1.size() && i != str2.size(); ++i) {
		if (str1[i] == str2[i])
			substr += str1[i];
		else
			break;
	}
	return substr;
}

hint_t HammerTools::CorrectAllReads() {
	// Now for the reconstruction step; we still have the reads in rv, correcting them in place.
	size_t changedReads = 0;
	size_t changedNucleotides = 0;

	int correct_nthreads = min( cfg::get().correct_nthreads, cfg::get().general_max_nthreads );

	TIMEDLN("Starting read correction in " << correct_nthreads << " threads.");

	// correcting paired files
	bool single_created = false;
	if (Globals::input_filenames.size() >= 2) {
		int iFile = 0;
		if (Globals::input_filename_bases.size() != 3) {
			Globals::input_filename_bases.push_back(
					getLargestPrefix(Globals::input_filename_bases[0], Globals::input_filename_bases[1]) + "unpaired");
			single_created = true;
		}
		std::string input_single_filename_base = Globals::input_filename_bases[2];

		ofstream ofcorl(HammerTools::getReadsFilename(cfg::get().input_working_dir, iFile,   Globals::iteration_no, "cor").c_str());
		ofstream ofbadl(HammerTools::getReadsFilename(cfg::get().input_working_dir, iFile,   Globals::iteration_no, "bad").c_str());
		ofstream ofcorr(HammerTools::getReadsFilename(cfg::get().input_working_dir, iFile+1, Globals::iteration_no, "cor").c_str());
		ofstream ofbadr(HammerTools::getReadsFilename(cfg::get().input_working_dir, iFile+1, Globals::iteration_no, "bad").c_str());
		ofstream ofunp (HammerTools::getReadsFilename(cfg::get().input_working_dir, 2,       Globals::iteration_no, "cor").c_str());

    HammerTools::CorrectPairedReadFiles(*Globals::kmers,
                                        changedReads, changedNucleotides,
                                        Globals::input_file_blob_positions[iFile], Globals::input_file_blob_positions[iFile+1], Globals::input_file_sizes[iFile],
                                        &ofbadl, &ofcorl, &ofbadr, &ofcorr, &ofunp);
		TIMEDLN("  " << Globals::input_filenames[iFile].c_str() << " and " << Globals::input_filenames[iFile+1].c_str() << " corrected as a pair.");
		// makes sense to change the input filenames for the next iteration immediately
		Globals::input_filenames[iFile] = HammerTools::getReadsFilename(cfg::get().input_working_dir, iFile, Globals::iteration_no, "cor");
		Globals::input_filenames[iFile+1] = HammerTools::getReadsFilename(cfg::get().input_working_dir, iFile+1, Globals::iteration_no, "cor");
		// and single file
		if (single_created) {
			Globals::input_filenames.push_back(
					HammerTools::getReadsFilename(cfg::get().input_working_dir, 2, Globals::iteration_no, "cor"));
		}
		++iFile;
	}

	// correcting single file
	if (!single_created && (Globals::input_filenames.size() == 3 || Globals::input_filenames.size() == 1)) {
		int iFile = Globals::input_filenames.size() - 1;
		ofstream ofgood(HammerTools::getReadsFilename(cfg::get().input_working_dir, iFile, Globals::iteration_no, "cor").c_str(), fstream::app);
		ofstream ofbad( HammerTools::getReadsFilename(cfg::get().input_working_dir, iFile, Globals::iteration_no, "bad").c_str());
    HammerTools::CorrectReadFile(*Globals::kmers,
                                 changedReads, changedNucleotides,
                                 Globals::input_file_blob_positions[iFile], Globals::input_file_sizes[iFile],
                                 &ofgood, &ofbad);
		TIMEDLN("  " << Globals::input_filenames[iFile].c_str() << " corrected.");
		// makes sense to change the input filenames for the next iteration immediately
		Globals::input_filenames[iFile] = HammerTools::getReadsFilename(cfg::get().input_working_dir, iFile, Globals::iteration_no, "cor");
		// delete output files from previous iteration
		//if (Globals::iteration_no > 0) {
		//	HammerTools::RemoveFile(HammerTools::getReadsFilename(cfg::get().input_working_dir, iFile, Globals::iteration_no - 1, "cor"));
		//	HammerTools::RemoveFile(HammerTools::getReadsFilename(cfg::get().input_working_dir, iFile, Globals::iteration_no - 1, "bad"));
		//}
	}

	TIMEDLN("Correction done. Changed " << changedNucleotides << " bases in " << changedReads << " reads.");
	return changedReads;
}


void HammerTools::ReadKmerNosFromFile( const string & fname, vector<hint_t> *kmernos ) {
	kmernos->clear();
	boost::shared_ptr<FIStream> fis = FIStream::init_buf(fname, 1 << cfg::get().general_file_buffer_exp);
	hint_t pos;
	hint_t prev_pos = -1;
	string buf;
	while (fis->fs.good()) {
		std::getline(fis->fs, buf);
		sscanf(buf.c_str(), "%lu", &pos);
		if (pos == prev_pos) break;
		kmernos->push_back(pos);
		prev_pos = pos;
	}
}

void HammerTools::RemoveFile(const string & fname) {
	if (cfg::get().general_remove_temp_files) {
		if (boost::filesystem::exists(boost::filesystem::path(fname))) {
			if(remove(fname.c_str()) != 0) {
				TIMEDLN("Error deleting file " + fname);
			}
		}
	}
}
