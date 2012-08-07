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
#include <boost/format.hpp>

#include <time.h>
#include <sys/resource.h>
#include <iomanip>

#include "read/ireadstream.hpp"
#include "mathfunctions.hpp"
#include "valid_kmer_generator.hpp"
#include "position_kmer.hpp"
#include "globals.hpp"
#include "config_struct_hammer.hpp"
#include "hammer_tools.hpp"
#include "mmapped_writer.hpp"
#include "kmer_index.hpp"

#include <sys/types.h>
#include <sys/wait.h>

#include <unordered_map>

using namespace std;

void HammerTools::ChangeNtoAinReadFiles() {
	std::vector<pid_t> pids;
	for (size_t iFile=0; iFile < Globals::input_filenames.size(); ++iFile) {
		string cur_filename = HammerTools::getFilename(cfg::get().input_working_dir, Globals::input_filename_bases[iFile]);
		pid_t cur_pid = vfork();
		if (cur_pid == 0) {
			INFO("  [" << getpid() << "] Child process for substituting Ns in " << Globals::input_filenames[iFile] << " starting.");
			string cmd = string("sed \'n;s/\\([ACGT]\\)N\\([ACGT]\\)/\\1A\\2/g;n;n\' ") + Globals::input_filenames[iFile].c_str() + " > " + cur_filename.c_str();
			int exitcode = system( cmd.c_str() );
			if (exitcode != 0) {
				INFO("  [" << getpid() << "] ERROR: finished with non-zero exit code " << exitcode);
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
        if (path::extension(Globals::input_filenames[iFile]) != ".gz")
			continue;

		stat(Globals::input_filenames[iFile].c_str(), &st);
		std::ifstream ifs(Globals::input_filenames[iFile], std::ios::in | std::ios::binary );
		ifs.read(f_cont, 2);
		ifs.close();
		if ( (f_cont[0] == (char)0x1f) && (f_cont[1] == (char)0x8b) ) {
			string newFilename = HammerTools::getFilename(cfg::get().input_working_dir, Globals::input_filename_bases[iFile]);
			string oldFilename = Globals::input_filenames[iFile];
			Globals::input_filenames[iFile] = newFilename;
            Globals::input_filename_bases[iFile] = path::basename(Globals::input_filename_bases[iFile]);

			pIDs[iFile] = vfork();
			if (pIDs[iFile] == 0) {
				string systemcall = string("gunzip -c ") + oldFilename + string(" > ") + newFilename;
				INFO("  [" << getpid() << "] " << systemcall);
				if (system(systemcall.c_str())) {
					INFO("  [" << getpid() << "] System error with unzipping input files!");
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

namespace hammer_tools {
size_t EstimateTotalReadSize(const std::vector<std::string> &fnames) {
  struct stat st;
  size_t totalReadSize = 0;
  for (auto I = fnames.begin(), E = fnames.end(); I != E; ++I) {
    stat(I->c_str(), &st);
    totalReadSize += st.st_size;
  }
  totalReadSize = totalReadSize / (2.5);
  return totalReadSize;
}
};

void HammerTools::InitializeSubKMerPositions() {
	ostringstream log_sstream;
	log_sstream.str("");
	Globals::subKMerPositions = new std::vector<uint32_t>(cfg::get().general_tau + 2);
	for (uint32_t i=0; i < (uint32_t)(cfg::get().general_tau + 1); ++i) {
		Globals::subKMerPositions->at(i) = (i * K / (cfg::get().general_tau + 1) );
		log_sstream << Globals::subKMerPositions->at(i) << " ";
	}
	Globals::subKMerPositions->at(cfg::get().general_tau + 1) = K;
	INFO("Hamming graph threshold tau=" << cfg::get().general_tau << ", k=" << K << ", subkmer positions = [ " << log_sstream.str() << "]" );
}

size_t HammerTools::ReadFileIntoBlob(const string & readsFilename, hint_t & curpos, hint_t & cur_read, bool reverse_complement) {
  INFO("Reading input file " << readsFilename);
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
      const char* qdata = q.data();
      // Verify user-provided character offset
      if (Globals::char_offset_user) {
        for (size_t i = 0; i < read_size; ++i)
          if (qdata[i] <= 0) {
            INFO(" Invalid quality value, probably phred offset specified was wrong");
            exit(-1);
          }
      }

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
  INFO("All files were read. Used " << curpos << " bytes out of " << Globals::blob_max_size << " allocated.");
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

bool HammerTools::internalCorrectReadProcedure(const std::string &seq, const KMerIndex &index,
                                               const PositionKMer &kmer, size_t pos, const KMerStat & stat,
                                               std::vector<std::vector<int> > & v,
                                               int & left, int & right, bool & isGood,
                                               ofstream * ofs,
                                               bool revcomp, bool correct_threshold, bool discard_singletons) {
  bool res = false;
  if (stat.isGoodForIterative() || (correct_threshold && stat.isGood())) {
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
    if (stat.change() &&
        (discard_singletons || index[stat.changeto].second.isGoodForIterative() ||
         (correct_threshold && stat.isGood()))) {
      if (ofs != NULL) *ofs << "\tchange to\n";
      //cout << "  kmer " << kmer.start() << " " << kmer.str() << " wants to change to " << stat.changeto << " " << km[stat.changeto].first.str() << endl;
      isGood = true;
      if ((int) pos < left)
        left = pos;
      if ((int) pos > right)
        right = pos;
      const PositionKMer & newkmer = index[stat.changeto].first;

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

size_t HammerTools::IterativeExpansionStep(int expand_iter_no, int nthreads, KMerIndex &index) {
  size_t res = 0;

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
    std::vector<unsigned> covered_by_solid(read_size, false);
    std::vector<hint_t> kmer_indices(read_size, -1);

    ValidKMerGenerator<K> gen(Globals::blob + pr.start(),
                              /* quality is not necessary */ NULL,
                              read_size);
    while (gen.HasMore()) {
      const Seq<K> &kmer = gen.kmer();
      auto it = index.seq_find(kmer);
      if (it != index.seq_end()) {
        size_t pos = it->second;
        size_t read_pos = gen.pos() - 1;

        kmer_indices[read_pos] = pos;
        if (index[pos].second.isGoodForIterative()) {
          for (size_t j = read_pos; j < read_pos + K; ++j)
            covered_by_solid[j] = true;
        }
      }
      gen.Next();
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
      if (kmer_indices[j] == (hint_t)-1 ) continue;
      if (!index[kmer_indices[j]].second.isGoodForIterative() &&
          !index[kmer_indices[j]].second.isMarkedGoodForIterative() ) {
#       pragma omp critical
        {
          ++res;
          index[kmer_indices[j]].second.makeGoodForIterative();
        }
      }
    }
  }

  if (cfg::get().expand_write_each_iteration) {
    ofstream oftmp( getFilename(cfg::get().input_working_dir, Globals::iteration_no, "goodkmers", expand_iter_no ).data());
    for (size_t n = 0; n < index.size(); ++n ) {
      if (index[n].second.isGoodForIterative() ) {
        oftmp << index[n].first.str() << "\n>" << index[n].first.start()
              << "  cnt=" << index[n].second.count << "  tql=" << (1-index[n].second.totalQual) << "\n";
      }
    }
  }

  return res;
}

void HammerTools::PrintKMerResult(std::ostream& outf, const vector<KMerCount> & kmers ) {
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

bool HammerTools::CorrectOneRead(const KMerIndex &index,
                                 size_t & changedReads, size_t & changedNucleotides,
                                 Read & r, bool correct_threshold, bool discard_singletons, bool discard_bad ) {
  bool isGood = false;
  std::string seq = r.getSequenceString();
  size_t read_size = seq.size();

  VERIFY(read_size >= K);

  // create auxiliary structures for consensus
  std::vector<int> vA(read_size, 0), vC(read_size, 0), vG(read_size, 0), vT(read_size, 0);
  std::vector<std::vector<int> > v;  // A=0, C=1, G=2, T=3
  v.push_back(vA); v.push_back(vC); v.push_back(vG); v.push_back(vT);
  isGood = false;

  // getting the leftmost and rightmost positions of a solid kmer
  int left = read_size; int right = -1;
  bool changedRead = false;
  ValidKMerGenerator<K> gen(seq.data(),
                            /* quality is not necessary */ NULL,
                            read_size);
  while (gen.HasMore()) {
    const Seq<K> &kmer = gen.kmer();
    auto it = index.seq_find(kmer);
    if (it != index.seq_end()) {
      size_t pos = it->second;
      size_t read_pos = gen.pos() - 1;

      const PositionKMer &kmer = index[pos].first;
      const KMerStat &stat = index[pos].second;
      changedRead = changedRead ||
                    internalCorrectReadProcedure(seq, index, kmer, read_pos, stat, v,
                                                 left, right, isGood, NULL, false, correct_threshold, discard_singletons);
    }
    gen.Next();
  }

  int left_rev = 0; int right_rev = read_size-(int)K;

  if (left <= right && left_rev <= right_rev) {
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
 

void HammerTools::CorrectReadsBatch(std::vector<bool> &res,
                                    std::vector<Read> &reads, size_t buf_size,
                                    size_t &changedReads, size_t &changedNucleotides,
                                    const KMerIndex &index) {
  unsigned correct_nthreads = min(cfg::get().correct_nthreads, cfg::get().general_max_nthreads);
  bool discard_singletons = cfg::get().bayes_discard_only_singletons;
  bool correct_threshold = cfg::get().correct_use_threshold;
  bool discard_bad = cfg::get().correct_discard_bad && !cfg::get().correct_notrim;
  if (HammerTools::doingMinimizers()) discard_bad = false;
 
  std::vector<size_t> changedReadBuf(correct_nthreads, 0);
  std::vector<size_t> changedNuclBuf(correct_nthreads, 0);
 
# pragma omp parallel for shared(reads, res, index) num_threads(correct_nthreads)
  for (size_t i = 0; i < buf_size; ++i) {
    if (reads[i].size() >= K) {
      res[i] =
          HammerTools::CorrectOneRead(index,
                                      changedReadBuf[omp_get_thread_num()], changedNuclBuf[omp_get_thread_num()],
                                      reads[i],
                                      correct_threshold, discard_singletons, discard_bad);
    } else
      res[i] = false;
  }
 
  for (unsigned i = 0; i < correct_nthreads; ++i ) {
    changedReads += changedReadBuf[i];
    changedNucleotides += changedNuclBuf[i];
  }
}
 
void HammerTools::CorrectReadFile(const KMerIndex &index,
                                  size_t & changedReads, size_t & changedNucleotides,
                                  const std::string &fname,
                                  ofstream *outf_good, ofstream *outf_bad ) {
  int qvoffset = cfg::get().input_qvoffset;
  int trim_quality = cfg::get().input_trim_quality;

  unsigned correct_nthreads = min(cfg::get().correct_nthreads, cfg::get().general_max_nthreads);
  size_t read_buffer_size = correct_nthreads * cfg::get().correct_readbuffer;
  std::vector<Read> reads(read_buffer_size);
  std::vector<bool> res(read_buffer_size, false);

  ireadstream irs(fname, qvoffset);
  VERIFY(irs.is_open());

  unsigned buffer_no = 0;
  while (!irs.eof()) {
    size_t buf_size = 0;
    for (; buf_size < read_buffer_size && !irs.eof(); ++buf_size) {
      irs >> reads[buf_size];
      reads[buf_size].trimNsAndBadQuality(trim_quality);
    }
    INFO("Prepared batch " << buffer_no << " of " << buf_size << " reads.");

    HammerTools::CorrectReadsBatch(res, reads, buf_size,
                                   changedReads, changedNucleotides,
                                   index);

    INFO("Processed batch " << buffer_no);
    for (size_t i = 0; i < buf_size; ++i) {
      reads[i].print(*(res[i] ? outf_good : outf_bad), qvoffset);
    }
    INFO("Written batch " << buffer_no);
    ++buffer_no;
  }
}

void HammerTools::CorrectPairedReadFiles(const KMerIndex &index,
                                         size_t &changedReads, size_t &changedNucleotides,
                                         const std::string &fnamel, const std::string &fnamer,
                                         ofstream * ofbadl, ofstream * ofcorl, ofstream * ofbadr, ofstream * ofcorr, ofstream * ofunp) {
  int qvoffset = cfg::get().input_qvoffset;
  int trim_quality = cfg::get().input_trim_quality;

  unsigned correct_nthreads = min(cfg::get().correct_nthreads, cfg::get().general_max_nthreads);
  size_t read_buffer_size = correct_nthreads * cfg::get().correct_readbuffer;
  std::vector<Read> l(read_buffer_size);
  std::vector<Read> r(read_buffer_size);
  std::vector<bool> left_res(read_buffer_size, false);
  std::vector<bool> right_res(read_buffer_size, false);
  

  unsigned buffer_no = 0;

  ireadstream irsl(fnamel, qvoffset), irsr(fnamer, qvoffset);
  VERIFY(irsl.is_open()); VERIFY(irsr.is_open());

  while (!irsl.eof() && !irsr.eof()) {
    size_t buf_size = 0;
    for (; buf_size < read_buffer_size && !irsl.eof() && !irsr.eof(); ++buf_size) {
      irsl >> l[buf_size]; irsr >> r[buf_size];
      l[buf_size].trimNsAndBadQuality(trim_quality);
      r[buf_size].trimNsAndBadQuality(trim_quality);
    }
    INFO("Prepared batch " << buffer_no << " of " << buf_size << " reads.");
  
    HammerTools::CorrectReadsBatch(left_res, l, buf_size,
                                   changedReads, changedNucleotides,
                                   index);
    HammerTools::CorrectReadsBatch(right_res, r, buf_size,
                                   changedReads, changedNucleotides,
                                   index);
 
    INFO("Processed batch " << buffer_no);
    for (size_t i = 0; i < buf_size; ++i) {
      if (left_res[i] && right_res[i]) {
        l[i].print(*ofcorl, qvoffset);
        r[i].print(*ofcorr, qvoffset);
      } else {
        l[i].print(*(left_res[i] ? ofunp : ofbadl), qvoffset);
        r[i].print(*(right_res[i] ? ofunp : ofbadr), qvoffset);
      }
    }
    INFO("Written batch " << buffer_no);
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

	INFO("Starting read correction in " << correct_nthreads << " threads.");

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

    HammerTools::CorrectPairedReadFiles(*Globals::kmer_index,
                                        changedReads, changedNucleotides,
                                        Globals::input_filenames[iFile], Globals::input_filenames[iFile+1],
                                        &ofbadl, &ofcorl, &ofbadr, &ofcorr, &ofunp);
		INFO("  " << Globals::input_filenames[iFile].c_str() << " and " << Globals::input_filenames[iFile+1].c_str() << " corrected as a pair.");
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
    HammerTools::CorrectReadFile(*Globals::kmer_index,
                                 changedReads, changedNucleotides,
                                 Globals::input_filenames[iFile],
                                 &ofgood, &ofbad);
		INFO("  " << Globals::input_filenames[iFile].c_str() << " corrected.");
		// makes sense to change the input filenames for the next iteration immediately
		Globals::input_filenames[iFile] = HammerTools::getReadsFilename(cfg::get().input_working_dir, iFile, Globals::iteration_no, "cor");
		// delete output files from previous iteration
		//if (Globals::iteration_no > 0) {
		//	HammerTools::RemoveFile(HammerTools::getReadsFilename(cfg::get().input_working_dir, iFile, Globals::iteration_no - 1, "cor"));
		//	HammerTools::RemoveFile(HammerTools::getReadsFilename(cfg::get().input_working_dir, iFile, Globals::iteration_no - 1, "bad"));
		//}
	}

	INFO("Correction done. Changed " << changedNucleotides << " bases in " << changedReads << " reads.");
	return changedReads;
}
