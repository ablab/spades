//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/*
 * hammer_tools.hpp
 *
 *  Created on: 08.07.2011
 *      Author: snikolenko
 */

#ifndef HAMMER_TOOLS_HPP
#define HAMMER_TOOLS_HPP

#include <unordered_map>
#include <algorithm>
#include <stdexcept>
#include <iomanip>
#include <fstream>
#include "read/read.hpp"
#include "read/ireadstream.hpp"
#include "sequence/seq.hpp"
#include "globals.hpp"
#include "kmer_stat.hpp"
#include "mmapped_reader.hpp"

using namespace std;

namespace hammer_tools {
/// estimate total read size in input read files
size_t EstimateTotalReadSize(const std::vector<std::string> &fnames);
};

/**
 * a container class for all general procedures in BayesHammer
 */
class HammerTools {
public:

	/// decompress gzipped input files if needed
	static void DecompressIfNeeded();

	/// change single Ns to As in input read files
	static void ChangeNtoAinReadFiles();

	/// initialize subkmer positions and log about it
	static void InitializeSubKMerPositions();

	/// read one input file into the blob and output current position and read number
	static size_t ReadFileIntoBlob(const string & readsFilename, hint_t & curpos, hint_t & cur_read);
	/// read all input files into the blob
	static void ReadAllFilesIntoBlob();

	/// leave only minimizers
	static void findMinimizers(vector< pair<hint_t, pair< double, size_t > > > & v, int num_minimizers,
                             vector< hint_t > & mmers, int which_first = 0 );
	/// check whether this is a minimizer iteration
	static bool doingMinimizers();

	/// do one step of iterative expansion, return the number of new solid k-mers
	static size_t IterativeExpansionStep(int expand_iter_no, int nthreads, KMerData &data);

	/// print out the resulting set of k-mers
	static void PrintKMerResult(std::ostream & outf, const vector<KMerStat> & kmers );

	/// correct one read
	static bool CorrectOneRead(const KMerData &data,
                             size_t & changedReads, size_t & changedNucleotides,
                             Read & r, bool correct_threshold, bool discard_singletons, bool discard_bad);
  /// parallel correction of batch of reads
	static void CorrectReadsBatch(std::vector<bool> &res, std::vector<Read> &reads, size_t buf_size,
                                size_t &changedReads, size_t &changedNucleotides,
                                const KMerData &data);
	/// correct reads in a given file
	static void CorrectReadFile(const KMerData &data,
                              size_t & changedReads, size_t & changedNucleotides,
                              const std::string &fname,
                              ofstream *outf_good, ofstream *outf_bad);
	/// correct reads in a given pair of files
	static void CorrectPairedReadFiles(const KMerData &data,
                                     size_t & changedReads, size_t & changedNucleotides,
                                     const std::string &fnamel, const std::string &fnamer,
                                     ofstream * ofbadl, ofstream * ofcorl, ofstream * ofbadr, ofstream * ofcorr, ofstream * ofunp);
	/// correct all reads
	static hint_t CorrectAllReads();

	static string getFilename( const string & dirprefix, const string & suffix );
	static string getFilename( const string & dirprefix, int iter_count, const string & suffix );
	static string getFilename( const string & dirprefix, int iter_count, const string & suffix, int suffix_num );
	static string getFilename( const string & dirprefix, int iter_count, const string & suffix, int suffix_num, const string & suffix2 );
	static string getFilename( const string & dirprefix, const string & suffix, int suffix_num );
	static string getReadsFilename( const string & dirprefix, int read_file_no, int iter_no, const string & suffix );
};



#endif

