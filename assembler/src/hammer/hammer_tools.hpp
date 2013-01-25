//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
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
#include "io/mmapped_reader.hpp"

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
	/// initialize subkmer positions and log about it
	static void InitializeSubKMerPositions();

	/// read one input file into the blob and output current position and read number
	static std::pair<size_t, size_t> ReadFileIntoBlob(const string & readsFilename, hint_t & curpos, hint_t & cur_read);
	/// read all input files into the blob
	static void ReadAllFilesIntoBlob();

	/// print out the resulting set of k-mers
	static void PrintKMerResult(std::ostream & outf, const vector<KMerStat> & kmers );

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

