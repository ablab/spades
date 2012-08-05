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
#include <zlib.h>
#include "read/read.hpp"
#include "read/ireadstream.hpp"
#include "union.hpp"
#include "sequence/seq.hpp"
#include "globals.hpp"
#include "kmer_stat.hpp"
#include "position_kmer.hpp"
#include "mmapped_reader.hpp"

using namespace std;

typedef Seq<K> Kmer;
typedef std::map<Kmer, KMerCount, Kmer::less2 > KMerMap;


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
	static size_t ReadFileIntoBlob(const string & readsFilename, hint_t & curpos, hint_t & cur_read, bool reverse_complement);
	/// read all input files into the blob
	static void ReadAllFilesIntoBlob();

	/// process a k-mer hash file
	static void ProcessKmerHashFile(const std::string &fname, std::vector<KMerCount> & kmcvec );
	/// count k-mers in input files
	static void CountKMersBySplitAndMerge();
	/// split kmers into files
	static void SplitKMers();

	/// leave only minimizers
	static void findMinimizers(vector< pair<hint_t, pair< double, size_t > > > & v, int num_minimizers,
                             vector< hint_t > & mmers, int which_first = 0 );
	/// check whether this is a minimizer iteration
	static bool doingMinimizers();
	/// fill map
	static void FillMapWithMinimizers( KMerMap & map );

	/// do one step of iterative expansion, return the number of new solid k-mers
	static size_t IterativeExpansionStep(int expand_iter_no, int nthreads, KMerIndex &index);

	/// print out the resulting set of k-mers
	static void PrintKMerResult(std::ostream & outf, const vector<KMerCount> & kmers );

	/// internal procedure
	static bool internalCorrectReadProcedure(const std::string & seq,
                                           const KMerIndex &index, const PositionKMer & kmer, size_t pos, const KMerStat & stat,
                                           std::vector<std::vector<int> > & v, int & left, int & right, bool & isGood,
                                           ofstream * ofs,
                                           bool revcomp, bool correct_threshold, bool discard_singletons);

	/// correct one read
	static bool CorrectOneRead(const KMerIndex &index,
                             size_t & changedReads, size_t & changedNucleotides,
                             Read & r, bool correct_threshold, bool discard_singletons, bool discard_bad);
  /// parallel correction of batch of reads
	static void CorrectReadsBatch(std::vector<bool> &res, std::vector<Read> &reads, size_t buf_size,
                                size_t &changedReads, size_t &changedNucleotides,
                                const KMerIndex &index);
	/// correct reads in a given file
	static void CorrectReadFile(const KMerIndex &index,
                              size_t & changedReads, size_t & changedNucleotides,
                              const std::string &fname,
                              ofstream *outf_good, ofstream *outf_bad);
	/// correct reads in a given pair of files
	static void CorrectPairedReadFiles(const KMerIndex &index,
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

