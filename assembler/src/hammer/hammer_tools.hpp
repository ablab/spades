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

#include <tr1/unordered_map>
#include <algorithm>
#include <stdexcept>
#include <iomanip>
#include <fstream>
#include <zlib.h>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/filesystem.hpp>
#include "read/read.hpp"
#include "read/ireadstream.hpp"
#include "union.hpp"
#include "sequence/seq.hpp"
#include "globals.hpp"
#include "kmer_stat.hpp"
#include "position_kmer.hpp"

using namespace std;

#define MAX_INT_64 1000000000000000000

#define TIMEDLN(a) print_full_stats(); cout << a << endl

typedef Seq<K> Kmer;
//typedef std::tr1::unordered_map<Kmer, KMerCount, typename Kmer::hash, typename Kmer::equal_to> KMerMap;
typedef std::map<Kmer, KMerCount, Kmer::less2 > KMerMap;

double oct2phred(string qoct, int qvoffset);
string encode3toabyte (const string & s);
void print_time();
void print_mem_usage();
void print_stats();
void print_full_stats();

/// structure for boost istreams
struct FIStream {
	std::string fn;
	boost::iostreams::filtering_istream fs;
	std::ifstream stdstream;
	std::vector<char> buffer;
	bool remove_it;
	FIStream(const string & fname);
	FIStream(const string & fname, bool input_output);
	FIStream(const string & fname, bool input_output, uint64_t bufsize);
	~FIStream();

	static boost::shared_ptr<FIStream> init(const string & fname, bool input_output = false);
	static boost::shared_ptr<FIStream> init_buf(const string & fname, uint64_t bufsize);
};

/// structure for boost ostreams
struct FOStream {
	std::string fn;
	boost::iostreams::filtering_ostream fs;
	std::ofstream stdstream;
	std::vector<char> buffer;
	bool remove_it;
	FOStream(const string & fname);
	FOStream(const string & fname, bool input_output);
	FOStream(const string & fname, bool input_output, uint64_t bufsize);
	~FOStream();

	static boost::shared_ptr<FOStream> init(const string & fname, bool input_output = false);
	static boost::shared_ptr<FOStream> init_buf(const string & fname, uint64_t bufsize);
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

	/// estimate total read size in input read files
	static hint_t EstimateTotalReadSize();

	/// initialize subkmer positions and log about it
	static void InitializeSubKMerPositions();

	/// read one input file into the blob and output current position and read number
	static void ReadFileIntoBlob(const string & readsFilename, hint_t & curpos, hint_t & cur_read, bool reverse_complement);
	/// read all input files into the blob
	static void ReadAllFilesIntoBlob();

	/// process a k-mer hash file
	static void ProcessKmerHashFile( boost::iostreams::filtering_istream & inStream, KMerNoHashMap & km, std::vector<KMerCount> & kmcvec );
	static void ProcessKmerHashVector( const vector< pair<hint_t, double> > & sv, KMerNoHashMap & km, std::vector<KMerCount> & vkmc );
	/// print a processed k-mer hash file
	static void PrintProcessedKmerHashFile( boost::iostreams::filtering_ostream & outStream, hint_t & kmer_num, KMerNoHashMap & km );
	static void PrintProcessedKmerHashFile( boost::iostreams::filtering_ostream & outf, hint_t & kmer_num, std::vector<KMerCount> & km );
	static void KmerHashUnique(const std::vector<KMerNo> & vec, std::vector<KMerCount> & vkmc);
	/// count k-mers in input files
	static void CountKMersBySplitAndMerge();
	/// split kmers into files
	static void SplitKMers();

	/// leave only minimizers
	static void findMinimizers( vector< pair<hint_t, pair< double, size_t > > > & v, int num_minimizers, int which_first = 0 );
	/// check whether this is a minimizer iteration
	static bool doingMinimizers();
	/// fill map
	static void FillMapWithMinimizers( KMerMap & map );

	/// do one step of iterative expansion, return the number of new solid k-mers
	static hint_t IterativeExpansionStep(int expand_iter_no, int nthreads, vector<KMerCount> & kmers);

	/// print out the resulting set of k-mers
	static void PrintKMerResult( boost::iostreams::filtering_ostream & outf, const vector<KMerCount> & kmers );

	/// internal procedure
	static bool internalCorrectReadProcedure( const Read & r, const hint_t readno, const string & seq,
			const vector<KMerCount> & km, const PositionKMer & kmer, const uint32_t pos, const KMerStat & stat,
			vector< vector<int> > & v, int & left, int & right, bool & isGood, ofstream * ofs, bool revcomp,
			bool correct_threshold, bool discard_singletons);

	/// correct one read
	static bool CorrectOneRead( const vector<KMerCount> & kmers, hint_t & changedReads, hint_t & changedNucleotides,
			hint_t readno, Read & r, size_t i, bool correct_threshold, bool discard_singletons, bool discard_bad );
	/// correct reads in a given file
	static void CorrectReadFile( const string & readsFilename, const vector<KMerCount> & kmers, hint_t & changedReads, hint_t & changedNucleotides, hint_t readno_start, ofstream *outf_good, ofstream *outf_bad );
	/// correct reads in a given pair of files
	static void CorrectPairedReadFiles( const string & readsFilenameLeft, const string & readsFilenameRight,
			const vector<KMerCount> & kmers, hint_t & changedReads, hint_t & changedNucleotides, hint_t readno_left_start, hint_t readno_right_start,
			ofstream * ofbadl, ofstream * ofcorl, ofstream * ofbadr, ofstream * ofcorr, ofstream * ofunp );
	/// correct all reads
	static hint_t CorrectAllReads();

	/// read k-mer indices from file
	static void ReadKmerNosFromFile( const string & fname, vector<hint_t> *kmernos );
	/// read a k-mer totals file
	//static void ReadKmersAndNosFromFile( const string & fname, vector<KMerCount> *kmers, vector<hint_t> *kmernos );
	/// read a k-mer result file
	static void ReadKmersWithChangeToFromFile( const string & fname, vector<KMerCount> *kmers, vector<hint_t> *kmernos );

	static string getFilename( const string & dirprefix, const string & suffix );
	static string getFilename( const string & dirprefix, int iter_count, const string & suffix );
	static string getFilename( const string & dirprefix, int iter_count, const string & suffix, int suffix_num );
	static string getFilename( const string & dirprefix, int iter_count, const string & suffix, int suffix_num, const string & suffix2 );
	static string getFilename( const string & dirprefix, const string & suffix, int suffix_num );
	static string getReadsFilename( const string & dirprefix, int read_file_no, int iter_no, const string & suffix );
	static void RemoveFile(const string & fname);
};



#endif

