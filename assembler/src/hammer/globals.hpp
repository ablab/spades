#ifndef HAMMER_GLOBALS_HPP_
#define HAMMER_GLOBALS_HPP_

#include "kmer_stat.hpp"
#include "kmerno.hpp"

struct Globals {
	static int qvoffset;
	static double error_rate;
	static int blocksize_quadratic_threshold;
	static double good_cluster_threshold;
	static double blob_margin;
	static bool paired_reads;
	static int trim_quality;


	static std::vector<PositionRead> * pr;
	static std::vector<Read> * rv;
	static std::vector<bool> * rv_bad;
	static std::vector<Read> * rvLeft;
	static std::vector<Read> * rvRight;
	static std::vector<Read> * rvLeft_bad;
	static std::vector<Read> * rvRight_bad;
	static hint_t revNo;
	static hint_t lastLeftNo;

	static char* blob;
	static char* blobquality;
	static hint_t blob_max_size;
	static hint_t blob_size;

	static uint64_t* blobhash;
	static KMerNoHashMap hm;

	static std::vector<uint32_t> * subKMerPositions;

	static void writeBlob( const char * fname );
	static void readBlob( const char * fname );
	static void writeBlobKMers( const char * fname );
	static void readBlobKMers( const char * fname );
	static void writeKMerCounts( const char * fname, const std::vector<KMerCount*> & kmers );
	static void readKMerCounts( const char * fname, std::vector<KMerCount*> * kmers );
};

#endif //  HAMMER_GLOBALS_HPP_

