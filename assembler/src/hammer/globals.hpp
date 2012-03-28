#ifndef HAMMER_GLOBALS_HPP_
#define HAMMER_GLOBALS_HPP_

#include "kmer_stat.hpp"
#include "kmerno.hpp"

struct Globals {
	static int iteration_no;
	static std::vector<std::string> input_filenames;
	static std::vector<std::string> input_filename_bases;
	static std::vector<hint_t> input_file_blob_positions;
	static std::vector<uint32_t> * subKMerPositions;
	static char* blob;
	static char* blobquality;
	static std::vector<PositionRead> * pr;
	static std::vector<hint_t> * kmernos;
	static std::vector<KMerCount> * kmers;
	static hint_t blob_max_size;
	static hint_t blob_size;
	static hint_t revNo;
	static char char_offset;

	static bool use_common_quality;
	static char common_quality;
	static double common_kmer_errprob;

	static void writeBlob( const char * fname );
	static void readBlob( const char * fname );
};

#endif //  HAMMER_GLOBALS_HPP_

