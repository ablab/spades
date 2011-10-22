#ifndef HAMMER_GLOBALS_HPP_
#define HAMMER_GLOBALS_HPP_

#include "kmer_stat.hpp"
#include "kmerno.hpp"

struct Globals {
	static int qvoffset;
	static double error_rate;
	static std::string working_dir;
	static int blocksize_quadratic_threshold;
	static double good_cluster_threshold;
	static double blob_margin;
	static bool paired_reads;
	static int trim_quality;
	static bool trim_left_right;
	static bool use_iterative_reconstruction;
	static bool reconstruction_in_full_iterations;
	static double iterative_reconstruction_threshold;
	static bool write_each_iteration_kmers;
	static int max_reconstruction_iterations;
	static bool read_kmers_after_clustering;
	static bool write_kmers_after_clustering;
	static bool regular_threshold_for_correction;
	static double special_nonsingleton_threshold;
	static bool discard_only_singletons;
	static string kmers_after_clustering;
	static bool use_true_likelihood;

	static bool conserve_memory;
	static int num_of_tmp_files;
	static int iteration_no;
	static bool skip_to_clustering;
	static bool skip_to_subvectors;
	static bool unload_blob_before_merge;

	static bool likelihood_e_step;
	static bool subtract_simplex_volume;

	static bool debug_output_clustering;
	static bool debug_output_likelihood;

	static std::vector<PositionRead> * pr;
	static std::vector<Read> * rv;
	static std::vector<bool> * rv_bad;
	static std::vector<Read> * rvLeft;
	static std::vector<Read> * rvRight;
	static std::vector<Read> * rvLeft_bad;
	static std::vector<Read> * rvRight_bad;
	static hint_t revNo;
	static hint_t lastLeftNo;
	static std::vector<hint_t> * kmernos;

	static char* blob;
	static char* blobquality;
	static char* totalquality;
	static hint_t blob_max_size;
	static hint_t blob_size;

	// static uint64_t* blobhash;
	static KMerNoHashMap hm;

	static std::vector<uint32_t> * subKMerPositions;

	static void writeBlob( const char * fname );
	static void readBlob( const char * fname );
	static void writeBlobKMers( const char * fname );
	static void readBlobKMers( const char * fname );
	static void writeKMerCounts( const char * fname, const std::vector<KMerCount*> & kmers );
	static void writeKMerHashMap( const char * fname, const KMerNoHashMap & hm );
	static void readKMerCounts( const char * fname, std::vector<KMerCount*> * kmers );
	static void readKMerHashMap( const char * fname, KMerNoHashMap * hm, std::vector<KMerCount*> * kmers );

	static void readKMerFiles( const char * f_num, const char * f_solid, const char * f_bad, std::vector<KMerCount*> * kmers );

};

#endif //  HAMMER_GLOBALS_HPP_

