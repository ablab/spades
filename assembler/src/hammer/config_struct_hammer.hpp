/*
 * config_struct_hammer.hpp
 *
 *  Created on: Aug 15, 2011
 *      Author: snikolenko
 */

#ifndef CONFIG_STRUCT_HAMMER_HPP_
#define CONFIG_STRUCT_HAMMER_HPP_

#include "config_common.hpp"

#define CONFIG_FILENAME "/home/snikolenko/algorithmic-biology/assembler/src/hammer/config.inp"

// struct for debruijn project's configuration file
struct hammer_config
{
	std::string working_dir;
	std::string reads;
	int tau;
	int num_threads;
	int num_iterations;
	int quality_offset;

	bool read_blob_and_kmers;
	bool write_blob_and_kmers;
	bool exit_after_writing_blob_and_kmers;
	std::string blob;
	std::string kmers;

	bool paired_reads;
	std::string reads_left;
	std::string reads_right;

	bool read_kmers_after_clustering;
	bool write_kmers_after_clustering;
	std::string kmers_after_clustering;

	double error_rate;
	int blocksize_quadratic_threshold;
	double good_cluster_threshold;
	double blob_margin;
	int trim_quality;

	bool trim_left_right;
	bool use_iterative_reconstruction;
	bool reconstruction_in_full_iterations;
	double iterative_reconstruction_threshold;
	int max_reconstruction_iterations;
	bool write_each_iteration_kmers;
	bool regular_threshold_for_correction;
	bool discard_only_singletons;
	double special_nonsingleton_threshold;
	bool use_true_likelihood;

	bool conserve_memory;
	int num_of_tmp_files;
	bool skip_to_clustering;
	bool skip_to_subvectors;
	bool skip_sorting_subvectors;
	bool skip_cluster_merging;
	int skip_iterative;
	int hamming_class_buffer;
	bool unload_blob_before_merge;

	bool likelihood_e_step;
	bool subtract_simplex_volume;
	bool change_n_to_random;

	bool debug_output_clustering;
	bool debug_output_likelihood;


	// new options
	// TODO: remove old options when these work
	int input_numfiles;
	bool input_paired;
	string input_file_0;
	string input_file_1;
	string input_file_2;
	string input_file_3;
	string input_file_4;
	string input_solid_kmers;
	string input_working_dir;
	int input_trim_quality;
	int input_qvoffset;
	bool input_read_solid_kmers;

	bool general_do_everything_after_first_iteration;
	bool general_reconstruct_only;
	bool general_change_n_to_a;
	int general_hard_memory_limit;
	int general_max_nthreads;
	int general_tau;
	int general_max_iterations;
	double general_blob_margin;

	bool count_do;
	int count_numfiles;
	int count_merge_nthreads;

	bool sort_do;

	bool subvectors_do;
	int subvectors_blocksize_quadratic_threshold;

	bool hamming_do;
	bool hamming_write_solid_kmers;
	bool hamming_write_bad_kmers;

	bool bayes_do;
	int bayes_nthreads;
	double bayes_quality_threshold;
	double bayes_singleton_threshold;
	double bayes_nonsingleton_threshold;
	bool bayes_discard_only_singletons;
	bool bayes_debug_output;
	bool bayes_use_hamming_dist;

	bool expand_do;
	int expand_max_iterations;
	int expand_nthreads;
	bool expand_write_each_iteration;
	bool expand_write_kmers_result;

	bool correct_do;
	bool correct_use_threshold;
	double correct_threshold;
	int correct_nthreads;
};


// main debruijn config load function
void load(hammer_config& cfg, boost::property_tree::ptree const& pt);

typedef config_common::config<hammer_config> cfg;

#endif

