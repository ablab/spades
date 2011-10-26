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
	bool unload_blob_before_merge;

	bool likelihood_e_step;
	bool subtract_simplex_volume;

	bool debug_output_clustering;
	bool debug_output_likelihood;
};


// main debruijn config load function
void load(boost::property_tree::ptree const& pt, hammer_config& cfg);

typedef config_common::config<hammer_config> cfg;

#endif

