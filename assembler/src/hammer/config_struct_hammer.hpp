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
	double iterative_reconstruction_threshold;
	int max_reconstruction_iterations;
};


// main debruijn config load function
void load(boost::property_tree::ptree const& pt, hammer_config& cfg)
{
	using config_common::load;
	// input options:
	load(pt, "working_dir", cfg.working_dir);
	load(pt, "reads", cfg.reads);
	load(pt, "tau", cfg.tau);
	load(pt, "num_threads", cfg.num_threads);
	load(pt, "num_iterations", cfg.num_iterations);
	load(pt, "quality_offset", cfg.quality_offset);

	load(pt, "paired_reads", cfg.paired_reads);
	if ( cfg.paired_reads ) {
		load(pt, "reads_left", cfg.reads_left);
		load(pt, "reads_right", cfg.reads_right);
	}

	load(pt, "read_blob_and_kmers", cfg.read_blob_and_kmers);
	load(pt, "write_blob_and_kmers", cfg.write_blob_and_kmers);
	load(pt, "exit_after_writing_blob_and_kmers", cfg.exit_after_writing_blob_and_kmers);
	if ( cfg.read_blob_and_kmers || cfg.write_blob_and_kmers ) {
		load(pt, "blob", cfg.blob);
		load(pt, "kmers", cfg.kmers);
	}

	load(pt, "read_kmers_after_clustering", cfg.read_kmers_after_clustering);
	load(pt, "write_kmers_after_clustering", cfg.write_kmers_after_clustering);
	load(pt, "kmers_after_clustering", cfg.kmers_after_clustering);

	load(pt, "error_rate", cfg.error_rate);
	load(pt, "blocksize_quadratic_threshold", cfg.blocksize_quadratic_threshold);
	load(pt, "good_cluster_threshold", cfg.good_cluster_threshold);
	load(pt, "blob_margin", cfg.blob_margin);
	load(pt, "trim_quality", cfg.trim_quality);

	load(pt, "trim_left_right", cfg.trim_left_right);
	load(pt, "use_iterative_reconstruction", cfg.use_iterative_reconstruction);
	load(pt, "iterative_reconstruction_threshold", cfg.iterative_reconstruction_threshold);
	load(pt, "max_reconstruction_iterations", cfg.max_reconstruction_iterations);
}

typedef config_common::config<hammer_config> cfg;

#endif

