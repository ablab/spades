/*
 * config_struct_hammer.hpp
 *
 *  Created on: Aug 15, 2011
 *      Author: snikolenko
 */

#ifndef CONFIG_STRUCT_HAMMER_HPP_
#define CONFIG_STRUCT_HAMMER_HPP_

#include "../debruijn/config_common.hpp"

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
};


// main debruijn config load function
void load(boost::property_tree::ptree const& pt, hammer_config& cfg)
{
	// input options:
	load(pt, "working_dir", cfg.working_dir);
	load(pt, "reads", cfg.reads);
	load(pt, "tau", cfg.tau);
	load(pt, "num_threads", cfg.num_threads);
	load(pt, "num_iterations", cfg.num_iterations);
	load(pt, "quality_offset", cfg.quality_offset);

	load(pt, "read_blob_and_kmers", cfg.read_blob_and_kmers);
	load(pt, "write_blob_and_kmers", cfg.write_blob_and_kmers);
	load(pt, "exit_after_writing_blob_and_kmers", cfg.exit_after_writing_blob_and_kmers);
	if ( cfg.read_blob_and_kmers || cfg.write_blob_and_kmers ) {
		load(pt, "blob", cfg.blob);
		load(pt, "kmers", cfg.kmers);
	}
}

typedef config<hammer_config> cfg;

#endif

