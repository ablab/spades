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
}

typedef config<hammer_config> cfg;

#endif

