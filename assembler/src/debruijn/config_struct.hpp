/*
 * config_struct.hpp
 *
 *  Created on: Aug 9, 2011
 *      Author: Alexey.Gurevich
 */

#ifndef CONFIG_STRUCT_HPP_
#define CONFIG_STRUCT_HPP_

#include "config_common.hpp"

const char* const CONFIG_FILENAME = "./src/debruijn/config.inp";
const size_t K = 55; // must be odd (so there is no k-mer which is equal to it's reverse-complimentary k-mer)

// struct for debruijn project's configuration file
struct debruijn_config
{
	struct tip_clipper
	{
	   size_t max_tip_length;
	   size_t max_coverage;
	   double max_relative_coverage;
	};

	struct bulge_remover
	{
		size_t max_length_div_K;
		double max_coverage;
		double max_relative_coverage;
		double max_delta;
		double max_relative_delta;
	};

	struct erroneous_connections_remover
	{
		double max_coverage;
		int max_length_div_K;
	};

	struct distance_estimator
	{
		size_t delta;
		size_t linkage_distance;
		size_t max_distance;
	};

	struct dataset
	{
		std::string first;
		std::string second;
		size_t IS;
		int LEN;
	};

	std::string input_dir;
	std::string output_dir;
	std::string dataset_name;
	std::string reference_genome;

	bool paired_mode;
	bool rectangle_mode;
	bool etalon_info_mode;
	bool from_saved_graph;

	std::string uncorrected_reads;

	tip_clipper tc;
	bulge_remover br;
	erroneous_connections_remover ec;
	distance_estimator de;

	dataset ds;
};


// specific load functions

void load(boost::property_tree::ptree const& pt, debruijn_config::tip_clipper& tc)
{
	load(pt, "max_tip_length", tc.max_tip_length);
	load(pt, "max_coverage", tc.max_coverage);
	load(pt, "max_relative_coverage", tc.max_relative_coverage);
}

void load(boost::property_tree::ptree const& pt, debruijn_config::bulge_remover& br)
{
	load(pt, "max_length_div_K", br.max_length_div_K);
	load(pt, "max_coverage", br.max_coverage);
	load(pt, "max_relative_coverage", br.max_relative_coverage);
	load(pt, "max_delta", br.max_delta);
	load(pt, "max_relative_delta", br.max_relative_delta);
}

void load(boost::property_tree::ptree const& pt, debruijn_config::erroneous_connections_remover& ec)
{
	load(pt, "max_coverage", ec.max_coverage);
	load(pt, "max_length_div_K", ec.max_length_div_K);
}

void load(boost::property_tree::ptree const& pt, debruijn_config::distance_estimator& de)
{
	load(pt, "delta", de.delta);
	load(pt, "linkage_distance", de.linkage_distance);
	load(pt, "max_distance", de.max_distance);
}

void load(boost::property_tree::ptree const& pt, debruijn_config::dataset& ds)
{
	load(pt, "first", ds.first);
	load(pt, "second", ds.second);
	load(pt, "IS", ds.IS);
	load(pt, "LEN", ds.LEN);
}

// main debruijn config load function
void load(boost::property_tree::ptree const& pt, debruijn_config& cfg)
{
	// input options:
	load(pt, "input_dir", cfg.input_dir);
	load(pt, "output_dir", cfg.output_dir);
	load(pt, "dataset", cfg.dataset_name);
	load(pt, "reference_genome", cfg.reference_genome);

	load(pt, "paired_mode", cfg.paired_mode);
	load(pt, "rectangle_mode", cfg.rectangle_mode);
	load(pt, "etalon_info_mode", cfg.etalon_info_mode);
	load(pt, "from_saved_graph", cfg.from_saved_graph);

	load(pt, "tc", cfg.tc); // tip clipper:
	load(pt, "br", cfg.br); // bulge remover:
	load(pt, "ec", cfg.ec); // erroneous connections remover:
	load(pt, "de", cfg.de); // distance estimator:
	load(pt, "uncorrected_reads", cfg.uncorrected_reads);
	load(pt, cfg.dataset_name, cfg.ds);
}

typedef config<debruijn_config> cfg;

#endif




