/*
 * config_struct.hpp
 *
 *  Created on: Aug 9, 2011
 *      Author: Alexey.Gurevich
 */

#ifndef CONFIG_STRUCT_HPP_
#define CONFIG_STRUCT_HPP_

#include "config_common.hpp"
#include <boost/bimap.hpp>
#include <sys/types.h>
#include <sys/stat.h>

namespace debruijn_graph
{
    enum working_stage
    {
        ws_construction,
        ws_paired_info_count,
        ws_simplification,
        ws_distance_estimation,
        ws_repeats_resolving,
        ws_n50_enlargement
    };

	const char* const cfg_filename = "./src/debruijn/config.info";
	const size_t K = 55; // must be odd (so there is no k-mer which is equal to it's reverse-complimentary k-mer)

	inline std::string MakeLaunchTimeDirName() {
		time_t rawtime;
		struct tm * timeinfo;
		char buffer[80];

		time(&rawtime);
		timeinfo = localtime(&rawtime);

		strftime(buffer, 80, "%m.%d_%H_%M", timeinfo);
		return std::string(buffer);
	}


	// struct for debruijn project's configuration file
	struct debruijn_config
	{
		typedef boost::bimap<std::string, working_stage> name_id_mapping;

		static const name_id_mapping FillStageInfo() {
			name_id_mapping working_stages_info;
			working_stages_info.insert(name_id_mapping::value_type("construction"       , ws_construction       ));
			working_stages_info.insert(name_id_mapping::value_type("paired_info_count"  , ws_paired_info_count  ));
			working_stages_info.insert(name_id_mapping::value_type("simplification"	    , ws_simplification     ));
			working_stages_info.insert(name_id_mapping::value_type("distance_estimation", ws_distance_estimation));
			working_stages_info.insert(name_id_mapping::value_type("repeats_resolving"	, ws_repeats_resolving  ));
			working_stages_info.insert(name_id_mapping::value_type("n50_enlargement"	, ws_n50_enlargement    ));
			return working_stages_info;
		}

		static const name_id_mapping& working_stages_info() {
			static name_id_mapping working_stages_info = FillStageInfo();
			return working_stages_info;
		}

		static const std::string& working_stage_name(working_stage stage_id) {
			name_id_mapping::right_const_iterator it = working_stages_info().right.find(stage_id);
			FATAL_ASSERT(it != working_stages_info().right.end(), "No name for working stage id = " << stage_id);
			return it->second;
		}

		static working_stage working_stage_id(std::string name) {
			name_id_mapping::left_const_iterator it = working_stages_info().left.find(name);
			FATAL_ASSERT(it != working_stages_info().left.end(), "There is no working stage with name = " << name);
			return it->second;
		}

		struct tip_clipper
		{
		   double max_tip_length_div_K;
		   double max_coverage;
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

		struct cheating_erroneous_connections_remover
		{
			size_t max_length;
			double coverage_gap;
			size_t sufficient_neighbour_length;
		};

		struct repeat_resolver
		{
			int mode;
			int near_vertex;
		};
		struct distance_estimator
		{
			size_t delta;
			size_t linkage_distance;
			size_t max_distance;
		};
		struct advanced_distance_estimator
		{
            size_t threshold;
            double range_coeff;
            double delta_coeff;
            double percentage;
			size_t cutoff;
			size_t minpeakpoints;
			double inv_density;
            double derivative_threshold;
            
		};

		struct dataset
		{
			std::string first;
			std::string second;
			size_t RL;
			size_t IS;
			int LEN;
		};

        std::string dataset_name;

        std::string input_dir;

		std::string output_root;
		std::string output_dir;
		std::string output_suffix;
		std::string output_saves;

		std::string reference_genome;

		std::string load_from;

		working_stage entry_point;

		bool paired_mode;
		bool rectangle_mode;
		bool etalon_info_mode;
		bool late_paired_info;
		bool advanced_estimator_mode;

		std::string uncorrected_reads;
		bool need_consensus;
		tip_clipper tc;
		bulge_remover br;
		erroneous_connections_remover ec;
		cheating_erroneous_connections_remover cec;
		distance_estimator de;
		advanced_distance_estimator ade;
		repeat_resolver rr;
		dataset ds;
	};


	// specific load functions

	inline void load(boost::property_tree::ptree const& pt, debruijn_config::tip_clipper& tc)
	{
		using config_common::load;
		load(pt, "max_tip_length_div_K", tc.max_tip_length_div_K);
		load(pt, "max_coverage", tc.max_coverage);
		load(pt, "max_relative_coverage", tc.max_relative_coverage);
	}

	inline void load(boost::property_tree::ptree const& pt, std::string const& key, working_stage& entry_point)
	{
		std::string ep = pt.get<std::string>(key);
		entry_point = debruijn_config::working_stage_id(ep);
	}

	inline void load(boost::property_tree::ptree const& pt, debruijn_config::bulge_remover& br)
	{
		using config_common::load;
		load(pt, "max_length_div_K", br.max_length_div_K);
		load(pt, "max_coverage", br.max_coverage);
		load(pt, "max_relative_coverage", br.max_relative_coverage);
		load(pt, "max_delta", br.max_delta);
		load(pt, "max_relative_delta", br.max_relative_delta);
	}

	inline void load(boost::property_tree::ptree const& pt, debruijn_config::erroneous_connections_remover& ec)
	{
		using config_common::load;
		load(pt, "max_coverage", ec.max_coverage);
		load(pt, "max_length_div_K", ec.max_length_div_K);
	}

	inline void load(boost::property_tree::ptree const& pt, debruijn_config::cheating_erroneous_connections_remover& cec)
	{
		using config_common::load;
		load(pt, "max_length", cec.max_length);
		load(pt, "coverage_gap", cec.coverage_gap);
		load(pt, "sufficient_neighbour_length", cec.sufficient_neighbour_length);
	}

	inline void load(boost::property_tree::ptree const& pt, debruijn_config::distance_estimator& de)
	{
		using config_common::load;
		load(pt, "delta", de.delta);
		load(pt, "linkage_distance", de.linkage_distance);
		load(pt, "max_distance", de.max_distance);
	}

	inline void load(boost::property_tree::ptree const& pt, debruijn_config::advanced_distance_estimator& ade)
	{
		using config_common::load;
		load(pt, "threshold", ade.threshold);
		load(pt, "range_coeff", ade.range_coeff);
		load(pt, "delta_coeff", ade.delta_coeff);
		load(pt, "percentage", ade.percentage);
		load(pt, "cutoff", ade.cutoff);
		load(pt, "minpeakpoints", ade.minpeakpoints);
		load(pt, "inv_density", ade.inv_density);
		load(pt, "derivative_threshold", ade.derivative_threshold);
	}

	inline void load(boost::property_tree::ptree const& pt, debruijn_config::repeat_resolver& rr)
	{
		using config_common::load;
		load(pt, "mode", rr.mode);
		load(pt, "near_vertex", rr.near_vertex);
	}

	inline void load(boost::property_tree::ptree const& pt, debruijn_config::dataset& ds)
	{
		using config_common::load;
		load(pt, "first", ds.first);
		load(pt, "second", ds.second);
		load(pt, "RL", ds.RL);
		load(pt, "IS", ds.IS);
		load(pt, "LEN", ds.LEN);
	}

	// main debruijn config load function
	inline void load(boost::property_tree::ptree const& pt, debruijn_config& cfg)
	{
		using config_common::load;
		// input options:
        load(pt, "dataset", cfg.dataset_name);

        load(pt, "input_dir"  , cfg.input_dir);

		std::string output_base;
		load(pt, "output_base", output_base);

		cfg.output_root  = output_base + cfg.dataset_name + "/K" + ToString(K) + "/";
		cfg.output_suffix= MakeLaunchTimeDirName() + "/";
		cfg.output_dir   = cfg.output_root + cfg.output_suffix;
		cfg.output_saves = cfg.output_dir + "saves/";

        load(pt, "load_from", cfg.load_from);
        cfg.load_from = cfg.output_root + cfg.load_from;

		load(pt, "entry_point", cfg.entry_point);

		load(pt, "reference_genome", cfg.reference_genome);
		//load(pt, "start_from", cfg.start_from);

		load(pt, "paired_mode", cfg.paired_mode);
		load(pt, "rectangle_mode", cfg.rectangle_mode);
		load(pt, "etalon_info_mode", cfg.etalon_info_mode);
		load(pt, "late_paired_info", cfg.late_paired_info);
		load(pt, "advanced_estimator_mode", cfg.advanced_estimator_mode);

		load(pt, "tc", cfg.tc); // tip clipper:
		load(pt, "br", cfg.br); // bulge remover:
		load(pt, "ec", cfg.ec); // erroneous connections remover:
		load(pt, "cec", cfg.cec); // cheating erroneous connections remover:
		load(pt, "de", cfg.de); // distance estimator:
		load(pt, "ade", cfg.ade); // advanced distance estimator:
		load(pt, "rr", cfg.rr); // repeat resolver:
		load(pt, "need_consensus", cfg.need_consensus);
		load(pt, "uncorrected_reads", cfg.uncorrected_reads);
		load(pt, cfg.dataset_name, cfg.ds);
	}

} // debruijn_graph

typedef config_common::config<debruijn_graph::debruijn_config> cfg;

#endif




