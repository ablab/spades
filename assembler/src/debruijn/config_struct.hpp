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


namespace debruijn
{
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


		enum working_stage {
			start	,
			after_construction		,
			after_pair_info_counting	,
			after_simplification		,
			after_distance_estimation	,
			after_repeat_resolving	,
			after_consensus
		};

		typedef boost::bimap<std::string, working_stage> name_id_mapping;

		static const name_id_mapping FillStageInfo() {
			name_id_mapping working_stages_info;
			working_stages_info.insert(name_id_mapping::value_type("start", start));
			working_stages_info.insert(name_id_mapping::value_type("after_construction", after_construction));
			working_stages_info.insert(name_id_mapping::value_type("after_pair_info_counting"	, after_pair_info_counting));
			working_stages_info.insert(name_id_mapping::value_type("after_simplification"		, after_simplification));
			working_stages_info.insert(name_id_mapping::value_type("after_distance_estimation"	, after_distance_estimation));
			working_stages_info.insert(name_id_mapping::value_type("after_repeat_resolving"		, after_repeat_resolving));
			working_stages_info.insert(name_id_mapping::value_type("after_consensus"			, after_consensus));
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
			size_t RL;
			size_t IS;
			int LEN;
		};

		std::string input_dir;
		std::string output_root;
		std::string output_dir;
		std::string output_dir_suffix;

		std::string previous_run_dir;
		std::string dataset_name;
		std::string reference_genome;
		std::string start_from;


//		working_stage entry_point;
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

	inline void load(boost::property_tree::ptree const& pt, debruijn_config::tip_clipper& tc)
	{
		using config_common::load;
		load(pt, "max_tip_length", tc.max_tip_length);
		load(pt, "max_coverage", tc.max_coverage);
		load(pt, "max_relative_coverage", tc.max_relative_coverage);
	}

	inline void load(boost::property_tree::ptree const& pt, std::string const& key, debruijn_config::working_stage& entry_point)
	{
		std::string ep = pt.get<std::string>(key);
	//	std::map<std::string, debruijn_config::working_stage> stages =
	//	{
	//			{"construction"			, debruijn_config::construction},
	//			{"pair_info_counting"	, debruijn_config::pair_info_counting},
	//			{"simplification"		, debruijn_config::simplification},
	//			{"distance_estimation"	, debruijn_config::distance_estimation},
	//			{"repeat_resolving"		, debruijn_config::repeat_resolving},
	//			{"consensus"			, debruijn_config::consensus}
	//	};
	//
	//	auto it = stages.find(ep);
	//	assert(it != stages.end());
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

	inline void load(boost::property_tree::ptree const& pt, debruijn_config::distance_estimator& de)
	{
		using config_common::load;
		load(pt, "delta", de.delta);
		load(pt, "linkage_distance", de.linkage_distance);
		load(pt, "max_distance", de.max_distance);
	}

	inline void load(boost::property_tree::ptree const& pt, debruijn_config::dataset& ds)
	{
		using config_common::load;
		load(pt, "first", ds.first);
		load(pt, "second", ds.second);
		load(pt, "IS", ds.IS);
		load(pt, "LEN", ds.LEN);
	}

	// main debruijn config load function
	inline void load(boost::property_tree::ptree const& pt, debruijn_config& cfg)
	{
		using config_common::load;
		// input options:
//		temporarily disabled
//		load(pt, "entry_point", cfg.entry_point);
		load(pt, "input_dir", cfg.input_dir);
	//	= cfg::get().output_dir
		load(pt, "previous_run_dir", cfg.previous_run_dir);
		load(pt, "dataset", cfg.dataset_name);

		load(pt, "output_dir", cfg.output_root);
		cfg.output_dir_suffix = MakeLaunchTimeDirName() + "." + cfg.dataset_name + "/";
		cfg.output_dir = cfg.output_root + cfg.output_dir_suffix;

		load(pt, "reference_genome", cfg.reference_genome);
		load(pt, "start_from", cfg.start_from);

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

}

typedef config_common::config<debruijn::debruijn_config> cfg;

#endif




