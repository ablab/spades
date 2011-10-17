/*
 * config_struct.hpp
 *
 *  Created on: Aug 9, 2011
 *      Author: Alexey.Gurevich
 */

#ifndef CONFIG_STRUCT_HPP_
#define CONFIG_STRUCT_HPP_

#include "config_common.hpp"
#include "k.hpp"
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
        ws_late_pair_info_count,
        ws_distance_estimation,
        ws_repeats_resolving,
        ws_n50_enlargement
    };

    enum simplification_mode
    {
    	sm_normal,
    	sm_cheating,
    	sm_chimeric,
    	sm_pair_info_aware
    };

	const char* const cfg_filename = "./src/debruijn/config.info";
//	const size_t K = 55; // must be odd (so there is no k-mer which is equal to it's reverse-complimentary k-mer)

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
		typedef boost::bimap<std::string, working_stage> stage_name_id_mapping;
		typedef boost::bimap<std::string, simplification_mode> simpl_mode_id_mapping;

		static const stage_name_id_mapping FillStageInfo() {
			stage_name_id_mapping working_stages_info;
			working_stages_info.insert(stage_name_id_mapping::value_type("construction"       , ws_construction       ));
			working_stages_info.insert(stage_name_id_mapping::value_type("paired_info_count"  , ws_paired_info_count  ));
			working_stages_info.insert(stage_name_id_mapping::value_type("simplification"	    , ws_simplification     ));
			working_stages_info.insert(stage_name_id_mapping::value_type("late_pair_info_count", ws_late_pair_info_count));
			working_stages_info.insert(stage_name_id_mapping::value_type("distance_estimation", ws_distance_estimation));
			working_stages_info.insert(stage_name_id_mapping::value_type("repeats_resolving"	, ws_repeats_resolving  ));
			working_stages_info.insert(stage_name_id_mapping::value_type("n50_enlargement"	, ws_n50_enlargement    ));
			return working_stages_info;
		}

		static const simpl_mode_id_mapping FillSimplifModeInfo() {
			simpl_mode_id_mapping simpl_mode_info;
			simpl_mode_info.insert(simpl_mode_id_mapping::value_type("normal", sm_normal));
			simpl_mode_info.insert(simpl_mode_id_mapping::value_type("pair_info_aware", sm_pair_info_aware));
			simpl_mode_info.insert(simpl_mode_id_mapping::value_type("cheating", sm_cheating));
			simpl_mode_info.insert(simpl_mode_id_mapping::value_type("chimeric", sm_chimeric));
			return simpl_mode_info;
		}

		static const simpl_mode_id_mapping& simpl_mode_info() {
			static simpl_mode_id_mapping simpl_mode_info = FillSimplifModeInfo();
			return simpl_mode_info;
		}

		static const stage_name_id_mapping& working_stages_info() {
			static stage_name_id_mapping working_stages_info = FillStageInfo();
			return working_stages_info;
		}

		static const std::string& simpl_mode_name(simplification_mode mode_id) {
			simpl_mode_id_mapping::right_const_iterator it = simpl_mode_info().right.find(mode_id);
			VERIFY_MSG(it != simpl_mode_info().right.end(), "No name for simplification mode id = " << mode_id);
			return it->second;
		}

		static simplification_mode simpl_mode_id(std::string name) {
			simpl_mode_id_mapping::left_const_iterator it = simpl_mode_info().left.find(name);
			VERIFY_MSG(it != simpl_mode_info().left.end(), "There is no simplification mode with name = " << name);
			return it->second;
		}

		static const std::string& working_stage_name(working_stage stage_id) {
			stage_name_id_mapping::right_const_iterator it = working_stages_info().right.find(stage_id);
			VERIFY_MSG(it != working_stages_info().right.end(), "No name for working stage id = " << stage_id);
			return it->second;
		}

		static working_stage working_stage_id(std::string name) {
			stage_name_id_mapping::left_const_iterator it = working_stages_info().left.find(name);
			VERIFY_MSG(it != working_stages_info().left.end(), "There is no working stage with name = " << name);
			return it->second;
		}

		struct simplification
		{
			struct tip_clipper
			{
			   double max_tip_length_div_K;
			   size_t max_tip_length;
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
				double  max_coverage;
				int     max_length_div_K;
			};

			struct cheating_erroneous_connections_remover
			{
				size_t max_length;
				double coverage_gap;
				size_t sufficient_neighbour_length;
			};

			struct pair_info_ec_remover
			{
				size_t max_length;
				size_t min_neighbour_length;
			};

			tip_clipper                            tc;
			bulge_remover                          br;
			erroneous_connections_remover          ec;
			cheating_erroneous_connections_remover cec;
			pair_info_ec_remover                   piec;
		};

		struct repeat_resolver
		{
			bool symmetric_resolve;
			int mode;
			int near_vertex;
		};
		struct distance_estimator
		{
			size_t delta;
			size_t linkage_distance;
			size_t max_distance;
			double filter_threshold;
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
			boost::optional<std::string> single_first;
			boost::optional<std::string> single_second;
			size_t RL;
			size_t IS;
			bool single_cell;
			std::string reference_genome;
			int LEN;
		};

        struct position_handler{
			int max_single_gap;
			std::string contigs_for_threading;
		};

        std::string dataset_name;
        std::string input_dir;
		std::string output_root;
		std::string output_dir;
		std::string output_suffix;
		std::string output_saves;

		bool use_single_reads;
		bool use_additional_contigs;
		std::string additional_contigs;

		std::string load_from;

		working_stage entry_point;
		simplification_mode simpl_mode;

		bool paired_mode;
		bool paired_info_statistics;
//		bool rectangle_mode;
		bool etalon_info_mode;
		bool late_paired_info;
		bool advanced_estimator_mode;
		bool componential_resolve;


		std::string uncorrected_reads;
		bool need_consensus;

		simplification simp;

		distance_estimator de;
		advanced_distance_estimator ade;
		repeat_resolver rr;
		dataset ds;
		position_handler pos;
	};


	// specific load functions

	inline void load(boost::property_tree::ptree const& pt, debruijn_config::simplification::tip_clipper& tc)
	{
		using config_common::load;
		load(pt, "max_tip_length_div_K" , tc.max_tip_length_div_K);
		load(pt, "max_tip_length" , tc.max_tip_length);
		load(pt, "max_coverage"		    , tc.max_coverage);
		load(pt, "max_relative_coverage", tc.max_relative_coverage);
	}

	inline void load(boost::property_tree::ptree const& pt, std::string const& key, working_stage& entry_point)
	{
		std::string ep = pt.get<std::string>(key);
		entry_point = debruijn_config::working_stage_id(ep);
	}

	inline void load(boost::property_tree::ptree const& pt, std::string const& key, simplification_mode& simp_mode)
	{
		std::string ep = pt.get<std::string>(key);
		simp_mode = debruijn_config::simpl_mode_id(ep);
	}

	inline void load(boost::property_tree::ptree const& pt, debruijn_config::simplification::bulge_remover& br)
	{
		using config_common::load;
		load(pt, "max_length_div_K", br.max_length_div_K);
		load(pt, "max_coverage", br.max_coverage);
		load(pt, "max_relative_coverage", br.max_relative_coverage);
		load(pt, "max_delta", br.max_delta);
		load(pt, "max_relative_delta", br.max_relative_delta);
	}

	inline void load(
	    boost::property_tree::ptree const& pt,
	    debruijn_config::simplification::pair_info_ec_remover& ec)
	{
		using config_common::load;
		load(pt, "max_length", ec.max_length);
		load(pt, "min_neighbour_length", ec.min_neighbour_length);
	}

	inline void load(
	    boost::property_tree::ptree const& pt,
	    debruijn_config::simplification::erroneous_connections_remover& ec)
	{
		using config_common::load;
		load(pt, "max_coverage", ec.max_coverage);
		load(pt, "max_length_div_K", ec.max_length_div_K);
	}

	inline void load(
	    boost::property_tree::ptree const& pt,
	    debruijn_config::simplification::cheating_erroneous_connections_remover& cec)
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
		load(pt, "filter_threshold", de.filter_threshold);
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
		load(pt, "symmetric_resolve", rr.symmetric_resolve);
		load(pt, "mode", rr.mode);
		load(pt, "near_vertex", rr.near_vertex);
	}

	inline void load(boost::property_tree::ptree const& pt, debruijn_config::position_handler& pos)
	{
		using config_common::load;
		load(pt, "max_single_gap", pos.max_single_gap);
		load(pt, "contigs_for_threading", pos.contigs_for_threading);
	}

	inline void load(boost::property_tree::ptree const& pt, debruijn_config::dataset& ds)
	{
		using config_common::load;
		load(pt, "first", ds.first);
		load(pt, "second", ds.second);
		ds.single_first = pt.get_optional<std::string>("single_first");
		ds.single_second = pt.get_optional<std::string>("single_second");
		load(pt, "RL", ds.RL);
		load(pt, "IS", ds.IS);
		load(pt, "single_cell", ds.single_cell);
		load(pt, "reference_genome", ds.reference_genome);
		load(pt, "LEN", ds.LEN);
	}

	inline void load(boost::property_tree::ptree const& pt, debruijn_config::simplification& simp)
	{
	    using config_common::load;

        load(pt, "tc" , simp.tc ); // tip clipper:
        load(pt, "br" , simp.br ); // bulge remover:
        load(pt, "ec" , simp.ec ); // erroneous connections remover:
        load(pt, "cec", simp.cec); // cheating erroneous connections remover:
        load(pt, "piec", simp.piec); // pair info aware erroneous connections remover:
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

		load(pt, "simpl_mode", cfg.simpl_mode);

		load(pt, "use_additional_contigs", cfg.use_additional_contigs);
		load(pt, "use_single_reads", cfg.use_single_reads);

		load(pt, "additional_contigs", cfg.additional_contigs);

		//load(pt, "reference_genome", cfg.reference_genome); moved to dataset
		//load(pt, "start_from", cfg.start_from);

		load(pt, "paired_mode", cfg.paired_mode);
		load(pt, "paired_info_statistics", cfg.paired_info_statistics);
//		load(pt, "rectangle_mode", cfg.rectangle_mode);
		load(pt, "etalon_info_mode", cfg.etalon_info_mode);
		load(pt, "late_paired_info", cfg.late_paired_info);
		load(pt, "componential_resolve", cfg.componential_resolve);
		load(pt, "advanced_estimator_mode", cfg.advanced_estimator_mode);
		load(pt, cfg.dataset_name, cfg.ds);
		load(pt, cfg.ds.single_cell ? "sc_de" : "usual_de", cfg.de);

		load(pt, "ade", cfg.ade); // advanced distance estimator:
		load(pt, "rr", cfg.rr); // repeat resolver:
		load(pt, "pos", cfg.pos); // position handler:
		load(pt, "need_consensus", cfg.need_consensus);
		load(pt, "uncorrected_reads", cfg.uncorrected_reads);

		load(pt, cfg.ds.single_cell ? "sc_simplification" : "usual_simplification", cfg.simp);
	}

} // debruijn_graph

typedef config_common::config<debruijn_graph::debruijn_config> cfg;

#endif




