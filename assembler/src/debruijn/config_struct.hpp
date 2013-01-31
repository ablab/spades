//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/*
 * config_struct.hpp
 *
 *  Created on: Aug 9, 2011
 *      Author: Alexey.Gurevich
 */

#ifndef CONFIG_STRUCT_HDIP_
#define CONFIG_STRUCT_HDIP_

#include "openmp_wrapper.h"

#include "config_common.hpp"
#include "path_extend/pe_config_struct.hpp"

#include <io/reader.hpp>
#include <io/ireader.hpp>
#include <io/easy_reader.hpp>
#include <io/converting_reader_wrapper.hpp>

#include <boost/bimap.hpp>

#include <sys/types.h>
#include <sys/stat.h>

namespace debruijn_graph {

enum working_stage {
	ws_construction,
	ws_simplification,
	ws_late_pair_info_count,
	ws_distance_estimation,
	ws_repeats_resolving,
	ws_n50_enlargement
};

enum simplification_mode {
	sm_normal,
	sm_cheating,
	sm_topology,
	sm_chimeric,
	sm_max_flow,
	sm_pair_info_aware
};

enum paired_metrics {
    pm_read_count,
    pm_product
};

enum estimation_mode {
    em_simple,
    em_weighted,
    em_extensive,
    em_smoothing
};

enum resolving_mode
{
    rm_none           ,
    rm_split          ,
    rm_path_extend    ,
    rm_combined       ,
    rm_split_scaff    ,
    rm_jump           ,
    rm_rectangles
};

enum info_printer_pos {
    ipp_default = 0,
    ipp_before_first_gap_closer,
    ipp_before_simplification,
    ipp_tip_clipping,
    ipp_bulge_removal,
    ipp_err_con_removal,
    ipp_before_final_err_con_removal,
    ipp_final_err_con_removal,
    ipp_final_tip_clipping,
    ipp_final_bulge_removal,
    ipp_removing_isolated_edges,
    ipp_final_simplified,
    ipp_before_repeat_resolution,

    ipp_total
};

namespace details {

inline const char* info_printer_pos_name(size_t pos) {
    const char* names[] = { "default", "before_first_gap_closer", "before_simplification", "tip_clipping",
        "bulge_removal", "err_con_removal", "before_final_err_con_removal",
        "final_err_con_removal", "final_tip_clipping",
        "final_bulge_removal", "removing_isolated_edges",
        "final_simplified", "before_repeat_resolution" };

    utils::check_array_size<ipp_total>(names);
    return names[pos];
}

} // namespace details

inline std::string MakeLaunchTimeDirName() {
    time_t rawtime;
    struct tm * timeinfo;
    char buffer[80];

    time(&rawtime);
    timeinfo = localtime(&rawtime);

    strftime(buffer, 80, "%m.%d_%H.%M.%S", timeinfo);
    return std::string(buffer);
}

// struct for debruijn project's configuration file
struct debruijn_config {
    typedef boost::bimap<string, working_stage> stage_name_id_mapping;
    typedef boost::bimap<string, simplification_mode> simpl_mode_id_mapping;
    typedef boost::bimap<string, estimation_mode> estimation_mode_id_mapping;
    typedef boost::bimap<string, paired_metrics> paired_metrics_id_mapping;
    typedef boost::bimap<string, resolving_mode> resolve_mode_id_mapping;

//  damn shit fix, it is to be removed! To determine is it started from run.sh or from spades.py
    bool run_mode;

	bool developer_mode;

	static const stage_name_id_mapping FillStageInfo() {
		stage_name_id_mapping::value_type info[] = {
      stage_name_id_mapping::value_type("construction", ws_construction),
      stage_name_id_mapping::value_type("simplification", ws_simplification),
      stage_name_id_mapping::value_type("late_pair_info_count", ws_late_pair_info_count),
      stage_name_id_mapping::value_type("distance_estimation", ws_distance_estimation),
      stage_name_id_mapping::value_type("repeats_resolving", ws_repeats_resolving),
      stage_name_id_mapping::value_type("n50_enlargement", ws_n50_enlargement)
    };

		return stage_name_id_mapping(info, utils::array_end(info));
	}

	static const simpl_mode_id_mapping FillSimplifModeInfo() {
		simpl_mode_id_mapping::value_type info[] = {
      simpl_mode_id_mapping::value_type("normal", sm_normal),
      simpl_mode_id_mapping::value_type("pair_info_aware", sm_pair_info_aware),
      simpl_mode_id_mapping::value_type("cheating", sm_cheating),
      simpl_mode_id_mapping::value_type("topology", sm_topology),
      simpl_mode_id_mapping::value_type("chimeric", sm_chimeric ),
      simpl_mode_id_mapping::value_type("max_flow", sm_max_flow)
    };

		return simpl_mode_id_mapping(info, utils::array_end(info));
	}

    static const estimation_mode_id_mapping FillEstimationModeInfo() {
        estimation_mode_id_mapping::value_type info[] = {
          estimation_mode_id_mapping::value_type("simple",      em_simple      ),
          estimation_mode_id_mapping::value_type("weighted",    em_weighted    ),
          estimation_mode_id_mapping::value_type("extensive",   em_extensive   ),
          estimation_mode_id_mapping::value_type("smoothing",    em_smoothing  ),
        };
        return estimation_mode_id_mapping(info, utils::array_end(info));
    }

    static const paired_metrics_id_mapping FillPairedMetricsInfo() {
        paired_metrics_id_mapping::value_type info[] = { 
          paired_metrics_id_mapping::value_type("read_count" ,      pm_read_count),
          paired_metrics_id_mapping::value_type("product"    ,      pm_product)
        };

        return paired_metrics_id_mapping(info, utils::array_end(info));
    }

	static const resolve_mode_id_mapping FillResolveModeInfo()
	{
		resolve_mode_id_mapping::value_type info [] = {
        resolve_mode_id_mapping::value_type("none"             , rm_none           ),
        resolve_mode_id_mapping::value_type("split"  		   , rm_split		   ),
        resolve_mode_id_mapping::value_type("path_extend"      , rm_path_extend    ),
        resolve_mode_id_mapping::value_type("combined"         , rm_combined       ),
        resolve_mode_id_mapping::value_type("split_scaff"      , rm_split_scaff    ),
        resolve_mode_id_mapping::value_type("jump"             , rm_jump           ),
        resolve_mode_id_mapping::value_type("rectangles"       , rm_rectangles     ),
    };

		return resolve_mode_id_mapping(info, utils::array_end(info));
	}


	static const stage_name_id_mapping& working_stages_info() {
		static stage_name_id_mapping working_stages_info = FillStageInfo();
		return working_stages_info;
	}

	static const simpl_mode_id_mapping& simpl_mode_info() {
		static simpl_mode_id_mapping simpl_mode_info = FillSimplifModeInfo();
		return simpl_mode_info;
	}

    static const estimation_mode_id_mapping& estimation_mode_info() {
        static estimation_mode_id_mapping est_mode_info = FillEstimationModeInfo();
        return est_mode_info;
    }

    static const paired_metrics_id_mapping& paired_metrics_info() {
        static paired_metrics_id_mapping paired_metrics_info = FillPairedMetricsInfo();
        return paired_metrics_info;
    }

	static const resolve_mode_id_mapping& resolve_mode_info() {
		static resolve_mode_id_mapping info = FillResolveModeInfo();
		return info;
	}

	static const std::string& simpl_mode_name(simplification_mode mode_id) {
		auto it = simpl_mode_info().right.find(mode_id);

		VERIFY_MSG(it != simpl_mode_info().right.end(),
				"No name for simplification mode id = " << mode_id);
		return it->second;
	}

	static simplification_mode simpl_mode_id(std::string name) {
		auto it = simpl_mode_info().left.find(name);
		VERIFY_MSG(it != simpl_mode_info().left.end(),
				"There is no simplification mode with name = " << name);

		return it->second;
	}

	static const std::string& working_stage_name(working_stage stage_id) {
		auto it = working_stages_info().right.find(stage_id);
		VERIFY_MSG(it != working_stages_info().right.end(),
				"No name for working stage id = " << stage_id);

		return it->second;
	}

	static working_stage working_stage_id(std::string name) {
		auto it = working_stages_info().left.find(name);
		VERIFY_MSG(it != working_stages_info().left.end(),
				"There is no working stage with name = " << name);

		return it->second;
	}

	static const std::string& paired_metrics_name(paired_metrics metrics_id) {
		auto it = paired_metrics_info().right.find(metrics_id);
		VERIFY_MSG(it != paired_metrics_info().right.end(),
				"No name for paired metrics id = " << metrics_id);

		return it->second;
	}

	static paired_metrics paired_metrics_id(std::string name) {
		auto it = paired_metrics_info().left.find(name);
		VERIFY_MSG(it != paired_metrics_info().left.end(),
				"There is no paired metrics with name = " << name);

		return it->second;
	}

    static const std::string& estimation_mode_name(estimation_mode est_id) {
        auto it = estimation_mode_info().right.find(est_id);
        VERIFY_MSG(it != estimation_mode_info().right.end(),
                "No name for estimation mode id = " << est_id);
        return it->second;
    }

    static estimation_mode estimation_mode_id(std::string name) { 
        auto it = estimation_mode_info().left.find(name);
        VERIFY_MSG(it != estimation_mode_info().left.end(),
                "There is no estimation mode with name = " << name);

        return it->second;
    }

	static const std::string& resolving_mode_name(resolving_mode mode_id) {
		auto it = resolve_mode_info().right.find(mode_id);
		VERIFY_MSG(it != resolve_mode_info().right.end(),
				"No name for resolving mode id = " << mode_id);

		return it->second;
	}

	static resolving_mode resolving_mode_id(std::string name) {
		auto it = resolve_mode_info().left.find(name);
		VERIFY_MSG(it != resolve_mode_info().left.end(),
				"There is no resolving mode with name = " << name);

		return it->second;
	}

	struct simplification {
		struct tip_clipper {
			string condition;
		};

		struct bulge_remover {
			double max_bulge_length_coefficient;
            size_t max_additive_length_coefficient;
			double max_coverage;
			double max_relative_coverage;
			double max_delta;
			double max_relative_delta;
		};

		struct erroneous_connections_remover {
			double	max_ec_length_coefficient;
            bool	estimate_max_coverage;
			double	max_coverage;
		};

		struct cheating_erroneous_connections_remover {
			double max_ec_length_coefficient;
			size_t sufficient_neighbour_length;
			double coverage_gap;
		};

		struct topology_based_ec_remover {
			size_t max_ec_length_coefficient;
			size_t uniqueness_length;
			size_t plausibility_length;
		};

		struct tr_based_ec_remover {
			size_t max_ec_length_coefficient;
			size_t uniqueness_length;
			double unreliable_coverage;
		};

		struct max_flow_ec_remover {
			double max_ec_length_coefficient;
			size_t uniqueness_length;
			size_t plausibility_length;
		};

		struct pair_info_ec_remover {
			double max_ec_length_coefficient;
			size_t min_neighbour_length;
		};

		struct isolated_edges_remover {
			size_t max_length;
			double max_coverage;
			size_t max_length_any_cov;
		};

		struct complex_bulge_remover {
			bool enabled;
			bool pics_enabled;
            std::string folder;
			double max_relative_length;
			size_t max_length_difference;
		};

		simplification_mode simpl_mode;
		tip_clipper tc;
		bulge_remover br;
		erroneous_connections_remover ec;
		cheating_erroneous_connections_remover cec;
		topology_based_ec_remover tec;
		tr_based_ec_remover trec;
		max_flow_ec_remover mfec;
		pair_info_ec_remover piec;
		isolated_edges_remover ier;
		complex_bulge_remover cbr;

		bool removal_checks_enabled;

		//typedef map<>
	};

	std::string uncorrected_reads;
	bool need_consensus;
	bool path_set_graph;
	double mismatch_ratio;
	simplification simp;

	struct repeat_resolver {
		bool symmetric_resolve;
		int mode;
		double inresolve_cutoff_proportion;
		int near_vertex;
		int max_distance;
		size_t max_repeat_length;
		bool kill_loops;
	};

	struct distance_estimator {
		double linkage_distance_coeff;
		double max_distance_coeff;
		double filter_threshold;
	};

	struct smoothing_distance_estimator {
		size_t threshold;
		double range_coeff;
		double delta_coeff;
		double percentage;
		size_t cutoff;
		size_t min_peak_points;
		double inv_density;
		double derivative_threshold;

	};

	struct dataset {
		vector<vector<std::string> > paired_reads;
		vector<std::string> single_reads;
		vector<vector<std::string> > original_paired_reads;
		vector<std::string> original_single_reads;
		boost::optional<std::string> jumping_first;
		boost::optional<std::string> jumping_second;
		boost::optional<size_t> jump_is;
		boost::optional<size_t> jump_rl;
		boost::optional<size_t> RL;
		boost::optional<double> IS;
		boost::optional<double> is_var;
		boost::optional<double> median;
		boost::optional<double> mad;
		map<int, size_t> hist;

		map<size_t, size_t> percentiles;
		optional<double> avg_coverage;
		bool single_cell;
		std::string reference_genome_filename;

		Sequence reference_genome;
	};

	struct position_handler {
		int max_single_gap;
		std::string contigs_for_threading;
		std::string contigs_to_analyze;
		bool late_threading;
		bool careful_labeling;

	};

	struct gap_closer {
		int     minimal_intersection;
		bool    before_simplify;
		bool    in_simplify;
		bool    after_simplify;
		double  weight_threshold;
	};

	struct info_printer {
        bool print_stats;
        bool write_components;
        string components_for_kmer;
        string components_for_genome_pos;
        bool write_components_along_genome;
        bool write_components_along_contigs;
        bool save_full_graph;
        bool write_error_loc;
        bool write_full_graph;
        bool write_full_nc_graph;
	};

	struct jump_cfg {
		bool load;
		double weight_threshold;
	};

	struct SAM_writer {
		bool produce_align_files;
		bool output_map_format;
		bool align_before_RR;
		bool align_after_RR;
		bool adjust_align;
		bool align_only_paired;
		bool output_broken_pairs;
		bool align_original_reads;
		optional<bool> print_quality;
	};

	struct graph_read_corr_cfg {
		bool enable;
		string output_dir;
		bool binary;
	};

	typedef map<info_printer_pos, info_printer> info_printers_t;

public:

	std::string dataset_file;
	std::string project_name;
	std::string input_dir;
	std::string output_base;
	std::string output_root;
	std::string output_dir;
	std::string output_suffix;
	std::string output_saves;
	std::string final_contigs_file;
	std::string log_filename;

	bool make_saves;
	bool output_pictures;
	bool output_nonfinal_contigs;
	bool compute_paths_number;

	bool use_additional_contigs;
	bool use_unipaths;
	std::string additional_contigs;

	std::string load_from;

	working_stage entry_point;

	bool paired_mode;
	bool additional_ec_removing;
	bool divide_clusters;;
	bool correct_mismatches;
	bool paired_info_statistics;
	bool paired_info_scaffolder;
	bool cut_bad_connections;
	bool componential_resolve;
	bool gap_closer_enable;
	bool SAM_writer_enable;

	//Convertion options
	size_t buffer_size;
	std::string temp_bin_reads_dir;
	std::string temp_bin_reads_path;
	std::string temp_bin_reads_info;
	std::string paired_read_prefix;
	std::string single_read_prefix;

	size_t K;

    bool use_multithreading;
	size_t max_threads;
	size_t max_memory;

//	size_t is_infinity;

    paired_metrics paired_metr;

    estimation_mode est_mode;

	resolving_mode rm;
	path_extend::pe_config::MainPEParamsT pe_params;

	distance_estimator de;
	smoothing_distance_estimator ade;
	repeat_resolver rr;
	bool use_scaffolder;
	bool mask_all;
	dataset ds;
	position_handler pos;
	gap_closer gc;
	jump_cfg jump;
	SAM_writer sw;
	graph_read_corr_cfg graph_read_corr;
	info_printers_t info_printers;
};

// specific load functions

inline void load(debruijn_config::simplification::tip_clipper& tc,
		boost::property_tree::ptree const& pt, bool complete) {
	using config_common::load;
	load(tc.condition					, pt, "condition"		 			);
}

inline void load(working_stage& entry_point,
		boost::property_tree::ptree const& pt, std::string const& key,
		bool complete) {
	std::string ep = pt.get<std::string>(key);
	entry_point = debruijn_config::working_stage_id(ep);
}

inline void load(resolving_mode& rm, boost::property_tree::ptree const& pt,
		std::string const& key, bool complete) {
	std::string ep = pt.get<std::string>(key);
	rm = debruijn_config::resolving_mode_id(ep);
}

inline void load(simplification_mode& simp_mode,
		boost::property_tree::ptree const& pt, std::string const& key,
		bool complete) {
	std::string ep = pt.get<std::string>(key);
	simp_mode = debruijn_config::simpl_mode_id(ep);
}

inline void load(estimation_mode& est_mode,
        boost::property_tree::ptree const& pt, std::string const& key,
        bool complete) {
    std::string ep = pt.get<std::string>(key);
    est_mode = debruijn_config::estimation_mode_id(ep);
}

inline void load(paired_metrics& paired_metr,
        boost::property_tree::ptree const& pt, std::string const& key,
        bool complete) {
    std::string ep = pt.get<std::string>(key);
    paired_metr = debruijn_config::paired_metrics_id(ep);
}

inline void load(debruijn_config::simplification::bulge_remover& br,
		boost::property_tree::ptree const& pt, bool complete) {
	using config_common::load;

	load(br.max_bulge_length_coefficient	, pt, "max_bulge_length_coefficient"	);
	load(br.max_additive_length_coefficient	, pt, "max_additive_length_coefficient"	);
	load(br.max_coverage         			, pt, "max_coverage"					);
	load(br.max_relative_coverage			, pt, "max_relative_coverage"			);
	load(br.max_delta		     			, pt, "max_delta"			    		);
	load(br.max_relative_delta   			, pt, "max_relative_delta"	    		);
}

inline void load(debruijn_config::simplification::pair_info_ec_remover& ec,
		boost::property_tree::ptree const& pt, bool complete) {
	using config_common::load;

	load(ec.max_ec_length_coefficient	, pt, "max_ec_length_coefficient"	);
	load(ec.min_neighbour_length		, pt, "min_neighbour_length"		);
}

inline void load(debruijn_config::simplification::isolated_edges_remover& ier,
		boost::property_tree::ptree const& pt, bool complete) {
	using config_common::load;

	load(ier.max_length					, pt, "max_length"					);
	load(ier.max_coverage				, pt, "max_coverage"				);
	load(ier.max_length_any_cov			, pt, "max_length_any_cov"			);
}

inline void load(debruijn_config::simplification::complex_bulge_remover& cbr,
		boost::property_tree::ptree const& pt, bool complete) {
	using config_common::load;

	load(cbr.enabled                    , pt, "enabled"                     );
	load(cbr.pics_enabled               , pt, "pics_enabled"                );
	load(cbr.folder                     , pt, "folder"                      );
	load(cbr.max_relative_length		, pt, "max_relative_length"			);
	load(cbr.max_length_difference		, pt, "max_length_difference"		);
}

inline void load(
		debruijn_config::simplification::erroneous_connections_remover& ec,
		boost::property_tree::ptree const& pt, bool complete) {
	using config_common::load;

	load(ec.max_ec_length_coefficient	, pt, "max_ec_length_coefficient"	);
	load(ec.estimate_max_coverage		, pt, "estimate_max_coverage"		);
	load(ec.max_coverage    			, pt, "max_coverage"    			);
}

inline void load(
		debruijn_config::simplification::cheating_erroneous_connections_remover& cec,
		boost::property_tree::ptree const& pt, bool complete) {
	using config_common::load;

	load(cec.max_ec_length_coefficient	, pt, "max_ec_length_coefficient"	);
	load(cec.coverage_gap               , pt, "coverage_gap"				);
	load(cec.sufficient_neighbour_length, pt, "sufficient_neighbour_length"	);
}

inline void load(
		debruijn_config::simplification::topology_based_ec_remover& tec,
		boost::property_tree::ptree const& pt, bool complete) {
	using config_common::load;

	load(tec.max_ec_length_coefficient	, pt, "max_ec_length_coefficient"	);
	load(tec.plausibility_length		, pt, "plausibility_length"			);
	load(tec.uniqueness_length			, pt, "uniqueness_length"			);
}

inline void load(debruijn_config::simplification::tr_based_ec_remover &trec,
		boost::property_tree::ptree const &pt, bool complaete) {
	using config_common::load;
	load(trec.max_ec_length_coefficient	, pt, "max_ec_length_coefficient"	);
	load(trec.unreliable_coverage		, pt, "unreliable_coverage"			);
	load(trec.uniqueness_length			, pt, "uniqueness_length"			);
}

inline void load(debruijn_config::simplification::max_flow_ec_remover& mfec,
		boost::property_tree::ptree const& pt, bool complete) {
	using config_common::load;

	load(mfec.max_ec_length_coefficient	, pt, "max_ec_length_coefficient"	);
	load(mfec.plausibility_length		, pt, "plausibility_length"			);
	load(mfec.uniqueness_length			, pt, "uniqueness_length"  			);
}

inline void load(debruijn_config::distance_estimator& de,
		boost::property_tree::ptree const& pt, bool complete) {
	using config_common::load;

	load(de.linkage_distance_coeff, pt, "linkage_distance_coeff");
	load(de.max_distance_coeff, pt, "max_distance_coeff");
	load(de.filter_threshold, pt, "filter_threshold");
}

inline void load(debruijn_config::smoothing_distance_estimator& ade,
		boost::property_tree::ptree const& pt, bool complete) {
	using config_common::load;

	load(ade.threshold, pt, "threshold");
	load(ade.range_coeff, pt, "range_coeff");
	load(ade.delta_coeff, pt, "delta_coeff");
	load(ade.percentage, pt, "percentage");
	load(ade.cutoff, pt, "cutoff");
	load(ade.min_peak_points, pt, "min_peak_points");
	load(ade.inv_density, pt, "inv_density");
	load(ade.derivative_threshold, pt, "derivative_threshold");
}

inline void load(debruijn_config::repeat_resolver& rr,
		boost::property_tree::ptree const& pt, bool complete) {
	using config_common::load;

	load(rr.symmetric_resolve, pt, "symmetric_resolve");
	load(rr.mode, pt, "mode");
	load(rr.inresolve_cutoff_proportion, pt, "inresolve_cutoff_proportion");
	load(rr.near_vertex, pt, "near_vertex");
	load(rr.max_distance, pt, "max_distance");
	load(rr.max_repeat_length, pt, "max_repeat_length");
	load(rr.kill_loops, pt, "kill_loops");
}

inline void load(debruijn_config::position_handler& pos,
		boost::property_tree::ptree const& pt, bool complete) {
	using config_common::load;
	load(pos.max_single_gap         , pt, "max_single_gap"		 );
	load(pos.contigs_for_threading  , pt, "contigs_for_threading");
	load(pos.contigs_to_analyze     , pt, "contigs_to_analyze"	 );
	load(pos.late_threading         , pt, "late_threading"		 );
	load(pos.careful_labeling       , pt, "careful_labeling"	 );
}

inline void load(debruijn_config::gap_closer& gc,
		boost::property_tree::ptree const& pt, bool complete) {
	using config_common::load;
	load(gc.minimal_intersection, pt, "minimal_intersection");
	load(gc.before_simplify     , pt, "before_simplify"     );
	load(gc.in_simplify         , pt, "in_simplify"         );
	load(gc.after_simplify      , pt, "after_simplify"      );
	load(gc.weight_threshold    , pt, "weight_threshold"    );
}

inline void load(debruijn_config::SAM_writer& sw,
		boost::property_tree::ptree const& pt, bool complete) {
	using config_common::load;
	load(sw.output_map_format   , pt, "output_map_format"   );
	load(sw.align_before_RR     , pt, "align_before_RR"     );
	load(sw.align_after_RR      , pt, "align_after_RR"      );
	load(sw.adjust_align        , pt, "adjust_align"        );
	load(sw.align_only_paired   , pt, "align_only_paired"   );
	load(sw.output_broken_pairs , pt, "output_broken_pairs" );
	load(sw.align_original_reads, pt, "align_original_reads");
	sw.print_quality = pt.get_optional<bool>("print_quality");
}

inline void load(debruijn_config::graph_read_corr_cfg& graph_read_corr,
		boost::property_tree::ptree const& pt, bool complete) {
	using config_common::load;
	load(graph_read_corr.enable    , pt, "enable"    );
	load(graph_read_corr.output_dir, pt, "output_dir");
	load(graph_read_corr.binary    , pt, "binary"    );
}

inline void load_paired_reads(vector<vector<std::string> >& vec, boost::property_tree::ptree const& pt, string const& key) {
	vector<std::string> strings;
	config_common::load(strings, pt, key);
	for (auto it = strings.begin(); it != strings.end(); ++it) {
		vector<std::string> paired_library;
		config_common::split(paired_library, *it);
		if (paired_library.size() < 1 || paired_library.size() > 2) {
			VERIFY_MSG(false, "Invalid library of " << paired_library.size() << " input files with paired reads: \"" << *it << "\"");
		}
		vec.push_back(paired_library);
	}
}

inline void load_single_reads(vector<std::string>& vec, boost::property_tree::ptree const& pt, string const& key) {
	vector<std::string> strings;
	config_common::load(strings, pt, key);
	for (auto it = strings.begin(); it != strings.end(); ++it) {
		config_common::split(vec, *it);
	}
}

inline void load(debruijn_config::dataset& ds,
		boost::property_tree::ptree const& pt, bool complete) {
	using config_common::load;

	load_paired_reads(ds.paired_reads, pt, "paired_reads");
	load_single_reads(ds.single_reads, pt, "single_reads");
	load_paired_reads(ds.original_paired_reads, pt, "original_paired_reads");
	load_single_reads(ds.original_single_reads, pt, "original_single_reads");
	load(ds.single_cell, pt, "single_cell");

	ds.jumping_first = pt.get_optional<std::string>("jumping_first");
	ds.jumping_second = pt.get_optional<std::string>("jumping_second");
	ds.jump_is = pt.get_optional<size_t>("jump_is");
	ds.jump_rl = pt.get_optional<size_t>("jump_rl");

	ds.RL = pt.get_optional<size_t>("RL");
	ds.is_var = pt.get_optional<size_t>("is_var");
	ds.IS = pt.get_optional<size_t>("IS");

	ds.reference_genome_filename = "";
	optional<std::string> refgen = pt.get_optional<std::string>(
			"reference_genome");
	if (refgen && *refgen != "N/A") {
		ds.reference_genome_filename = *refgen;
	}
}

inline void load_reference_genome(debruijn_config::dataset& ds,
		std::string input_dir) {
	if (ds.reference_genome_filename == "") {
		ds.reference_genome = Sequence();
		return;
	}
	if (ds.reference_genome_filename[0] != '/')
		ds.reference_genome_filename = input_dir + ds.reference_genome_filename;
	checkFileExistenceFATAL(ds.reference_genome_filename);
	io::Reader genome_stream(ds.reference_genome_filename);
	io::SingleRead genome;
	genome_stream >> genome;
	VERIFY(genome.IsValid());
//	if (VERIFY(genome.IsValid())) 
   // {
	ds.reference_genome = genome.sequence();
        //INFO("Reference genome loaded. Length " << ds.reference_genome.size());
//        cout << "Reference genome loaded. Length " << ds.reference_genome.size() << endl;
//	} 
//    else {
//		//INFO("Reference genome (" + ds.reference_genome_filename + ") has non-ACGT characters. Skipping it");
//		cout << "Reference genome (" + ds.reference_genome_filename + ") has non-ACGT characters. Skipping it" << endl;
//		ds.reference_genome = Sequence();
//		ds.reference_genome_filename = "";
//	}
}

inline void load(debruijn_config::simplification& simp,
		boost::property_tree::ptree const& pt, bool complete) {
	using config_common::load;

	load(simp.simpl_mode, pt, "simpl_mode");

	load(simp.tc, pt, "tc"); // tip clipper:
	load(simp.br, pt, "br"); // bulge remover:
	load(simp.ec, pt, "ec"); // erroneous connections remover:
	load(simp.cec, pt, "cec"); // cheating erroneous connections remover:
	load(simp.tec, pt, "tec"); // topology aware erroneous connections remover:
	load(simp.trec, pt, "trec"); // topology and reliability based erroneous connections remover:
	// need fix in config file
    load(simp.mfec, pt, "mfec"); // max flow erroneous connections remover:
	load(simp.piec, pt, "piec"); // pair info aware erroneous connections remover:
	load(simp.ier, pt, "ier"); // isolated edges remover
	load(simp.cbr, pt, "cbr"); // complex bulge remover

	load(simp.removal_checks_enabled, pt, "removal_checks_enabled");
}

inline void load(debruijn_config::info_printer& printer, boost::property_tree::ptree const& pt, bool complete)
{
    using config_common::load;

    load(printer.print_stats				  ,	pt, "print_stats"                   , complete);
    load(printer.write_components			  , pt, "write_components"              , complete);
    load(printer.components_for_kmer		  , pt, "components_for_kmer"           , complete);
    load(printer.components_for_genome_pos	  , pt, "components_for_genome_pos"     , complete);
    load(printer.write_components_along_genome,	pt, "write_components_along_genome" , complete);
    load(printer.write_components_along_contigs,pt, "write_components_along_contigs" , complete);
    load(printer.save_full_graph			  ,	pt, "save_full_graph"			  	, complete);
    load(printer.write_full_graph			  , pt, "write_full_graph"            , complete);
    load(printer.write_full_nc_graph		  , pt, "write_full_nc_graph"            , complete);
    load(printer.write_error_loc			  , pt, "write_error_loc"            , complete);
}

inline void load(debruijn_config::info_printers_t& printers,
		boost::property_tree::ptree const& pt, bool complete) {
	using config_common::load;
	using details::info_printer_pos_name;

	debruijn_config::info_printer def;
	load(def, pt, info_printer_pos_name(ipp_default), true);

	for (size_t pos = ipp_default + 1; pos != ipp_total; ++pos) {
		debruijn_config::info_printer printer(def);
		load(printer, pt, info_printer_pos_name(pos), false);

		printers[info_printer_pos(pos)] = printer;
	}
}

inline void load(debruijn_config::jump_cfg& jump,
		boost::property_tree::ptree const& pt, bool complete) {
	using config_common::load;

	load(jump.load, pt, "load");
	load(jump.weight_threshold, pt, "weight_threshold");
}

// main debruijn config load function
inline void load(debruijn_config& cfg, boost::property_tree::ptree const& pt, bool complete) {
	using config_common::load;

    load(cfg.K, pt, "K");

	// input options:
	load(cfg.dataset_file, pt, "dataset");
	// input dir is based on dataset file location (all pathes in datasets are relative to its location)
    cfg.input_dir = path::parent_path(cfg.dataset_file);
	if (cfg.input_dir[cfg.input_dir.length() - 1] != '/') {
		cfg.input_dir += '/';
	}

	load(cfg.output_base, pt, "output_base");
	if (cfg.output_base[cfg.output_base.length() - 1] != '/') {
		cfg.output_base += '/';
	}

	// TODO: remove this shit
	load(cfg.run_mode, pt, "run_mode");
    
	if (cfg.run_mode) {
	    load(cfg.project_name, pt, "project_name");
        cfg.output_root = cfg.project_name.empty()
			? (cfg.output_base + "/K" + ToString(cfg.K) + "/")
			: (cfg.output_base + cfg.project_name + "/K" + ToString(cfg.K) + "/");
        cfg.output_suffix = MakeLaunchTimeDirName() + "/";
        cfg.output_dir = cfg.output_root + cfg.output_suffix;
    } 
    else {
        cfg.output_root = cfg.output_base + "/K" + ToString(cfg.K) + "/";
        
        cfg.output_dir = cfg.output_root;
    }

	cfg.output_saves = cfg.output_dir + "saves/";

	load(cfg.log_filename, pt, "log_filename");

	load(cfg.developer_mode, pt, "developer_mode");
	if (cfg.developer_mode) {
		load(cfg.make_saves, pt, "make_saves");
		load(cfg.output_pictures, pt, "output_pictures");
		load(cfg.output_nonfinal_contigs, pt, "output_nonfinal_contigs");
		load(cfg.compute_paths_number, pt, "compute_paths_number");
	} else {
		cfg.make_saves = false;
		cfg.output_pictures = false;
		cfg.output_nonfinal_contigs = false;
		cfg.compute_paths_number = false;
	}
	if (!cfg.make_saves) {
		load(cfg.make_saves, pt, "force_make_saves");
	}

	load(cfg.load_from, pt, "load_from");
	if (cfg.load_from[0] != '/') { // relative path
		if (cfg.run_mode) 
            cfg.load_from = cfg.output_root + cfg.load_from;
        else 
            cfg.load_from = cfg.output_dir + cfg.load_from;
	}

	load(cfg.entry_point, pt, "entry_point");

	load(cfg.use_additional_contigs, pt, "use_additional_contigs");
	load(cfg.use_unipaths, pt, "use_unipaths");

	load(cfg.additional_contigs, pt, "additional_contigs");

	load(cfg.paired_mode, pt, "paired_mode");
	load(cfg.additional_ec_removing, pt, "additional_ec_removing");
	load(cfg.divide_clusters, pt, "divide_clusters");
	load(cfg.correct_mismatches, pt, "correct_mismatches");
	load(cfg.paired_info_statistics, pt, "paired_info_statistics");
	load(cfg.paired_info_scaffolder, pt, "paired_info_scaffolder");
	load(cfg.cut_bad_connections, pt, "cut_bad_connections");
	load(cfg.componential_resolve, pt, "componential_resolve");
	load(cfg.gap_closer_enable, pt, "gap_closer_enable");
	load(cfg.SAM_writer_enable, pt, "SAM_writer_enable");


    load(cfg.buffer_size, pt, "buffer_size");
    cfg.buffer_size <<= 20; //turn MB to bytes

	load(cfg.temp_bin_reads_dir, pt, "temp_bin_reads_dir");
    if (cfg.temp_bin_reads_dir[cfg.temp_bin_reads_dir.length() - 1] != '/') {
        cfg.temp_bin_reads_dir += '/';
    }
    cfg.temp_bin_reads_path = cfg.project_name.empty()
            ? (cfg.output_base + "/" + cfg.temp_bin_reads_dir)
            : (cfg.output_base + cfg.project_name + "/" + cfg.temp_bin_reads_dir);

    cfg.temp_bin_reads_info = cfg.temp_bin_reads_path + "INFO";

	cfg.paired_read_prefix = cfg.temp_bin_reads_path + "_paired";
	cfg.single_read_prefix =cfg.temp_bin_reads_path +  "_single";

	load(cfg.use_multithreading, pt, "use_multithreading");
	load(cfg.max_threads, pt, "max_threads");
  // Fix number of threads according to OMP capabilities.
  cfg.max_threads = std::min(cfg.max_threads, (size_t)omp_get_max_threads());
  // Inform OpenMP runtime about this :)
  omp_set_num_threads(cfg.max_threads);

  load(cfg.max_memory, pt, "max_memory");

	checkFileExistenceFATAL(cfg.dataset_file);
	boost::property_tree::ptree ds_pt;
	boost::property_tree::read_info(cfg.dataset_file, ds_pt);
	load(cfg.ds, ds_pt, true);

	load(cfg.de, pt, (cfg.ds.single_cell ? "sc_de" : "usual_de"));

	load(cfg.ade, pt, (cfg.ds.single_cell ? "sc_ade" : "usual_ade")); // advanced distance estimator:
	load(cfg.rr, pt, (cfg.ds.single_cell ? "sc_rr" : "usual_rr")); // repeat resolver:
	load(cfg.pos, pt, "pos"); // position handler:

    load(cfg.paired_metr, pt, "paired_metrics");

    load(cfg.est_mode, pt, "estimation_mode");

	load(cfg.rm               , pt, "resolving_mode"   );
    cfg.pe_params.name = cfg.ds.single_cell ? "singlecell" : "multicell";
    load(cfg.pe_params, pt, "andrey_params"    );
    load(cfg.use_scaffolder, pt, "use_scaffolder");
    load(cfg.mask_all, pt, "mask_all");


	load(cfg.gc, pt, "gap_closer");
	load(cfg.sw, pt, "SAM_writer");
	load(cfg.graph_read_corr, pt, "graph_read_corr");
	load(cfg.need_consensus, pt, "need_consensus");
	load(cfg.uncorrected_reads, pt, "uncorrected_reads");
	load(cfg.path_set_graph, pt, "path_set_graph");
	load(cfg.mismatch_ratio, pt, "mismatch_ratio");

	load(
			cfg.simp,
			pt,
			(cfg.ds.single_cell ? "sc_simplification" : "usual_simplification"));

    cfg.simp.cbr.folder = cfg.output_dir + cfg.simp.cbr.folder + "/";

	load(cfg.info_printers, pt, "info_printers");
	load(cfg.jump, pt, "jump");
    
	load_reference_genome(cfg.ds, cfg.input_dir);

//	cfg.is_infinity = 100000000;
}

} // debruijn_graph

typedef config_common::config<debruijn_graph::debruijn_config> cfg;

namespace debruijn_graph {

inline string input_file(string filename) {
	if (filename[0] == '/')
		return filename;
	return cfg::get().input_dir + filename;
}

}

#endif
