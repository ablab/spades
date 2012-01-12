/*
 * config_struct.hpp
 *
 *  Created on: Aug 9, 2011
 *      Author: Alexey.Gurevich
 */

#ifndef CONFIG_STRUCT_HDIP_
#define CONFIG_STRUCT_HDIP_

#include "config_common.hpp"
#include "long_contigs/lc_config_struct.hpp"
#include "k.hpp"
#include <sys/types.h>
#include <sys/stat.h>
#include <io/reader.hpp>

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
	sm_normal           ,
	sm_cheating         ,
	sm_topology         ,
	sm_chimeric         ,
	sm_pair_info_aware
};

enum resolving_mode
{
	rm_none           ,
	rm_dima           ,
	rm_andrew         ,
	rm_combined       ,
	rm_jump
};

enum info_printer_pos
{
    ipp_default                   = 0,
    ipp_before_simplification        ,
    ipp_tip_clipping                 ,
    ipp_bulge_removal                ,
    ipp_err_con_removal              ,
    ipp_before_final_err_con_removal ,
    ipp_final_err_con_removal        ,
    ipp_final_tip_clipping           ,
    ipp_final_bulge_removal          ,
    ipp_removing_isolated_edges      ,
    ipp_final_simplified             ,
    ipp_before_repeat_resolution     ,

    ipp_total
};

namespace details
{

inline const char* info_printer_pos_name(size_t pos)
{
    const char* names[] =
    {
       "default"                     ,
       "before_simplification"       ,
       "tip_clipping"                ,
       "bulge_removal"               ,
       "err_con_removal"             ,
       "before_final_err_con_removal",
       "final_err_con_removal"       ,
       "final_tip_clipping"          ,
       "final_bulge_removal"         ,
       "removing_isolated_edges"     ,
       "final_simplified"            ,
       "before_repeat_resolution"
    };

    utils::check_array_size<ipp_total>(names);
    return names[pos];
}

} // namespace details

const char* const cfg_filename = "./src/debruijn/config.info";

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
struct debruijn_config
{
	typedef bimap<string, working_stage      > stage_name_id_mapping;
	typedef bimap<string, simplification_mode> simpl_mode_id_mapping;
	typedef bimap<string, resolving_mode     > resolve_mode_id_mapping;

	static const stage_name_id_mapping FillStageInfo()
	{
	    stage_name_id_mapping::value_type info [] =
		{
		        {"construction"         , ws_construction           },
		        {"paired_info_count"    , ws_paired_info_count      },
		        {"simplification"       , ws_simplification         },
		        {"late_pair_info_count" , ws_late_pair_info_count   },
		        {"distance_estimation"  , ws_distance_estimation    },
		        {"repeats_resolving"    , ws_repeats_resolving      },
		        {"repeats_resolving"    , ws_repeats_resolving      },
		        {"n50_enlargement"      , ws_n50_enlargement        }
		};

		return stage_name_id_mapping(info, utils::array_end(info));
	}

	static const simpl_mode_id_mapping FillSimplifModeInfo()
	{
		simpl_mode_id_mapping::value_type info [] =
        {
                {"normal"           , sm_normal         },
                {"pair_info_aware"  , sm_pair_info_aware},
                {"cheating"         , sm_cheating       },
                {"topology"         , sm_topology       },
                {"chimeric"         , sm_chimeric       }
        };

		return simpl_mode_id_mapping(info, utils::array_end(info));
	}

	static const resolve_mode_id_mapping FillResolveModeInfo()
	{
		resolve_mode_id_mapping::value_type info [] =
        {
                {"none"             , rm_none           },
                {"dima"  			, rm_dima			},
                {"andrew"           , rm_andrew         },
                {"combined"         , rm_combined       },
                {"jump"             , rm_jump           },
        };

		return resolve_mode_id_mapping(info, utils::array_end(info));
	}

	static const simpl_mode_id_mapping& simpl_mode_info() {
		static simpl_mode_id_mapping simpl_mode_info = FillSimplifModeInfo();
		return simpl_mode_info;
	}

	static const stage_name_id_mapping& working_stages_info() {
		static stage_name_id_mapping working_stages_info = FillStageInfo();
		return working_stages_info;
	}

	static const resolve_mode_id_mapping& resolve_mode_info() {
		static resolve_mode_id_mapping info = FillResolveModeInfo();
		return info;
	}

	static const std::string& simpl_mode_name(simplification_mode mode_id) {
		auto it = simpl_mode_info().right.find(mode_id);

		VERIFY_MSG(it != simpl_mode_info().right.end(), "No name for simplification mode id = " << mode_id);
		return it->second;
	}

	static simplification_mode simpl_mode_id(std::string name) {
		auto it = simpl_mode_info().left.find(name);
		VERIFY_MSG(it != simpl_mode_info().left.end(), "There is no simplification mode with name = " << name);

		return it->second;
	}

	static const std::string& working_stage_name(working_stage stage_id) {
		auto it = working_stages_info().right.find(stage_id);
		VERIFY_MSG(it != working_stages_info().right.end(), "No name for working stage id = " << stage_id);

		return it->second;
	}

	static working_stage working_stage_id(std::string name) {
		auto it = working_stages_info().left.find(name);
		VERIFY_MSG(it != working_stages_info().left.end(), "There is no working stage with name = " << name);

		return it->second;
	}

	static const std::string& resolving_mode_name(resolving_mode mode_id) {
		auto it = resolve_mode_info().right.find(mode_id);
		VERIFY_MSG(it != resolve_mode_info().right.end(), "No name for resolving mode id = " << mode_id);

		return it->second;
	}

	static resolving_mode resolving_mode_id(std::string name) {
		auto it = resolve_mode_info().left.find(name);
		VERIFY_MSG(it != resolve_mode_info().left.end(), "There is no resolving mode with name = " << name);

		return it->second;
	}

	struct simplification {
		struct tip_clipper {
			size_t max_tip_length;
			double max_coverage;
			double max_relative_coverage;
		};

		struct bulge_remover {
			double max_length_div_K;
			double max_coverage;
			double max_relative_coverage;
			double max_delta;
			double max_relative_delta;
		};

		struct erroneous_connections_remover {
            optional<double> 	threshold_percentile;
			double 				max_coverage;
			int 				max_length_div_K;
		};

		struct cheating_erroneous_connections_remover {
			size_t max_length;
			double coverage_gap;
			size_t sufficient_neighbour_length;
		};

		struct topology_based_ec_remover {
			size_t max_length;
			size_t uniqueness_length;
			size_t plausibility_length;
		};

		struct pair_info_ec_remover {
			size_t max_length;
			size_t min_neighbour_length;
		};

		simplification_mode simpl_mode;
		tip_clipper tc;
		bulge_remover br;
		erroneous_connections_remover ec;
		cheating_erroneous_connections_remover cec;
		topology_based_ec_remover tec;
		pair_info_ec_remover piec;

		double isolated_min_len;
		bool   removal_checks_enabled;

		//typedef map<>
	};


		std::string uncorrected_reads;
		bool need_consensus;
		bool path_set_graph;
		simplification simp;

	struct repeat_resolver {
		bool symmetric_resolve;
		int mode;
		double inresolve_cutoff_proportion;
		int near_vertex;
	};
	struct distance_estimator {
		size_t delta;
		size_t linkage_distance;
		size_t max_distance;
		double filter_threshold;
	};
	struct advanced_distance_estimator {
		size_t threshold;
		double range_coeff;
		double delta_coeff;
		double percentage;
		size_t cutoff;
		size_t minpeakpoints;
		double inv_density;
		double derivative_threshold;

	};

	struct dataset {
		std::string first;
		std::string second;
		boost::optional<std::string> single_first;
		boost::optional<std::string> single_second;
		boost::optional<std::string> jumping_first;
		boost::optional<std::string> jumping_second;
		boost::optional<size_t> jump_is;
		boost::optional<size_t> jump_rl;
		size_t RL;
		size_t IS;
		bool single_cell;
		std::string reference_genome_filename;
		Sequence reference_genome;
	};

	struct position_handler {
		int max_single_gap;
		std::string contigs_for_threading;
		bool late_threading;
	};


	struct gap_closer {
		int minimal_intersection;
	};

	struct info_printer
	{
       bool     print_stats;
       bool     detailed_dot_write;
       bool     write_components;
       string   components_for_kmer;
       bool     write_components_along_genome;
       bool		save_full_graph;
	};

	typedef map<info_printer_pos, info_printer> info_printers_t;


public:

	std::string dataset_name;
	std::string input_dir;
    std::string output_base;
	std::string output_root;
	std::string output_dir;
	std::string output_suffix;
	std::string output_saves;

	bool use_single_reads;
	bool use_additional_contigs;
	bool etalon_graph_mode;
	std::string additional_contigs;

	std::string load_from;

	working_stage entry_point;

	bool paired_mode;
	bool paired_info_statistics;
//	bool rectangle_mode;
	bool etalon_info_mode;
	bool late_paired_info;
	bool advanced_estimator_mode;
	bool componential_resolve;

	resolving_mode rm;
	long_contigs::lc_config::lc_params andrey_params;

	distance_estimator          de;
	advanced_distance_estimator ade;
	repeat_resolver             rr;
	dataset                     ds;
	position_handler            pos;
	gap_closer                  gc;

	info_printers_t info_printers;
};

// specific load functions

inline void load(debruijn_config::simplification::tip_clipper& tc, boost::property_tree::ptree const& pt, bool complete)
{
	using config_common::load;
	load(tc.max_tip_length			, pt, "max_tip_length"		 );
	load(tc.max_coverage			, pt, "max_coverage"		 );
	load(tc.max_relative_coverage	, pt, "max_relative_coverage");
}

inline void load(working_stage& entry_point, boost::property_tree::ptree const& pt, std::string const& key, bool complete)
{
	std::string ep = pt.get<std::string>(key);
	entry_point = debruijn_config::working_stage_id(ep);
}

inline void load(resolving_mode& rm, boost::property_tree::ptree const& pt, std::string const& key, bool complete)
{
	std::string ep = pt.get<std::string>(key);
	rm = debruijn_config::resolving_mode_id(ep);
}

inline void load(simplification_mode& simp_mode, boost::property_tree::ptree const& pt, std::string const& key, bool complete)
{
	std::string ep = pt.get<std::string>(key);
	simp_mode = debruijn_config::simpl_mode_id(ep);
}

inline void load(debruijn_config::simplification::bulge_remover& br, boost::property_tree::ptree const& pt, bool complete)
{
	using config_common::load;

	load(br.max_length_div_K     , pt, "max_length_div_K"		);
	load(br.max_coverage         , pt, "max_coverage"			);
	load(br.max_relative_coverage, pt, "max_relative_coverage"  );
	load(br.max_delta		     , pt, "max_delta"			    );
	load(br.max_relative_delta   , pt, "max_relative_delta"	    );
}

inline void load(debruijn_config::simplification::pair_info_ec_remover& ec, boost::property_tree::ptree const& pt, bool complete)
{
	using config_common::load;

	load(ec.max_length,           pt, "max_length"			);
	load(ec.min_neighbour_length, pt, "min_neighbour_length");
}

inline void load(debruijn_config::simplification::erroneous_connections_remover& ec, boost::property_tree::ptree const& pt, bool complete)
{
	using config_common::load;

	ec.threshold_percentile = pt.get_optional<double>("threshold_percentile");

	load(ec.max_coverage    	, pt, "max_coverage"    );
	load(ec.max_length_div_K	, pt, "max_length_div_K");
}

inline void load(debruijn_config::simplification::cheating_erroneous_connections_remover& cec, boost::property_tree::ptree const& pt, bool complete)
{
	using config_common::load;

	load(cec.max_length                 , pt, "max_length"  );
	load(cec.coverage_gap               , pt, "coverage_gap");
	load(cec.sufficient_neighbour_length, pt, "sufficient_neighbour_length");
}

inline void load(debruijn_config::simplification::topology_based_ec_remover& tec, boost::property_tree::ptree const& pt, bool complete)
{
	using config_common::load;

	load(tec.max_length         , pt, "max_length"         );
	load(tec.plausibility_length, pt, "plausibility_length");
	load(tec.uniqueness_length  , pt, "uniqueness_length"  );
}

inline void load(debruijn_config::distance_estimator& de, boost::property_tree::ptree const& pt, bool complete)
{
	using config_common::load;

	load(de.linkage_distance, pt, "linkage_distance");
    load(de.delta           , pt, "delta"           );
	load(de.max_distance    , pt, "max_distance"	);
	load(de.filter_threshold, pt, "filter_threshold");
}

inline void load(debruijn_config::advanced_distance_estimator& ade, boost::property_tree::ptree const& pt, bool complete)
{
	using config_common::load;

	load(ade.threshold      , pt, "threshold"       );
	load(ade.range_coeff    , pt, "range_coeff"     );
	load(ade.delta_coeff    , pt, "delta_coeff"		);
	load(ade.percentage     , pt, "percentage"		);
	load(ade.cutoff         , pt, "cutoff"			);
	load(ade.minpeakpoints  , pt, "minpeakpoints"	);
	load(ade.inv_density    , pt, "inv_density"		);
	load(ade.derivative_threshold, pt, "derivative_threshold");
}

inline void load(debruijn_config::repeat_resolver& rr, boost::property_tree::ptree const& pt, bool complete)
{
	using config_common::load;

	load(rr.symmetric_resolve, pt, "symmetric_resolve");
	load(rr.mode             , pt, "mode"			  );
	load(rr.inresolve_cutoff_proportion, pt, "inresolve_cutoff_proportion");
	load(rr.near_vertex      , pt, "near_vertex"	  );
}

inline void load(debruijn_config::position_handler& pos, boost::property_tree::ptree const& pt, bool complete)
{
	using config_common::load;

	load(pos.max_single_gap         , pt, "max_single_gap"		 );
	load(pos.contigs_for_threading  , pt, "contigs_for_threading");
	load(pos.late_threading         , pt, "late_threading"		 );
}

inline void load(debruijn_config::gap_closer& gc, boost::property_tree::ptree const& pt, bool complete)
{
	using config_common::load;
	load(gc.minimal_intersection, pt, "minimal_intersection");
}

inline void load(debruijn_config::dataset& ds, boost::property_tree::ptree const& pt, bool complete)
{
	using config_common::load;

	load(ds.first , pt, "first" );
	load(ds.second, pt, "second");

	ds.single_first  = pt.get_optional<std::string>("single_first");
	ds.single_second = pt.get_optional<std::string>("single_second");

	ds.jumping_first  = pt.get_optional<std::string>("jumping_first");
	ds.jumping_second = pt.get_optional<std::string>("jumping_second");
	ds.jump_is = pt.get_optional<size_t>("jump_is");
	ds.jump_rl = pt.get_optional<size_t>("jump_rl");

	load(ds.RL, pt, "RL");
	load(ds.IS, pt, "IS");
	load(ds.single_cell, pt, "single_cell");

	ds.reference_genome_filename = "";
	boost::optional<std::string> refgen = pt.get_optional<std::string>("reference_genome");
	if (refgen && *refgen != "N/A") {
		ds.reference_genome_filename = *refgen;
	}
}

inline void load_reference_genome(debruijn_config::dataset& ds, std::string input_dir) {
	if (ds.reference_genome_filename == "") {
		ds.reference_genome = Sequence();
		return;
	}
	std::string genome_filename = input_dir + ds.reference_genome_filename;
	checkFileExistenceFATAL(genome_filename);
	io::Reader<io::SingleRead> genome_stream(genome_filename);
	io::SingleRead genome;
	genome_stream >> genome;
	ds.reference_genome = genome.sequence();
//		std::string genome;
//		genome = full_genome.GetSequenceString().substr(0, cfg::get().ds.LEN); // cropped
//		return Sequence(genome);
}

inline void load(debruijn_config::simplification& simp, boost::property_tree::ptree const& pt, bool complete)
{
	using config_common::load;

	load(simp.simpl_mode, pt, "simpl_mode");

	load(simp.tc  , pt, "tc"  ); // tip clipper:
	load(simp.br  , pt, "br"  ); // bulge remover:
	load(simp.ec  , pt, "ec"  ); // erroneous connections remover:
	load(simp.cec , pt, "cec" ); // cheating erroneous connections remover:
	load(simp.tec , pt, "tec" ); // topology aware erroneous connections remover:
	load(simp.piec, pt, "piec"); // pair info aware erroneous connections remover:

	load(simp.isolated_min_len      , pt, "isolated_min_len"      );
	load(simp.removal_checks_enabled, pt, "removal_checks_enabled");
}

inline void load(debruijn_config::info_printer& printer, boost::property_tree::ptree const& pt, bool complete)
{
    using config_common::load;

    load(printer.print_stats				  ,	pt, "print_stats"                   , complete);
    load(printer.detailed_dot_write			  , pt, "detailed_dot_write"            , complete);
    load(printer.write_components			  , pt, "write_components"              , complete);
    load(printer.components_for_kmer		  , pt, "components_for_kmer"           , complete);
    load(printer.write_components_along_genome,	pt, "write_components_along_genome" , complete);
    load(printer.save_full_graph			  ,	pt, "save_full_graph"			  	, complete);
}

inline void load(debruijn_config::info_printers_t& printers, boost::property_tree::ptree const& pt, bool complete)
{
    using config_common::load;
    using details::info_printer_pos_name;

    debruijn_config::info_printer def;
    load(def, pt, info_printer_pos_name(ipp_default), true);

    for (size_t pos = ipp_default + 1; pos != ipp_total; ++pos)
    {
        debruijn_config::info_printer printer(def);
        load(printer, pt, info_printer_pos_name(pos), false);

        printers[info_printer_pos(pos)] = printer;
    }
}

// main debruijn config load function
inline void load(debruijn_config& cfg, boost::property_tree::ptree const& pt, bool complete)
{
	using config_common::load;
	// input options:
	load(cfg.dataset_name, pt, "dataset"  );
	load(cfg.input_dir   , pt, "input_dir");

	load(cfg.output_base, pt, "output_base");

	cfg.output_root     = cfg.output_base + cfg.dataset_name + "/K" + ToString(K) + "/";
	cfg.output_suffix   = MakeLaunchTimeDirName() + "/";
	cfg.output_dir      = cfg.output_root + cfg.output_suffix;
	cfg.output_saves    = cfg.output_dir + "saves/";

	load(cfg.load_from, pt, "load_from");
	cfg.load_from = cfg.output_root + cfg.load_from;

	load(cfg.entry_point            , pt, "entry_point"            );

	load(cfg.etalon_graph_mode      , pt, "etalon_graph_mode"      );
	load(cfg.use_single_reads       , pt, "use_single_reads"       );
	load(cfg.use_additional_contigs , pt, "use_additional_contigs" );

	load(cfg.additional_contigs     , pt, "additional_contigs"     );

	load(cfg.paired_mode            , pt, "paired_mode"            );
	load(cfg.paired_info_statistics , pt, "paired_info_statistics" );
	load(cfg.etalon_info_mode       , pt, "etalon_info_mode"       );
	load(cfg.late_paired_info       , pt, "late_paired_info"       );
	load(cfg.componential_resolve   , pt, "componential_resolve"   );
	load(cfg.advanced_estimator_mode, pt, "advanced_estimator_mode");
	load(cfg.ds                     , pt, cfg.dataset_name         );

	load(cfg.de, pt, (cfg.ds.single_cell ? "sc_de" : "usual_de"));


	load(cfg.ade              , pt, "ade"              ); // advanced distance estimator:
	load(cfg.rr               , pt, (cfg.ds.single_cell ? "sc_rr" : "usual_rr")               ); // repeat resolver:
	load(cfg.pos              , pt, "pos"              ); // position handler:

	load(cfg.rm               , pt, "resolving_mode"   );
	if (cfg.rm == rm_andrew || cfg.rm == rm_combined || cfg.rm == rm_jump) {
	    cfg.andrey_params.param_set_name = cfg.ds.single_cell ? "singlecell" : "multicell";
	    load(cfg.andrey_params, pt, "andrey_params"    );
	}

	load(cfg.gc               , pt, "gap_closer"       );
	load(cfg.need_consensus   , pt, "need_consensus"   );
	load(cfg.uncorrected_reads, pt, "uncorrected_reads");
	load(cfg.path_set_graph, pt, "path_set_graph");

	load(cfg.simp, pt, (cfg.ds.single_cell ? "sc_simplification" : "usual_simplification"));
	load(cfg.info_printers, pt, "info_printers");

	load_reference_genome(cfg.ds, cfg.input_dir);
}


} // debruijn_graph

typedef config_common::config<debruijn_graph::debruijn_config> cfg;

#endif
