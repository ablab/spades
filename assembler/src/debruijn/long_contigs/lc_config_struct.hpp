/*
 * lc_config_struct.hpp
 *
 *  Created on: Aug 16, 2011
 *      Author: Alexey.Gurevich
 */

#ifndef LC_CONFIG_STRUCT_HPP_
#define LC_CONFIG_STRUCT_HPP_

#include "config_common.hpp"
#include <vector>

namespace long_contigs {

const char* const lc_cfg_filename = "./src/debruijn/long_contigs/lc_config.info";

// struct for long_contigs subproject's configuration file
struct lc_config
{

	struct dataset
	{
		struct paired_lib_t {
			size_t read_size;
			size_t insert_size;

			size_t de_delta;
			size_t is_delta;
			size_t var;

	        bool precounted;
	        std::string precounted_path;

	        std::string first;
	        std::string second;

	        bool has_advanced;
	        std::string advanced;

			bool has_raw;
			std::string raw;
		};

		struct etalon_paired_lib_t {
			size_t read_size;
			size_t insert_size;

			size_t de_delta;
			size_t is_delta;
			size_t var;
		};

        std::string reference_genome;
        int LEN;

		std::string graph_file;

	    std::string param_set;

//	    std::vector<paired_lib_t> paired_lib;
//	    std::vector<etalon_paired_lib_t> etalon_paired_lib;
	};



	struct param_set {
		struct symmetrization {
			bool cut_tips;
			size_t min_conjugate_len;
		};

		struct seed_selection
		{
			double min_coverage;
			double start_egde_coverage;
			bool   glue_seeds;
			bool   check_trusted;
			bool   remove_untrusted;
			size_t max_cycles;
			double trusted_threshold;

			size_t short_single;
			size_t chimeric_len;
			size_t chimeric_delta;

			symmetrization sym;
		} ss;

		struct extension_selection
		{
			bool   use_weight_function_first;
			double weight_fun_threshold;
			double weight_threshold;

			bool use_delta_first;
			int  etalon_distance_dev;
			int max_iter;
			int max_depth;

			double priority_coeff;

			bool fix_weight;
			bool use_advanced;
			double advanced_coeff;
		} es;

		struct loops_removal
		{
			bool investigation;
			size_t loop_to_investigate;
			size_t max_exits;
			size_t max_loop_len;

			size_t max_loops;
			bool full_loop_removal;

			bool stop_on_long;

			bool exlude_cycle;
		} lr;

		struct stop_criteria
		{
			bool all_seeds;
			double edge_coverage;
			double len_coverage;
		} sc;

		struct filter_options
		{
			bool remove_duplicates;
			bool remove_subpaths;
			bool remove_overlaps;

			double length_percent;
			double conjugate_percent;

			bool remove_single;
			bool remove_similar;
			double similar_edges;
			double similar_length;

			bool remove_sefl_conjugate;
			double conj_len_percent;
			bool break_sc;
			size_t chimeric_delta;

			double agreed_coeff;

			bool write_uncovered_edges;

			symmetrization sym;
		} fo;

	};

	struct research {
		bool research_mode;
		bool detailed_output;

		bool fiter_by_edge;
		size_t edge_length;

		bool force_to_cycle;
		size_t cycle_priority_edge;
	};

	struct utils {
		int mode;
		std::string file1;
		std::string file2;

		std::string clustered;
		std::string advanced;
		size_t insert_size;
		size_t read_size;
		size_t dev;

	};

	bool from_file;
	bool syminfo;
	bool paired_info_only;
	bool paired_info_for_resolved;
	std::string resolved_graph;
	bool cluster_paired_info;
	bool etalon_mode;

    bool write_real_paired_info;
    std::string paired_info_file_prefix;

	bool use_new_metrics;
	std::string dataset_name;
    utils u;

    dataset ds;

    struct lc_params {
        bool write_seeds;
        bool write_overlaped_paths;
        bool write_paths;
        bool write_path_loc;
        bool write_contigs;
        bool write_raw_paired_info;

        bool write_graph;
        bool print_stats;


        bool total_symmetric_mode;
        bool first_grow_forward;

        std::string param_set_name;

        param_set ps;
        research rs;
    } params;
};

inline void load(lc_config::utils& u, boost::property_tree::ptree const& pt, bool complete)
{
	using config_common::load;
	load(u.mode, pt, "mode");
	load(u.file1, pt, "file1");
	load(u.file2, pt, "file2");

	load(u.advanced, pt, "advanced");
	load(u.clustered, pt, "clustered");
	load(u.insert_size, pt, "insert_size");
	load(u.read_size, pt, "read_size");
	load(u.dev, pt, "dev");
}


// specific load functions
inline void load(lc_config::research& rs, boost::property_tree::ptree const& pt, bool complete)
{
	using config_common::load;
	load(rs.research_mode  , pt, "research_mode"  );
	load(rs.detailed_output, pt, "detailed_output");
	load(rs.fiter_by_edge  , pt, "fiter_by_edge"  );
	load(rs.force_to_cycle , pt, "force_to_cycle" );
	load(rs.edge_length    , pt, "edge_length"	  );
	load(rs.cycle_priority_edge, pt, "cycle_priority_edge");
}

inline void load(lc_config::dataset::paired_lib_t& rl, boost::property_tree::ptree const& pt, bool complete)
{
	using config_common::load;
	load(rl.read_size      , pt, "read_size"    );
	load(rl.insert_size    , pt, "insert_size"  );

	load(rl.var            , pt, "var"          );
	load(rl.de_delta       , pt, "de_delta"     );
	load(rl.is_delta       , pt, "is_delta"   );

	load(rl.precounted     , pt, "precounted"   );
	load(rl.precounted_path, pt, "precounted_path");

	load(rl.first          , pt, "first"        );
	load(rl.second         , pt, "second"       );

	load(rl.has_advanced   , pt, "has_advanced" );
	load(rl.advanced       , pt, "advanced"     );

	load(rl.has_raw            , pt, "has_raw"      );
	load(rl.raw        , pt, "raw"          );
}

inline void load(lc_config::dataset::etalon_paired_lib_t& el, boost::property_tree::ptree const& pt, bool complete)
{
	using config_common::load;
	load(el.read_size  , pt, "read_size"  );
	load(el.insert_size, pt, "insert_size");

	load(el.var        , pt, "var"        );
	load(el.de_delta   , pt, "de_delta"   );
	load(el.is_delta   , pt, "is_delta" );
}

inline void load(lc_config::dataset& ds, boost::property_tree::ptree const& pt, bool complete)
{
	using config_common::load;
    boost::optional<std::string> rg = pt.get_optional<std::string>("reference_genome");
    ds.reference_genome = (rg) ? *rg : "";

    if (ds.reference_genome == "N/A")
        ds.reference_genome = "";

    load(ds.LEN, pt, "length");

	load(ds.graph_file, pt, "graph_file");
	load(ds.param_set, pt, "param_set");
//	load(ds.paired_lib, pt, "real_libs");
//	load(ds.etalon_paired_lib, pt, "etalon_libs");
}


inline void load(lc_config::param_set::seed_selection& ss, boost::property_tree::ptree const& pt, bool complete)
{
	using config_common::load;

	load(ss.min_coverage       , pt, "min_coverage"     );
	load(ss.start_egde_coverage, pt, "start_egde_coverage");
	load(ss.glue_seeds         , pt, "glue_seeds"       );
	load(ss.max_cycles         , pt, "max_cycles"       );
	load(ss.check_trusted      , pt, "check_trusted"    );
    load(ss.remove_untrusted   , pt, "remove_untrusted" );
    load(ss.trusted_threshold  , pt, "trusted_threshold");
    load(ss.sym                , pt, "sym"              );
    load(ss.short_single       , pt, "short_single"     );
    load(ss.chimeric_delta     , pt, "chimeric_delta"   );
    load(ss.chimeric_len       , pt, "chimeric_len"     );
}

inline void load(lc_config::param_set::extension_selection& es, boost::property_tree::ptree const& pt, bool complete)
{
	using config_common::load;

	load(es.use_weight_function_first , pt, "use_weight_function_first");
	load(es.weight_fun_threshold      , pt, "weight_fun_threshold"     );
	load(es.weight_threshold          , pt, "weight_threshold"         );
	load(es.use_delta_first           , pt, "use_delta_first"          );
	load(es.etalon_distance_dev       , pt, "etalon_distance_dev"      );
	load(es.max_iter                  , pt, "max_iter"                 );
	load(es.priority_coeff            , pt, "priority_coeff"           );
	load(es.max_depth                 , pt, "max_depth"                );
	load(es.use_advanced              , pt, "use_advanced"             );
	load(es.fix_weight                , pt, "fix_weight"               );
	load(es.advanced_coeff            , pt, "advanced_coeff"           );
}

inline void load(lc_config::param_set::loops_removal& lr, boost::property_tree::ptree const& pt, bool complete)
{
	using config_common::load;
	load(lr.investigation      , pt,"investigation"      );
	load(lr.loop_to_investigate, pt,"loop_to_investigate");
	load(lr.max_exits          , pt,"max_exits"          );
	load(lr.max_loop_len       , pt,"max_loop_len"       );

	load(lr.max_loops          , pt,"max_loops"          );
	load(lr.full_loop_removal  , pt,"full_loop_removal"  );
	load(lr.exlude_cycle       , pt,"exlude_cycle"       );

	load(lr.stop_on_long       , pt,"stop_on_long"       );
}

inline void load(lc_config::param_set::stop_criteria& sc, boost::property_tree::ptree const& pt, bool complete)
{
	using config_common::load;
	load(sc.all_seeds    ,  pt, "all_seeds"    );
	load(sc.edge_coverage,  pt, "edge_coverage");
	load(sc.len_coverage ,  pt, "len_coverage" );
}

inline void load(lc_config::param_set::filter_options& fo, boost::property_tree::ptree const& pt, bool complete)
{
	using config_common::load;
	load(fo.remove_overlaps      , pt, "remove_overlaps"      );
	load(fo.remove_subpaths      , pt, "remove_subpaths"      );
	load(fo.remove_duplicates    , pt, "remove_duplicates"    );
	load(fo.conjugate_percent    , pt, "conjugate_percent"    );
	load(fo.length_percent       , pt, "length_percent"       );
	load(fo.remove_single        , pt, "remove_single"        );

	load(fo.remove_similar       , pt, "remove_similar"       );
    load(fo.similar_edges        , pt, "similar_edges"        );
    load(fo.similar_length       , pt, "similar_length"       );

    load(fo.remove_sefl_conjugate, pt, "remove_sefl_conjugate");
    load(fo.conj_len_percent     , pt, "conj_len_percent"     );
    load(fo.break_sc             , pt, "break_sc"             );
    load(fo.chimeric_delta       , pt, "chimeric_delta"       );

    load(fo.agreed_coeff         , pt, "agreed_coeff"         );

    load(fo.sym                  , pt, "sym"                  );

    load(fo.write_uncovered_edges, pt, "write_uncovered_edges");
}

inline void load(lc_config::param_set::symmetrization& sym, boost::property_tree::ptree const& pt, bool complete)
{
	using config_common::load;
	load(sym.min_conjugate_len, pt, "min_conjugate_len");
	load(sym.cut_tips         , pt, "cut_tips"         );

}

inline void load(lc_config::param_set& p, boost::property_tree::ptree const& pt, bool complete) {

    using config_common::load;

    load(p.ss, pt, "ss");
    load(p.es, pt, "es");
    load(p.fo, pt, "fo");
    load(p.lr, pt, "lr");
    load(p.sc, pt, "sc");
}

inline void load(lc_config::lc_params& p, boost::property_tree::ptree const& pt, bool complete) {
    using config_common::load;
    load(p.write_seeds            , pt, "write_seeds"           );
    load(p.write_overlaped_paths  , pt, "write_overlaped_paths" );
    load(p.write_paths            , pt, "write_paths"           );
    load(p.write_path_loc         , pt, "write_path_loc"        );
    load(p.write_contigs          , pt, "write_contigs"         );
    load(p.write_raw_paired_info  , pt, "write_raw_paired_info" );
    load(p.write_graph            , pt, "write_graph"           );
    load(p.print_stats            , pt, "print_stats"           );


    load(p.total_symmetric_mode   , pt, "total_symmetric_mode"  );
    load(p.first_grow_forward     , pt, "first_grow_forward"    );

    load(p.rs                     , pt, "research"              );

    load(p.ps                     , pt,  p.param_set_name       );
}


// main long contigs config load function
inline void load(lc_config& lc_cfg, boost::property_tree::ptree const& pt, bool complete)
{
	using config_common::load;
	load(lc_cfg.from_file              , pt, "from_file"             );
	load(lc_cfg.syminfo                , pt, "syminfo"               );
	load(lc_cfg.paired_info_only       , pt, "paired_info_only"      );
	load(lc_cfg.paired_info_for_resolved, pt, "paired_info_for_resolved");
	load(lc_cfg.resolved_graph         , pt, "resolved_graph"       );

	load(lc_cfg.cluster_paired_info    , pt, "cluster_paired_info"   );
	load(lc_cfg.use_new_metrics        , pt, "use_new_metrics"       );
    load(lc_cfg.write_real_paired_info , pt, "write_real_paired_info");
    load(lc_cfg.paired_info_file_prefix, pt, "paired_info_file_prefix");
    load(lc_cfg.etalon_mode            , pt, "etalon_mode"           );

    load(lc_cfg.dataset_name           , pt, "dataset"               );
    load(lc_cfg.u                      , pt, "utils"                 );

    load(lc_cfg.ds                     , pt,  lc_cfg.dataset_name    );

    lc_cfg.params.param_set_name = lc_cfg.ds.param_set;
    load(lc_cfg.params                 , pt, "lc_params"             );
}

} // namespace lc

typedef config_common::config<long_contigs::lc_config> lc_cfg;

long_contigs::lc_config::lc_params params;

#endif /* CONFIG_STRUCT_HPP_ */
