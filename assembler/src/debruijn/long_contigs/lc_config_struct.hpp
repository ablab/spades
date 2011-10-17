/*
 * lc_config_struct.hpp
 *
 *  Created on: Aug 16, 2011
 *      Author: Alexey.Gurevich
 */

#ifndef LC_CONFIG_STRUCT_HPP_
#define LC_CONFIG_STRUCT_HPP_

#include "config_common.hpp"
#include "../config_struct.hpp"
#include <vector>

namespace long_contigs {

const char* const lc_cfg_filename = "./src/debruijn/long_contigs/lc_config.info";

// struct for long_contigs subproject's configuration file
struct lc_config
{
	struct dataset
	{
		std::string graph_file;
	};

	struct rl_dataset
	{
		bool precounted;
	    std::string precounted_path;
		std::string first;
		std::string second;
		bool has_advanced;
		std::string advanced;
		std::string raw;
	};

	struct real_lib
	{
		size_t read_size;
		size_t insert_size;
		size_t var;
		size_t is_delta;

		rl_dataset ds;
	};

	struct etalon_lib
	{
		size_t read_size;
		size_t insert_size;
	};


	struct symmetrization {
		bool cut_tips;
		size_t min_conjugate_len;
	};

	struct seed_selection
	{
	    double min_coverage;
	    bool   glue_seeds;
	    bool   check_trusted;
	    bool   remove_untrusted;
	    size_t max_cycles;
	    double trusted_threshold;

	    symmetrization sym;
	};

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
	};

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
	};

	struct stop_criteria
	{
		bool all_seeds;
		double edge_coverage;
		double len_coverage;
	};


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

	    symmetrization sym;
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
	//size_t real_libs_count;
	//size_t etalon_libs_count;

	bool use_new_metrics;

	bool write_seeds;
	bool write_overlaped_paths;
	bool write_paths;
	bool write_contigs;
	bool write_real_paired_info;
	bool write_raw_paired_info;
	bool cluster_paired_info;
	bool write_graph;
	bool print_stats;
	std::string paired_info_file_prefix;

	bool total_symmetric_mode;
	bool first_grow_forward;

	dataset ds;

	std::vector<real_lib> real_libs;
	std::vector<etalon_lib> etalon_libs;

	seed_selection ss;
	extension_selection es;
	loops_removal lr;
	stop_criteria sc;
	filter_options fo;
	research rs;
	utils u;
	symmetrization sym;
};

void load(boost::property_tree::ptree const& pt, lc_config::utils& u)
{
	using config_common::load;
	load(pt, "mode", u.mode);
	load(pt, "file1", u.file1);
	load(pt, "file2", u.file2);

	load(pt, "advanced", u.advanced);
	load(pt, "clustered", u.clustered);
	load(pt, "insert_size", u.insert_size);
	load(pt, "read_size", u.read_size);
	load(pt, "dev", u.dev);
}


// specific load functions
void load(boost::property_tree::ptree const& pt, lc_config::research& rs)
{
	using config_common::load;
	load(pt, "research_mode", rs.research_mode);
	load(pt, "detailed_output", rs.detailed_output);
	load(pt, "fiter_by_edge", rs.fiter_by_edge);
	load(pt, "force_to_cycle", rs.force_to_cycle);
	load(pt, "edge_length", rs.edge_length);
	load(pt, "cycle_priority_edge", rs.cycle_priority_edge);
}

void load(boost::property_tree::ptree const& pt, lc_config::rl_dataset& ds)
{
	using config_common::load;
	load(pt, "precounted", ds.precounted);
	load(pt, "precounted_path", ds.precounted_path);
	load(pt, "first", ds.first);
	load(pt, "second", ds.second);
	load(pt, "has_advanced", ds.has_advanced);
	load(pt, "advanced", ds.advanced);
	load(pt, "raw", ds.raw);
}

void load(boost::property_tree::ptree const& pt, lc_config::real_lib& rl)
{
	using config_common::load;
	load(pt, "read_size", rl.read_size);
	load(pt, "insert_size", rl.insert_size);
	load(pt, "var", rl.var);
	load(pt, "is_delta", rl.is_delta);
	load(pt, cfg::get().dataset_name, rl.ds);
}

void load(boost::property_tree::ptree const& pt, lc_config::etalon_lib& el)
{
	using config_common::load;
	load(pt, "insert_size", el.insert_size);
	load(pt, "read_size", el.read_size);
}

void load(boost::property_tree::ptree const& pt, lc_config::dataset& ds)
{
	using config_common::load;
	load(pt, "graph_file", ds.graph_file);
}

void load(boost::property_tree::ptree const& pt, lc_config::seed_selection& ss)
{
	using config_common::load;
	load(pt, "min_coverage", ss.min_coverage);
	load(pt, "glue_seeds", ss.glue_seeds);
	load(pt, "max_cycles", ss.max_cycles);
	load(pt, "check_trusted", ss.check_trusted);
    load(pt, "remove_untrusted", ss.remove_untrusted);
    load(pt, "trusted_threshold", ss.trusted_threshold);
    load(pt, "sym", ss.sym);
}

void load(boost::property_tree::ptree const& pt, lc_config::extension_selection& es)
{
	using config_common::load;
	load(pt, "use_weight_function_first", es.use_weight_function_first);
	load(pt, "weight_fun_threshold", es.weight_fun_threshold);
	load(pt, "weight_threshold", es.weight_threshold);
	load(pt, "use_delta_first", es.use_delta_first);
	load(pt, "etalon_distance_dev", es.etalon_distance_dev);
	load(pt, "max_iter", es.max_iter);
	load(pt, "priority_coeff", es.priority_coeff);
	load(pt, "max_depth", es.max_depth);
	load(pt, "use_advanced", es.use_advanced);
	load(pt, "fix_weight", es.fix_weight);
	load(pt, "advanced_coeff", es.advanced_coeff);
}

void load(boost::property_tree::ptree const& pt, lc_config::loops_removal& lr)
{
	using config_common::load;
	load(pt, "investigation", lr.investigation);
	load(pt, "loop_to_investigate", lr.loop_to_investigate);
	load(pt, "max_exits", lr.max_exits);
	load(pt, "max_loop_len", lr.max_loop_len);

	load(pt, "max_loops", lr.max_loops);
	load(pt, "full_loop_removal", lr.full_loop_removal);
	load(pt, "exlude_cycle", lr.exlude_cycle);

	load(pt, "stop_on_long", lr.stop_on_long);
}

void load(boost::property_tree::ptree const& pt, lc_config::stop_criteria& sc)
{
	using config_common::load;
	load(pt, "all_seeds", sc.all_seeds);
	load(pt, "edge_coverage", sc.edge_coverage);
	load(pt, "len_coverage", sc.len_coverage);
}

void load(boost::property_tree::ptree const& pt, lc_config::filter_options& fo)
{
	using config_common::load;
	load(pt, "remove_overlaps", fo.remove_overlaps);
	load(pt, "remove_subpaths", fo.remove_subpaths);
	load(pt, "remove_duplicates", fo.remove_duplicates);
	load(pt, "conjugate_percent", fo.conjugate_percent);
	load(pt, "length_percent", fo.length_percent);
	load(pt, "remove_single", fo.remove_single);

	load(pt, "remove_similar", fo.remove_similar);
    load(pt, "similar_edges", fo.similar_edges);
    load(pt, "similar_length", fo.similar_length);

    load(pt, "remove_sefl_conjugate", fo.remove_sefl_conjugate);
    load(pt, "conj_len_percent", fo.conj_len_percent);
    load(pt, "break_sc", fo.break_sc);
    load(pt, "chimeric_delta", fo.chimeric_delta);

    load(pt, "agreed_coeff", fo.agreed_coeff);

    load(pt, "sym", fo.sym);
}

void load(boost::property_tree::ptree const& pt, lc_config::symmetrization& sym)
{
	using config_common::load;
	load(pt, "min_conjugate_len", sym.min_conjugate_len);
	load(pt, "cut_tips", sym.cut_tips);

}


// main long contigs config load function
void load(boost::property_tree::ptree const& pt, lc_config& lc_cfg)
{
	using config_common::load;
	load(pt, "from_file", lc_cfg.from_file);
	load(pt, "syminfo", lc_cfg.syminfo);
	load(pt, "paired_info_only", lc_cfg.paired_info_only);
	//load(pt, "real_libs_count", cfg.real_libs_count);
	//load(pt, "etalon_libs_count", cfg.etalon_libs_count);

	load(pt, "use_new_metrics", lc_cfg.use_new_metrics);
	load(pt, "write_seeds", lc_cfg.write_seeds);
	load(pt, "write_overlaped_paths", lc_cfg.write_overlaped_paths);
	load(pt, "write_paths", lc_cfg.write_paths);
	load(pt, "write_contigs", lc_cfg.write_contigs);
	load(pt, "write_real_paired_info", lc_cfg.write_real_paired_info);
	load(pt, "write_raw_paired_info", lc_cfg.write_raw_paired_info);
	load(pt, "write_graph", lc_cfg.write_graph);
	load(pt, "print_stats", lc_cfg.print_stats);
	load(pt, "cluster_paired_info", lc_cfg.cluster_paired_info);
	load(pt, "paired_info_file_prefix", lc_cfg.paired_info_file_prefix);

	load(pt, "total_symmetric_mode", lc_cfg.total_symmetric_mode);
	load(pt, "first_grow_forward", lc_cfg.first_grow_forward);

	load(pt, cfg::get().dataset_name, lc_cfg.ds);
	load(pt, "real_libs", lc_cfg.real_libs);
	load(pt, "etalon_libs", lc_cfg.etalon_libs);

	load(pt, "ss", lc_cfg.ss);
	load(pt, "es", lc_cfg.es);
	load(pt, "lr", lc_cfg.lr);
	load(pt, "sc", lc_cfg.sc);
	load(pt, "fo", lc_cfg.fo);
	load(pt, "research", lc_cfg.rs);
	load(pt, "utils", lc_cfg.u);

}

} // namespace lc

typedef config_common::config<long_contigs::lc_config> lc_cfg;

#endif /* CONFIG_STRUCT_HPP_ */
