//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

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

namespace path_extend {

const char * const pe_cfg_filename = "./config/debruijn/path_extend/lc_config.info";

// struct for long_contigs subproject's configuration file
struct pe_config
{

	struct DatasetT
	{
		struct PairedLibT {
			size_t read_size;
			size_t insert_size;
			size_t var;

	        std::string path;
		};

        std::string param_set;
		std::string graph_file;

		std::vector<PairedLibT> libs;

		//boost::optional<std::string> reference_genome;
	};


	struct OutputParamsT {
	    bool write_seeds;
	    bool write_overlaped_paths;
	    bool write_paths;
	    bool write_path_loc;

	    void DisableAll() {
	        write_seeds = false;
	        write_overlaped_paths = false;
	        write_paths = false;
	        write_path_loc = false;
	    }
	};

    struct VisualizeParamsT {
        bool print_seeds;
        bool print_overlaped_paths;
        bool print_paths;

        void DisableAll() {
            print_seeds = false;
            print_overlaped_paths = false;
            print_paths = false;
        }
    };

    struct ResearchT {
        bool on;

        bool count_seed_weight;
        bool count_path_weight;

        bool fiter_seeds_by_id;
        std::vector<size_t> seed_ids;
    };

	struct ParamSetT {

	    std::string metric;
	    bool normalize_weight;
	    bool normalize_by_coverage;

	    bool improve_paired_info;

		struct SeedSelectionT
		{
	        std::string metric;

			double min_coverage;
			double start_egde_coverage;
			size_t max_cycles;

	        bool exclude_chimeric;
	        int chimeric_delta;

			bool check_trusted;
			double threshold;

		} seed_selection;


		struct ExtensionOptionsT
		{
		    std::string metric;

			bool try_deep_search;

			struct SelectOptionsT {
			    boost::optional<double> single_threshold;
			    double weight_threshold;
			    double priority_coeff;

			    SelectOptionsT() {

			    }

			    SelectOptionsT(const SelectOptionsT& so): single_threshold(so.single_threshold), weight_threshold(so.weight_threshold), priority_coeff(so.priority_coeff) {

			    }
			}
			select_options;

		} extension_options;


		struct ScaffolderOptionsT {
		    bool on;
	        int cutoff;
	        double rel_cutoff;
	        double sum_threshold;

	        bool cluster_info;
	        double cl_threshold;

	        bool fix_gaps;
	        double min_gap_score;
	        double max_must_overlap;
	        double max_can_overlap;
	        int short_overlap;
	        int artificial_gap;
		} scaffolder_options;


		struct LoopRemovalT
		{
			bool inspect_short_loops;

			size_t max_loops;
			bool full_loop_removal;

		} loop_removal;


		struct FilterOptionsT
		{
			bool remove_overlaps;
		} filter_options;
	};

	struct UtilsT {
		int mode;
		std::string file1;
		std::string file2;

		std::string clustered;
		std::string advanced;
		size_t insert_size;
		size_t read_size;
		size_t dev;
	};

    struct MainPEParamsT {
	    std::string name;

	    bool debug_output;

	    OutputParamsT output;
	    VisualizeParamsT viz;
	    ParamSetT param_set;
    } params;

	std::string dataset_name;
	DatasetT dataset;

};

inline void load(pe_config::OutputParamsT& o, boost::property_tree::ptree const& pt, bool complete) {
    using config_common::load;

    load(o.write_seeds,             pt, "write_seeds"           );
    load(o.write_overlaped_paths,   pt, "write_overlaped_paths" );
    load(o.write_paths,             pt, "write_paths"           );
    load(o.write_path_loc,          pt, "write_path_loc"        );
}

inline void load(pe_config::VisualizeParamsT& o, boost::property_tree::ptree const& pt, bool complete) {
    using config_common::load;
    load(o.print_seeds,             pt, "print_seeds"           );
    load(o.print_overlaped_paths,   pt, "print_overlaped_paths" );
    load(o.print_paths,             pt, "print_paths"           );
}

inline void load(pe_config::UtilsT& u, boost::property_tree::ptree const& pt, bool complete)
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


inline void load(pe_config::DatasetT::PairedLibT& pl, boost::property_tree::ptree const& pt, bool complete)
{
	using config_common::load;
	load(pl.read_size  , pt, "read_size"  );
	load(pl.insert_size, pt, "insert_size");
	load(pl.var        , pt, "var"        );
	load(pl.path       , pt, "path"       );
}

inline void load(pe_config::DatasetT& ds, boost::property_tree::ptree const& pt, bool complete)
{
	using config_common::load;

	//ds.reference_genom = pt.get_optional<std::string>("reference_genome");

	load(ds.graph_file, pt, "graph_file");

	load(ds.param_set,  pt, "param_set");
	//TODO
	//load(ds.libs,       pt,     "libs" );
}


inline void load(pe_config::ParamSetT::SeedSelectionT& ss, boost::property_tree::ptree const& pt, bool complete)
{
	using config_common::load;

	load(ss.min_coverage       , pt, "min_coverage"     );
	load(ss.start_egde_coverage, pt, "start_egde_coverage");
	load(ss.max_cycles         , pt, "max_cycles"       );
    load(ss.exclude_chimeric   , pt, "exclude_chimeric");
    load(ss.chimeric_delta     , pt, "chimeric_delta"    );

	load(ss.check_trusted      , pt, "check_trusted"    );
	load(ss.threshold          , pt, (ss.metric + "_trusted_threshold").c_str());
}

inline void load(pe_config::ParamSetT::ExtensionOptionsT::SelectOptionsT& so, boost::property_tree::ptree const& pt, bool complete) {
    using config_common::load;

    load(so.priority_coeff,     pt, "priority_coeff");
    load(so.weight_threshold,   pt, "weight_threshold");
    so.single_threshold = pt.get_optional<double>("single_threshold");
}

inline void load(pe_config::ParamSetT::ExtensionOptionsT& es, boost::property_tree::ptree const& pt, bool complete)
{
	using config_common::load;

	load(es.try_deep_search , pt, "try_deep_search");
	load(es.select_options  , pt, es.metric.c_str()          );
}

inline void load(pe_config::ParamSetT::LoopRemovalT& lr, boost::property_tree::ptree const& pt, bool complete)
{
	using config_common::load;
	load(lr.inspect_short_loops, pt,"inspect_short_loops"      );

	load(lr.max_loops          , pt,"max_loops"          );
	load(lr.full_loop_removal  , pt,"full_loop_removal"  );
}


inline void load(pe_config::ParamSetT::FilterOptionsT& fo, boost::property_tree::ptree const& pt, bool complete)
{
	using config_common::load;
	load(fo.remove_overlaps      , pt, "remove_overlaps"      );
}


inline void load(pe_config::ParamSetT::ScaffolderOptionsT& so, boost::property_tree::ptree const& pt, bool complete)
{
    using config_common::load;
    load(so.on      , pt, "on"      );
    load(so.cutoff      , pt, "cutoff"      );
    load(so.rel_cutoff      , pt, "rel_cutoff"      );
    load(so.sum_threshold      , pt, "sum_threshold"      );

    load(so.cluster_info      , pt, "cluster_info"      );
    load(so.cl_threshold      , pt, "cl_threshold"      );

    load(so.fix_gaps      , pt, "fix_gaps"      );
    load(so.min_gap_score      , pt, "min_gap_score"      );
    load(so.max_must_overlap      , pt, "max_must_overlap"      );
    load(so.max_can_overlap      , pt, "max_can_overlap"      );
    load(so.short_overlap      , pt, "short_overlap"      );
    load(so.artificial_gap      , pt, "artificial_gap"      );
}

inline void load(pe_config::ParamSetT& p, boost::property_tree::ptree const& pt, bool complete) {

    using config_common::load;

    load(p.metric, pt,  "metric");
    load(p.normalize_weight, pt,  "normalize_weight");
    load(p.normalize_by_coverage, pt,  "normalize_by_coverage");

    load(p.improve_paired_info, pt,  "improve_paired_info");


    p.seed_selection.metric = p.metric;
    p.extension_options.metric = p.metric;

    load(p.seed_selection,    pt, "seed_selection");
    load(p.extension_options, pt, "extension_options");
    load(p.scaffolder_options, pt, "scaffolder");
    load(p.loop_removal,      pt, "loop_removal");
    load(p.filter_options,    pt, "filter_options");
}

inline void load(pe_config::MainPEParamsT& p, boost::property_tree::ptree const& pt, bool complete) {
    using config_common::load;

    load(p.debug_output      , pt,  "debug_output"   );

    load(p.output      , pt,  "output"   );
    load(p.viz         , pt,  "visualize");
    load(p.param_set   , pt,  p.name.c_str()   );

    if (!p.debug_output) {
        p.output.DisableAll();
        p.viz.DisableAll();
    }
}


// main long contigs config load function
inline void load(pe_config& pe_cfg, boost::property_tree::ptree const& pt, bool complete)
{
	using config_common::load;

    load(pe_cfg.dataset_name           , pt, "dataset"               );
    load(pe_cfg.dataset                , pt,  pe_cfg.dataset_name.c_str()  );

    pe_cfg.params.name = pe_cfg.dataset.param_set;
    load(pe_cfg.params                 , pt, "pe_params"             );
}

}

typedef config_common::config<path_extend::pe_config> pe_cfg;

path_extend::pe_config::MainPEParamsT params;

#endif /* CONFIG_STRUCT_HPP_ */
