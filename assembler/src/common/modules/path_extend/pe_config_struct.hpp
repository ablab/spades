//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

/*
 * lc_config_struct.hpp
 *
 *  Created on: Aug 16, 2011
 *      Author: Alexey.Gurevich
 */

#ifndef LC_CONFIG_STRUCT_HPP_
#define LC_CONFIG_STRUCT_HPP_

#include "pipeline/config_singl.hpp"
#include "utils/cpp_utils.hpp"

#include <boost/optional.hpp>
#include <boost/property_tree/ptree_fwd.hpp>
#include <boost/bimap.hpp>

#include <string>
#include <vector>

namespace path_extend {

enum scaffolding_mode {
    sm_old,
    sm_2015,
    sm_combined,
    sm_old_pe_2015
};

//Both this functions return always true, right?
//still necessary?
inline bool IsScaffolder2015Enabled(const scaffolding_mode mode) {
    return (mode == sm_old_pe_2015 || mode == sm_2015 || mode == sm_combined);
}

inline bool IsOldPEEnabled(const scaffolding_mode mode) {
    return (mode == sm_old_pe_2015 || mode == sm_old || mode == sm_combined);
}

// struct for path extend subproject's configuration file
struct pe_config {
    typedef boost::bimap<std::string, scaffolding_mode> scaffolding_mode_id_mapping;

    static const scaffolding_mode_id_mapping FillSMInfo() {
        scaffolding_mode_id_mapping::value_type info[] = {
            scaffolding_mode_id_mapping::value_type("old", sm_old),
            scaffolding_mode_id_mapping::value_type("2015", sm_2015),
            scaffolding_mode_id_mapping::value_type("combined", sm_combined),
            scaffolding_mode_id_mapping::value_type("old_pe_2015", sm_old_pe_2015)
        };

        return scaffolding_mode_id_mapping(info, utils::array_end(info));
    }

    static const scaffolding_mode_id_mapping &scaffolding_mode_info() {
        static scaffolding_mode_id_mapping scaffolding_mode_info = FillSMInfo();
        return scaffolding_mode_info;
    }

    static const std::string &scaffolding_mode_name(scaffolding_mode sm) {
        auto it = scaffolding_mode_info().right.find(sm);
        VERIFY_MSG(it != scaffolding_mode_info().right.end(),
                   "No name for scaffolding mode id = " << sm);

        return it->second;
    }

    static scaffolding_mode scaffolding_mode_id(std::string name) {
        auto it = scaffolding_mode_info().left.find(name);
        VERIFY_MSG(it != scaffolding_mode_info().left.end(),
                   "There is no scaffolding mode with name = " << name);

        return it->second;
    }

    struct OutputParamsT {
        bool write_overlaped_paths;
        bool write_paths;

        void DisableAll() {
            write_overlaped_paths = false;
            write_paths = false;
        }
    };


    struct VisualizeParamsT {
        bool print_overlaped_paths;
        bool print_paths;

        void DisableAll() {
            print_overlaped_paths = false;
            print_paths = false;
        }
    };

    struct ParamSetT {
        scaffolding_mode sm;

        bool normalize_weight;
        size_t split_edge_length;

        bool multi_path_extend;

        struct OverlapRemovalOptionsT {
            bool enabled;
            bool end_start_only;
            bool cut_all;
        } overlap_removal;

        struct SimpleCoverageResolver {
            bool enabled;
            double coverage_delta;
            double min_upper_coverage;
        } simple_coverage_resolver;

        struct ExtensionOptionsT {
            double single_threshold;
            double weight_threshold;
            double priority_coeff;
            size_t max_repeat_length;
        } extension_options;

        ExtensionOptionsT mate_pair_options;


        struct ScaffolderOptionsT {
            bool enabled;
            int cutoff;
            int hard_cutoff;
            double rel_cutoff;
            double sum_threshold;

            bool cluster_info;
            double cl_threshold;

            bool use_la_gap_joiner;
            double min_gap_score;
            double max_can_overlap;
            int short_overlap;
            size_t artificial_gap;

            double var_coeff;
            double basic_overlap_coeff;

            size_t min_overlap_length;
            double flank_multiplication_coefficient;
            int flank_addition_coefficient;

            boost::optional<int> min_overlap_for_rna_scaffolding;
        } scaffolder_options;

        struct PathFiltrationT {
            bool enabled;
            size_t min_length;
            double rel_cutoff;
            size_t isolated_min_length;
            double isolated_min_cov;
            double rel_isolated_cutoff;
            size_t min_length_for_low_covered;
            double rel_low_covered_cutoff;
            double min_coverage;
        };

        std::map<std::string, PathFiltrationT> path_filtration;

        bool use_coordinated_coverage;

        struct CoordinatedCoverageT {
            size_t max_edge_length_in_repeat;
            double delta;
            size_t min_path_len;
        } coordinated_coverage;

        struct Scaffolding2015 {
            double relative_weight_cutoff;

            size_t unique_length_upper_bound;
            size_t unique_length_lower_bound;
            size_t unique_length_step;

            size_t graph_connectivity_max_edges;
        } scaffolding2015;

        struct ScaffoldGraphParamsT {
            bool construct;
            bool output;
            size_t always_add;
            size_t never_add;
            double relative_threshold;
            bool use_graph_connectivity;
            size_t max_path_length;
        } scaffold_graph_params;

        struct GenomeConsistencyCheckerParamsT {
            size_t max_gap;
            double relative_max_gap;
            bool use_main_storage;
            size_t unresolvable_jump;
            size_t unique_length;
        } genome_consistency_checker;

        struct LoopTraversalParamsT {
            size_t min_edge_length ;
            size_t max_component_size;
            size_t max_path_length;
        } loop_traversal;

        struct UniquenessAnalyserParamsT  {
            bool enabled;
            double unique_coverage_variation;

            double nonuniform_coverage_variation;
            double uniformity_fraction_threshold;
        } uniqueness_analyser;

    };

    struct LongReads {
        double filtering;
        double weight_priority;
        double unique_edge_priority;
        size_t min_significant_overlap;
    };

    struct AllLongReads {
        LongReads single_reads;
        LongReads pacbio_reads;
        LongReads contigs;
        LongReads meta_contigs;
    };


    struct MainPEParamsT {
        bool debug_output;
        std::string etc_dir;

        OutputParamsT output;
        VisualizeParamsT viz;
        ParamSetT param_set;
        AllLongReads long_reads;
    }; // params;

};

void load(pe_config::ParamSetT &p, boost::property_tree::ptree const &pt, bool complete = true);
void load(pe_config::MainPEParamsT &p, boost::property_tree::ptree const &pt, bool complete = true);

}

//typedef config_common::config<path_extend::pe_config> pe_cfg;

#endif /* CONFIG_STRUCT_HPP_ */
