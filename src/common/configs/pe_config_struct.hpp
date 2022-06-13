//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2015-2022 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef LC_CONFIG_STRUCT_HPP_
#define LC_CONFIG_STRUCT_HPP_

#include "configs/config_singl.hpp"
#include "utils/cpp_utils.hpp"

#include <boost/property_tree/ptree_fwd.hpp>

#include <filesystem>
#include <map>
#include <optional>
#include <string>
#include <vector>

namespace path_extend {

enum class scaffolding_mode {
    sm_old,
    sm_2015,
    sm_combined,
    sm_old_pe_2015,
    undefined
};

//Both this functions return always true, right?
//still necessary?
inline bool IsScaffolder2015Enabled(const scaffolding_mode mode) {
    return (mode == scaffolding_mode::sm_old_pe_2015 || mode == scaffolding_mode::sm_2015 || mode == scaffolding_mode::sm_combined);
}

inline bool IsOldPEEnabled(const scaffolding_mode mode) {
    return (mode == scaffolding_mode::sm_old_pe_2015 || mode == scaffolding_mode::sm_old || mode == scaffolding_mode::sm_combined);
}

// struct for path extend subproject's configuration file
struct pe_config {
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
            double coverage_margin;
            double min_upper_coverage;
            double max_coverage_variation;
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
            double rel_cov_cutoff;
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

            std::optional<int> min_overlap_for_rna_scaffolding;
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
        LongReads rna_long_reads;
    };

  struct ReadCloud {
    std::string tslr_dataset;

    size_t edge_tail_len;
    size_t frame_size;

    bool read_cloud_gap_closer_on;
    bool read_cloud_resolution_on;

    size_t long_edge_length_min_upper_bound;
    size_t long_edge_length_max_upper_bound;
    size_t long_edge_length_lower_bound;

    double gap_closer_connection_score_threshold;
    double gap_closer_relative_coverage_threshold;
    size_t gap_closer_connection_length_threshold;
    size_t gap_closer_scan_bound;

    size_t min_training_edges;
    size_t min_training_total_length;
    size_t optimal_training_total_length;

    double relative_score_threshold;

    struct scaffold_graph_construction {
      double score_percentile;
      double cluster_length_percentile;
      size_t count_threshold;
      double relative_coverage_threshold;
      size_t connection_length_threshold;
      size_t connection_count_threshold;
      double split_procedure_strictness;
      size_t transitive_distance_threshold;
      size_t path_scaffolder_tail_threshold;
      size_t path_scaffolder_count_threshold;
      size_t min_edge_length_for_barcode_collection;
      double path_scaffolding_score;
      double ultralong_edge_length_percentile;
      double short_edge_threshold;
    };

    scaffold_graph_construction scaff_con;

    struct scaffold_polisher {
      size_t max_scaffold_dijkstra_distance;
      size_t read_count_threshold;
      size_t path_cluster_linkage_distance;
      size_t path_cluster_min_reads;
      double share_threshold;
      double path_cluster_relative_threshold;
    };

    scaffold_polisher scaff_pol;

    struct read_cloud_extender {
      size_t reliable_edge_length;
      size_t tail_threshold;
      size_t distance_bound;
      size_t seed_edge_length;
      double extender_score_threshold;
      double relative_coverage_threshold;
      size_t barcode_threshold;
      size_t score_function_tail_threshold;
    };

    read_cloud_extender read_ext;

    struct path_searching {
      size_t max_path_growing_iterations;
      size_t max_paths_to_process;
      size_t max_edge_visits;
    };

    path_searching path_search;

    struct stats {
      std::filesystem::path genome_path;
      std::filesystem::path base_contigs_path;
      std::filesystem::path cloud_contigs_path;
      std::filesystem::path scaffold_graph_statistics;
    };

    stats statistics;

    bool path_scaffolding_on;
    bool debug_mode;
  };

    struct MainPEParamsT {
        bool debug_output;
        std::filesystem::path etc_dir;

        OutputParamsT output;
        VisualizeParamsT viz;
        ParamSetT param_set;
        AllLongReads long_reads;
        ReadCloud read_cloud;
    }; // params;

};

void load(pe_config::ParamSetT &p, boost::property_tree::ptree const &pt, bool complete = true);
void load(pe_config::MainPEParamsT &p, boost::property_tree::ptree const &pt, bool complete = true);

}

//typedef config_common::config<path_extend::pe_config> pe_cfg;

#endif /* CONFIG_STRUCT_HPP_ */
