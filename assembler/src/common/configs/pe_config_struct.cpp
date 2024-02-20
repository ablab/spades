//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "configs/pe_config_struct.hpp"
#include "configs/config_common.hpp"
#include "utils/logger/logger.hpp"

#include <llvm/ADT/StringSwitch.h>

namespace path_extend {

//convert string to vector of words separated by space
std::vector<std::string> StringToVector(const std::string& s) {
    std::string word;
    std::vector<std::string> res;
    for (size_t i = 0; i < s.length(); ++i) {
        if (s[i] == ' ') {
            if (word != "") {
                res.push_back(word);
                word = "";
            }
        }
        else {
            word += s[i];

        }
    }
    if (word != "") {
        res.push_back(word);
    }
    return res;
}

void load(scaffolding_mode &sm, boost::property_tree::ptree const& pt, std::string const& key, bool complete) {
    if (complete || pt.find(key) != pt.not_found()) {
        std::string ep = pt.get<std::string>(key);
        sm = llvm::StringSwitch<scaffolding_mode>(ep)
             .Case("old", scaffolding_mode::sm_old)
             .Case("2015", scaffolding_mode::sm_2015)
             .Case("combined", scaffolding_mode::sm_combined)
             .Case("old_pe_2015", scaffolding_mode::sm_old_pe_2015)
             .Default(scaffolding_mode::undefined);
        CHECK_FATAL_ERROR(sm != scaffolding_mode::undefined, "Invalid scaffolding mode");
    }
}

void load(pe_config::ParamSetT::ScaffoldGraphParamsT& sg, boost::property_tree::ptree const& pt, bool /*complete*/) {
    using config_common::load;
    load(sg.construct,          pt, "construct"         );
    load(sg.output,             pt, "output"            );
    load(sg.always_add,         pt, "always_add"        );
    load(sg.never_add,          pt, "never_add"         );
    load(sg.relative_threshold,  pt, "relative_threshold" );
    load(sg.use_graph_connectivity, pt, "use_graph_connectivity");
    load(sg.max_path_length,    pt, "max_path_length"   );
}

void load(pe_config::OutputParamsT& o, boost::property_tree::ptree const& pt, bool complete) {
  using config_common::load;

  load(o.write_overlaped_paths,   pt, "write_overlaped_paths" , complete);
  load(o.write_paths,             pt, "write_paths"           , complete);
}

void load(pe_config::VisualizeParamsT& o, boost::property_tree::ptree const& pt, bool complete) {
  using config_common::load;
  load(o.print_overlaped_paths,   pt, "print_overlaped_paths" , complete);
  load(o.print_paths,             pt, "print_paths"           , complete);
}

void load(pe_config::ParamSetT::ExtensionOptionsT& es,
          boost::property_tree::ptree const& pt, bool complete) {
    using config_common::load;
    load(es.priority_coeff, pt, "priority_coeff", complete);
    load(es.weight_threshold, pt, "weight_threshold", complete);
    load(es.single_threshold, pt, "single_threshold", complete);
    load(es.max_repeat_length, pt, "max_repeat_length", complete);
}


void load(pe_config::ParamSetT::CoordinatedCoverageT& coord_cov,
          boost::property_tree::ptree const& pt, bool complete) {
    using config_common::load;
    load(coord_cov.max_edge_length_in_repeat, pt, "max_edge_length_repeat", complete);
    load(coord_cov.delta, pt, "delta", complete);
    load(coord_cov.min_path_len, pt, "min_path_len", complete);
}

void load(pe_config::ParamSetT::ScaffolderOptionsT& so, 
            boost::property_tree::ptree const& pt, bool complete) {
    using config_common::load;
    load(so.enabled, pt, "enabled"      , complete);
    load(so.cutoff      , pt, "cutoff", complete);
    load(so.hard_cutoff      , pt, "hard_cutoff", complete);
    load(so.rel_cov_cutoff      , pt, "rel_cov_cutoff", complete);
    load(so.sum_threshold      , pt, "sum_threshold", complete);

    load(so.cluster_info      , pt, "cluster_info", complete);
    load(so.cl_threshold      , pt, "cl_threshold", complete);

    load(so.use_la_gap_joiner      , pt, "use_la_gap_joiner", complete);
    load(so.min_gap_score      , pt, "min_gap_score", complete);
    load(so.max_can_overlap      , pt, "max_can_overlap", complete);
    load(so.short_overlap      , pt, "short_overlap", complete);
    load(so.artificial_gap      , pt, "artificial_gap", complete);
    load(so.min_overlap_length, pt, "min_overlap_length", complete);
    load(so.flank_multiplication_coefficient, pt, "flank_multiplication_coefficient", complete);
    load(so.flank_addition_coefficient, pt, "flank_addition_coefficient", complete);

    load(so.var_coeff          , pt, "var_coeff", complete);
    load(so.basic_overlap_coeff, pt, "basic_overlap_coeff", complete);

    if (pt.count("min_overlap_for_rna_scaffolding")) {
        VERIFY_MSG(!so.min_overlap_for_rna_scaffolding, "Option can be loaded only once");
        so.min_overlap_for_rna_scaffolding.reset();
        load(so.min_overlap_for_rna_scaffolding, pt, "min_overlap_for_rna_scaffolding");
    }
}


void load(pe_config::ParamSetT::PathFiltrationT& pf,
          boost::property_tree::ptree const& pt, bool complete)
{
    using config_common::load;
    load(pf.enabled      , pt, "enabled"      , complete);
    if (pf.enabled) {
        load(pf.min_length      , pt, "min_length"      , complete);
        load(pf.isolated_min_length      , pt, "isolated_min_length"      , complete);
        load(pf.isolated_min_cov      , pt, "isolated_min_cov"      , complete);
        load(pf.min_length_for_low_covered      , pt, "min_length_for_low_covered"      , complete);
        load(pf.rel_cutoff      , pt, "rel_cutoff"      , complete);
        load(pf.rel_isolated_cutoff      , pt, "rel_isolated_cutoff"      , complete);
        load(pf.rel_low_covered_cutoff      , pt, "rel_low_covered_cutoff"      , complete);
        load(pf.min_coverage      , pt, "min_coverage"      , complete);
    }
}

void load(pe_config::ParamSetT::GenomeConsistencyCheckerParamsT& gcc,
          boost::property_tree::ptree const& pt, bool complete)
{
    using config_common::load;
    load(gcc.max_gap      , pt, "max_gap"      , complete);
    load(gcc.relative_max_gap      , pt, "relative_max_gap"      , complete);
    load(gcc.use_main_storage      , pt, "use_main_storage"      , complete);
    load(gcc.unresolvable_jump      , pt, "unresolvable_jump"      , complete);
    load(gcc.unique_length      , pt, "unique_length"      , complete);

}

void load(pe_config::ParamSetT::OverlapRemovalOptionsT& ors,
          boost::property_tree::ptree const& pt, bool complete) {
    using config_common::load;
    load(ors.enabled, pt, "enabled"      , complete);
    load(ors.end_start_only, pt, "end_start_only"      , complete);
    load(ors.cut_all, pt, "cut_all"      , complete);
}

void load(pe_config::ParamSetT::SimpleCoverageResolver& scr,
          boost::property_tree::ptree const& pt, bool complete)
{
    using config_common::load;
    load(scr.enabled      , pt, "enabled"      , complete);
    load(scr.coverage_margin      , pt, "coverage_margin"      , complete);
    load(scr.min_upper_coverage      , pt, "min_upper_coverage"      , complete);
    load(scr.max_coverage_variation      , pt, "max_coverage_variation"      , complete);

}

void load(pe_config::ParamSetT& p, boost::property_tree::ptree const& pt, bool complete) {
    using config_common::load;
    load(p.sm, pt, "scaffolding_mode", complete);
    load(p.normalize_weight, pt,  "normalize_weight", complete);
    load(p.overlap_removal, pt, "overlap_removal", complete);
    load(p.multi_path_extend, pt, "multi_path_extend", complete);
    load(p.extension_options, pt, "extension_options", complete);
    load(p.mate_pair_options, pt, "mate_pair_options", complete);
    load(p.scaffolder_options, pt, "scaffolder", complete);
    load(p.coordinated_coverage, pt, "coordinated_coverage", complete);
    load(p.use_coordinated_coverage, pt, "use_coordinated_coverage", complete);
    load(p.scaffolding2015, pt, "scaffolding2015", complete);
    load(p.scaffold_graph_params, pt, "scaffold_graph", complete);

    std::string path_cleaning_presets;
    load(path_cleaning_presets, pt, "path_cleaning_presets", complete);
    auto presets = StringToVector(path_cleaning_presets);
    for (auto &key : presets) {
        pe_config::ParamSetT::PathFiltrationT path_filtration;
        std::string config_key = key == "default" ? "path_cleaning" : key + "_path_cleaning";
        load(path_filtration, pt, config_key, complete);
        p.path_filtration[key] = path_filtration;
    }
    load(p.genome_consistency_checker, pt, "genome_consistency_checker", complete);
    load(p.uniqueness_analyser, pt, "uniqueness_analyser", complete);
    load(p.loop_traversal, pt, "loop_traversal", complete);
    load(p.simple_coverage_resolver, pt, "simple_coverage_resolver", complete);
}

void load(pe_config::LongReads& p, boost::property_tree::ptree const& pt,
          bool complete) {
    using config_common::load;
    load(p.filtering, pt, "filtering", complete);
    load(p.weight_priority, pt, "weight_priority", complete);
    load(p.unique_edge_priority, pt, "unique_edge_priority", complete);
    load(p.min_significant_overlap, pt, "min_significant_overlap", complete);

}

void load(pe_config::ParamSetT::LoopTraversalParamsT& p, boost::property_tree::ptree const& pt,
          bool complete) {
    using config_common::load;
    load(p.min_edge_length, pt, "min_edge_length", complete);
    load(p.max_component_size, pt, "max_component_size", complete);
    load(p.max_path_length, pt, "max_path_length", complete);
}

void load(pe_config::ParamSetT::UniquenessAnalyserParamsT& p, boost::property_tree::ptree const& pt,
          bool complete) {
    using config_common::load;
    load(p.enabled, pt, "enabled", complete);
    load(p.nonuniform_coverage_variation, pt, "nonuniform_coverage_variation", complete);
    load(p.uniformity_fraction_threshold, pt, "uniformity_fraction_threshold", complete);
    load(p.unique_coverage_variation, pt, "unique_coverage_variation", complete);
}

void load(pe_config::ParamSetT::Scaffolding2015& p, boost::property_tree::ptree const& pt,
          bool complete) {
    using config_common::load;
    load(p.unique_length_lower_bound, pt, "unique_length_lower_bound", complete);
    load(p.unique_length_upper_bound, pt, "unique_length_upper_bound", complete);
    load(p.unique_length_step, pt, "unique_length_step", complete);
    load(p.graph_connectivity_max_edges, pt, "graph_connectivity_max_edges", complete);
    load(p.relative_weight_cutoff, pt, "relative_weight_cutoff", complete);
}

void load(pe_config::AllLongReads& p, boost::property_tree::ptree const& pt,
          bool complete) {
    using config_common::load;
    load(p.pacbio_reads, pt, "pacbio_reads", complete);
    load(p.single_reads, pt, "single_reads", complete);
    load(p.contigs, pt, "contigs", complete);
    load(p.meta_contigs, pt, "meta_untrusted_contigs", complete);
    load(p.rna_long_reads, pt, "rna_long_reads", complete);
}

void load(pe_config::ReadCloud::stats& statistics,
          boost::property_tree::ptree const& pt, bool /*complete*/) {
    using config_common::load;
    load(statistics.genome_path, pt, "genome_path");
    load(statistics.base_contigs_path, pt, "base_contigs_path");
    load(statistics.cloud_contigs_path, pt, "cloud_contigs_path");
    load(statistics.scaffold_graph_statistics, pt, "scaffold_graph_statistics");
}

void load(pe_config::ReadCloud::scaffold_polisher& scaff_pol,
          boost::property_tree::ptree const& pt, bool /*complete*/) {
    using config_common::load;
    load(scaff_pol.share_threshold, pt, "share_threshold");
    load(scaff_pol.read_count_threshold, pt, "read_count_threshold");
    load(scaff_pol.max_scaffold_dijkstra_distance, pt, "max_scaffold_distance");
    load(scaff_pol.path_cluster_linkage_distance, pt, "path_cluster_linkage_distance");
    load(scaff_pol.path_cluster_min_reads, pt, "path_cluster_min_reads");
    load(scaff_pol.path_cluster_relative_threshold, pt, "path_cluster_score_threshold");
}

void load(pe_config::ReadCloud::scaffold_graph_construction& scaff_con,
          boost::property_tree::ptree const &pt, bool /*complete*/) {
    using config_common::load;
    load(scaff_con.score_percentile, pt, "score_percentile");
    load(scaff_con.cluster_length_percentile, pt, "cluster_length_percentile");
    load(scaff_con.count_threshold, pt, "count_threshold");
    load(scaff_con.initial_distance, pt, "initial_distance");
    load(scaff_con.relative_coverage_threshold, pt, "relative_coverage_threshold");
    load(scaff_con.connection_length_threshold, pt, "connection_length_threshold");
    load(scaff_con.connection_count_threshold, pt, "connection_count_threshold");
    load(scaff_con.split_procedure_strictness, pt, "split_strictness");
    load(scaff_con.transitive_distance_threshold, pt, "transitive_distance_threshold");
    load(scaff_con.path_scaffolder_tail_threshold, pt, "path_scaffolder_tail_threshold");
    load(scaff_con.path_scaffolder_count_threshold, pt, "path_scaffolder_count_threshold");
    load(scaff_con.min_edge_length_for_barcode_collection, pt, "min_edge_length_for_barcode_collection");
    load(scaff_con.path_scaffolding_score, pt, "path_scaffolding_score");
    load(scaff_con.ultralong_edge_length_percentile, pt, "ultralong_edge_length_percentile");
    load(scaff_con.short_edge_threshold, pt, "short_edge_threshold");
}

void load(pe_config::ReadCloud::read_cloud_extender& read_ext,
          boost::property_tree::ptree const &pt, bool /*complete*/) {
    using config_common::load;
    load(read_ext.reliable_edge_length, pt, "reliable_edge_length");
    load(read_ext.tail_threshold, pt, "tail_threshold");
    load(read_ext.distance_bound, pt, "distance_bound");
    load(read_ext.seed_edge_length, pt, "seed_edge_length");
    load(read_ext.extender_score_threshold, pt, "extender_score_threshold");
    load(read_ext.relative_coverage_threshold, pt, "relative_coverage_threshold");
    load(read_ext.barcode_threshold, pt, "barcode_threshold");
    load(read_ext.score_function_tail_threshold, pt, "score_function_tail_threshold");
}

void load(pe_config::ReadCloud::path_searching &path_search,
          boost::property_tree::ptree const &pt, bool /*complete*/) {
    using config_common::load;
    load(path_search.max_path_growing_iterations, pt, "max_path_growing_iterations");
    load(path_search.max_paths_to_process, pt, "max_paths_to_process");
    load(path_search.max_edge_visits, pt, "max_edge_visits");
}

void load(pe_config::ReadCloud& read_cloud,
          boost::property_tree::ptree const& pt, bool /*complete*/) {
    using config_common::load;
    load(read_cloud.tslr_dataset, pt, "tslr_dataset");
    load(read_cloud.edge_tail_len, pt, "edge_tail_len");
    load(read_cloud.frame_size, pt, "frame_size");
    load(read_cloud.read_cloud_gap_closer_on, pt, "read_cloud_gap_closer_on");
    load(read_cloud.read_cloud_resolution_on, pt, "read_cloud_resolution_on");
    load(read_cloud.scaff_pol, pt, "scaffold_polisher");
    load(read_cloud.scaff_con, pt, "scaffold_graph_construction");
    load(read_cloud.read_ext, pt, "read_cloud_extender");
    load(read_cloud.long_edge_length_min_upper_bound, pt, "long_edge_length_min_upper_bound");
    load(read_cloud.long_edge_length_max_upper_bound, pt, "long_edge_length_max_upper_bound");
    load(read_cloud.long_edge_length_lower_bound, pt, "long_edge_length_lower_bound");
    load(read_cloud.min_training_edges, pt, "min_training_edges");
    load(read_cloud.min_training_total_length, pt, "min_training_total_length");
    load(read_cloud.optimal_training_total_length, pt, "optimal_training_total_length");
    load(read_cloud.path_search, pt, "path_searching");
    load(read_cloud.statistics, pt, "statistics");
    load(read_cloud.path_scaffolding_on, pt, "path_scaffolding_on");
    load(read_cloud.debug_mode, pt, "debug_mode");
    load(read_cloud.gap_closer_connection_score_threshold, pt, "gap_closer_connection_score_threshold");
    load(read_cloud.gap_closer_relative_coverage_threshold, pt, "gap_closer_relative_coverage_threshold");
    load(read_cloud.gap_closer_connection_length_threshold, pt, "gap_closer_connection_length_threshold");
    load(read_cloud.gap_closer_scan_bound, pt, "gap_closer_scan_bound");
    load(read_cloud.relative_score_threshold, pt, "relative_score_threshold");
}

void load(pe_config::MainPEParamsT& p, boost::property_tree::ptree const& pt,
          bool complete) {
    using config_common::load;
    load(p.debug_output, pt, "debug_output", complete);
    load(p.output, pt, "output", complete);
    load(p.viz, pt, "visualize", complete);
    load(p.param_set, pt, "params", complete);
    load(p.long_reads, pt, "long_reads", complete);
    load(p.read_cloud, pt, "read_cloud", complete);
    if (!p.debug_output) {
        p.output.DisableAll();
        p.viz.DisableAll();
    }
    p.etc_dir = "path_extend";
}

//// main long contigs config load function
//void load(pe_config& pe_cfg, boost::property_tree::ptree const& pt, bool complete) {
//  using config_common::load;
//
//  load(pe_cfg.dataset_name           , pt, "dataset", complete);
//  load(pe_cfg.params                 , pt, "pe_params", complete);
//}

};

