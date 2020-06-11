//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************
#pragma once

#include "config_singl.hpp"
#include "library.hpp"
#include "library_data.hpp"
#include "utils/verify.hpp"
#include "modules/path_extend/pe_config_struct.hpp"
#include "configs/aligner_config.hpp"
#include "common/utils/logger/logger.hpp"

#include <boost/optional.hpp>
#include "math/xmath.h"

namespace debruijn_graph {
namespace config {

enum class info_printer_pos : char {
    default_pos = 0,
    before_raw_simplification,
    before_first_gap_closer,
    before_simplification,
    before_post_simplification,
    final_simplified,
    final_gap_closed,
    before_repeat_resolution,

    total
};

std::vector<std::string> InfoPrinterPosNames();

enum class pipeline_type : char {
    base = 0,
    isolate,
    mda,
    meta,
    moleculo,
    rna,
    plasmid,
    large_genome,
    metaextrachromosomal,
    rnaviral,
    total
};

class PipelineHelper {
public:
    static bool IsPlasmidPipeline(const pipeline_type pipeline);
    static bool IsMetagenomicPipeline(const pipeline_type pipeline);
    static bool IsRNAPipeline(const pipeline_type pipeline);
};


std::vector<std::string> PipelineTypeNames();

enum class resolving_mode : char {
    none = 0,
    path_extend,

    total
};

std::vector<std::string> ResolveModeNames();

enum class single_read_resolving_mode : char {
    none = 0,
    only_single_libs,
    all,

    total
};

enum class output_broken_scaffolds: char {
    none = 0,
    break_gaps,
    break_all,

    total
};

enum class Checkpoints : char {
    None = 0,
    Last,
    All
};

std::vector<std::string> SingleReadResolveModeNames();

template<typename mode_t>
mode_t ModeByName(const std::string& name, const std::vector<std::string>& names) {
    auto it = std::find(names.begin(), names.end(), name);
    CHECK_FATAL_ERROR(it != names.end(), "Unrecognized mode name: " << name);
    return mode_t(it - names.begin());
}

template<typename mode_t>
std::string ModeName(const mode_t& mode, const std::vector<std::string>& names) {
    auto mode_id = static_cast<size_t>(mode);
    CHECK_FATAL_ERROR(mode_id < names.size(), "Unrecognized mode id: " << mode_id);
    return names[mode_id];
}

struct dataset {
    typedef io::DataSet<LibraryData>::Library Library;

    dataset() :
        RL(0),
        no_merge_RL(0),
        aRL(0.0),
        average_coverage(0.0) {}

    size_t RL;
    size_t no_merge_RL;
    double aRL;
    double average_coverage;
    io::DataSet<LibraryData> reads;
    std::vector<std::string> reference_genome;
};

// struct for debruijn project's configuration file
struct debruijn_config {
    pipeline_type mode;
    bool uneven_depth;

    bool developer_mode;

    bool preserve_raw_paired_index;

    struct simplification {
        struct tip_clipper {
            std::string condition;
            tip_clipper() {}
            tip_clipper(std::string condition_) : condition(condition_) {}
        };

        struct dead_end_clipper {
            std::string condition;
            bool enabled;
        };

        struct topology_tip_clipper {
            double length_coeff;
            size_t uniqueness_length;
            size_t plausibility_length;
        };

        struct complex_tip_clipper {
            bool enabled;
            double max_relative_coverage;
            size_t max_edge_len;
            std::string condition;
        };

        struct bulge_remover {
            bool enabled;
            bool main_iteration_only;
            double max_bulge_length_coefficient;
            size_t max_additive_length_coefficient;
            double max_coverage;
            double max_relative_coverage;
            size_t max_delta;
            double max_relative_delta;
            size_t max_number_edges;
            size_t dijkstra_vertex_limit;
            bool parallel;
            size_t buff_size;
            double buff_cov_diff;
            double buff_cov_rel_diff;
            double min_identity;
        };

        struct erroneous_connections_remover {
            std::string condition;
            erroneous_connections_remover() {}
            erroneous_connections_remover(std::string condition_) : condition(condition_) {}
        };

        struct relative_coverage_ec_remover {
            bool enabled;
            size_t max_ec_length;
            double rcec_ratio;
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

        struct interstrand_ec_remover {
            size_t max_ec_length_coefficient;
            size_t uniqueness_length;
            size_t span_distance;
        };

        struct max_flow_ec_remover {
            bool enabled;
            double max_ec_length_coefficient;
            size_t uniqueness_length;
            size_t plausibility_length;
        };

        struct low_covered_edge_remover {
            bool enabled;
            double coverage_threshold;
        };

        struct isolated_edge_remover {
            bool enabled;
            bool use_rl_for_max_length; // max_length will be taken max with read_length
            size_t max_length;
            double max_coverage;
            bool use_rl_for_max_length_any_cov;
            size_t max_length_any_cov;
            size_t rl_threshold_increase; // add this value to read length if used, i.e. flags above are set
        };

        struct complex_bulge_remover {
            bool enabled;
            double max_relative_length;
            size_t max_length_difference;
        };

        struct hidden_ec_remover {
            bool enabled;
            size_t uniqueness_length;
            double unreliability_threshold;
            double relative_threshold;
        };

        struct relative_coverage_edge_disconnector {
            bool enabled;
            double diff_mult;
            size_t edge_sum;
            double unconditional_diff_mult;
        };

        struct relative_coverage_comp_remover {
            bool enabled;
            double coverage_gap;
            double length_coeff;
            double tip_allowing_length_coeff;
            size_t max_ec_length_coefficient;
            double max_coverage_coeff;
            size_t vertex_count_limit;
        };

        struct init_cleaning {
            std::string self_conj_condition;

            bool early_it_only;
            double activation_cov;
            isolated_edge_remover ier;
            std::string tip_condition;
            std::string ec_condition;
            double disconnect_flank_cov;
        };

        size_t cycle_iter_count;

        bool topology_simplif_enabled;
        tip_clipper tc;
        dead_end_clipper dead_end;
        complex_tip_clipper complex_tc;
        topology_tip_clipper ttc;
        bulge_remover br;
        erroneous_connections_remover ec;
        relative_coverage_ec_remover rcec;
        relative_coverage_comp_remover rcc;
        relative_coverage_edge_disconnector red;
        topology_based_ec_remover tec;
        tr_based_ec_remover trec;
        interstrand_ec_remover isec;
        max_flow_ec_remover mfec;
        isolated_edge_remover ier;
        complex_bulge_remover cbr;
        hidden_ec_remover her;
        low_covered_edge_remover lcer;
        tip_clipper final_tc;
        bulge_remover final_br;
        bulge_remover subspecies_br;
        init_cleaning init_clean;
    };

    struct construction {
        struct early_tip_clipper {
            bool enable;
            boost::optional<size_t> length_bound;
            early_tip_clipper() : enable(false) {}
        };

        early_tip_clipper early_tc;
        bool keep_perfect_loops;
        unsigned read_cov_threshold;
        size_t read_buffer_size;
        construction() :
                keep_perfect_loops(true),
                read_cov_threshold(0),
                read_buffer_size(0) {}
    };

    simplification simp;
    boost::optional<simplification> preliminary_simp;

    struct distance_estimator {
        double linkage_distance_coeff;
        double max_distance_coeff;
        double max_distance_coeff_scaff;
        double clustered_filter_threshold;
        unsigned raw_filter_threshold;
        double rounding_thr;
        double rounding_coeff;
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

    struct ambiguous_distance_estimator {
        bool enabled;
        double haplom_threshold;
        double relative_length_threshold;
        double relative_seq_threshold;
    };

    struct plasmid {
        size_t long_edge_length;
        size_t edge_length_for_median;
        double relative_coverage;
        size_t small_component_size;
        double small_component_relative_coverage;
        size_t min_component_length;
        size_t min_isolated_length;
        std::string reference_removal;
//settings for iterative removal;
        bool iterative_coverage_elimination;
        size_t additive_step; //5
        double relative_step; //1.3
        size_t max_length; //1000000
        bool output_linear; //false
        size_t min_circular_length; //500
        size_t min_linear_length; //1000
    };

    struct position_handler {
        size_t max_mapping_gap;
        size_t max_gap_diff;
        std::string contigs_for_threading;
        std::string contigs_to_analyze;
        bool late_threading;
        bool careful_labeling;
    };

    struct gap_closer {
        int minimal_intersection;
        bool before_raw_simplify;
        bool before_simplify;
        bool after_simplify;
        double weight_threshold;
        int max_dist_to_tip;
    };

    struct ss_coverage_splitter_t {
        bool enabled;
        size_t bin_size;
        size_t min_edge_len;
        double min_edge_coverage;
        double coverage_margin;
        double min_flanking_coverage;
    };

    struct info_printer {
        bool basic_stats;
        bool lib_info;
        bool extended_stats;
        bool write_components;
        std::string components_for_kmer;
        std::string components_for_genome_pos;
        bool write_components_along_genome;
        bool write_components_along_contigs;
        bool save_full_graph;
        bool save_all;
        bool save_graph_pack;
        bool write_error_loc;
        bool write_full_graph;
        bool write_full_nc_graph;
    };

    struct graph_read_corr_cfg {
        bool enable;
        std::string output_dir;
        bool binary;
    };

    struct kmer_coverage_model {
        double probability_threshold;
        double strong_probability_threshold;
        double coverage_threshold;
        bool use_coverage_threshold;
    };

    struct time_tracing {
        bool enable;
        unsigned granularity;
    };
    
    typedef std::map<info_printer_pos, info_printer> info_printers_t;

    std::string dataset_file;
    std::string input_dir;
    std::string output_base;
    std::string output_dir;
    std::string tmp_dir;
    Checkpoints checkpoints;
    std::string output_saves;
    std::string log_filename;
    std::string series_analysis;

    bool output_pictures;
    bool output_nonfinal_contigs;
    bool compute_paths_number;

    bool use_additional_contigs;
    std::string additional_contigs;

    struct scaffold_correction {
        std::string scaffolds_file;
        bool output_unfilled;
        size_t max_insert;
        size_t max_cut_length;
    };

    struct truseq_analysis {
        std::string scaffolds_file;
        std::string genome_file;
    };

    struct contig_output {
        std::string contigs_name;
        std::string scaffolds_name;
        output_broken_scaffolds obs_mode;
    };

    struct strand_specificity {
        bool ss_enabled;
        bool antisense;
    };

    struct hmm_matching {
        std::string hmm_set;
        size_t component_size_part;
        bool start_only_from_tips;
        bool set_copynumber;
    };
    
    contig_output co;

    boost::optional<scaffold_correction> sc_cor;
    truseq_analysis tsa;

    std::string load_from;
    std::string entry_point;

    bool rr_enable;
    bool two_step_rr;
    bool use_intermediate_contigs;
    size_t min_edge_length_for_is_count;

    single_read_resolving_mode single_reads_rr;
    bool use_single_reads;

    bool correct_mismatches;
    bool paired_info_statistics;
    bool paired_info_scaffolder;
    bool gap_closer_enable;

    size_t max_repeat_length;

    // Conversion options
    std::string temp_bin_reads_dir;
    std::string temp_bin_reads_path;
    std::string paired_read_prefix;
    std::string single_read_prefix;

    size_t K;

    bool main_iteration;

    unsigned max_threads;
    size_t max_memory;

    resolving_mode rm;
    path_extend::pe_config::MainPEParamsT pe_params;
    boost::optional<path_extend::pe_config::MainPEParamsT> prelim_pe_params;
    bool avoid_rc_connections;

    construction con;
    distance_estimator de;
    smoothing_distance_estimator ade;
    ambiguous_distance_estimator amb_de;
    pacbio_processor pb;
    bool use_scaffolder;
    dataset ds;
    position_handler pos;
    gap_closer gc;
    ss_coverage_splitter_t ss_coverage_splitter;
    graph_read_corr_cfg graph_read_corr;
    info_printers_t info_printers;
    kmer_coverage_model kcm;
    bwa_aligner bwa;
    boost::optional<plasmid> pd;
    boost::optional<hmm_matching> hm;
    size_t flanking_range;

    bool calculate_coverage_for_each_lib;
    strand_specificity ss;
    time_tracing tt;

    bool need_mapping;

    debruijn_config() :
            use_single_reads(false) {

    }
};

void init_libs(io::DataSet<LibraryData> &dataset, size_t max_threads,
               const std::string &temp_bin_reads_path);
void load(debruijn_config& cfg, const std::vector<std::string> &filenames);
void load(debruijn_config& cfg, const std::string &filename);
void load_lib_data(const std::string& prefix);
void load_lib_data(std::istream& is);
void write_lib_data(const std::string& prefix);
void write_lib_data(std::ostream& os);

} // config
} // debruijn_graph


typedef config_common::config<debruijn_graph::config::debruijn_config> cfg;
