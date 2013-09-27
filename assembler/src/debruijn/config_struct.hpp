//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef CONFIG_STRUCT_HDIP_
#define CONFIG_STRUCT_HDIP_

#include "config_singl.hpp"
#include "cpp_utils.hpp"
#include "sequence/sequence.hpp"
#include "path_extend/pe_config_struct.hpp"

#include "io/library.hpp"
#include "io/binary_streams.hpp"
#include "io/rc_reader_wrapper.hpp"
#include "io/read_stream_vector.hpp"

#include <boost/bimap.hpp>
#include "xmath.h"

namespace debruijn_graph {

enum working_stage {
  ws_construction,
  ws_simplification,
  ws_late_pair_info_count,
  ws_distance_estimation,
  ws_pacbio_aligning,
  ws_repeats_resolving
};

enum construction_mode {
	con_old, con_extention
};

enum estimation_mode {
  em_simple, em_weighted, em_extensive, em_smoothing
};

enum resolving_mode {
  rm_none,
  rm_path_extend,
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
  const char* names[] = { "default", "before_first_gap_closer",
                          "before_simplification", "tip_clipping", "bulge_removal",
                          "err_con_removal", "before_final_err_con_removal",
                          "final_err_con_removal", "final_tip_clipping",
                          "final_bulge_removal", "removing_isolated_edges",
                          "final_simplified", "before_repeat_resolution" };

  utils::check_array_size < ipp_total > (names);
  return names[pos];
}

} // namespace details

// struct for debruijn project's configuration file
struct debruijn_config {

  typedef boost::bimap<std::string, working_stage> stage_name_id_mapping;
  typedef boost::bimap<std::string, construction_mode> construction_mode_id_mapping;
  typedef boost::bimap<std::string, estimation_mode> estimation_mode_id_mapping;
  typedef boost::bimap<std::string, resolving_mode> resolve_mode_id_mapping;

  //  bad fix, it is to be removed! To determine is it started from run.sh or from spades.py
  bool run_mode;

  bool developer_mode;

  static const stage_name_id_mapping FillStageInfo() {
    stage_name_id_mapping::value_type info[] = {
      stage_name_id_mapping::value_type("construction", ws_construction),
      stage_name_id_mapping::value_type("simplification", ws_simplification),
      stage_name_id_mapping::value_type("late_pair_info_count", ws_late_pair_info_count),
      stage_name_id_mapping::value_type("distance_estimation",  ws_distance_estimation),
      stage_name_id_mapping::value_type("pacbio_aligning",  ws_pacbio_aligning),
      stage_name_id_mapping::value_type("repeats_resolving", ws_repeats_resolving)
    };

    return stage_name_id_mapping(info, utils::array_end(info));
  }

  static const construction_mode_id_mapping FillConstructionModeInfo() {
     construction_mode_id_mapping::value_type info[] = {
         construction_mode_id_mapping::value_type("old", con_old),
         construction_mode_id_mapping::value_type("extension", con_extention), };
     return construction_mode_id_mapping(info, utils::array_end(info));
  }

  static const estimation_mode_id_mapping FillEstimationModeInfo() {
    estimation_mode_id_mapping::value_type info[] = {
      estimation_mode_id_mapping::value_type("simple", em_simple),
      estimation_mode_id_mapping::value_type("weighted", em_weighted),
      estimation_mode_id_mapping::value_type("extensive", em_extensive),
      estimation_mode_id_mapping::value_type("smoothing", em_smoothing)
    };
    return estimation_mode_id_mapping(info, utils::array_end(info));
  }

  static const resolve_mode_id_mapping FillResolveModeInfo() {
    resolve_mode_id_mapping::value_type info[] = {
      resolve_mode_id_mapping::value_type("none", rm_none),
      resolve_mode_id_mapping::value_type("path_extend", rm_path_extend),
      resolve_mode_id_mapping::value_type("rectangles", rm_rectangles)
    };

    return resolve_mode_id_mapping(info, utils::array_end(info));
  }

  static const stage_name_id_mapping& working_stages_info() {
    static stage_name_id_mapping working_stages_info = FillStageInfo();
    return working_stages_info;
  }

  static const construction_mode_id_mapping& construction_mode_info() {
    static construction_mode_id_mapping con_mode_info =
    FillConstructionModeInfo();
    return con_mode_info;
  }

  static const estimation_mode_id_mapping& estimation_mode_info() {
    static estimation_mode_id_mapping est_mode_info = FillEstimationModeInfo();
    return est_mode_info;
  }

  static const resolve_mode_id_mapping& resolve_mode_info() {
    static resolve_mode_id_mapping info = FillResolveModeInfo();
    return info;
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

  static const std::string& construction_mode_name(construction_mode con_id) {
    auto it = construction_mode_info().right.find(con_id);
    VERIFY_MSG(it != construction_mode_info().right.end(),
        "No name for construction mode id = " << con_id);
    return it->second;
  }

  static construction_mode construction_mode_id(std::string name) {
    auto it = construction_mode_info().left.find(name);
    VERIFY_MSG(it != construction_mode_info().left.end(),
        "There is no construction mode with name = " << name);

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
      std::string condition;
    };

    struct topology_tip_clipper {
      double length_coeff;
      size_t uniqueness_length;
      size_t plausibility_length;
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
      std::string condition;
    };

    struct relative_coverage_ec_remover {
      size_t max_ec_length_coefficient;
      double max_coverage_coeff;
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

    tip_clipper tc;
    topology_tip_clipper ttc;
    bulge_remover br;
    erroneous_connections_remover ec;
    relative_coverage_ec_remover rec;
    topology_based_ec_remover tec;
    tr_based_ec_remover trec;
    interstrand_ec_remover isec;
    max_flow_ec_remover mfec;
    isolated_edges_remover ier;
    complex_bulge_remover cbr;
  };

  struct construction {
      struct early_tip_clipper {
          bool enable;
          boost::optional<size_t> length_bound;
      };

      construction_mode con_mode;
      early_tip_clipper early_tc;
      bool keep_perfect_loops;
  };

  std::string uncorrected_reads;
  bool need_consensus;
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

  struct pacbio_processor {
//align and traverse.
	std::string pacbio_reads;

	size_t  pacbio_k; //13
	bool additional_debug_info; //false
	bool pacbio_optimized_sw; //false
	double compression_cutoff;// 0.6
	double domination_cutoff; //1.5
	double path_limit_stretching; //1.3
	double path_limit_pressing;//0.7
//	double gap_closing_relative_iterations; // 20.0
	int gap_closing_iterations; //5000;

//gap_closer
	size_t long_seq_limit; //400
	int split_cutoff; //100
	int match_value; // 1
	int mismatch_penalty; //1
	int insertion_penalty; //2
	int deletion_penalty; //2
  };

  struct DataSetData {
    size_t read_length;
    double mean_insert_size;
    double insert_size_deviation;
    double median_insert_size;
    double insert_size_mad;
    std::map<int, size_t> insert_size_distribution;

    uint64_t total_nucls;
    double average_coverage;
    double pi_threshold;

    std::string paired_read_prefix;
    std::string single_read_prefix;
    size_t thread_num;

    typedef io::IReader<io::SingleReadSeq> SequenceSingleReadStream;
    typedef io::IReader<io::PairedReadSeq> SequencePairedReadStream;

    DataSetData()
                : read_length(0),
                  mean_insert_size(0.0),
                  insert_size_deviation(0.0),
                  median_insert_size(0.0),
                  insert_size_mad(0.0),
                  total_nucls(0),
                  average_coverage(0.0),
                  pi_threshold(0.0) {
        }
    };

  struct dataset {
    io::DataSet<DataSetData> reads;

    size_t max_read_length;
    double average_coverage;

    size_t RL() const { return max_read_length; }
    void set_RL(size_t RL) {
        max_read_length = RL;
    }

    double avg_coverage() const { return average_coverage; }
    void set_avg_coverage(double avg_coverage) {
        average_coverage = avg_coverage;
        for (size_t i = 0; i < reads.lib_count(); ++i) {
            reads[i].data().average_coverage = avg_coverage;
        }
    }

    bool single_cell;
    std::string reference_genome_filename;
    std::string reads_filename;

    Sequence reference_genome;

    dataset(): max_read_length(0), average_coverage(0.0) {
    }
  };

  struct position_handler {
    int max_single_gap;
    std::string contigs_for_threading;
    std::string contigs_to_analyze;
    bool late_threading;
    bool careful_labeling;
  };

  struct gap_closer {
    int minimal_intersection;
    bool before_simplify;
    bool in_simplify;
    bool after_simplify;
    double weight_threshold;
  };

  struct info_printer {
    bool print_stats;
    bool write_components;
    std::string components_for_kmer;
    std::string components_for_genome_pos;
    bool write_components_along_genome;
    bool write_components_along_contigs;
    bool save_full_graph;
    bool write_error_loc;
    bool write_full_graph;
    bool write_full_nc_graph;
  };

  struct graph_read_corr_cfg {
    bool enable;
    std::string output_dir;
    bool binary;
  };

  typedef std::map<info_printer_pos, info_printer> info_printers_t;

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
  bool topology_simplif_enabled;
  bool use_unipaths;
  std::string additional_contigs;


  bool pacbio_test_on;
  bool coverage_based_rr_on;
  struct coverage_based_rr {
  double coverage_threshold_one_list;
  double coverage_threshold_match;
  double coverage_threshold_global;
  double tandem_ratio_lower_threshold;
  double tandem_ratio_upper_threshold;
  double repeat_length_upper_threshold;
 };

  coverage_based_rr cbrr;

  std::string load_from;

  working_stage entry_point;

  bool paired_mode;
  bool long_single_mode;
  bool divide_clusters;

  bool mismatch_careful;
  bool correct_mismatches;
  bool paired_info_statistics;
  bool paired_info_scaffolder;
  bool cut_bad_connections;
  bool componential_resolve;
  bool gap_closer_enable;

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

  estimation_mode est_mode;

  resolving_mode rm;
  path_extend::pe_config::MainPEParamsT pe_params;
  bool avoid_rc_connections;

  construction con;
  distance_estimator de;
  smoothing_distance_estimator ade;
  repeat_resolver rr;
  pacbio_processor pb;
  bool use_scaffolder;
  bool mask_all;
  dataset ds;
  position_handler pos;
  gap_closer gc;
  graph_read_corr_cfg graph_read_corr;
  info_printers_t info_printers;
};

void load(debruijn_config& cfg, const std::string &filename);
void load_lib_data(const std::string& prefix);
void write_lib_data(const std::string& prefix);
} // debruijn_graph

typedef config_common::config<debruijn_graph::debruijn_config> cfg;

#endif
