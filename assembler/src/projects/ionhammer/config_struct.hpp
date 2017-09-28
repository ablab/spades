//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef __HAMMER_IT_CONFIG_HPP__
#define __HAMMER_IT_CONFIG_HPP__

#include "pipeline/config_singl.hpp"
#include "pipeline/library.hpp"

namespace hammer_config {
enum class HammerStage {
  KMerCounting = 1,
  HammingClustering = 2,
  SubClustering = 3,
  ReadCorrection = 4
};

enum class CenterType { COUNT_ARGMAX, CONSENSUS, BY_POSTERIOR_QUALITY };


  struct hammer_config {
  io::DataSet<> dataset;

  std::string working_dir;
  std::string output_dir;

  unsigned max_nthreads;
  unsigned tau;
  unsigned hard_memory_limit;

  size_t count_split_buffer;

  double kmer_qual_threshold;
  double center_qual_threshold;
  double delta_score_threshold;
  bool keep_uncorrected_ends;

  bool debug_mode;
  HammerStage start_stage;

  double sample_rate = 1.0;
  unsigned max_full_del = 1;
  unsigned max_indel = 3;
  unsigned max_second_indel = 1;
  unsigned max_from_zero_insertion = 1;
  std::string oracle_path = "";

  unsigned subcluster_min_count = 15;
  double good_threshold = -0.69;
  double skip_threshold = -0.01;
  double subcluster_threshold = -0.001;
  bool subcluster_filter_by_count_enabled = true;
  int queue_limit_multiplier = 200;
  double dist_one_subcluster_alpha = 0.6;
  double subcluster_qual_mult = 1.0;
  double subcluster_count_mult = 0.4;

  double correction_penalty = -7;
  double bad_kmer_penalty = -20;
  double count_dist_eps = 1e-3;
  double count_dist_skip_quantile = 0.05;

  unsigned noise_filter_count_threshold = 3;
  CenterType center_type = CenterType::COUNT_ARGMAX;
};

void load(hammer_config& cfg, const std::string& filename);
}  // namespace hammer_config

typedef config_common::config<hammer_config::hammer_config> cfg;

#endif  // __HAMMER_IT_CONFIG_HPP__
