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
};

void load(hammer_config& cfg, const std::string &filename);
}

typedef config_common::config<hammer_config::hammer_config> cfg;

#endif // __HAMMER_IT_CONFIG_HPP__
