#ifndef __HAMMER_IT_CONFIG_HPP__
#define __HAMMER_IT_CONFIG_HPP__

#include "config_singl.hpp"

#include "io/library.hpp"

namespace hammer_config {
enum class HammerStage {
  KMerCounting = 1,
  HammingClustering = 2,
  SubClustering = 3,
  ReadCorrection = 4
};

struct hammer_config {
  io::DataSet dataset;

  unsigned max_nthreads;
  unsigned tau;
  unsigned hard_memory_limit;
  double kmer_qual_threshold;
  bool keep_uncorrected_ends;

  bool debug_mode;
  HammerStage start_stage;
};

void load(hammer_config& cfg, const std::string &filename);
}

typedef config_common::config<hammer_config::hammer_config> cfg;

#endif // __HAMMER_IT_CONFIG_HPP__
