//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "config_struct.hpp"

#include "utils/parallel/openmp_wrapper.h"

#include "llvm/Support/YAMLParser.h"
#include "llvm/Support/YAMLTraits.h"

#include <string>

using namespace llvm;

namespace llvm {
namespace yaml {
template <>
struct ScalarEnumerationTraits<hammer_config::HammerStage> {
  static void enumeration(yaml::IO &io, hammer_config::HammerStage &value) {
    io.enumCase(value, "count", hammer_config::HammerStage::KMerCounting);
    io.enumCase(value, "hamcluster",
                hammer_config::HammerStage::HammingClustering);
    io.enumCase(value, "subcluster", hammer_config::HammerStage::SubClustering);
    io.enumCase(value, "correct", hammer_config::HammerStage::ReadCorrection);
  }
};

  template <>
  struct ScalarEnumerationTraits<hammer_config::CenterType> {
    static void enumeration(yaml::IO &io, hammer_config::CenterType &value) {
      io.enumCase(value, "count_max", hammer_config::CenterType::COUNT_ARGMAX);
      io.enumCase(value, "consensus", hammer_config::CenterType::CONSENSUS);
      io.enumCase(value, "posterior_consensus", hammer_config::CenterType::BY_POSTERIOR_QUALITY);
    }
  };
}  // namespace yaml
}  // namespace llvm

// FIXME: This is temporary
class DataSetReader {
 public:
  DataSetReader(yaml::IO &) {}
  DataSetReader(yaml::IO &, io::DataSet<> &) {}

  io::DataSet<> denormalize(yaml::IO &) { return io::DataSet<>(path); }

  std::string path;
};

namespace llvm {
namespace yaml {
template <>
struct MappingTraits<hammer_config::hammer_config> {
  static void mapping(yaml::IO &io, hammer_config::hammer_config &cfg) {
    yaml::MappingNormalization<DataSetReader, io::DataSet<>> dataset(
        io, cfg.dataset);

    io.mapRequired("dataset", dataset->path);
    io.mapOptional("working_dir", cfg.working_dir, std::string("."));
    io.mapOptional("output_dir", cfg.output_dir, std::string("."));
    io.mapRequired("hard_memory_limit", cfg.hard_memory_limit);
    io.mapOptional("count_split_buffer", cfg.count_split_buffer, 0ul);
    io.mapOptional("max_nthreads", cfg.max_nthreads, 1u);

    io.mapOptional("oracle_path", cfg.oracle_path, std::string(""));
    io.mapOptional("max_full_del", cfg.max_full_del, 1u);
    io.mapOptional("max_second_indel", cfg.max_second_indel, 1u);
    io.mapOptional("max_indel", cfg.max_indel, 3u);
    io.mapOptional("max_from_zero_insertion", cfg.max_from_zero_insertion, 1u);

    io.mapOptional("sample_rate", cfg.sample_rate, 1.0);
    io.mapOptional("subcluster_min_count", cfg.subcluster_min_count, 15u);
    io.mapOptional("good_threshold", cfg.good_threshold, -0.69);
    io.mapOptional("skip_threshold", cfg.skip_threshold, -0.01);
    io.mapOptional("subcluster_threshold", cfg.subcluster_threshold, -0.001);
    io.mapOptional("subcluster_filter_by_count", cfg.subcluster_filter_by_count_enabled, true);
    io.mapOptional("queue_limit_multiplier", cfg.queue_limit_multiplier, 500);
    io.mapOptional("dist_one_subcluster_alpha", cfg.dist_one_subcluster_alpha, 0.51);
    io.mapOptional("subcluster_qual_mult", cfg.subcluster_qual_mult, 1.0);
    io.mapOptional("subcluster_count_mult", cfg.subcluster_count_mult, 0.3);
    io.mapOptional("correction_penalty", cfg.correction_penalty, -7.0);
    io.mapOptional("bad_kmer_penalty", cfg.bad_kmer_penalty, -20.0);
    io.mapOptional("count_dist_eps", cfg.count_dist_eps, 1e-3);
    io.mapOptional("count_dist_skip_quantile", cfg.count_dist_skip_quantile, 0.05);
    io.mapOptional("noise_filter_count_threshold", cfg.noise_filter_count_threshold, 3u);
    io.mapOptional("center_type", cfg.center_type, hammer_config::CenterType::COUNT_ARGMAX);


    io.mapRequired("kmer_qual_threshold", cfg.kmer_qual_threshold);
    io.mapRequired("center_qual_threshold", cfg.center_qual_threshold);
    io.mapRequired("delta_score_threshold", cfg.delta_score_threshold);
    io.mapRequired("keep_uncorrected_ends", cfg.keep_uncorrected_ends);
    io.mapRequired("tau", cfg.tau);
    io.mapOptional("debug_mode", cfg.debug_mode, false);
    io.mapOptional("start_stage", cfg.start_stage,
                   hammer_config::HammerStage::KMerCounting);
  }
};
}  // namespace yaml
}  // namespace llvm

namespace hammer_config {
void load(hammer_config &cfg, const std::string &filename) {
  ErrorOr<std::unique_ptr<MemoryBuffer>> Buf = MemoryBuffer::getFile(filename);
  if (!Buf) throw(std::string("Failed to load config file ") + filename);

  yaml::Input yin(*Buf.get());
  yin >> cfg;

  if (yin.error()) throw(std::string("Failed to load config file ") + filename);

  cfg.max_nthreads = spades_set_omp_threads(cfg.max_nthreads);
}
}  // namespace hammer_config
