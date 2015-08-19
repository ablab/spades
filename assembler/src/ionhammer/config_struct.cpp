#include "config_struct.hpp"

#include "openmp_wrapper.h"

#include <yaml-cpp/yaml.h>
#include <string>

namespace YAML {
template<>
struct convert<hammer_config::HammerStage> {
  static bool decode(const YAML::Node &node, hammer_config::HammerStage &rhs) {
    std::string val = node.as<std::string>();

    if (val == "count") {
      rhs = hammer_config::HammerStage::KMerCounting;
      return true;
    } else if (val == "hamcluster") {
      rhs = hammer_config::HammerStage::HammingClustering;
      return true;
    } else if (val == "subcluster") {
      rhs = hammer_config::HammerStage::SubClustering;
      return true;
    } else if (val == "correct") {
      rhs = hammer_config::HammerStage::ReadCorrection;
      return true;
    }

    return false;
  }
};
}


namespace hammer_config {
void load(hammer_config& cfg, const std::string &filename) {
  YAML::Node config = YAML::LoadFile(filename);

  cfg.dataset.load(config["dataset"].as<std::string>());

  cfg.working_dir = config["working_dir"].as<std::string>(".");
  cfg.output_dir = config["output_dir"].as<std::string>(".");

  // FIXME: Make trivial deserialization trivial
  cfg.hard_memory_limit = config["hard_memory_limit"].as<unsigned>();

  cfg.count_split_buffer = config["count_split_buffer"].as<size_t>(0);

  cfg.max_nthreads = config["max_nthreads"].as<unsigned>();
  // Fix number of threads according to OMP capabilities.
  cfg.max_nthreads = std::min(cfg.max_nthreads, (unsigned)omp_get_max_threads());
  // Inform OpenMP runtime about this :)
  omp_set_num_threads(cfg.max_nthreads);

  cfg.kmer_qual_threshold = config["kmer_qual_threshold"].as<double>();
  cfg.center_qual_threshold = config["center_qual_threshold"].as<double>();
  cfg.delta_score_threshold = config["delta_score_threshold"].as<double>();
  cfg.keep_uncorrected_ends = config["keep_uncorrected_ends"].as<bool>();
  cfg.tau = config["tau"].as<unsigned>();

  cfg.debug_mode = config["debug_mode"].as<bool>();
  cfg.start_stage = config["start_stage"].as<HammerStage>();
}
}
