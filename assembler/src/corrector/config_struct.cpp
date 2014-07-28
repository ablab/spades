#include "config_struct.hpp"

#include "openmp_wrapper.h"

#include <yaml-cpp/yaml.h>
#include <string>


namespace corrector {
void load(corrector_config& cfg, const std::string &filename) {
//  INFO("loading config from " << filename);
  YAML::Node config = YAML::LoadFile(filename);
  INFO("loaded");
//  INFO(config["dataset"].as<std::string>("."));
//  std::string dataset = config["dataset"].as<std::string>(".");
//  YAML::Node tmp = YAML::LoadFile(config["dataset"].as<std::string>("."));
//  INFO("yamlnode ok");
//  cfg.dataset.load(dataset);
//  cfg.dataset = config["dataset"].as<std::string>(".");
  INFO("A");
  cfg.dataset.load(config["dataset"].as<std::string>());

  cfg.work_dir = config["work_dir"].as<std::string>(".");
  INFO("B");
  cfg.output_dir = config["output_dir"].as<std::string>(".");
  INFO("C");
  cfg.max_nthreads = config["max_nthreads"].as<unsigned>();
  INFO("D");
  // Fix number of threads according to OMP capabilities.
  cfg.max_nthreads = std::min(cfg.max_nthreads, (unsigned)omp_get_max_threads());
  INFO("E");
  cfg.strategy = config["strategy"].as<std::string>(".");
  INFO("F");
// Inform OpenMP runtime about this :)
  omp_set_num_threads(cfg.max_nthreads);
}
}
