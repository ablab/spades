#pragma once

#include "config_singl.hpp"

#include "io/library.hpp"

namespace corrector {

struct corrector_config {
  io::DataSet<> dataset;

  std::string work_dir;
  std::string output_dir;
  unsigned max_nthreads;
  std::string strategy;
};

void load(corrector::corrector_config& cfg, const std::string &filename);
}

typedef config_common::config<corrector::corrector_config> corr_cfg;
