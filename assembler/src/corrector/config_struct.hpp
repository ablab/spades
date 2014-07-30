#pragma once

#include "config_singl.hpp"
#include "../debruijn/config_struct.hpp"

#include "io/library.hpp"

namespace corrector {

struct corrector_config {
//  io::DataSet<debruijn_graph::debruijn_config::DataSetData> dataset;
//std::string dataset;
  io::DataSet<> dataset;
  std::string work_dir;
  std::string output_dir;
  unsigned max_nthreads;
  std::string strategy;
  std::string bwa;
};

void load(corrector::corrector_config& cfg, const std::string &filename);
}


typedef config_common::config<corrector::corrector_config> corr_cfg;
