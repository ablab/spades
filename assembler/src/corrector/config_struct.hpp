#pragma once

#include "config_singl.hpp"

#include "io/library.hpp"

namespace corrector {
// WTF: Enum naming style
enum strategy{all_reads, mapped_squared, not_started, majority_only};
struct corrector_config {
    io::DataSet<> dataset;
    std::string work_dir;
    std::string output_dir;
    unsigned max_nthreads;
    strategy strat;
    std::string bwa;
};

void load(corrector::corrector_config& cfg, const std::string &filename);
}

typedef config_common::config<corrector::corrector_config> corr_cfg;
