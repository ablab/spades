//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "pipeline/config_singl.hpp"

#include "pipeline/library.hpp"

namespace corrector {
enum class Strategy {
    AllReads = 1,
    MappedSquared = 2,
    AllExceptJustStarted = 3,
    MajorityOnly = 4
};
struct corrector_config {
    io::DataSet<> dataset;
    std::string work_dir;
    std::string output_dir;
    unsigned max_nthreads;
    Strategy strat;
    std::string bwa;
    std::string log_filename;
};

void load(corrector::corrector_config& cfg, const std::string &filename);
}

typedef config_common::config<corrector::corrector_config> corr_cfg;
