//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2015-2022 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "configs/config_singl.hpp"

#include "library/library.hpp"

namespace corrector {
enum class Strategy {
    AllReads = 1,
    MappedSquared = 2,
    AllExceptJustStarted = 3,
    MajorityOnly = 4
};
struct corrector_config {
    io::DataSet<> dataset;
    std::filesystem::path work_dir;
    std::filesystem::path output_dir;
    unsigned max_nthreads;
    Strategy strat;
    std::string bwa;
    std::filesystem::path log_filename;
};

void load(corrector::corrector_config& cfg, const std::filesystem::path &filename);
}

typedef config_common::config<corrector::corrector_config> corr_cfg;
