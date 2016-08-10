//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "dataset_processor.hpp"
#include "pipeline/config_struct.hpp"

#include "utils/logger/log_writers.hpp"
#include "config_struct.hpp"
#include "utils/segfault_handler.hpp"

#include "version.hpp"

#include <sys/types.h>
#include <sys/stat.h>
#include <string>

using namespace std;
void create_console_logger() {
    using namespace logging;

    logger *lg = create_logger("");
    lg->add_writer(std::make_shared<console_writer>());
    attach_logger(lg);
}

int main(int argc, char** argv) {
    perf_counter pc;

    srand(42);
    srandom(42);
    try {
        create_console_logger();

        if (argc != 3) {
            WARN("Wrong argument number");
            return 1;
        }
        string contig_name(argv[2]);
        string cfg_file(argv[1]);
        corr_cfg::create_instance(cfg_file);
        string work_dir = corr_cfg::get().work_dir;
        if (!path::check_existence(corr_cfg::get().output_dir))
            path::make_dir(corr_cfg::get().output_dir);
        if (!path::check_existence(corr_cfg::get().work_dir))
            path::make_dir(corr_cfg::get().work_dir);

        INFO("Starting MismatchCorrector, built from " SPADES_GIT_REFSPEC ", git revision " SPADES_GIT_SHA1);

        corrector::DatasetProcessor dp(contig_name, corr_cfg::get().work_dir, corr_cfg::get().output_dir, corr_cfg::get().max_nthreads);
        dp.ProcessDataset();
    } catch (std::string const &s) {
        std::cerr << s;
        return EINTR;
    }
    unsigned ms = (unsigned) pc.time_ms();
    unsigned secs = (ms / 1000) % 60;
    unsigned mins = (ms / 1000 / 60) % 60;
    unsigned hours = (ms / 1000 / 60 / 60);

    INFO("Correcting time: " << hours << " hours " << mins << " minutes " << secs << " seconds");

    return 0;
}
