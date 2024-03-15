//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "dataset_processor.hpp"

#include "utils/logger/log_writers.hpp"
#include "utils/segfault_handler.hpp"

#include "version.hpp"

#include "config_struct.hpp"

#include <sys/types.h>
#include <sys/stat.h>
#include <string>

using namespace std;

void create_console_logger(const filesystem::path& dir) {
    using namespace logging;

    filesystem::path log_props_file = corr_cfg::get().log_filename;

    if (!exists(log_props_file))
        log_props_file = dir / corr_cfg::get().log_filename;
    cout << log_props_file;
    logger *lg = create_logger(exists(log_props_file) ? log_props_file : "");
    lg->add_writer(std::make_shared<console_writer>());
    //lg->add_writer(std::make_shared<mutex_writer>(std::make_shared<console_writer>()));
    attach_logger(lg);
}

int main(int argc, char** argv) {
    utils::perf_counter pc;

    srand(42);
    srandom(42);
    try {
        if (argc != 3) {
            std::cerr << "Wrong argument number\n";
            return 1;
        }
        string contig_name(argv[2]);
        filesystem::path cfg_file(argv[1]);
        corr_cfg::create_instance(cfg_file);
        create_console_logger(cfg_file.parent_path());
        if (!exists(corr_cfg::get().output_dir))
            create_directory(corr_cfg::get().output_dir);
        if (!exists(corr_cfg::get().work_dir))
            create_directory(corr_cfg::get().work_dir);

        START_BANNER("mismatch corrector");
        INFO("Maximum # of threads to use (adjusted due to OMP capabilities): " << corr_cfg::get().max_nthreads);

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
