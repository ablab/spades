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

void create_console_logger(const string& dir) {
    using namespace logging;

    string log_props_file = corr_cfg::get().log_filename;

    if (!fs::FileExists(log_props_file))
        log_props_file = fs::append_path(dir, corr_cfg::get().log_filename);
    cout << log_props_file;
    logger *lg = create_logger(fs::FileExists(log_props_file) ? log_props_file : "");
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
            WARN("Wrong argument number");
            return 1;
        }
        string contig_name(argv[2]);
        string cfg_file(argv[1]);
        corr_cfg::create_instance(cfg_file);
        string cfg_dir = fs::parent_path(cfg_file);
        create_console_logger(cfg_dir);
        string work_dir = corr_cfg::get().work_dir;
        if (!fs::check_existence(corr_cfg::get().output_dir))
            fs::make_dir(corr_cfg::get().output_dir);
        if (!fs::check_existence(corr_cfg::get().work_dir))
            fs::make_dir(corr_cfg::get().work_dir);

        INFO("Starting MismatchCorrector, built from " SPADES_GIT_REFSPEC ", git revision " SPADES_GIT_SHA1);
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
