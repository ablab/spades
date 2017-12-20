//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

/*
 * TruSeq Analysis Main
 */
#include "utils/logger/log_writers.hpp"
#include "utils/segfault_handler.hpp"
#include "utils/memory_limit.hpp"
#include "utils/filesystem/copy_file.hpp"
#include "pipeline/config_struct.hpp"
#include "analysis_pipeline.hpp"

void load_config(string cfg_filename) {
    fs::CheckFileExistenceFATAL(cfg_filename);

    cfg::create_instance(cfg_filename);

    if (!cfg::get().project_name.empty()) {
        make_dir(cfg::get().output_base + cfg::get().project_name);
    }

    make_dir(cfg::get().output_dir);
    make_dir(cfg::get().tmp_dir);

    if (cfg::get().developer_mode)
        make_dir(cfg::get().output_saves);

    make_dir(cfg::get().temp_bin_reads_path);
}

void create_console_logger(string cfg_filename) {
    using namespace logging;

    string log_props_file = cfg::get().log_filename;

    if (!fs::FileExists(log_props_file))
        log_props_file = fs::append_path(fs::parent_path(cfg_filename), cfg::get().log_filename);

    logger *lg = create_logger(fs::FileExists(log_props_file) ? log_props_file : "");
    lg->add_writer(std::make_shared<console_writer>());
    attach_logger(lg);
}

int main(int /*argc*/, char** argv) {
    utils::perf_counter pc;

    const size_t GB = 1 << 30;

    srand(42);
    srandom(42);

    try {
        using namespace debruijn_graph;

        string cfg_filename = argv[1];

        load_config          (cfg_filename);
        create_console_logger(cfg_filename);

        VERIFY(cfg::get().K >= runtime_k::MIN_K && cfg::get().K < runtime_k::MAX_K);
        VERIFY(cfg::get().K % 2 != 0);

        // read configuration file (dataset path etc.)

        utils::limit_memory(cfg::get().max_memory * GB);

        // assemble it!
        INFO("Assembling dataset (" << cfg::get().dataset_file << ") with K=" << cfg::get().K);

        spades::run_truseq_analysis();
    } catch (std::bad_alloc const& e) {
        std::cerr << "Not enough memory to run SPAdes. " << e.what() << std::endl;
        return EINTR;
    } catch (std::exception const& e) {
        std::cerr << "Exception caught " << e.what() << std::endl;
        return EINTR;
    } catch (...) {
        std::cerr << "Unknown exception caught " << std::endl;
        return EINTR;
    }

    unsigned ms = (unsigned)pc.time_ms();
    unsigned secs = (ms / 1000) % 60;
    unsigned mins = (ms / 1000 / 60) % 60;
    unsigned hours = (ms / 1000 / 60 / 60);
    INFO("Assembling time: " << hours << " hours " << mins << " minutes " << secs << " seconds");

    // OK
    return 0;
}
