//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

/*
 * Assembler Main
 */
#include "utils/logger/log_writers.hpp"

#include "utils/memory_limit.hpp"
#include "utils/segfault_handler.hpp"
#include "launch.hpp"
#include "utils/copy_file.hpp"
#include "version.hpp"

void load_config(const vector<string>& cfg_fns) {
    for (const auto& s : cfg_fns) {
        path::CheckFileExistenceFATAL(s);
    }

    cfg::create_instance(cfg_fns);

    if (!cfg::get().project_name.empty()) {
        make_dir(cfg::get().output_base + cfg::get().project_name);
    }

    make_dir(cfg::get().output_dir);
    make_dir(cfg::get().tmp_dir);

    if (cfg::get().developer_mode)
        make_dir(cfg::get().output_saves);

    make_dir(cfg::get().temp_bin_reads_path);
}

void create_console_logger(const string& dir) {
    using namespace logging;

    string log_props_file = cfg::get().log_filename;

    if (!path::FileExists(log_props_file))
        log_props_file = path::append_path(dir, cfg::get().log_filename);

    logger *lg = create_logger(path::FileExists(log_props_file) ? log_props_file : "");
    lg->add_writer(std::make_shared<console_writer>());
    attach_logger(lg);
}

int main(int argc, char **argv) {
    perf_counter pc;

    const size_t GB = 1 << 30;

    srand(42);
    srandom(42);

    try {
        using namespace debruijn_graph;

        string cfg_dir = path::parent_path(argv[1]);

        vector<string> cfg_fns;
        for (int i = 1; i < argc; ++i) {
           cfg_fns.push_back(argv[i]);
        }

        load_config(cfg_fns);

        create_console_logger(cfg_dir);

        for (const auto& cfg_fn : cfg_fns)
            INFO("Loading config from " << cfg_fn);

        VERIFY(cfg::get().K >= runtime_k::MIN_K && cfg::get().K < runtime_k::MAX_K);
        VERIFY(cfg::get().K % 2 != 0);

        // read configuration file (dataset path etc.)

        limit_memory(cfg::get().max_memory * GB);

        // assemble it!
        INFO("Starting SPAdes, built from "
                     SPADES_GIT_REFSPEC
                     ", git revision "
                     SPADES_GIT_SHA1);
        INFO("Maximum k-mer length: " << runtime_k::MAX_K);
        INFO("Assembling dataset (" << cfg::get().dataset_file << ") with K=" << cfg::get().K);

        spades::assemble_genome();

    } catch (std::bad_alloc const &e) {
        std::cerr << "Not enough memory to run SPAdes. " << e.what() << std::endl;
        return EINTR;
    } catch (std::exception const &e) {
        std::cerr << "Exception caught " << e.what() << std::endl;
        return EINTR;
    } catch (...) {
        std::cerr << "Unknown exception caught " << std::endl;
        return EINTR;
    }

    unsigned ms = (unsigned) pc.time_ms();
    unsigned secs = (ms / 1000) % 60;
    unsigned mins = (ms / 1000 / 60) % 60;
    unsigned hours = (ms / 1000 / 60 / 60);
    INFO("Assembling time: " << hours << " hours " << mins << " minutes " << secs << " seconds");

    // OK
    return 0;
}
