//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

// just to check that headers from include and debruijn folders are correctly included
#include "vis_logger.hpp"

#include "standard_vis.hpp"
#include "utils/segfault_handler.hpp"
#include "utils/stacktrace.hpp"
#include "pipeline/config_struct.hpp"
#include "io/reads/io_helper.hpp"
#include "utils/stl_utils.hpp"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include "utils/memory_limit.hpp"
#include "io/dataset_support/read_converter.hpp"

#include "debruijn_online_visualizer.hpp"

void create_console_logger(string const& cfg_filename) {
    using namespace logging;

    string log_props_file = cfg::get().log_filename;

    if (!fs::FileExists(log_props_file))
        log_props_file = fs::append_path(fs::parent_path(cfg_filename), cfg::get().log_filename);

    logger *lg = create_logger(fs::FileExists(log_props_file) ? log_props_file : "");
    lg->add_writer(std::make_shared<console_writer>());

    attach_logger(lg);
}

int main(int argc, char** argv) {
    const size_t GB = 1 << 30;

    try {
        VERIFY(argc > 1)
        using namespace online_visualization;
        string cfg_filename = argv[1];
        fs::CheckFileExistenceFATAL(cfg_filename);

        cfg::create_instance(cfg_filename);

        VERIFY(cfg::get().K >= runtime_k::MIN_K && cfg::get().K < runtime_k::MAX_K);
        VERIFY(cfg::get().K % 2 != 0);

        create_console_logger(cfg_filename);
        cout << "\nGAF (Graph Analysis Framework) started" << endl;
        cout << "Print help to see readme file" << endl;
        utils::limit_memory(cfg::get().max_memory * GB);

        DebruijnOnlineVisualizer online_vis;
        online_vis.init();
        online_vis.run();
    }
    catch (std::exception const& e) {
        std::cerr << "Exception caught " << e.what() << std::endl;
        return EINTR;
    }
    catch (...) {
        std::cerr << "Unknown exception caught " << std::endl;
        return EINTR;
    }
    return 0;
}
