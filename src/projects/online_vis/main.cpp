//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2015-2022 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

// just to check that headers from include and debruijn folders are correctly included
#include "vis_logger.hpp"

#include "standard_vis.hpp"
#include "utils/segfault_handler.hpp"
#include "utils/stacktrace.hpp"
#include "configs/config_struct.hpp"
#include "io/reads/io_helper.hpp"
#include "utils/stl_utils.hpp"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include "utils/memory_limit.hpp"
#include "io/dataset_support/read_converter.hpp"

#include "debruijn_online_visualizer.hpp"

void create_console_logger(filesystem::path const& cfg_filename) {
    using namespace logging;

    filesystem::path log_props_file = cfg::get().log_filename;

    if (!exists(log_props_file))
        log_props_file = cfg_filename.parent_path() / cfg::get().log_filename;

    logger *lg = create_logger(exists(log_props_file) ? log_props_file : "");
    lg->add_writer(std::make_shared<console_writer>());

    attach_logger(lg);
}

int main(int argc, char** argv) {
    const size_t GB = 1 << 30;

    try {
        VERIFY(argc > 1)
        using namespace online_visualization;
        filesystem::path cfg_filename = argv[1];
        CHECK_FATAL_ERROR(exists(cfg_filename), "File " << cfg_filename << " doesn't exist or can't be read!");

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
