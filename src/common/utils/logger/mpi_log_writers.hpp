//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "utils/filesystem/path_helper.hpp"
#include "logger.hpp"

#include <iostream>

#include "config.hpp"

#include <iostream>
#include <mutex>
#include <mpi.h>

namespace logging {

struct mpi_console_writer : public writer {

    void write_msg(double time, size_t cmem, size_t max_rss, level l, const char *file, size_t line_num,
                   const char *source, const char *msg) {
        const std::string node_info = nodeinfo();
        if (cmem != -1ull)
            std::cout << fmt::format("NODE {:s} | {:14s} {:>5s} / {:<5s} {:6.6s} {:24.24s} ({:26.26s}:{:4d})   {:s}",
                                     node_info,
                                     utils::human_readable_time(time), utils::human_readable_memory(cmem),
                                     utils::human_readable_memory(max_rss), logging::level_name(l),
                                     source, fs::filename(file), int(line_num), msg)
            << std::endl;
        else
            std::cout << fmt::format("NODE {:s} | {:14s} {:^5s} {:6.6s} {:24.24s} ({:26.26s}:{:4d})   {:s}",
                                     node_info,
                                     utils::human_readable_time(time), utils::human_readable_memory(max_rss),
                                     logging::level_name(l), source, fs::filename(file), int(line_num), msg)
                      << std::endl;
    }

private:
    std::string nodeinfo() const {
        int initialized, finalized;
        MPI_Initialized(&initialized);
        MPI_Finalized(&finalized);
        if (initialized && !finalized) {
            int world_rank, world_size;
            MPI_Comm_size(MPI_COMM_WORLD, &world_size);
            MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
            return fmt::format("{:>2d}/{:<2d}", world_rank, world_size);
        } else {
            return fmt::format("{:^5}", "N/A");
        }

    }
};

} // namespace logging
