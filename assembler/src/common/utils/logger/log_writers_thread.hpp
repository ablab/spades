//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "utils/filesystem/path_helper.hpp"
#include "logger.hpp"

#include <iostream>
#include <fstream>

#include "config.hpp"

#include <iostream>
#include <mutex>
#include "utils/parallel/openmp_wrapper.h"

namespace logging {

struct console_writer_thread : public writer {

    void write_msg(double time, size_t cmem, size_t max_rss, level l, const char *file, size_t line_num,
                   const char *source, const char *msg) {
        int thread = omp_get_thread_num();
        if (cmem != -1ull)
            std::cout << fmt::format("thread #{:<2d} {:14s} {:>5s} / {:<5s} {:6.6s} {:24.24s} ({:26.26s}:{:4d})   {:s}",
                                     thread,
                                     utils::human_readable_time(time), utils::human_readable_memory(cmem),
                                     utils::human_readable_memory(max_rss), logging::level_name(l),
                                     source, fs::filename(file), int(line_num), msg)
            << std::endl;
        else
            std::cout << fmt::format("thread #{:<2d} {:14s} {:^5s} {:6.6s} {:24.24s} ({:26.26s}:{:4d})   {:s}",
                                     thread,
                                     utils::human_readable_time(time), utils::human_readable_memory(max_rss),
                                     logging::level_name(l), source, fs::filename(file), int(line_num), msg)
                      << std::endl;
    }

};

class file_writer_thread : public writer {
public:
    file_writer_thread(const std::string &filename) : fout(filename) {}

    void write_msg(double time, size_t cmem, size_t max_rss, level l, const char *file, size_t line_num,
                   const char *source, const char *msg) {
        int thread = omp_get_thread_num();
        if (cmem != -1ull)
            fout << fmt::format("thread #{:<2d} {:14s} {:>5s} / {:<5s} {:6.6s} {:24.24s} ({:26.26s}:{:4d})   {:s}",
                                thread,
                                utils::human_readable_time(time), utils::human_readable_memory(cmem),
                                utils::human_readable_memory(max_rss), logging::level_name(l),
                                source, fs::filename(file), int(line_num), msg)
            << std::endl;
        else
            fout << fmt::format("thread #{:<2d} {:14s} {:^5s} {:6.6s} {:24.24s} ({:26.26s}:{:4d})   {:s}",
                                thread,
                                utils::human_readable_time(time), utils::human_readable_memory(max_rss),
                                logging::level_name(l), source, fs::filename(file), int(line_num), msg)
                << std::endl;
    }

private:
    std::ofstream fout;
};

} // logging
