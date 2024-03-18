//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "log_writers_thread.hpp"
#include "utils/parallel/openmp_wrapper.h"

namespace logging {

void console_writer_thread::write_msg(double time, size_t cmem, size_t max_rss, level l, const std::filesystem::path& file, size_t line_num,
                                      const char *source, const char *msg) {
    int thread = omp_get_thread_num();
    if (cmem != -1ull)
        std::cout << fmt::format("thread #{:<2d} {:14s} {:>5s} / {:<5s} {:6.6s} {:24.24s} ({:26.26s}:{:4d})   {:s}",
                                 thread,
                                 utils::human_readable_time(time), utils::human_readable_memory(cmem),
                                 utils::human_readable_memory(max_rss), logging::level_name(l),
                                 source, file.filename().c_str(), int(line_num), msg)
                  << std::endl;
    else
        std::cout << fmt::format("thread #{:<2d} {:14s} {:^5s} {:6.6s} {:24.24s} ({:26.26s}:{:4d})   {:s}",
                                 thread,
                                 utils::human_readable_time(time), utils::human_readable_memory(max_rss),
                                 logging::level_name(l), source, file.filename().c_str(), int(line_num), msg)
                  << std::endl;
}

void file_writer_thread::write_msg(double time, size_t cmem, size_t max_rss, level l, const std::filesystem::path& file, size_t line_num,
                                   const char *source, const char *msg) {
    int thread = omp_get_thread_num();
    if (cmem != -1ull)
        fout << fmt::format("thread #{:<2d} {:14s} {:>5s} / {:<5s} {:6.6s} {:24.24s} ({:26.26s}:{:4d})   {:s}",
                            thread,
                            utils::human_readable_time(time), utils::human_readable_memory(cmem),
                            utils::human_readable_memory(max_rss), logging::level_name(l),
                            source, file.filename().c_str(), int(line_num), msg)
             << std::endl;
    else
        fout << fmt::format("thread #{:<2d} {:14s} {:^5s} {:6.6s} {:24.24s} ({:26.26s}:{:4d})   {:s}",
                            thread,
                            utils::human_readable_time(time), utils::human_readable_memory(max_rss),
                            logging::level_name(l), source, file.filename().c_str(), int(line_num), msg)
             << std::endl;
}

}
