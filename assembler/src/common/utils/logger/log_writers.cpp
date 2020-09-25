//***************************************************************************
//* Copyright (c) 2020 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "log_writers.hpp"
#include "config.hpp"

#include <filesystem>
#include <iostream>

namespace logging {

void console_writer::write_msg(double time, size_t cmem, size_t max_rss, level l, const std::filesystem::path& file, size_t line_num,
                               const char *source, const char *msg) {
    if (cmem != -1ull)
        std::cout << fmt::format("{:14s} {:>5s} / {:<5s} {:6.6s} {:24.24s} ({:26.26s}:{:4d})   {:s}",
                                 utils::human_readable_time(time), utils::human_readable_memory(cmem),
                                 utils::human_readable_memory(max_rss), logging::level_name(l),
                                 source, file.filename().c_str(), int(line_num), msg)
                  << std::endl;
    else
        std::cout << fmt::format("{:14s} {:^5s} {:6.6s} {:24.24s} ({:26.26s}:{:4d})   {:s}",
                                 utils::human_readable_time(time), utils::human_readable_memory(max_rss),
                                 logging::level_name(l), source,
                                 file.filename().c_str(),
                                 int(line_num), msg)
                  << std::endl;
}

void mutex_writer::write_msg(double time, size_t cmem, size_t max_rss, level l, const std::filesystem::path& file, size_t line_num,
                             const char *source, const char *msg) {
    std::lock_guard<std::mutex> guard(writer_mutex_);
    writer_->write_msg(time, cmem, max_rss, l, file, line_num, source, msg);
}

} // logging

