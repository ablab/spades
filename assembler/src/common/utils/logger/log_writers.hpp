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

#include "config.hpp"

#include <iostream>
#include <mutex>

namespace logging {

struct console_writer : public writer {

    void write_msg(double time, size_t cmem, size_t max_rss, level l, const char *file, size_t line_num,
                   const char *source, const char *msg) {
        if (cmem != -1ull)
            std::cout << fmt::format("{:14s} {:>5s} / {:<5s} {:6.6s} {:24.24s} ({:26.26s}:{:4d})   {:s}",
                                     utils::human_readable_time(time), utils::human_readable_memory(cmem),
                                     utils::human_readable_memory(max_rss), logging::level_name(l),
                                     source, fs::filename(file), int(line_num), msg)
            << std::endl;
        else
            std::cout << fmt::format("{:14s} {:^5s} {:6.6s} {:24.24s} ({:26.26s}:{:4d})   {:s}",
                                     utils::human_readable_time(time), utils::human_readable_memory(max_rss),
                                     logging::level_name(l), source, fs::filename(file), int(line_num), msg)
                      << std::endl;
    }

};

class mutex_writer : public writer {
    std::mutex writer_mutex_;
    std::shared_ptr<writer> writer_;

public:

    mutex_writer(std::shared_ptr<writer> writer) : writer_(writer) {}

    void write_msg(double time, size_t cmem, size_t max_rss, level l, const char *file, size_t line_num,
                   const char *source, const char *msg) override {
        std::lock_guard<std::mutex> guard(writer_mutex_);
        writer_->write_msg(time, cmem, max_rss, l, file, line_num, source, msg);
    }
};

} // logging
