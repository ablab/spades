//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "utils/path_helper.hpp"
#include "logger.hpp"

#include <iostream>

#include "config.hpp"

#include <iostream>

namespace logging {

struct console_writer : public writer {
#ifdef SPADES_USE_JEMALLOC

    void write_msg(double time, size_t cmem, size_t max_rss, level l, const char *file, size_t line_num,
                   const char *source, const char *msg) {
        std::cout << fmt::format("{:14s} {:>5s} / {:<5s} {:6.6s} {:24.24s} ({:26.26s}:{:4d})   {:s}",
                                 human_readable_time(time), human_readable_memory(cmem),
                                 human_readable_memory(max_rss), logging::level_name(l),
                                 source, path::filename(file), int(line_num), msg)
        << std::endl;
    }

#else
void write_msg(double time, size_t max_rss, level l, const char* file, size_t line_num, const char* source, const char* msg) {
  std::cout << fmt::format("{:14s} {:^5s} {:6.6s} {:24.24s} ({:26.26s}:{:4d})   {:s}",
                           human_readable_time(time), human_readable_memory(max_rss), logging::level_name(l),
                           source, path::filename(file), int(line_num), msg)
            << std::endl;
}
#endif
};

} // logging
