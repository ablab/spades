//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "logger.hpp"
#include <mutex>

namespace logging {

struct console_writer : public writer {
    void write_msg(double time, size_t cmem, size_t max_rss, level l, const std::filesystem::path& file, size_t line_num,
                   const char *source, const char *msg);
};

class mutex_writer : public writer {
    std::mutex writer_mutex_;
    std::shared_ptr<writer> writer_;
public:
    mutex_writer(std::shared_ptr<writer> writer) : writer_(writer) {}

    void write_msg(double time, size_t cmem, size_t max_rss, level l, const std::filesystem::path& file, size_t line_num,
                   const char *source, const char *msg);
};

} // logging
