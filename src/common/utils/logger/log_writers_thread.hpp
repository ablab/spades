//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2015-2022 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "logger.hpp"

namespace logging {

struct console_writer_thread : public writer {
    void write_msg(double time, size_t cmem, size_t max_rss, level l, const std::filesystem::path& file, size_t line_num,
                   const char *source, const char *msg);
};

class file_writer_thread : public writer {
public:
    file_writer_thread(const std::string &filename) : fout(filename) {}

    void write_msg(double time, size_t cmem, size_t max_rss, level l, const std::filesystem::path& file, size_t line_num,
                   const char *source, const char *msg);

private:
    std::ofstream fout;
};

} // logging
