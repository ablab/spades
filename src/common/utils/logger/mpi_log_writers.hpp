//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "logger.hpp"

namespace logging {

struct mpi_console_writer : public writer {

    void write_msg(double time, size_t cmem, size_t max_rss, level l, const std::filesystem::path& file, size_t line_num,
                   const char *source, const char *msg);
private:
    std::string nodeinfo() const;
};

} // namespace logging
