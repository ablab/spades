//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once
#include "logger.hpp"

namespace logging
{

struct console_writer
    : public writer
{
    void write_msg(double time, level l, const char* file, size_t line_num, const char* source, const char* msg)
    {
        fs::path filepath(file);

        std::cout
            << str(boost::format("%14s %6s %-24.24s (%-26.26s:%4d)   %s")
                        % human_readable_time(time) % logging::level_name(l) % source % filepath.filename().c_str() % int(line_num) % msg)
            << std::endl;
    }
};

} // logging
