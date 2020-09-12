//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "utils/logger/log_writers.hpp"
#include <string>

#pragma once

namespace io {
template<class T>
class IdMapper;
}

namespace debruijn_graph {
class GraphPack;
}

namespace toolchain {

void create_console_logger(logging::level log_level = logging::L_INFO);
io::IdMapper<std::string> *LoadGraphPack(debruijn_graph::GraphPack &gp, const std::string &filename);

}

