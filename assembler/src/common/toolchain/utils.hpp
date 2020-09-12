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
// FIXME: Make fwd
class DeBruijnGraph;
typedef DeBruijnGraph ConjugateDeBruijnGraph;
typedef ConjugateDeBruijnGraph Graph;

class GraphPack;
}

namespace toolchain {

void create_console_logger(logging::level log_level = logging::L_INFO);
io::IdMapper<std::string> *LoadGraphFromGFA(debruijn_graph::Graph &graph,
                                            const std::string &filename);
io::IdMapper<std::string> *LoadGraphPack(debruijn_graph::GraphPack &gp, const std::string &filename);
io::IdMapper<std::string> *LoadBaseGraph(debruijn_graph::Graph &g, const std::string &filename);
io::IdMapper<std::string> *LoadBaseGraph(debruijn_graph::GraphPack &gp, const std::string &filename);


}

