//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <string>

namespace graph_pack {
class GraphPack;
}

// FIXME: Add fwd header
namespace debruijn_graph {
class DeBruijnGraph;
typedef DeBruijnGraph ConjugateDeBruijnGraph;
typedef ConjugateDeBruijnGraph Graph;

//Legacy wrappers
namespace graphio {

bool ScanBasicGraph(const std::string &file_name, Graph &g);
bool ScanGraphPack(const std::string &file_name, graph_pack::GraphPack &gp);

};

} // namespace debruijn_graph
