//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <string>

// FIXME: Add fwd header
namespace debruijn_graph {
class DeBruijnGraph;
typedef DeBruijnGraph ConjugateDeBruijnGraph;
typedef ConjugateDeBruijnGraph Graph;
class GraphPack;

//Legacy wrappers
namespace graphio {

bool ScanBasicGraph(const std::string &file_name, Graph &g);
bool ScanGraphPack(const std::string &file_name, GraphPack &gp);

};

} // namespace debruijn_graph
