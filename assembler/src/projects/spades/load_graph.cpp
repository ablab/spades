//***************************************************************************
//* Copyright (c) 2020 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "load_graph.hpp"
#include "io/graph/gfa_reader.hpp"
#include "utils/logger/logger.hpp"

#include <iterator>

namespace debruijn_graph {

void LoadGraph::run(graph_pack::GraphPack &gp, const char*) {
    std::filesystem::path path;
    for (const auto &lib : cfg::get().ds.reads.libraries()) {
        if (lib.is_assembly_graph()) {
            path = *lib.single_begin();
            break;
        }
    }
    CHECK_FATAL_ERROR(exists(path), "File " << path << " doesn't exist or can't be read!");
    
    gfa::GFAReader gfa(path);
    unsigned k = gfa.to_graph(gp.get_mutable<Graph>());
    INFO("GFA segments: " << gfa.num_edges() << ", links: " << gfa.num_links());
    if (k == -1U)
        FATAL_ERROR("Failed to determine GFA k-mer length");
    if (k % 2 != 1)
        FATAL_ERROR("GFA used k-mer length must be odd (k=" << k << ")");
    if (k != gp.k())
        FATAL_ERROR("GFA used k-mer length (k=" << k << ") must match the command line settings (k=" << gp.k() <<")");
}

}
