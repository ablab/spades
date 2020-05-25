//***************************************************************************
//* Copyright (c) 2020 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "load_graph.hpp"
#include "io/graph/gfa_reader.hpp"
#include "utils/filesystem/path_helper.hpp"
#include "utils/logger/logger.hpp"

#include <iterator>

namespace debruijn_graph {

void LoadGraph::run(GraphPack &gp, const char*) {
    std::string path;
    for (const auto &lib : cfg::get().ds.reads.libraries()) {
        if (lib.is_assembly_graph()) {
            path = *lib.single_begin();
            break;
        }
    }
    fs::CheckFileExistenceFATAL(path);
    
    gfa::GFAReader gfa(path);
    INFO("GFA segments: " << gfa.num_edges() << ", links: " << gfa.num_links());
    if (gfa.k() == -1U)
        FATAL_ERROR("Failed to determine GFA k-mer length");
    if (gfa.k() % 2 != 1)
        FATAL_ERROR("GFA used k-mer length must be odd (k=" << gfa.k() << ")");
    if (gfa.k() != gp.k())
        FATAL_ERROR("GFA used k-mer length (k=" << gfa.k() << ") must match the command line settings (k=" << gp.k() <<")");
    gfa.to_graph(gp.get_mutable<Graph>());
}

}
