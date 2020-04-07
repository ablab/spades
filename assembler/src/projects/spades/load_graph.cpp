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
    VERIFY_MSG(gfa.k() != -1U, "Failed to determine k-mer length");
    VERIFY_MSG(gfa.k() % 2 == 1, "k-mer length must be odd");
    VERIFY_MSG(gfa.k() == gp.k(), "k-mer length must match the command line settings");
    gfa.to_graph(gp.get_mutable<Graph>());
}

}
