//***************************************************************************
//* Copyright (c) 2020 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "utils.hpp"
#include "utils/logger/log_writers.hpp"
#include "io/graph/gfa_reader.hpp"
#include "io/binary/graph_pack.hpp"

namespace toolchain {

void create_console_logger(logging::level log_level) {
    using namespace logging;

    logger *lg = create_logger("", log_level);
    lg->add_writer(std::make_shared<console_writer>());
    attach_logger(lg);
}

io::IdMapper<std::string> *LoadGraphPack(debruijn_graph::GraphPack &gp, const std::string &filename) {
    auto &graph = gp.get_mutable<debruijn_graph::Graph>();
    io::IdMapper<std::string> *id_mapper = nullptr;
    if (fs::extension(filename) == ".gfa") {
        id_mapper = new io::IdMapper<std::string>();
        gfa::GFAReader gfa(filename);
        INFO("GFA segments: " << gfa.num_edges() << ", links: " << gfa.num_links());

        if (gfa.k() == -1U)
            FATAL_ERROR("Failed to determine GFA k-mer length");
        if (gfa.k() % 2 != 1)
            FATAL_ERROR("GFA used k-mer length must be odd (k=" << gfa.k() << ")");
        if (gfa.k() != gp.k())
            FATAL_ERROR("GFA used k-mer length (k=" << gfa.k() << ") must match the command line settings (k=" << gp.k() <<")");

        gfa.to_graph(graph, id_mapper);
    } else {
        io::binary::BasePackIO().Load(filename, gp);
    }
    INFO("Graph loaded. Total vertices: " << graph.size() << " Total edges: " << graph.e_size());
    return id_mapper;
}

}

