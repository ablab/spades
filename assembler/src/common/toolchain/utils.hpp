//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "utils/logger/log_writers.hpp"
#include "io/graph/gfa_reader.hpp"
#include "io/binary/graph_pack.hpp"

#pragma once

namespace toolchain {

static void create_console_logger() {
    using namespace logging;

    logger *lg = create_logger("");
    lg->add_writer(std::make_shared<console_writer>());
    attach_logger(lg);
}

static bool ends_with(const std::string &s, const std::string &p) {
    if (s.size() < p.size())
        return false;

    return (s.compare(s.size() - p.size(), p.size(), p) == 0);
}

static void PrintGraphInfo(debruijn_graph::Graph &g) {
    size_t sz = 0;
    for (auto it = g.ConstEdgeBegin(); !it.IsEnd(); ++it)
        sz += 1;

    INFO("Graph loaded. Total vertices: " << g.size() << " Total edges: " << sz);
}

static io::IdMapper<std::string> *LoadGraph(debruijn_graph::conj_graph_pack &gp, const std::string &filename) {
    io::IdMapper<std::string> *id_mapper = nullptr;
    if (ends_with(filename, ".gfa")) {
        id_mapper = new io::IdMapper<std::string>();
        gfa::GFAReader gfa(filename);
        INFO("GFA segments: " << gfa.num_edges() << ", links: " << gfa.num_links());
        gfa.to_graph(gp.g, id_mapper);
    } else {
        io::binary::BasePackIO<Graph>().Load(filename, gp);
    }
    PrintGraphInfo(gp.g);
    return id_mapper;
}

}

