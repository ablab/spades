//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2020-2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "utils.hpp"
#include "io/graph/gfa_reader.hpp"
#include "io/binary/graph_pack.hpp"

#include "utils/logger/logger.hpp"

namespace toolchain {

void create_console_logger(logging::level log_level) {
    using namespace logging;

    logger *lg = create_logger("", log_level);
    lg->add_writer(std::make_shared<console_writer>());
    attach_logger(lg);
}


io::IdMapper<std::string> *LoadGraphFromGFA(debruijn_graph::Graph &graph,
                                            const std::filesystem::path &filename) {
    io::IdMapper<std::string> *id_mapper = new io::IdMapper<std::string>();
    gfa::GFAReader gfa(filename);
    INFO("GFA segments: " << gfa.num_edges() << ", links: " << gfa.num_links());

    unsigned gk = graph.k();
    unsigned k = gfa.to_graph(graph, id_mapper);
    if (k == -1U)
        FATAL_ERROR("Failed to determine GFA k-mer length");
    if (k % 2 != 1)
        FATAL_ERROR("GFA used k-mer length must be odd (k=" << k << ")");
    if (k != gk)
        FATAL_ERROR("GFA used k-mer length (k=" << k << ") must match the command line settings (k=" << gk <<")");

    return id_mapper;
}


io::IdMapper<std::string> *LoadGraphPack(graph_pack::GraphPack &gp, const std::filesystem::path &filename) {
    auto &graph = gp.get_mutable<debruijn_graph::Graph>();
    io::IdMapper<std::string> *id_mapper = nullptr;
    if (filename.extension() == ".gfa") {
        id_mapper = LoadGraphFromGFA(graph, filename);
    } else {
        io::binary::BasePackIO().Load(filename, gp);
    }
    INFO("Graph loaded. Total vertices: " << graph.size() << " Total edges: " << graph.e_size());
    return id_mapper;
}

io::IdMapper<std::string> *LoadBaseGraph(debruijn_graph::Graph &graph,
                                         const std::filesystem::path &filename) {
    io::IdMapper<std::string> *id_mapper = nullptr;
    if (filename.extension() == ".gfa") {
        id_mapper = LoadGraphFromGFA(graph, filename);
    } else {
        io::binary::BasicGraphIO<debruijn_graph::Graph>().Load(filename, graph);
    }
    INFO("Graph loaded. Total vertices: " << graph.size() << " Total edges: " << graph.e_size());
    return id_mapper;
}

io::IdMapper<std::string> *LoadBaseGraph(graph_pack::GraphPack &gp, const std::filesystem::path &filename) {
    auto &graph = gp.get_mutable<debruijn_graph::Graph>();
    return LoadBaseGraph(graph, filename);
}


}

