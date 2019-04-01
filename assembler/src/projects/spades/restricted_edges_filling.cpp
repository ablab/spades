//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "restricted_edges_filling.hpp"
#include "io/reads/file_reader.hpp"

#include <unordered_set>

namespace debruijn_graph {
static void FillRestrictedEdges(conj_graph_pack &gp) {
    if (!fs::check_existence(cfg::get().output_dir + "temp_anti/restricted_edges.fasta"))
        return;

    auto mapper = MapperInstance(gp);
    auto reader = io::FixingWrapper(io::FileReadStream(cfg::get().output_dir + "temp_anti/restricted_edges.fasta"));

    std::unordered_set<EdgeId> domain_edges;
    while (!reader.eof()) {
        io::SingleRead read;
        reader >> read;
        std::vector<EdgeId> edges = mapper->MapSequence(read.sequence()).simple_path();
        for (auto edge : edges) {
            domain_edges.insert(edge);
            domain_edges.insert(gp.g.conjugate(edge));
        }
    }

    DEBUG("Number of restricted edges - " << domain_edges.size());
    gp.emplace<SmartEdgeSet<std::unordered_set<EdgeId>, Graph>>(gp.g, domain_edges);
}

void RestrictedEdgesFilling::run(conj_graph_pack &gp, const char*) {
    FillRestrictedEdges(gp);
}

}
