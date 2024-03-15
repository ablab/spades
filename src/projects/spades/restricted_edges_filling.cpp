//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "restricted_edges_filling.hpp"
#include "pipeline/graph_pack_helpers.h"
#include "alignment/sequence_mapper.hpp"
#include "assembly_graph/core/graph_iterators.hpp"
#include "io/reads/file_reader.hpp"
#include "io/reads/wrapper_collection.hpp"
#include "pipeline/sequence_mapper_gp_api.hpp"
#include "assembly_graph/core/graph_iterators.hpp"
#include "io/dataset_support/dataset_readers.hpp"

#include <unordered_set>

namespace debruijn_graph {

static void MapRestrictedEdgesFromFile(graph_pack::GraphPack &gp, const std::string &filename) {
    auto mapper = MapperInstance(gp);
    auto reader = io::FixingWrapper(io::FileReadStream(filename));

    const auto &graph = gp.get<Graph>();
    std::unordered_set<EdgeId> domain_edges;
    while (!reader.eof()) {
        io::SingleRead read;
        reader >> read;
        std::vector<EdgeId> edges = mapper->MapSequence(read.sequence()).simple_path();
        for (auto edge : edges) {
            domain_edges.insert(edge);
            domain_edges.insert(graph.conjugate(edge));
        }
    }

    DEBUG("Number of restricted edges - " << domain_edges.size());
    gp.add("restricted_edges", omnigraph::SmartEdgeSet<std::unordered_set<EdgeId>, Graph>(graph, domain_edges));
}

static void MapRestrictedEdgesFromTrustedContigs(graph_pack::GraphPack &gp) {
    if (!gp.get_mutable<KmerMapper<Graph>>().IsAttached())
        gp.get_mutable<KmerMapper<Graph>>().Attach();

    EnsureBasicMapping(gp);
    std::vector<size_t> trusted_contigs;
    for (size_t lib_id = 0; lib_id < cfg::get().ds.reads.lib_count(); ++lib_id) {
        if (cfg::get().ds.reads[lib_id].type() == io::LibraryType::TrustedContigs)
            trusted_contigs.push_back(lib_id);
    }

    std::unordered_set<EdgeId> restricted_edges;
    auto mapper = MapperInstance(gp);
    io::SingleRead reference;
    const auto &graph = gp.get<Graph>();

    for (auto lib_id : trusted_contigs) {
        auto stream = io::single_easy_reader(cfg::get().ds.reads[lib_id], false, false);
        while(!stream.eof()) {
            stream >> reference;
            std::vector<EdgeId> edges = mapper->MapSequence(reference.sequence()).simple_path();
            for (auto edge : edges) {
                restricted_edges.insert(edge);
                restricted_edges.insert(graph.conjugate(edge));
            }
        }
    }

    INFO("Number of restricted edges - " << restricted_edges.size());
    gp.add("restricted_edges", omnigraph::SmartEdgeSet<std::unordered_set<EdgeId>, Graph>(graph, restricted_edges));
}

static void FillRestrictedEdges(graph_pack::GraphPack &gp) {
    if (cfg::get().sewage)
        MapRestrictedEdgesFromTrustedContigs(gp);

    if (!exists(cfg::get().output_dir / "temp_anti/restricted_edges.fasta"))
        return;
    MapRestrictedEdgesFromFile(gp, cfg::get().output_dir / "temp_anti/restricted_edges.fasta");
}

void RestrictedEdgesFilling::run(graph_pack::GraphPack &gp, const char*) {
    FillRestrictedEdges(gp);
}

}
