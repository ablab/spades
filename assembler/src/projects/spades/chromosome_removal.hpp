//***************************************************************************
//* Copyright (c) 2015 Saint-Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include "pipeline/stage.hpp"
#include "assembly_graph/core/graph.hpp"
#include "assembly_graph/graph_support/coverage_uniformity_analyzer.hpp"

namespace debruijn_graph {

class ChromosomeRemoval : public spades::AssemblyStage {
public:
    ChromosomeRemoval()
            : AssemblyStage("Chromosome Removal", "chromosome_removal"), long_component_(), long_vertex_component_(),deadends_count_() { }

    void run(conj_graph_pack &gp, const char *);

private:
    std::unordered_map <EdgeId, size_t> long_component_;
    std::unordered_map <VertexId, size_t> long_vertex_component_;
    std::unordered_map <EdgeId, size_t> deadends_count_;

    size_t CalculateComponentSize(debruijn_graph::EdgeId e, Graph &g_);

    double RemoveLongGenomicEdges(conj_graph_pack &gp, size_t long_edge_bound, double coverage_limits,
                                  double external_chromosome_coverage = 0);
    void PlasmidSimplify(conj_graph_pack &gp, size_t long_edge_bound,
                                            std::function<void(typename Graph::EdgeId)> removal_handler = 0);
    void CompressAll(Graph &g);
    void DeleteAndCompress(EdgeId e, Graph &g);
};
}
