//***************************************************************************
//* Copyright (c) 2015-2019 Saint-Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include "assembly_graph/core/graph.hpp"
#include "pipeline/graph_pack.hpp"
#include "pipeline/config_struct.hpp"
#include <unordered_set>
#include <unordered_map>


namespace debruijn_graph {

class ChromosomeRemover {
public:
    ChromosomeRemover(GraphPack &gp, size_t ext_limit, config::debruijn_config::plasmid plasmid_config)
            : gp_(gp), ext_limit_(ext_limit), plasmid_config_(plasmid_config), chromosome_coverage_((double) ext_limit), long_component_(), long_vertex_component_(),
              deadends_count_(), component_list_(), full_name_(std::string("chromosome_removal") + (ext_limit == 0 ? std::string(""):std::to_string(ext_limit))) {
    }

    void run(GraphPack &gp, const char *);
    void RunMetaPipeline();
    void RunIsolatedPipeline();
    void FilterSmallComponents();
private:

    GraphPack& gp_;
    size_t ext_limit_;
    config::debruijn_config::plasmid plasmid_config_;
    double chromosome_coverage_;
    std::unordered_map <EdgeId, size_t> long_component_;
    std::unordered_map <VertexId, size_t> long_vertex_component_;
    std::unordered_map <EdgeId, size_t> deadends_count_;
    std::vector<std::vector<EdgeId>> component_list_;

    std::string full_name_;
    const size_t max_iteration_count = 30;

    size_t CalculateComponentSize(debruijn_graph::EdgeId e, const Graph &g);

    double RemoveLongGenomicEdges(size_t long_edge_bound, double coverage_limits,
                                  double external_chromosome_coverage = 0);
    void FillForbiddenSet(Graph &g, std::unordered_set<VertexId> &forbidden);
    void PlasmidSimplify(size_t long_edge_bound,
                         std::function<void(typename Graph::EdgeId)> removal_handler = 0);
    void CompressAll(Graph &g);

    void RemoveNearlyEverythingByCoverage(double limit);

    void CoverageFilter(double coverage_cutoff);
    std::vector<std::vector<EdgeId>> GetNineShapeComponents ();
    void OutputSuspiciousComponents ();

//  Not used, for debug purpose only
    void ReferenceBasedRemoveChromosomal();

    DECL_LOGGER("ChromosomeRemover");
};


} //debruijn_graph
