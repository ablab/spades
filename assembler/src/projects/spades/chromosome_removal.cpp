//****************************************************************************
//***************************************************************************
//* Copyright (c) 2015-2019 Saint-Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "chromosome_removal.hpp"
#include "assembly_graph/core/basic_graph_stats.hpp"
#include "modules/chromosome_remover.hpp"


namespace debruijn_graph {
using namespace std;

void ChromosomeRemoval::run(GraphPack &gp, const char*) {
    const auto &graph = gp.get<Graph>();
    ChromosomeRemover remover(gp, ext_limit_, *cfg::get().pd);
    INFO("Starting chromosomal removal procedure, external coverage limit set to " << ext_limit_ <<", there are " << graph.size() << " vertices in graph");
    if (cfg::get().mode == config::pipeline_type::metaextrachromosomal ) {
        remover.RunMetaPipeline();
    } else {
        remover.RunIsolatedPipeline();
    }
    INFO("Counting average coverage after genomic edge removal");
    omnigraph::AvgCovereageCounter<Graph> cov_counter(graph);
    cfg::get_writable().ds.average_coverage = cov_counter.Count();
    INFO("Average coverage = " << cfg::get().ds.average_coverage);
}


} //debruijn_graph
