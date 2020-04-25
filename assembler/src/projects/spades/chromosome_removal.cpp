//****************************************************************************
//***************************************************************************
//* Copyright (c) 2015-2019 Saint-Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "assembly_graph/core/basic_graph_stats.hpp"
#include "chromosome_removal.hpp"
#include "modules/chromosome_remover.hpp"


namespace debruijn_graph {
using namespace std;

void ChromosomeRemoval::run(GraphPack &gp, const char*) {
    ChromosomeRemover remover(gp, ext_limit_, *cfg::get().pd);
    INFO("Starting chromosomal removal procedure, external coverage limit set to " << ext_limit_ <<", there are " << gp.g.size() << " vertices in graph");
    if (cfg::get().mode == config::pipeline_type::metaplasmid) {
        remover.RunMetaPipeline();
    } else {
        remover.RunIsolatedPipeline();
    }
    remover.FilterSmallComponents();

    INFO("Counting average coverage after genomic edge removal");
    omnigraph::AvgCovereageCounter<Graph> cov_counter(gp.get<Graph>());
    cfg::get_writable().ds.average_coverage = cov_counter.Count();
    INFO("Average coverage = " << cfg::get().ds.average_coverage);
}


} //debruijn_graph
