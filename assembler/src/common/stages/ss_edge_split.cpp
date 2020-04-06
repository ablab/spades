//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "ss_edge_split.hpp"

#include <modules/alignment/sequence_mapper_notifier.hpp>
#include <modules/alignment/sequence_mapper.hpp>
#include <modules/alignment/rna/ss_coverage_filler.hpp>
#include <io/dataset_support/read_converter.hpp>

namespace debruijn_graph {

using namespace config;
void SSEdgeSplit::run(GraphPack& gp, const char *) {
    using namespace omnigraph;

    if (!cfg::get().ss.ss_enabled) {
        INFO("Dataset is not strand-specific, strand-specific edge splitter will not run");
        return;
    }

    gp.EnsureBasicMapping();

    for (size_t i = 0; i < cfg::get().ds.reads.lib_count(); ++i) {
        auto &lib = cfg::get_writable().ds.reads[i];
        if (!lib.is_graph_constructable()) {
            continue;
        }

        auto& reads = cfg::get_writable().ds.reads[i];
        if (reads.data().unmerged_read_length < gp.k()) {
            INFO("Reads are too short for SS coverage splitter");
            continue;
        }

        INFO("Running strand-specific edge splitter only for library # " << i);
        SequenceMapperNotifier notifier(gp, cfg::get().ds.reads.lib_count());
        const auto& params = cfg::get().ss_coverage_splitter;
        SSCoverageSplitter splitter(gp.get_mutable<Graph>(), params.bin_size, params.min_edge_len,
            params.min_edge_coverage, params.coverage_margin, params.min_flanking_coverage);
        SSBinCoverageFiller ss_coverage_filler(splitter);
        notifier.Subscribe(i, &ss_coverage_filler);

        INFO("Selecting usual mapper");
        auto mapper_ptr = MapperInstance(gp);
        auto single_streams = io::single_binary_readers(reads, /*followed_by_rc*/ false, /*map_paired*/true);
        notifier.ProcessLibrary(single_streams, i, *mapper_ptr);

        auto &index = gp.get_mutable<EdgeIndex<Graph>>();
        if (index.IsAttached())
            index.Detach();
        index.clear();

        splitter.SplitEdges();
        break;
    }
}

};
