//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "io/reads/vector_reader.hpp"
#include "io/dataset_support/dataset_readers.hpp"
#include "pipeline/graph_pack.hpp"
#include "io/dataset_support/read_converter.hpp"

#include "modules/graph_construction.hpp"
#include "assembly_graph/stats/picture_dump.hpp"
#include "construction.hpp"

namespace debruijn_graph {

template<class Read>
void construct_graph(io::ReadStreamList<Read>& streams,
                     conj_graph_pack& gp, io::SingleStreamPtr contigs_stream = io::SingleStreamPtr()) {
    config::debruijn_config::construction params = cfg::get().con;
    params.early_tc.enable &= !cfg::get().gap_closer_enable;

    ReadStatistics stats = ConstructGraphWithCoverage(params, streams, gp.g,
                                                      gp.index, gp.flanking_cov, contigs_stream);
    size_t rl = stats.max_read_length_;

    if (!cfg::get().ds.RL()) {
        INFO("Figured out: read length = " << rl);
        cfg::get_writable().ds.set_RL(rl);
        cfg::get_writable().ds.set_aRL((double) stats.bases_ / (double) stats.reads_);
    } else if (cfg::get().ds.RL() != rl)
        WARN("In datasets.info, wrong RL is specified: " << cfg::get().ds.RL() << ", not " << rl);
}

void Construction::run(conj_graph_pack &gp, const char*) {
    // Has to be separate stream for not counting it in coverage
    io::ReadStreamList<io::SingleRead> trusted_contigs;
    if (cfg::get().use_additional_contigs) {
        DEBUG("Contigs from previous K will be used: " << cfg::get().additional_contigs);
        trusted_contigs.push_back(io::EasyStream(cfg::get().additional_contigs, true));
    }

    bool trusted_contigs_exist = false;
    for (const auto& lib : cfg::get().ds.reads) {
        if (lib.type() != io::LibraryType::TrustedContigs)
            continue;

        for (const auto& read : lib.single_reads()) {
            trusted_contigs.push_back(io::EasyStream(read, true));
            trusted_contigs_exist = true;
        }
    }

    if (trusted_contigs_exist)
        INFO("Trusted contigs will be used in graph construction");
    auto contigs_stream = MultifileWrap(trusted_contigs);

    auto& dataset = cfg::get_writable().ds;
    std::vector<size_t> libs_for_construction;
    for (size_t i = 0; i < dataset.reads.lib_count(); ++i)
        if (dataset.reads[i].is_graph_contructable())
            libs_for_construction.push_back(i);

    auto streams = io::single_binary_readers_for_libs(dataset, libs_for_construction, true, true);
    construct_graph<io::SingleReadSeq>(streams, gp, contigs_stream);
}

} //namespace debruijn_graph
