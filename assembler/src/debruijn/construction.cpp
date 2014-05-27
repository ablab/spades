//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "standard.hpp"
#include "io/easy_reader.hpp"
#include "io/vector_reader.hpp"
#include "dataset_readers.hpp"
#include "graph_pack.hpp"
#include "read_converter.hpp"
#include "omni/visualization/graph_labeler.hpp"

#include "graph_construction.hpp"
#include "stats/debruijn_stats.hpp"
#include "positions.hpp"
#include "construction.hpp"

namespace debruijn_graph {

template<class Read>
void construct_graph(io::ReadStreamList<Read>& streams,
                     conj_graph_pack& gp, io::SingleStreamPtr contigs_stream = io::SingleStreamPtr()) {
    debruijn_config::construction params = cfg::get().con;
    params.early_tc.enable &= !cfg::get().gap_closer_enable;

    size_t rl = ConstructGraphWithCoverage(params, streams, gp.g,
                                           gp.index, gp.flanking_cov, contigs_stream);

    if (!cfg::get().ds.RL()) {
        INFO("Figured out: read length = " << rl);
        cfg::get_writable().ds.set_RL(rl);
    } else if (cfg::get().ds.RL() != rl)
        WARN("In datasets.info, wrong RL is specified: " << cfg::get().ds.RL() << ", not " << rl);
}

void Construction::run(conj_graph_pack &gp, const char*) {
    // Has to be separate stream for not counting it in coverage
    io::ReadStreamList<io::SingleRead> trusted_contigs;
    if (cfg::get().use_additional_contigs) {
        INFO("Contigs from previous K will be used");
        trusted_contigs.push_back(io::EasyStream(cfg::get().additional_contigs, true));
    }

    bool trusted_contigs_exist = false;
    for (size_t i = 0; i < cfg::get().ds.reads.lib_count(); ++i) {
        if (cfg::get().ds.reads[i].type() == io::LibraryType::TrustedContigs) {
            for (auto it = cfg::get().ds.reads[i].single_begin();
                    it != cfg::get().ds.reads[i].single_end();
                    ++it) {
                trusted_contigs.push_back(io::EasyStream(*it, true));
                trusted_contigs_exist = true;
            }
        }
    }
    if (trusted_contigs_exist)
        INFO("Trusted contigs will be used in graph construction");
    auto contigs_stream = MultifileWrap(trusted_contigs);

    std::vector<size_t> libs_for_construction;
    for (size_t i = 0; i < cfg::get().ds.reads.lib_count(); ++i)
        if (cfg::get().ds.reads[i].is_graph_contructable())
            libs_for_construction.push_back(i);

    if (cfg::get().use_multithreading) {
        auto streams = single_binary_readers_for_libs(libs_for_construction, true, true);
        construct_graph<io::SingleReadSeq>(streams, gp, contigs_stream);
    } else {
        io::ReadStreamList<io::SingleRead> streams(single_easy_reader_for_libs(libs_for_construction, true, true));
        construct_graph<io::SingleRead>(streams, gp, contigs_stream);
    }
}

} //namespace debruijn_graph
