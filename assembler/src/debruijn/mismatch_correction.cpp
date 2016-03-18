//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "mismatch_correction.hpp"

#include "mismatch_shall_not_pass.hpp"
#include "io/dataset_support/read_converter.hpp"

namespace debruijn_graph {

void MismatchCorrection::run(conj_graph_pack &gp, const char*) {
    gp.EnsureBasicMapping();
    std::vector<size_t> libs;
    for (size_t i = 0; i < cfg::get().ds.reads.lib_count(); ++i) {
        if (cfg::get().ds.reads[i].is_mismatch_correctable())
            libs.push_back(i);
    }
    auto streams = single_binary_readers_for_libs(libs, true,  true);
    size_t corrected = MismatchShallNotPass<conj_graph_pack, io::SingleReadSeq>(gp, 2).ParallelStopAllMismatches(streams, 1);
    INFO("Corrected " << corrected << " nucleotides");
}

}
