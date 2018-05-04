//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include <modules/mismatch_shall_not_pass.hpp>
#include "mismatch_correction.hpp"

#include "io/dataset_support/read_converter.hpp"

namespace debruijn_graph {

void MismatchCorrection::run(conj_graph_pack &gp, const char*) {
    gp.EnsureBasicMapping();

    auto& dataset = cfg::get_writable().ds;
    std::vector<size_t> libs;
    for (size_t i = 0; i < dataset.reads.lib_count(); ++i) {
        if (dataset.reads[i].is_mismatch_correctable())
            libs.push_back(i);
    }
    auto streams = io::single_binary_readers_for_libs(dataset.reads, libs);
    size_t corrected = mismatches::MismatchShallNotPass(gp, 2).
                       ParallelStopAllMismatches(streams, 1);
    INFO("Corrected " << corrected << " nucleotides");
}

}
