//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "mismatch_correction.hpp"

#include "mismatch_shall_not_pass.hpp"
#include "read_converter.hpp"

namespace debruijn_graph {

void MismatchCorrection::run(conj_graph_pack &gp, const char*) {
    gp.EnsureBasicMapping();
    auto streams = single_binary_readers(true,  true);
    size_t corrected = MismatchShallNotPass<conj_graph_pack, io::SingleReadSeq>(gp, 2).ParallelStopAllMismatches(streams, 1);
    INFO("Corrected " << corrected << " nucleotides");
}

}
