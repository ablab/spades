//***************************************************************************
//* Copyright (c) 2015-2021 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "projects/spades/pair_info_count.hpp"
#include "mpi/alignment/sequence_mapper_notifier_mpi.hpp"
#include "mpi/pipeline/mpi_stage.hpp"

namespace debruijn_graph {
class PairInfoCountMPI : public PairInfoCountBase, public spades_mpi::MPIAssemblyStage {
public:
    PairInfoCountMPI(bool preliminary = false)
            : MPIAssemblyStage(preliminary ? "Preliminary Paired Information Counting" : "Paired Information Counting",
                               preliminary ? "late_pair_info_count_preliminary" : "late_pair_info_count") {}

    void run(graph_pack::GraphPack &gp, const char* s) override {
        execute(gp, s, MapLibFuncMPI(), partask::overall_num_threads());
    }
};

}
