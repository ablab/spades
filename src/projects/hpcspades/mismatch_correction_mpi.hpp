//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "projects/spades/mismatch_correction.hpp"
#include "common/alignment/sequence_mapper_notifier_mpi.hpp"
#include "common/pipeline/mpi_stage.hpp"
#include "pipeline/graph_pack_helpers.h"

namespace debruijn_graph {
class MismatchCorrectionMPI : public spades_mpi::MPIAssemblyStage {
public:
    MismatchCorrectionMPI()
            : MPIAssemblyStage("Mismatch Correction", "mismatch_correction") {}

    void run(graph_pack::GraphPack &gp, const char *) override {
        EnsureBasicMapping(gp);
        size_t corrected = mismatches::MismatchShallNotPass(MapLibFuncMPI(), gp, 2, partask::overall_num_threads()).
                           ParallelStopAllMismatches(1);
        INFO("Corrected " << corrected << " nucleotides");
    }
};
}

