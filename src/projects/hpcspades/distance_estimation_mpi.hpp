//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "common/pipeline/mpi_stage.hpp"
#include "projects/spades/distance_estimation.hpp"

namespace debruijn_graph {
class DistanceEstimationMPI : public DistanceEstimationBase, public spades_mpi::MPIAssemblyStage {
public:
    DistanceEstimationMPI(bool preliminary = false)
            : MPIAssemblyStage(preliminary ? "Preliminary Distance Estimation" : "Distance Estimation",
                               preliminary ? "distance_estimation_preliminary" : "distance_estimation") {}

    void run(graph_pack::GraphPack &gp, const char *) override;
};
}

