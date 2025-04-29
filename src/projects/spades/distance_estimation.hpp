//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2015-2022 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <common/pipeline/mpi_stage.hpp>
#include "pipeline/stage.hpp"

namespace debruijn_graph {

class DistanceEstimation : public spades::MPIAssemblyStage {
 public:
    DistanceEstimation(bool preliminary = false)
        : MPIAssemblyStage(preliminary ? "Preliminary Distance Estimation" : "Distance Estimation",
                           preliminary ? "distance_estimation_preliminary" : "distance_estimation") {}

    void run(graph_pack::GraphPack &gp, const char*) override;
};
}

