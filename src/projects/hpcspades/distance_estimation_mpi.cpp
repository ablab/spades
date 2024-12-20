//***************************************************************************
//* Copyright (c) 2021 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "distance_estimation_mpi.hpp"

#include "mpi/paired_info/distance_estimation_utils.hpp"

namespace debruijn_graph {
void DistanceEstimationMPI::run(graph_pack::GraphPack &gp, const char* s) {
    DistanceEstimationBase::run(gp, s, distance_estimation::EstimatePairedDistancesMPI,
                                distance_estimation::EstimateScaffoldingDistancesMPI);
}
}
