//***************************************************************************
//* Copyright (c) 2021 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "distance_estimation.hpp"

namespace omnigraph {
namespace de {
using namespace debruijn_graph;

void DistanceEstimatorMPI::Estimate(OutPairedIndex &result, size_t nthreads) const {
    this->Init();
    const auto &index = this->index();

    DEBUG("Collecting edge infos");
    std::vector<EdgeId> edges;
    for (EdgeId e : this->graph().edges())
        edges.push_back(e);

    partask::TaskRegistry treg;
    auto dist_estimator_mpi = treg.add<DistanceEstimatorTask>(std::cref(index), std::cref(*dist_estimator_),
                                                              std::ref(result));
    treg.listen();

    if (partask::master()) {
        dist_estimator_mpi(edges, nthreads);
    }
    treg.stop_listening();
    partask::broadcast(result);
}
}
}
