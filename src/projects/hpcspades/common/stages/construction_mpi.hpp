//***************************************************************************
//* Copyright (c) 2021 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "projects/hpcspades/common/pipeline/mpi_stage.hpp"

namespace debruijn_graph {
struct ConstructionStorage;

class ConstructionMPI : public spades_mpi::MPICompositeStageDeferred<ConstructionStorage> {
public:
    ConstructionMPI();
    ~ConstructionMPI();

    void init(graph_pack::GraphPack &gp, const char *) override;
    void fini(graph_pack::GraphPack &gp) override;
};
}

