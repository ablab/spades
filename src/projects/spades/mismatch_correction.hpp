//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2015-2022 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "pipeline/stage.hpp"
#include "pipeline/mpi_stage.hpp"

namespace debruijn_graph {

class MismatchCorrection : public spades::MPIAssemblyStage {
public:
    MismatchCorrection()
            : MPIAssemblyStage("Mismatch Correction", "mismatch_correction") { }

    void run(graph_pack::GraphPack &gp, const char *) override;
};

}

