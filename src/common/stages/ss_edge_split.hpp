//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2019-2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "pipeline/stage.hpp"

namespace debruijn_graph {


class SSEdgeSplit: public spades::AssemblyStage {

public:
    SSEdgeSplit() : AssemblyStage("Strand-specific Edge Spliting", "ss_edge_splitting") {}

    void run(graph_pack::GraphPack& gp, const char *) override;
};
}


