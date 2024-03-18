//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2020-2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "pipeline/stage.hpp"

namespace debruijn_graph {

class LoadGraph : public spades::AssemblyStage {
  public:
    LoadGraph()
        : AssemblyStage("Load and prepare assembly graph", "load_graph") {}

    void run(graph_pack::GraphPack &gp, const char*) override;
};

}
