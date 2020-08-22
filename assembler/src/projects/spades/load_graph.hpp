//***************************************************************************
//* Copyright (c) 2020 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "common/pipeline/stage.hpp"

namespace debruijn_graph {

class LoadGraph : public spades::AssemblyStage {
  public:
    LoadGraph()
        : AssemblyStage("Load and prepare assembly graph", "load_graph") {}

    void run(GraphPack &gp, const char*) override;
};

}
