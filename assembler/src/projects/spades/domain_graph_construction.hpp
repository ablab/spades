//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "common/pipeline/stage.hpp"

namespace debruijn_graph {

//todo rename
class DomainGraphConstruction : public spades::AssemblyStage {
  public:
    DomainGraphConstruction()
        : AssemblyStage("Domain Graph Construction", "domain_graph_construction") {}

    void run(GraphPack &gp, const char*) override;
};

}
