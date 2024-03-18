//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2015-2022 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "pipeline/stage.hpp"

namespace debruijn_graph {

struct ConstructionStorage;

class Construction : public spades::CompositeStageDeferred<ConstructionStorage> {
public:
    Construction();
    ~Construction();

    void init(graph_pack::GraphPack &gp, const char *) override;
    void fini(graph_pack::GraphPack &gp) override;
};

}

