//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
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

    void init(debruijn_graph::conj_graph_pack &gp, const char *) override;
    void fini(debruijn_graph::conj_graph_pack &gp) override;
};

}

