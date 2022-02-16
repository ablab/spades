//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "pipeline/stage.hpp"

namespace debruijn_graph {

//todo rename
    class RestrictedEdgesFilling : public spades::AssemblyStage {
    public:
        RestrictedEdgesFilling()
                : AssemblyStage("Restricted Edges Filling", "restricted_edges_filling") {}

        void run(graph_pack::GraphPack &gp, const char*) override;
    };

}
