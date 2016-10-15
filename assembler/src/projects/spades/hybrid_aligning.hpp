//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "pipeline/stage.hpp"

namespace debruijn_graph {

class HybridLibrariesAligning : public spades::AssemblyStage {
public:
    HybridLibrariesAligning()
            : AssemblyStage("Hybrid Aligning", "hybrid_aligning") {
    }
    void run(conj_graph_pack &gp, const char*);
};

}

