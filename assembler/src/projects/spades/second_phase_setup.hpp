//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "pipeline/stage.hpp"

namespace debruijn_graph {

//todo rename
class SecondPhaseSetup : public spades::AssemblyStage {
public:
    SecondPhaseSetup()
            : AssemblyStage("Second Phase Setup", "second_phase_setup") { }

    void run(conj_graph_pack &gp, const char *);
};

}
