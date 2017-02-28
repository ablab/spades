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
private:
    string contig_name_prefix_;

public:
    SecondPhaseSetup(const string& contig_name_prefix = "")
            : AssemblyStage("Second Phase Setup", "second_phase_setup"),contig_name_prefix_(contig_name_prefix)  { }

    void run(conj_graph_pack &gp, const char *);
};

}
