//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "pipeline/stage.hpp"

namespace debruijn_graph {

class MismatchCorrection : public spades::AssemblyStage {
public:
    MismatchCorrection()
            : AssemblyStage("Mismatch Correction", "mismatch_correction") { }

    void run(conj_graph_pack &gp, const char *);
};

}

