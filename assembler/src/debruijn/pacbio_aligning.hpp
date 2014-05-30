//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include "stage.hpp"

namespace debruijn_graph {

class PacBioAligning : public spades::AssemblyStage {
public:
    PacBioAligning()
            : AssemblyStage("PacBio Aligning", "pacbio_aligning") {
    }
    void run(conj_graph_pack &gp, const char*);
};

}

