//***************************************************************************
//* Copyright (c) 2015-2019 Saint-Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include "pipeline/stage.hpp"

namespace debruijn_graph {

class ChromosomeRemoval : public spades::AssemblyStage {
public:
//TODO: multiple stages with same name, saves should be switched off here.
    ChromosomeRemoval(size_t ext_limit = 0)
            : AssemblyStage("Chromosome Removal", "chromosome_removal"), ext_limit_(ext_limit) {
    }

    void run(GraphPack &gp, const char *) override;

private:
    size_t ext_limit_;
    DECL_LOGGER("ChromosomeRemoval");
};
}
