//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include "stage.hpp"

namespace debruijn_graph {

class RepeatResolution : public spades::AssemblyStage {
  public:
    RepeatResolution()
        : AssemblyStage("Repeat Resolving", "repeat_resolving") {}

    void load(conj_graph_pack& gp) {}
    void save(const conj_graph_pack& gp) const {}

    void run(conj_graph_pack &gp);
};

class ContigOutput : public spades::AssemblyStage {
  public:
    ContigOutput()
        : AssemblyStage("Contig Output", "contig_output") {}

    void load(conj_graph_pack& gp) {}
    void save(const conj_graph_pack& gp) const {}

    void run(conj_graph_pack &gp);
};

}

