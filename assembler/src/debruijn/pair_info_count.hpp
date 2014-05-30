//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include "stage.hpp"

namespace debruijn_graph {

class PairInfoCount : public spades::AssemblyStage {
  public:
    PairInfoCount()
        : AssemblyStage("Paired Information Counting", "late_pair_info_count") {}

    void run(conj_graph_pack &gp, const char*);
};

}

