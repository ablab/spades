//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include "stage.hpp"

namespace debruijn_graph {

class Simplification : public spades::AssemblyStage {
  public:
    Simplification()
        : AssemblyStage("Simplification", "simplification") {}

    void run(conj_graph_pack &gp);
};

}

