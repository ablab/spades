//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include "stage.hpp"

namespace debruijn_graph {

//todo rename
class PairInfoRemover : public spades::AssemblyStage {
  public:
	PairInfoRemover()
        : AssemblyStage("Pair Info Remover", "pair_info_remover") {}

    void run(conj_graph_pack &gp, const char*);
};

}
