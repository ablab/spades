//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "pipeline/stage.hpp"

namespace debruijn_graph {

class PairInfoCount : public spades::AssemblyStage {
  public:
    PairInfoCount(bool preliminary = false)
        : AssemblyStage(preliminary ? "Preliminary Paired Information Counting" : "Paired Information Counting",
                        preliminary ? "late_pair_info_count_preliminary" : "late_pair_info_count") {}

    void run(conj_graph_pack &gp, const char*);
};

}

