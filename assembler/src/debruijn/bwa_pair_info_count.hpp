//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "stage.hpp"

namespace debruijn_graph {

class BWAPairInfoCount : public spades::AssemblyStage {
  public:
    BWAPairInfoCount(bool preliminary = false)
        : AssemblyStage(preliminary ? "Preliminary BWA Paired Information Counting" : "BWA Paired Information Counting",
                        preliminary ? "bwa_late_pair_info_count_preliminary" : "bwa_late_pair_info_count") {}

    void run(conj_graph_pack &gp, const char*);
};

}
