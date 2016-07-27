//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "pipeline/stage.hpp"
#include "assembly_graph/graph_alignment/sequence_mapper_notifier.hpp"

namespace debruijn_graph {

void ProcessSingleReads(conj_graph_pack &gp,
                        size_t ilib,
                        bool use_binary = true,
                        bool map_paired = false,
                        SequenceMapperListener* mapping_listener_ptr = nullptr);

class PairInfoCount : public spades::AssemblyStage {
  public:
    PairInfoCount(bool preliminary = false)
        : AssemblyStage(preliminary ? "Preliminary Paired Information Counting" : "Paired Information Counting",
                        preliminary ? "late_pair_info_count_preliminary" : "late_pair_info_count") {}

    void run(conj_graph_pack &gp, const char*);
};

}

