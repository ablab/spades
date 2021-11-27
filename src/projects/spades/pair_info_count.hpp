//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2015-2022 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "alignment/sequence_mapper_notifier.hpp"
#include "pipeline/stage.hpp"

namespace debruijn_graph {
class PairInfoCountBase {
public:
    void execute(graph_pack::GraphPack &gp, const char *, const MapLibBase&, size_t num_readers=0);
};

class PairInfoCount : public PairInfoCountBase, public spades::AssemblyStage {
 public:
    PairInfoCount(bool preliminary = false)
        : AssemblyStage(preliminary ? "Preliminary Paired Information Counting" : "Paired Information Counting",
                           preliminary ? "late_pair_info_count_preliminary" : "late_pair_info_count") {}

    void run(graph_pack::GraphPack &gp, const char*) override;
};

}

