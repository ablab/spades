//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2015-2022 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "pipeline/stage.hpp"

namespace debruijn_graph {

class RepeatResolution : public spades::AssemblyStage {
    const bool preliminary_;
public:
    RepeatResolution(bool preliminary = false)
            : AssemblyStage(preliminary ? "Preliminary Repeat Resolving" : "Repeat Resolving",
                            preliminary ? "repeat_resolving_preliminary" : "repeat_resolving"),
              preliminary_(preliminary) { }

    void load(graph_pack::GraphPack &, const std::filesystem::path &, const char *) override;
    void save(const graph_pack::GraphPack &, const std::filesystem::path &, const char *) const override;
    void run(graph_pack::GraphPack &gp, const char *) override;
};

}

