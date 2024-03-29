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

class RawSimplification : public spades::AssemblyStage {
    //TODO this is a dirty hack to get a different behavior for last meta iteration!
    //Has nothing to do with two-step concept...
    const bool preliminary_;
public:
    RawSimplification(bool preliminary = false)
            : AssemblyStage("Raw Simplification", "raw_simplification"),
              preliminary_(preliminary) { }

    void run(graph_pack::GraphPack &gp, const char *) override;
};

class Simplification : public spades::AssemblyStage {
    const bool preliminary_;
public:
    Simplification(bool preliminary = false)
            : AssemblyStage(preliminary ? "Preliminary Simplification" : "Simplification",
                            preliminary ? "simplification_preliminary" : "simplification"),
              preliminary_(preliminary) { }

    void run(graph_pack::GraphPack &gp, const char *) override;
};

class SimplificationCleanup : public spades::AssemblyStage {
public:
    SimplificationCleanup()
            : AssemblyStage("Simplification Cleanup", "simplification_cleanup") { }

    void run(graph_pack::GraphPack &gp, const char *) override;
};

}

