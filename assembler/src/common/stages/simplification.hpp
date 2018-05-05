//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
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

    void run(conj_graph_pack &gp, const char *);
};

class Simplification : public spades::AssemblyStage {
    const bool preliminary_;
public:
    Simplification(bool preliminary = false)
            : AssemblyStage(preliminary ? "Preliminary Simplification" : "Simplification",
                            preliminary ? "simplification_preliminary" : "simplification"),
              preliminary_(preliminary) { }

    void run(conj_graph_pack &gp, const char *);
};

class SimplificationCleanup : public spades::AssemblyStage {
public:
    SimplificationCleanup()
            : AssemblyStage("Simplification Cleanup", "simplification_cleanup") { }

    void run(conj_graph_pack &gp, const char *);
};

}

