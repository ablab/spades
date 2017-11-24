//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
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

    void load(conj_graph_pack &, const std::string &, const char *) { }

    void save(const conj_graph_pack &, const std::string &, const char *) const { }

    void run(conj_graph_pack &gp, const char *);
};

}

