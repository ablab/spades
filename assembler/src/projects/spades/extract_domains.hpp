//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "pipeline/stage.hpp"

namespace debruijn_graph {

class ExtractDomains : public spades::AssemblyStage {
public:
    ExtractDomains()
            : AssemblyStage("Extract Domains", "extract_domains") { }

    void run(GraphPack &gp, const char *) override;
};

}
