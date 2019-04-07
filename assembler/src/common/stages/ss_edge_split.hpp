//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "pipeline/stage.hpp"

namespace debruijn_graph {


class SSEdgeSplit: public spades::AssemblyStage {

public:
    SSEdgeSplit() : AssemblyStage("Strand-specific Edge Spliting", "ss_edge_splitting") {}

    void run(conj_graph_pack& gp, const char *);
};
}


