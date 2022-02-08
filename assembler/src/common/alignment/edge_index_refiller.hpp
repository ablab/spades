//***************************************************************************
//* Copyright (c) 2016 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "assembly_graph/core/graph.hpp"

#include <string>
#include <vector>

namespace debruijn_graph {

// The stuff is template here to provide interface w/o including any headers
// In our case both EdgeIndex and Graph are very complex template objects we
// do not want to pull the corresponding headers here until we untangle all
// the mess
struct EdgeIndexRefiller {
    std::string workdir_;

    EdgeIndexRefiller(const std::string &workdir);

    template<class EdgeIndex>
    void Refill(EdgeIndex &index, const Graph &g);

    template<class EdgeIndex>
    void Refill(EdgeIndex &index, const Graph &g,
                const std::vector<typename Graph::EdgeId> &edges);
};

}
