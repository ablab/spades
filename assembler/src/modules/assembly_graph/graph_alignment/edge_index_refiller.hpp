//***************************************************************************
//* Copyright (c) 2016 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

namespace debruijn_graph {

struct EdgeIndexRefiller {
    template<class EdgeIndex, class Graph>
    void Refill(EdgeIndex &index, const Graph &g);
};

}
