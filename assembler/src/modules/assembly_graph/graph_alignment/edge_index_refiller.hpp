//***************************************************************************
//* Copyright (c) 2016 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

namespace debruijn_graph {

// The stuff is template here to provide interface w/o including any headers
// In our case both EdgeIndex and Graph are very complex template objects we
// do not want to pull the corresponding headers here until we untangle all
// the mess
struct EdgeIndexRefiller {
    template<class EdgeIndex, class Graph>
    void Refill(EdgeIndex &index, const Graph &g);
};

}
