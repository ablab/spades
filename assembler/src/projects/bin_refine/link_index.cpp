//***************************************************************************
//* Copyright (c) 2021 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "link_index.hpp"

using namespace binning;

void GraphLinkIndex::Init(const debruijn_graph::Graph &g) {
    for (EdgeId e : g.canonical_edges()) {
        for (EdgeId o : g_.OutgoingEdges(g.EdgeEnd(e)))
            add(e, o);

        for (EdgeId i : g_.IncomingEdges(g.EdgeStart(e)))
            add(e, i);
    }
}

