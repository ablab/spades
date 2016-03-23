//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "compare_standard.hpp"
#include "coloring.hpp"

namespace cap {

template<class Graph>
struct cap_graph_pack {
    typedef Graph graph_t;
    typedef string contig_id_t;
    typedef typename Graph::EdgeId EdgeId;
    Graph g;
    omnigraph::GraphElementFinder<Graph> element_finder;
    ColorHandler<Graph> coloring;
//    map<contig_id_t, vector<EdgeId>> red_paths;
//    map<contig_id_t, vector<EdgeId>> blue_paths;
    EdgesPositionHandler<Graph> edge_pos;

    cap_graph_pack(size_t k) :
            g(k), element_finder(g), coloring(g), edge_pos(g) {

    }
};

}
