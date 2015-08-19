#pragma once

#include "graph_pack.hpp"
#include "standard_base.hpp"

#include <readline/readline.h>
#include <readline/history.h>

namespace online_visualization {
    typedef debruijn_graph::conj_graph_pack GraphPack;
    typedef GraphPack::graph_t Graph;
    typedef GraphPack::index_t Index;
    typedef EdgesPositionHandler<Graph> EdgePos;
    typedef Graph::VertexId VertexId;
    typedef Graph::EdgeId EdgeId;
}
