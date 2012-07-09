#pragma once

#include "standard.hpp"
#include "graph_pack.hpp"
#include "standard_base.hpp"

typedef debruijn_graph::conj_graph_pack GraphPack;
typedef GraphPack::graph_t Graph;
typedef EdgesPositionHandler<Graph> EdgePos;
typedef Graph::VertexId VertexId;
typedef Graph::EdgeId EdgeId;
