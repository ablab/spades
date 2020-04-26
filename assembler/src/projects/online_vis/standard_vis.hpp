//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "pipeline/graph_pack.hpp"
#include "assembly_graph/core/graph.hpp"
#include "assembly_graph/handlers/edges_position_handler.hpp"
#include "modules/alignment/edge_index.hpp"

#include <readline/readline.h>
#include <readline/history.h>

//TODO: remove this sometime
using namespace std;
using namespace fs;

namespace online_visualization {
    using debruijn_graph::GraphPack;
    using debruijn_graph::Graph;
    using Index = debruijn_graph::EdgeIndex<Graph>;
    using EdgePos = omnigraph::EdgesPositionHandler<Graph>;
    using VertexId = Graph::VertexId;
    using EdgeId = Graph::EdgeId;
}
