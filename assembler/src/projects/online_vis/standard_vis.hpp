//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "pipeline/graph_pack.hpp"

#include <readline/readline.h>
#include <readline/history.h>

//TODO: remove this sometime
using namespace std;
using namespace fs;

namespace online_visualization {
    typedef debruijn_graph::conj_graph_pack GraphPack;
    typedef GraphPack::graph_t Graph;
    typedef GraphPack::index_t Index;
    typedef EdgesPositionHandler<Graph> EdgePos;
    typedef Graph::VertexId VertexId;
    typedef Graph::EdgeId EdgeId;
}
