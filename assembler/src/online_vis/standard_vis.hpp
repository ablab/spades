//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include "graph_pack.hpp"
#include "standard_base.hpp"
#include "omni/total_labeler.hpp"
#include "omni/visualization/graph_labeler.hpp"

#include <readline/readline.h>
#include <readline/history.h>

namespace online_visualization {

typedef debruijn_graph::conj_graph_pack GraphPack;
typedef GraphPack::graph_t Graph;
typedef GraphPack::index_t Index;
typedef EdgesPositionHandler<Graph> EdgePos;
typedef Graph::VertexId VertexId;
typedef Graph::EdgeId EdgeId;

typedef omnigraph::TotalLabelerGraphStruct  <debruijn_graph::ConjugateDeBruijnGraph>	total_labeler_graph_struct;
typedef omnigraph::TotalLabeler             <debruijn_graph::ConjugateDeBruijnGraph>    total_labeler;
typedef omnigraph::GraphLabeler             <debruijn_graph::ConjugateDeBruijnGraph>    graph_labeler;

}
