#pragma once

#include "common/assembly_graph/components/graph_component.hpp"
#include "common/assembly_graph/core/graph.hpp"

namespace gfa_tools {

omnigraph::GraphComponent<debruijn_graph::Graph> FindAroundSequenceScope(debruijn_graph::Graph& g,  
                                                    const debruijn_graph::EdgeId edge_id, const size_t depth);

}
