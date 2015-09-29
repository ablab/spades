//
// Created by andrey on 17.09.15.
//

#include "scaffold_graph.hpp"


namespace path_extend { namespace scaffold_graph {

std::atomic<ScaffoldGraph::ScaffoldEdgeIdT> ScaffoldGraph::ScaffoldEdge::scaffold_edge_id_{0};

} //scaffold_graph
} //path_extend