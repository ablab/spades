#include "scaffold_graph_extractor.hpp"

namespace path_extend {
vector<ScaffoldGraphExtractor::ScaffoldEdge> ScaffoldGraphExtractor::ExtractUnivocalEdges(const ScaffoldGraph& scaffold_graph) {
    vector<ScaffoldEdge> result;
    for (const ScaffoldGraph::ScaffoldEdge& edge: scaffold_graph.edges()) {
        if (scaffold_graph.HasUniqueOutgoing(edge.getStart()) and scaffold_graph.HasUniqueIncoming(edge.getEnd())) {
            result.push_back(edge);
        }
    }
    return result;
}
}
