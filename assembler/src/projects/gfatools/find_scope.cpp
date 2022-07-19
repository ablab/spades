#include "find_scope.hpp"

#include "common/assembly_graph/dijkstra/dijkstra_algorithm.hpp"

namespace gfa_tools {

omnigraph::GraphComponent<debruijn_graph::Graph> FindAroundSequenceScope(debruijn_graph::Graph& g, 
                                                    const debruijn_graph::EdgeId edge_id, const size_t depth) {
    typedef omnigraph::ComposedDijkstraSettings<debruijn_graph::Graph,
                                              omnigraph::BoundedEdgeLenCalculator<debruijn_graph::Graph>,
                                              omnigraph::BoundProcessChecker<debruijn_graph::Graph>,
                                              omnigraph::BoundPutChecker<debruijn_graph::Graph>,
                                              omnigraph::UnorientedNeighbourIteratorFactory<debruijn_graph::Graph> > ComposedDijkstraSettings;
    typedef omnigraph::Dijkstra<debruijn_graph::Graph, ComposedDijkstraSettings> Dijkstra;
    
    Dijkstra dijkstra_algorithm(g, ComposedDijkstraSettings(
        omnigraph::BoundedEdgeLenCalculator<debruijn_graph::Graph>(g, 0),
        omnigraph::BoundProcessChecker<debruijn_graph::Graph>(depth),
        omnigraph::BoundPutChecker<debruijn_graph::Graph>(depth),
        omnigraph::UnorientedNeighbourIteratorFactory<debruijn_graph::Graph>(g)
    ));                                          
    dijkstra_algorithm.Run(g.EdgeStart(edge_id));
    std::vector<debruijn_graph::VertexId> reached_vertices = dijkstra_algorithm.ReachedVertices();
    
    return omnigraph::GraphComponent<debruijn_graph::Graph>::FromVertices(g, reached_vertices.begin(), reached_vertices.end());
}
}
