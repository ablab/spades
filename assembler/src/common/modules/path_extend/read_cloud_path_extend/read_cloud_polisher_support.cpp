#include "read_cloud_polisher_support.hpp"
set<EdgeId> path_extend::SupportedEdgesGraphExtractor::ExtractSupportedEdges(const EdgeId start,
                                                                             const EdgeId end) const {
    const size_t MAX_GAP_SIZE = 5000;
    const size_t min_short_edge_length = read_cloud_predicate_->GetParams().edge_length_threshold_;

    DEBUG("Extracting supported edges");
    DijkstraHelper<debruijn_graph::Graph>::BoundedDijkstra forward_dijkstra(
        DijkstraHelper<debruijn_graph::Graph>::CreateBoundedDijkstra(g_, MAX_GAP_SIZE));
    TRACE("Running dijkstra from " << start.int_id());
    forward_dijkstra.Run(g_.EdgeEnd(start));
    set<EdgeId> forward_edges;
    for (auto v: forward_dijkstra.ReachedVertices()) {
        if (forward_dijkstra.GetDistance(v) < MAX_GAP_SIZE) {
            for (auto connected: g_.OutgoingEdges(v)) {
                forward_edges.insert(connected);
            }
        }
    }
    DEBUG("Extracted " << forward_edges.size() << " forward edges");
    bool end_found = forward_edges.find(end) != forward_edges.end();
    if (not end_found) {
        DEBUG("Start and end are not connected!");
    }

    DijkstraHelper<debruijn_graph::Graph>::BackwardBoundedDijkstra backward_dijkstra(
        DijkstraHelper<debruijn_graph::Graph>::CreateBackwardBoundedDijkstra(g_, MAX_GAP_SIZE));
    backward_dijkstra.Run(g_.EdgeStart(end));
    set<EdgeId> backward_edges;
    for (auto v: backward_dijkstra.ReachedVertices()) {
        if (backward_dijkstra.GetDistance(v) < MAX_GAP_SIZE) {
            for (auto connected: g_.IncomingEdges(v)) {
                TRACE("Checking edge " << connected.int_id());
                backward_edges.insert(connected);
            }
        }
    }
    DEBUG("Extracted " << backward_edges.size() << " backward edges");

    set<EdgeId> intersection;
    std::set_intersection(forward_edges.begin(), forward_edges.end(), backward_edges.begin(), backward_edges.end(),
                          std::inserter(intersection, intersection.begin()));
    DEBUG("Intersection size: " << intersection.size());
    set<EdgeId> supported_edges;
    auto is_supported = [this, min_short_edge_length](const EdgeId& edge) {
      return g_.length(edge) > min_short_edge_length and read_cloud_predicate_->Check(edge);
    };
    std::copy_if(intersection.begin(), intersection.end(), std::inserter(supported_edges, supported_edges.end()), is_supported);
    DEBUG("Extracted " << supported_edges.size() << " supported edges");
    return supported_edges;
}
path_extend::SupportedEdgesGraphExtractor::SupportedEdgesGraphExtractor(
        const Graph &g_,
        const shared_ptr<path_extend::LongEdgePairGapCloserPredicate> &read_cloud_predicate_)
    : g_(g_), read_cloud_predicate_(read_cloud_predicate_) {}
