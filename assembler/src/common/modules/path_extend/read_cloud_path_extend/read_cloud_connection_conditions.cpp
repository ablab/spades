#include "read_cloud_connection_conditions.hpp"
#include "read_cloud_dijkstras.hpp"

namespace path_extend {

map<EdgeId, double> AssemblyGraphUniqueConnectionCondition::ConnectedWith(EdgeId e) const {
    VERIFY_MSG(interesting_edge_set_.find(e) != interesting_edge_set_.end(),
               " edge " << e.int_id() << " not applicable for connection condition");
    if (stored_distances_.find(e) != stored_distances_.end()) {
        return stored_distances_[e];
    }
    stored_distances_.insert(make_pair(e, map<debruijn_graph::EdgeId, double>()));
    for (auto connected: g_.OutgoingEdges(g_.EdgeEnd(e))) {
        if (interesting_edge_set_.find(connected) != interesting_edge_set_.end()) {
            stored_distances_[e].insert(make_pair(connected, 1));
        }
    }

    auto dij = omnigraph::CreateUniqueDijkstra(g_, max_connection_length_, unique_storage_);
    dij.Run(g_.EdgeEnd(e));
    for (auto v: dij.ReachedVertices()) {
        for (auto connected: g_.OutgoingEdges(v)) {
            if (interesting_edge_set_.find(connected) != interesting_edge_set_.end() &&
                dij.GetDistance(v) < max_connection_length_ && connected != g_.conjugate(e)) {
                stored_distances_[e].insert(make_pair(connected, 1));
            }
        }
    }
    return stored_distances_[e];
}
AssemblyGraphUniqueConnectionCondition::AssemblyGraphUniqueConnectionCondition(const Graph& g,
                                                                               size_t max_connection_length,
                                                                               const ScaffoldingUniqueEdgeStorage& unique_edges)
    : AssemblyGraphConnectionCondition(g, max_connection_length, unique_edges), unique_storage_(unique_edges) {}
bool AssemblyGraphUniqueConnectionCondition::IsLast() const {
    return false;
}

//fixme use total barcodes from index
BarcodeScoreFunction::BarcodeScoreFunction(const size_t read_count_threshold,
                                           const size_t tail_threshold,
                                           const barcode_index::FrameBarcodeIndexInfoExtractor& barcode_extractor_,
                                           const Graph& graph)
    : read_count_threshold_(read_count_threshold),
      tail_threshold_(tail_threshold),
      total_barcodes_(1000000),
      barcode_extractor_(barcode_extractor_),
      graph_(graph){}
double BarcodeScoreFunction::GetScore(const scaffold_graph::ScaffoldGraph::ScaffoldEdge& edge) const {
    EdgeId first = edge.getStart();
    EdgeId second = edge.getEnd();
    size_t first_length = graph_.length(first);
    VERIFY(graph_.length(first) >= tail_threshold_);
    VERIFY(graph_.length(second) >= tail_threshold_);
    auto first_barcodes = barcode_extractor_.GetBarcodesFromRange(first, read_count_threshold_,
                                                                  first_length - tail_threshold_, first_length);
    auto second_barcodes = barcode_extractor_.GetBarcodesFromRange(second, read_count_threshold_, 0, tail_threshold_);
    size_t shared_count = barcode_extractor_.CountSharedBarcodesWithFilter(first, second, read_count_threshold_, tail_threshold_);
//    double second_coverage = graph_.coverage(second);
    double first_barcodes_size = static_cast<double>(first_barcodes.size());
    double second_barcodes_size = static_cast<double>(second_barcodes.size());
    double total_barcodes = static_cast<double>(total_barcodes_);

    return static_cast<double>(shared_count) * total_barcodes / (first_barcodes_size * second_barcodes_size);
}
LongGapDijkstraPredicate::LongGapDijkstraPredicate(const Graph& g,
                                                   const ScaffoldingUniqueEdgeStorage& unique_storage_,
                                                   const barcode_index::FrameBarcodeIndexInfoExtractor& barcode_extractor_,
                                                   const LongGapDijkstraParams& params)
    : g(g),
      unique_storage_(unique_storage_),
      barcode_extractor_(barcode_extractor_),
      params_(params) {}

bool LongGapDijkstraPredicate::Check(const scaffold_graph::ScaffoldGraph::ScaffoldEdge& scaffold_edge) const {
    auto barcode_intersection = barcode_extractor_.GetSharedBarcodesWithFilter(scaffold_edge.getStart(), scaffold_edge.getEnd(),
                                                                               params_.count_threshold_, params_.tail_threshold_);
    auto long_gap_dijkstra = CreateLongGapCloserDijkstra(g, params_.distance_, unique_storage_,
                                                         barcode_extractor_, barcode_intersection,
                                                         params_.barcode_threshold_, params_.count_threshold_,
                                                         params_.tail_threshold_, params_.len_threshold_);
    long_gap_dijkstra.Run(g.EdgeEnd(scaffold_edge.getStart()));
    VertexId target = g.EdgeStart(scaffold_edge.getEnd());
    for (const auto& vertex: long_gap_dijkstra.ReachedVertices()) {
        if (vertex == target) {
            return true;
        }
    }
    return false;
}
LongGapDijkstraParams::LongGapDijkstraParams(const size_t barcode_threshold_,
                                             const size_t count_threshold_,
                                             const size_t tail_threshold_,
                                             const size_t len_threshold_,
                                             const size_t distance)
    : barcode_threshold_(barcode_threshold_),
      count_threshold_(count_threshold_),
      tail_threshold_(tail_threshold_),
      len_threshold_(len_threshold_),
      distance_(distance){}

EdgeSplitPredicate::EdgeSplitPredicate(const Graph& g_,
                                       const barcode_index::FrameBarcodeIndexInfoExtractor& barcode_extractor_,
                                       const size_t count_threshold_,
                                       double strictness)
    : g_(g_),
      barcode_extractor_(barcode_extractor_),
      count_threshold_(count_threshold_),
      strictness_(strictness) {}
bool EdgeSplitPredicate::CheckOrderingForThreeSegments(const vector<EdgeSplitPredicate::BarcodeId>& first,
                                                       const vector<EdgeSplitPredicate::BarcodeId>& second,
                                                       const vector<EdgeSplitPredicate::BarcodeId>& third,
                                                       double strictness) const {
    vector<BarcodeId> first_second_intersection;
    std::set_intersection(first.begin(), first.end(), second.begin(), second.end(),
                          std::back_inserter(first_second_intersection));
    vector<BarcodeId> first_third_intersection;
    std::set_intersection(first.begin(), first.end(), third.begin(), third.end(),
                          std::back_inserter(first_third_intersection));

    bool result = math::ge(static_cast<double>(first_second_intersection.size()),
                           strictness * static_cast<double>(first_third_intersection.size()));
    if (not result) {
        DEBUG("First second: " << first_second_intersection.size());
        DEBUG("First third: " << first_third_intersection.size());
        DEBUG("First: " << first.size());
        DEBUG("Second: " << second.size());
        DEBUG("Third: " << third.size());
    }
    return result;
}

bool EdgeSplitPredicate::CheckOrderingForFourSegments(const vector<EdgeSplitPredicate::BarcodeId>& first,
                                                      const vector<EdgeSplitPredicate::BarcodeId>& second,
                                                      const vector<EdgeSplitPredicate::BarcodeId>& third,
                                                      const vector<EdgeSplitPredicate::BarcodeId>& fourth) const {
    vector<BarcodeId> first_fourth_intersection;
    std::set_intersection(first.begin(), first.end(), fourth.begin(), fourth.end(),
                          std::back_inserter(first_fourth_intersection));
    vector<BarcodeId> second_third_intersection;
    std::set_intersection(second.begin(), second.end(), third.begin(), third.end(),
                          std::back_inserter(second_third_intersection));
    if (second_third_intersection.size() <= first_fourth_intersection.size()) {
        DEBUG("Second third: " << second_third_intersection.size());
        DEBUG("First fourth: " << first_fourth_intersection.size());
    }
    return second_third_intersection.size() > first_fourth_intersection.size();
}
bool EdgeSplitPredicate::Check(const ScaffoldEdge& scaffold_edge) const {
    EdgeId start = scaffold_edge.getStart();
    EdgeId end = scaffold_edge.getEnd();
    size_t start_length = g_.length(start);
    size_t end_length = g_.length(end);
    vector<BarcodeId> first_half_of_start = barcode_extractor_.GetBarcodesFromRange(start, count_threshold_,
                                                                                    0, start_length / 2);
    vector<BarcodeId> second_half_of_start = barcode_extractor_.GetBarcodesFromRange(start, count_threshold_,
                                                                                    start_length - start_length / 2,
                                                                                    start_length);
    vector<BarcodeId> first_half_of_end = barcode_extractor_.GetBarcodesFromRange(end, count_threshold_, 0, end_length / 2);
    vector<BarcodeId> second_half_of_end = barcode_extractor_.GetBarcodesFromRange(end, count_threshold_,
                                                                                   end_length - end_length / 2,
                                                                                   end_length);
    bool next_conjugate_check = CheckOrderingForThreeSegments(second_half_of_start, first_half_of_end,
                                                              second_half_of_end, strictness_);
    if (not next_conjugate_check) {
        DEBUG("Next conjugate check failed");
    }
    bool previous_conjugate_check = CheckOrderingForThreeSegments(first_half_of_end, second_half_of_start,
                                                                  first_half_of_start, strictness_);
    if (not previous_conjugate_check) {
        DEBUG("Previous conjugate check failed");
    }
    bool previous_check = CheckOrderingForFourSegments(first_half_of_start, second_half_of_start,
                                                       first_half_of_end, second_half_of_end);
    if (not previous_check) {
        DEBUG("Previous check failed.");
    }
    return next_conjugate_check and previous_conjugate_check and previous_check;

}
EdgeInTheMiddlePredicate::EdgeInTheMiddlePredicate(const Graph& g_,
                                                   const barcode_index::FrameBarcodeIndexInfoExtractor& barcode_extractor_,
                                                   size_t count_threshold,
                                                   double shared_fraction_threshold)
    : g_(g_), barcode_extractor_(barcode_extractor_), count_threshold_(count_threshold),
      shared_fraction_threshold_(shared_fraction_threshold) {}
bool EdgeInTheMiddlePredicate::IsCorrectOrdering(const EdgeId& first, const EdgeId& second, const EdgeId& third) {
        vector<BarcodeId> first_barcodes = barcode_extractor_.GetBarcodesFromRange(first, count_threshold_, 0, g_.length(first));
        vector<BarcodeId> second_barcodes = barcode_extractor_.GetBarcodesFromRange(second, count_threshold_, 0, g_.length(second));
        vector<BarcodeId> third_barcodes = barcode_extractor_.GetBarcodesFromRange(third, count_threshold_, 0, g_.length(third));
        vector<BarcodeId> first_third_intersection;
        std::set_intersection(first_barcodes.begin(), first_barcodes.end(), third_barcodes.begin(), third_barcodes.end(),
                              std::back_inserter(first_third_intersection));
        vector<BarcodeId> all_intersection;
        std::set_intersection(first_third_intersection.begin(), first_third_intersection.end(), second_barcodes.begin(),
                              second_barcodes.end(), std::back_inserter(all_intersection));
        double shared_fraction = static_cast<double>(all_intersection.size()) / static_cast<double>(first_third_intersection.size());
        DEBUG("Second barcodes: " << second_barcodes.size());
        DEBUG("First third: " << first_third_intersection.size());
        DEBUG("Shared fraction: " << shared_fraction);
        return math::ge(shared_fraction, shared_fraction_threshold_);
}
NextFarEdgesPredicate::NextFarEdgesPredicate(const Graph& g_,
                                             const barcode_index::FrameBarcodeIndexInfoExtractor& barcode_extractor_,
                                             const size_t count_threshold_,
                                             const double shared_fraction_threshold_,
                                             const std::function<vector<NextFarEdgesPredicate::ScaffoldVertex>
                                                                     (NextFarEdgesPredicate::ScaffoldVertex)>& candidates_getter_)
    : g_(g_),
      barcode_extractor_(barcode_extractor_),
      count_threshold_(count_threshold_),
      shared_fraction_threshold_(shared_fraction_threshold_),
      candidates_getter_(candidates_getter_) {}

bool NextFarEdgesPredicate::Check(const ScaffoldEdgePredicate::ScaffoldEdge& scaffold_edge) const {
    EdgeInTheMiddlePredicate predicate(g_, barcode_extractor_, count_threshold_, shared_fraction_threshold_);
    ScaffoldVertex candidate = scaffold_edge.getEnd();
    ScaffoldVertex current_vertex = scaffold_edge.getStart();
    vector<ScaffoldVertex> other_candidates = candidates_getter_(current_vertex);
    for (const auto& other: other_candidates) {
        if (other != candidate and predicate.IsCorrectOrdering(current_vertex, other, candidate)) {
            return false;
        }
    }
    return true;
}



TransitiveEdgesPredicate::TransitiveEdgesPredicate(const scaffold_graph::ScaffoldGraph& graph, const Graph& g, size_t distance_threshold) :
    scaffold_graph_(graph), g_(g), distance_threshold_(distance_threshold) {}
bool TransitiveEdgesPredicate::Check(const ScaffoldEdgePredicate::ScaffoldEdge& scaffold_edge) const {
    ScaffoldVertex current = scaffold_edge.getStart();
    ScaffoldVertex candidate = scaffold_edge.getEnd();
    //fixme replace with dijkstra and length threshold
    DEBUG("Checking edge (" << current.int_id() << ", " << candidate.int_id() << ")");
    DEBUG("Lengths: " << g_.length(current) << ", " << g_.length(candidate));
    SimpleSearcher simple_searcher(scaffold_graph_, g_, distance_threshold_);
    auto reachable_vertices = simple_searcher.GetReachableVertices(current, scaffold_edge);
    for (const auto& vertex: reachable_vertices) {
        if (candidate == vertex) {
            DEBUG("Found another path, false");
            return false;
        }
    }
    DEBUG("True";)
    return true;
}
SimpleSearcher::SimpleSearcher(const scaffold_graph::ScaffoldGraph& graph_, const Graph& g,  size_t distance_)
    : scaff_graph_(graph_), g_(g), distance_threshold_(distance_) {}
vector<SimpleSearcher::ScaffoldVertex> SimpleSearcher::GetReachableVertices(const SimpleSearcher::ScaffoldVertex& vertex,
                                                                            const ScaffoldGraph::ScaffoldEdge& restricted_edge) {
    vector<ScaffoldVertex> result;
    VertexWithDistance new_vertex(vertex, 0);
    std::queue<VertexWithDistance> vertex_queue;
    vertex_queue.push(new_vertex);
    unordered_set<ScaffoldVertex> visited;
    while (not vertex_queue.empty()) {
        auto current_vertex = vertex_queue.front();
        vertex_queue.pop();
        DEBUG("Popped vertex with length: " << g_.length(current_vertex.vertex));
        DEBUG("Id: " << current_vertex.vertex.int_id());
        DEBUG("Distance: " << current_vertex.distance);
        if (current_vertex.distance <= distance_threshold_) {
            DEBUG("Passed threshold. Processing")
            ProcessVertex(vertex_queue, current_vertex, visited, restricted_edge);
            DEBUG("Processing finished");
            result.push_back(current_vertex.vertex);
        }
    }
    return result;
}

void SimpleSearcher::ProcessVertex(std::queue<VertexWithDistance>& vertex_queue,
                                   const VertexWithDistance& vertex,
                                   std::unordered_set<ScaffoldVertex>& visited,
                                   const ScaffoldGraph::ScaffoldEdge& restricted_edge) {
    size_t current_distance = vertex.distance;
    size_t new_distance = current_distance + 1;
    for (const ScaffoldGraph::ScaffoldEdge& edge: scaff_graph_.OutgoingEdges(vertex.vertex)) {
        DEBUG("Checking vertex: " << edge.getEnd().int_id());
        DEBUG("With length: " << g_.length(edge.getEnd()));
        DEBUG("Visited: " << (visited.find(edge.getEnd()) != visited.end()));
        DEBUG("Edge restricted: " << AreEqual(edge, restricted_edge));
        if (visited.find(edge.getEnd()) == visited.end() and not AreEqual(edge, restricted_edge)) {
            DEBUG("Passed");
            vertex_queue.emplace(edge.getEnd(), new_distance);
            visited.insert(edge.getEnd());
        }
    }
}
bool SimpleSearcher::AreEqual(const scaffold_graph::ScaffoldGraph::ScaffoldEdge& first,
                              const scaffold_graph::ScaffoldGraph::ScaffoldEdge& second) {
    return first.getStart() == second.getStart() and first.getEnd() == second.getEnd();
}

SimpleSearcher::VertexWithDistance::VertexWithDistance(const SimpleSearcher::ScaffoldVertex& vertex, size_t distance)
    : vertex(vertex), distance(distance) {}
}
