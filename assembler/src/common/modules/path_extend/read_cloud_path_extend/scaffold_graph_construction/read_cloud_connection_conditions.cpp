//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "read_cloud_connection_conditions.hpp"

#include "common/assembly_graph/dijkstra/read_cloud_dijkstra/read_cloud_dijkstras.hpp"
#include "common/assembly_graph/dijkstra/read_cloud_dijkstra/path_extend_dijkstras.hpp"

namespace path_extend {
namespace read_cloud {

std::map<EdgeId, double> AssemblyGraphUniqueConnectionCondition::ConnectedWith(EdgeId e) const {
    VERIFY_MSG(interesting_edge_set_.find(e) != interesting_edge_set_.end(),
               " edge " << e.int_id() << " not applicable for connection condition");
    if (stored_distances_.find(e) != stored_distances_.end()) {
        return stored_distances_[e];
    }
    stored_distances_.insert(make_pair(e, std::map<debruijn_graph::EdgeId, double>()));
    for (auto connected: g_.OutgoingEdges(g_.EdgeEnd(e))) {
        if (interesting_edge_set_.find(connected) != interesting_edge_set_.end()) {
            stored_distances_[e].insert(std::make_pair(connected, 1));
        }
    }

    ReadCloudDijkstraHelper helper;
    auto dij = helper.CreateUniqueDijkstra(g_, max_connection_length_, unique_storage_);
    dij.Run(g_.EdgeEnd(e));
    for (auto v: dij.ReachedVertices()) {
        for (auto connected: g_.OutgoingEdges(v)) {
            if (interesting_edge_set_.find(connected) != interesting_edge_set_.end() &&
                dij.GetDistance(v) < max_connection_length_ && connected != g_.conjugate(e)) {
                stored_distances_[e].insert(std::make_pair(connected, static_cast<double>(dij.GetDistance(v))));
            }
        }
    }
    return stored_distances_[e];
}
AssemblyGraphUniqueConnectionCondition::AssemblyGraphUniqueConnectionCondition(const Graph &g,
                                                                               size_t max_connection_length,
                                                                               const ScaffoldingUniqueEdgeStorage &unique_edges)
    : AssemblyGraphConnectionCondition(g, max_connection_length, unique_edges), unique_storage_(unique_edges) {}
bool AssemblyGraphUniqueConnectionCondition::IsLast() const {
    return false;
}

double NormalizedBarcodeScoreFunction::GetScore(const scaffold_graph::ScaffoldGraph::ScaffoldEdge &edge) const {
    auto first = edge.getStart();
    auto second = edge.getEnd();
    DEBUG("Checking edge " << edge.getStart().int_id() << " -> " << edge.getEnd().int_id());
    size_t first_length = first.GetLengthFromGraph(graph_);
    size_t second_length = second.GetLengthFromGraph(graph_);
    size_t first_size = barcode_extractor_->GetTailSize(first);
    size_t second_size = barcode_extractor_->GetHeadSize(second);
    if (first_size == 0 or second_size == 0) {
        DEBUG("No barcodes on one of the long edges");
        return 0.0;
    }
    size_t shared_count = barcode_extractor_->GetIntersectionSize(first, second);
    size_t min_size = std::min(first_size, second_size);
    double containment_index = static_cast<double>(shared_count) / static_cast<double>(min_size);
    if (math::ge(containment_index, 0.05)) {
        DEBUG("First length: " << first_length);
        DEBUG("Second length: " << second_length);
        DEBUG("First size: " << first_size);
        DEBUG("Second size: " << second_size);
        DEBUG("Intersection: " << shared_count);
        DEBUG("Score: " << containment_index);
    }
    VERIFY(math::ge(1.0, containment_index));
//    double first_coverage = first.GetCoverageFromGraph(graph_);
//    double second_coverage = second.GetCoverageFromGraph(graph_);
    return containment_index;
}
NormalizedBarcodeScoreFunction::NormalizedBarcodeScoreFunction(
        const Graph &graph_,
        std::shared_ptr<barcode_index::ScaffoldVertexIndexInfoExtractor> barcode_extractor_) :
    AbstractBarcodeScoreFunction(graph_, barcode_extractor_) {}

ReadCloudMiddleDijkstraPredicate::ReadCloudMiddleDijkstraPredicate(
        const Graph &g,
        const ScaffoldingUniqueEdgeStorage &unique_storage_,
        std::shared_ptr<barcode_index::SimpleIntersectingScaffoldVertexExtractor> short_edge_extractor,
        const std::shared_ptr<barcode_index::SimpleScaffoldVertexIndexInfoExtractor> long_edge_extractor,
        const ReadCloudMiddleDijkstraParams &params)
    : g(g),
      unique_storage_(unique_storage_),
      short_edge_extractor_(short_edge_extractor),
      long_edge_extractor_(long_edge_extractor),
      params_(params) {}

bool ReadCloudMiddleDijkstraPredicate::Check(const scaffold_graph::ScaffoldGraph::ScaffoldEdge &scaffold_edge) const {
    DEBUG("Checking edge " << scaffold_edge.getStart().int_id() << " -> " << scaffold_edge.getEnd().int_id());
    ReadCloudDijkstraHelper helper;
    auto long_gap_dijkstra = helper.CreateLongGapCloserDijkstra(g, params_.distance_, unique_storage_,
                                                                short_edge_extractor_, long_edge_extractor_,
                                                                scaffold_edge.getStart(), scaffold_edge.getEnd(),
                                                                params_.edge_pair_gap_closer_params_);
    DEBUG("Created dijkstra");
    long_gap_dijkstra.Run(scaffold_edge.getStart().GetEndGraphVertex(g));
    DEBUG("Dijkstra finished");
    VertexId target = scaffold_edge.getEnd().GetStartGraphVertex(g);
    for (const auto &vertex: long_gap_dijkstra.ReachedVertices()) {
        if (vertex == target) {
            return true;
        }
    }
    return false;
}

EdgeSplitPredicate::EdgeSplitPredicate(
        const Graph &g_,
        std::shared_ptr<barcode_index::SimpleScaffoldVertexIndexInfoExtractor> barcode_extractor,
        const size_t count_threshold,
        double strictness)
    : g_(g_),
      barcode_extractor_(barcode_extractor),
      count_threshold_(count_threshold),
      strictness_(strictness) {}
bool EdgeSplitPredicate::CheckOrderingForThreeSegments(const barcode_index::SimpleVertexEntry &first,
                                                       const barcode_index::SimpleVertexEntry &second,
                                                       const barcode_index::SimpleVertexEntry &third,
                                                       double strictness) const {
    std::vector<BarcodeId> first_second_intersection;
    std::set_intersection(first.begin(), first.end(), second.begin(), second.end(),
                          std::back_inserter(first_second_intersection));
    std::vector<BarcodeId> first_third_intersection;
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

bool EdgeSplitPredicate::CheckOrderingForFourSegments(const barcode_index::SimpleVertexEntry &first,
                                                      const barcode_index::SimpleVertexEntry &second,
                                                      const barcode_index::SimpleVertexEntry &third,
                                                      const barcode_index::SimpleVertexEntry &fourth) const {
    std::vector<BarcodeId> first_fourth_intersection;
    std::set_intersection(first.begin(), first.end(), fourth.begin(), fourth.end(),
                          std::back_inserter(first_fourth_intersection));
    std::vector<BarcodeId> second_third_intersection;
    std::set_intersection(second.begin(), second.end(), third.begin(), third.end(),
                          std::back_inserter(second_third_intersection));
    if (second_third_intersection.size() <= first_fourth_intersection.size()) {
        DEBUG("Second third: " << second_third_intersection.size());
        DEBUG("First fourth: " << first_fourth_intersection.size());
    }
    return second_third_intersection.size() > first_fourth_intersection.size();
}
bool EdgeSplitPredicate::Check(const ScaffoldEdge &scaffold_edge) const {
    auto start = scaffold_edge.getStart();
    auto end = scaffold_edge.getEnd();
    auto first_half_of_start = barcode_extractor_->GetHeadEntry(start);
    auto second_half_of_start = barcode_extractor_->GetTailEntry(start);
    auto first_half_of_end = barcode_extractor_->GetHeadEntry(end);
    auto second_half_of_end = barcode_extractor_->GetTailEntry(end);
    DEBUG("First half of start: " << first_half_of_start.size());
    DEBUG("Second half of start: " << second_half_of_start.size());
    DEBUG("First half of end: " << first_half_of_end.size());
    DEBUG("Second half of end: " << second_half_of_end.size());
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
EdgeInTheMiddlePredicate::EdgeInTheMiddlePredicate(const Graph &g_,
                                                   const barcode_index::FrameBarcodeIndexInfoExtractor &barcode_extractor,
                                                   size_t count_threshold,
                                                   double shared_fraction_threshold)
    : g_(g_), barcode_extractor_(barcode_extractor), count_threshold_(count_threshold),
      shared_fraction_threshold_(shared_fraction_threshold) {}
bool EdgeInTheMiddlePredicate::IsCorrectOrdering(const EdgeId &first, const EdgeId &second, const EdgeId &third) {
    std::vector<BarcodeId>
        first_barcodes = barcode_extractor_.GetBarcodesFromRange(first, count_threshold_, 0, g_.length(first));
    std::vector<BarcodeId>
        second_barcodes = barcode_extractor_.GetBarcodesFromRange(second, count_threshold_, 0, g_.length(second));
    std::vector<BarcodeId>
        third_barcodes = barcode_extractor_.GetBarcodesFromRange(third, count_threshold_, 0, g_.length(third));
    std::vector<BarcodeId> first_third_intersection;
    std::set_intersection(first_barcodes.begin(), first_barcodes.end(), third_barcodes.begin(), third_barcodes.end(),
                          std::back_inserter(first_third_intersection));
    std::vector<BarcodeId> all_intersection;
    std::set_intersection(first_third_intersection.begin(), first_third_intersection.end(), second_barcodes.begin(),
                          second_barcodes.end(), std::back_inserter(all_intersection));
    double shared_fraction =
        static_cast<double>(all_intersection.size()) / static_cast<double>(first_third_intersection.size());
    DEBUG("Second barcodes: " << second_barcodes.size());
    DEBUG("First third: " << first_third_intersection.size());
    DEBUG("Shared fraction: " << shared_fraction);
    return math::ge(shared_fraction, shared_fraction_threshold_);
}

TransitiveEdgesPredicate::TransitiveEdgesPredicate(const scaffold_graph::ScaffoldGraph &graph,
                                                   const Graph &g,
                                                   size_t distance_threshold) :
    scaffold_graph_(graph), g_(g), distance_threshold_(distance_threshold) {}
bool TransitiveEdgesPredicate::Check(const ScaffoldEdgePredicate::ScaffoldEdge &scaffold_edge) const {
    ScaffoldVertex current = scaffold_edge.getStart();
    ScaffoldVertex candidate = scaffold_edge.getEnd();
    //fixme replace with dijkstra and length threshold
    DEBUG("Checking edge (" << current.int_id() << ", " << candidate.int_id() << ")");
    SimpleSearcher simple_searcher(scaffold_graph_, g_, distance_threshold_);
    auto reachable_vertices = simple_searcher.GetReachableVertices(current, scaffold_edge);
    for (const auto &vertex: reachable_vertices) {
        if (candidate == vertex) {
            DEBUG("Found another path, false");
            return false;
        }
    }
    DEBUG("True";)
    return true;
}
SimpleSearcher::SimpleSearcher(const scaffold_graph::ScaffoldGraph &graph, const Graph &g, size_t distance)
    : scaff_graph_(graph), g_(g), distance_threshold_(distance) {}
std::vector<SimpleSearcher::ScaffoldVertex> SimpleSearcher::GetReachableVertices(
        const SimpleSearcher::ScaffoldVertex &vertex,
        const ScaffoldGraph::ScaffoldEdge &restricted_edge) {
    std::vector<ScaffoldVertex> result;
    VertexWithDistance new_vertex(vertex, 0);
    std::queue<VertexWithDistance> vertex_queue;
    vertex_queue.push(new_vertex);
    std::unordered_set<ScaffoldVertex> visited;
    visited.insert(vertex);
    visited.insert(vertex.GetConjugateFromGraph(g_));
    visited.insert(restricted_edge.getEnd().GetConjugateFromGraph(g_));
    while (not vertex_queue.empty()) {
        auto current_vertex = vertex_queue.front();
        vertex_queue.pop();
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

void SimpleSearcher::ProcessVertex(std::queue<VertexWithDistance> &vertex_queue,
                                   const VertexWithDistance &vertex,
                                   std::unordered_set<ScaffoldVertex> &visited,
                                   const ScaffoldGraph::ScaffoldEdge &restricted_edge) {
    size_t current_distance = vertex.distance;
    size_t new_distance = current_distance + 1;
    for (const ScaffoldGraph::ScaffoldEdge &edge: scaff_graph_.OutgoingEdges(vertex.vertex)) {
        DEBUG("Checking vertex: " << edge.getEnd().int_id());
        DEBUG("Visited: " << (visited.find(edge.getEnd()) != visited.end()));
        DEBUG("Edge restricted: " << AreEqual(edge, restricted_edge));
        if (visited.find(edge.getEnd()) == visited.end() and not AreEqual(edge, restricted_edge)) {
            DEBUG("Passed");
            vertex_queue.emplace(edge.getEnd(), new_distance);
            visited.insert(edge.getEnd());
        }
    }
}
bool SimpleSearcher::AreEqual(const scaffold_graph::ScaffoldGraph::ScaffoldEdge &first,
                              const scaffold_graph::ScaffoldGraph::ScaffoldEdge &second) {
    return first.getStart() == second.getStart() and first.getEnd() == second.getEnd();
}

SimpleSearcher::VertexWithDistance::VertexWithDistance(const SimpleSearcher::ScaffoldVertex &vertex, size_t distance)
    : vertex(vertex), distance(distance) {}

AbstractBarcodeScoreFunction::AbstractBarcodeScoreFunction(
        const Graph &graph_,
        const std::shared_ptr<barcode_index::ScaffoldVertexIndexInfoExtractor> barcode_extractor)
    :
    graph_(graph_),
    barcode_extractor_(barcode_extractor) {}
TrivialBarcodeScoreFunction::TrivialBarcodeScoreFunction(
    const Graph &graph_,
    std::shared_ptr<barcode_index::ScaffoldVertexIndexInfoExtractor> barcode_extractor_,
    const size_t read_count_threshold_,
    const size_t tail_threshold_) : AbstractBarcodeScoreFunction(graph_,
                                                                 barcode_extractor_),
                                    read_count_threshold_(read_count_threshold_),
                                    tail_threshold_(tail_threshold_) {}
double TrivialBarcodeScoreFunction::GetScore(const scaffold_graph::ScaffoldGraph::ScaffoldEdge &edge) const {
    size_t shared_count = barcode_extractor_->GetIntersectionSize(edge.getStart(), edge.getEnd());

    return static_cast<double>(shared_count);
}
bool CompositeConnectionPredicate::Check(const scaffold_graph::ScaffoldGraph::ScaffoldEdge &scaffold_edge) const {
    DEBUG("Start composite check");
    auto start = scaffold_edge.getStart();
    auto end = scaffold_edge.getEnd();
    DEBUG("Checking edge " << start.int_id() << " -> " << end.int_id());
    const auto &start_entry = long_edge_extractor_->GetTailEntry(start);
    const auto &end_entry = long_edge_extractor_->GetHeadEntry(end);
    auto short_edge_score_function = std::make_shared<RepetitiveVertexEntryScoreFunction>(short_edge_extractor_);
    auto pair_entry_processor = std::make_shared<TwoSetsBasedPairEntryProcessor>(start_entry, end_entry,
                                                                                 short_edge_score_function);
    auto scaffold_vertex_predicate = ConstructScaffoldVertexPredicate(start, end, pair_entry_processor);
    auto pe_extension_chooser = search_parameter_pack_.extension_chooser;
    auto multi_chooser = std::make_shared<path_extend::PredicateExtensionChooser>(gp_.g, scaffold_vertex_predicate,
                                                                                  pe_extension_chooser);
    ExtenderSearcher extender_searcher(gp_, multi_chooser, search_parameter_pack_.search_params,
                                       search_parameter_pack_.searching_extender_params, length_bound_);

    BidirectionalPath *initial_path = start.ToPath(gp_.g);
    VertexId target_vertex = end.GetStartGraphVertex(gp_.g);
    return extender_searcher.IsReachable(target_vertex, initial_path);
}

CompositeConnectionPredicate::CompositeConnectionPredicate(
        const conj_graph_pack &gp,
        std::shared_ptr<barcode_index::SimpleIntersectingScaffoldVertexExtractor> short_edge_extractor,
        std::shared_ptr<barcode_index::SimpleScaffoldVertexIndexInfoExtractor> barcode_extractor,
        const ScaffoldingUniqueEdgeStorage &unique_storage,
        const size_t length_bound,
        const ReadCloudSearchParameterPack &search_parameter_pack,
        const LongEdgePairGapCloserParams &predicate_params, bool scaffolding_mode) :
    gp_(gp),
    short_edge_extractor_(short_edge_extractor),
    long_edge_extractor_(barcode_extractor),
    unique_storage_(unique_storage),
    length_bound_(length_bound),
    search_parameter_pack_(search_parameter_pack),
    predicate_params_(predicate_params),
    scaffolding_mode_(scaffolding_mode) {}
std::shared_ptr<ScaffoldVertexPredicate> CompositeConnectionPredicate::ConstructScaffoldVertexPredicate(
        const ScaffoldVertex &start, const ScaffoldVertex &end,
        std::shared_ptr<PairEntryProcessor> entry_processor) const {
    auto long_gap_cloud_predicate = std::make_shared<LongEdgePairGapCloserPredicate>(gp_.g, short_edge_extractor_,
                                                                                     predicate_params_,
                                                                                     start, end,
                                                                                     entry_processor);
    size_t length_threshold = unique_storage_.min_length();

    auto length_predicate = std::make_shared<LengthChecker>(length_threshold, gp_.g);
    auto scaffold_vertex_predicate = std::make_shared<AndPredicate>(length_predicate, long_gap_cloud_predicate);
    return scaffold_vertex_predicate;
}

ReadCloudMiddleDijkstraParams::ReadCloudMiddleDijkstraParams(const size_t count_threshold_,
                                                             const size_t tail_threshold_,
                                                             const size_t distance_,
                                                             const LongEdgePairGapCloserParams &edge_pair_gap_closer_params_)
    : count_threshold_(count_threshold_),
      tail_threshold_(tail_threshold_),
      distance_(distance_),
      edge_pair_gap_closer_params_(edge_pair_gap_closer_params_) {}
}
}