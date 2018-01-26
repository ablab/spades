#include "read_cloud_connection_conditions.hpp"
#include "read_cloud_dijkstras.hpp"
#include "path_extend_dijkstras.hpp"
#include "path_extender.hpp"
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
                stored_distances_[e].insert(make_pair(connected, static_cast<double>(dij.GetDistance(v))));
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
    size_t first_length = first.getLengthFromGraph(graph_);
    size_t second_length = second.getLengthFromGraph(graph_);
    DEBUG("First length: " << first_length);
    DEBUG("Second length: " << second_length);
    size_t first_size = barcode_extractor_->GetTailSize(first);
    size_t second_size = barcode_extractor_->GetHeadSize(second);
    if (first_size == 0 or second_size == 0) {
        DEBUG("No barcodes on one of the long edges");
        return 0.0;
    }
    size_t shared_count = barcode_extractor_->GetIntersectionSize(first, second);
    DEBUG("First size: " << first_size);
    DEBUG("Second size: " << second_size);
    DEBUG("Intersection: " << shared_count);
    size_t min_size = std::min(first_size, second_size);
    double containment_index = static_cast<double>(shared_count) / static_cast<double>(min_size);
    DEBUG("Score: " << containment_index);
    VERIFY(math::ge(1.0, containment_index));
//    double first_coverage = first.getCoverageFromGraph(graph_);
//    double second_coverage = second.getCoverageFromGraph(graph_);
    return containment_index;
}
NormalizedBarcodeScoreFunction::NormalizedBarcodeScoreFunction(
    const Graph &graph_,
    shared_ptr<barcode_index::ScaffoldVertexIndexInfoExtractor> barcode_extractor_,
    const size_t read_count_threshold_,
    const size_t tail_threshold_) : AbstractBarcodeScoreFunction(graph_, barcode_extractor_),
                                    read_count_threshold_(read_count_threshold_),
                                    tail_threshold_(tail_threshold_),
    //fixme use total barcodes from index
                                    total_barcodes_(1000000) {}

ReadCloudMiddleDijkstraPredicate::ReadCloudMiddleDijkstraPredicate(
        const Graph &g,
        const ScaffoldingUniqueEdgeStorage &unique_storage_,
        shared_ptr<barcode_index::SimpleIntersectingScaffoldVertexExtractor> short_edge_extractor,
        const shared_ptr<barcode_index::SimpleScaffoldVertexIndexInfoExtractor> long_edge_extractor,
        const ReadCloudMiddleDijkstraParams &params)
    : g(g),
      unique_storage_(unique_storage_),
      short_edge_extractor_(short_edge_extractor),
      long_edge_extractor_(long_edge_extractor),
      params_(params) {}

bool ReadCloudMiddleDijkstraPredicate::Check(const scaffold_graph::ScaffoldGraph::ScaffoldEdge &scaffold_edge) const {
    DEBUG("Checking edge " << scaffold_edge.getStart().int_id() << " -> " << scaffold_edge.getEnd().int_id());
    auto barcode_intersection =
        long_edge_extractor_->GetIntersection(scaffold_edge.getStart(), scaffold_edge.getEnd());
    DEBUG("Intersection size: " << barcode_intersection.size());
    auto long_gap_dijkstra = CreateLongGapCloserDijkstra(g, params_.distance_, unique_storage_,
                                                         short_edge_extractor_, long_edge_extractor_,
                                                         scaffold_edge.getStart(), scaffold_edge.getEnd(),
                                                         params_.edge_pair_gap_closer_params_);
    DEBUG("Created dijkstra");
    long_gap_dijkstra.Run(scaffold_edge.getStart().getEndGraphVertex(g));
    DEBUG("Dijkstra finished");
    VertexId target = scaffold_edge.getEnd().getStartGraphVertex(g);
    for (const auto &vertex: long_gap_dijkstra.ReachedVertices()) {
        if (vertex == target) {
            return true;
        }
    }
    return false;
}

EdgeSplitPredicate::EdgeSplitPredicate(const Graph &g_,
                                       const barcode_index::FrameBarcodeIndexInfoExtractor &barcode_extractor_,
                                       const size_t count_threshold_,
                                       double strictness)
    : g_(g_),
      barcode_extractor_(barcode_extractor_),
      count_threshold_(count_threshold_),
      strictness_(strictness) {}
bool EdgeSplitPredicate::CheckOrderingForThreeSegments(const vector<EdgeSplitPredicate::BarcodeId> &first,
                                                       const vector<EdgeSplitPredicate::BarcodeId> &second,
                                                       const vector<EdgeSplitPredicate::BarcodeId> &third,
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

bool EdgeSplitPredicate::CheckOrderingForFourSegments(const vector<EdgeSplitPredicate::BarcodeId> &first,
                                                      const vector<EdgeSplitPredicate::BarcodeId> &second,
                                                      const vector<EdgeSplitPredicate::BarcodeId> &third,
                                                      const vector<EdgeSplitPredicate::BarcodeId> &fourth) const {
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
bool EdgeSplitPredicate::Check(const ScaffoldEdge &scaffold_edge) const {
    path_extend::scaffold_graph::EdgeGetter getter;

    EdgeId start = getter.GetEdgeFromScaffoldVertex(scaffold_edge.getStart());
    EdgeId end = getter.GetEdgeFromScaffoldVertex(scaffold_edge.getEnd());
    size_t start_length = g_.length(start);
    size_t end_length = g_.length(end);
    vector<BarcodeId> first_half_of_start = barcode_extractor_.GetBarcodesFromRange(start, count_threshold_,
                                                                                    0, start_length / 2);
    vector<BarcodeId> second_half_of_start = barcode_extractor_.GetBarcodesFromRange(start, count_threshold_,
                                                                                     start_length - start_length / 2,
                                                                                     start_length);
    vector<BarcodeId>
        first_half_of_end = barcode_extractor_.GetBarcodesFromRange(end, count_threshold_, 0, end_length / 2);
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
EdgeInTheMiddlePredicate::EdgeInTheMiddlePredicate(const Graph &g_,
                                                   const barcode_index::FrameBarcodeIndexInfoExtractor &barcode_extractor_,
                                                   size_t count_threshold,
                                                   double shared_fraction_threshold)
    : g_(g_), barcode_extractor_(barcode_extractor_), count_threshold_(count_threshold),
      shared_fraction_threshold_(shared_fraction_threshold) {}
bool EdgeInTheMiddlePredicate::IsCorrectOrdering(const EdgeId &first, const EdgeId &second, const EdgeId &third) {
    vector<BarcodeId>
        first_barcodes = barcode_extractor_.GetBarcodesFromRange(first, count_threshold_, 0, g_.length(first));
    vector<BarcodeId>
        second_barcodes = barcode_extractor_.GetBarcodesFromRange(second, count_threshold_, 0, g_.length(second));
    vector<BarcodeId>
        third_barcodes = barcode_extractor_.GetBarcodesFromRange(third, count_threshold_, 0, g_.length(third));
    vector<BarcodeId> first_third_intersection;
    std::set_intersection(first_barcodes.begin(), first_barcodes.end(), third_barcodes.begin(), third_barcodes.end(),
                          std::back_inserter(first_third_intersection));
    vector<BarcodeId> all_intersection;
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
SimpleSearcher::SimpleSearcher(const scaffold_graph::ScaffoldGraph &graph_, const Graph &g, size_t distance_)
    : scaff_graph_(graph_), g_(g), distance_threshold_(distance_) {}
vector<SimpleSearcher::ScaffoldVertex> SimpleSearcher::GetReachableVertices(const SimpleSearcher::ScaffoldVertex &vertex,
                                                                            const ScaffoldGraph::ScaffoldEdge &restricted_edge) {
    vector<ScaffoldVertex> result;
    VertexWithDistance new_vertex(vertex, 0);
    std::queue<VertexWithDistance> vertex_queue;
    vertex_queue.push(new_vertex);
    unordered_set<ScaffoldVertex> visited;
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

AbstractBarcodeScoreFunction::AbstractBarcodeScoreFunction(const Graph &graph_,
                                                           const shared_ptr<barcode_index::ScaffoldVertexIndexInfoExtractor> barcode_extractor_)
    :
    graph_(graph_),
    barcode_extractor_(barcode_extractor_) {}
TrivialBarcodeScoreFunction::TrivialBarcodeScoreFunction(
    const Graph &graph_,
    shared_ptr<barcode_index::ScaffoldVertexIndexInfoExtractor> barcode_extractor_,
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
//    std::set<std::pair<size_t, size_t>> interesting_edges;
//    interesting_edges.emplace(419428075, 419445323);
//    interesting_edges.emplace(419430151, 419393178);
//    interesting_edges.emplace(419442929, 419440649);
//    interesting_edges.emplace(419442959, 419374051);
//
//    auto current_edge_pair = std::make_pair<size_t, size_t>(start_edge.int_id(), end_edge.int_id());
//    if (interesting_edges.find(current_edge_pair) == interesting_edges.end()) {
//        return false;
//    }

    auto extension_chooser = ConstructSimpleExtensionChooser();
    auto start = scaffold_edge.getStart();
    auto end = scaffold_edge.getEnd();
    auto pair_entry_extractor = make_shared<path_extend::TwoSetsBasedPairEntryProcessor>(
        long_edge_extractor_->GetTailEntry(start), long_edge_extractor_->GetHeadEntry(end), short_edge_extractor_);
    auto long_gap_cloud_predicate = make_shared<path_extend::LongEdgePairGapCloserPredicate>(gp_.g, short_edge_extractor_,
                                                                                             predicate_params_,
                                                                                             start, end,
                                                                                             pair_entry_extractor);
    size_t length_threshold = unique_storage_.min_length();

    auto length_predicate = make_shared<path_extend::LengthChecker>(length_threshold, gp_.g);
    auto scaffold_vertex_predicate =
        make_shared<path_extend::AndChecker>(length_predicate, long_gap_cloud_predicate);

    DEBUG("Checking edge " << start.int_id() << " -> " << end.int_id());

    auto multi_chooser =
        make_shared<path_extend::MultiExtensionChooser>(gp_.g, scaffold_vertex_predicate, extension_chooser);
    const auto &lib = params_.dataset_info.reads[params_.paired_lib_index_];
    shared_ptr<PairedInfoLibrary> paired_lib = MakeNewLib(gp_.g, lib, gp_.clustered_indices[params_.paired_lib_index_]);
    GraphCoverageMap cover_map(gp_.g);
    UsedUniqueStorage used_unique_storage(unique_storage_);
    const bool investigate_loops = true;
    const bool short_loop_resolver = false;
    auto opts = params_.pe_params_.pset.extension_options;

    path_extend::QueueContainer paths_container;
    BidirectionalPath* initial_path = start.toPath(gp_.g);
    paths_container.push(initial_path);
    path_extend::SearchingMultiExtender multi_extender(gp_, cover_map, used_unique_storage, multi_chooser,
                                                       paired_lib->GetISMax(), investigate_loops, short_loop_resolver,
                                                       opts.weight_threshold, length_bound_, paths_container);

    const size_t max_path_growing_iterations = 1000;
    const size_t max_paths_to_process = 200;
    const size_t max_edge_visits = 10;
    size_t path_processing_iterations = 0;

    VertexId start_vertex = start.getEndGraphVertex(gp_.g);
    VertexId target_vertex = end.getStartGraphVertex(gp_.g);
    DEBUG("Start: " << start_vertex.int_id());
    DEBUG("Target: " << target_vertex.int_id());

    if (start_vertex == target_vertex) {
        return true;
    }

    std::unordered_map<EdgeId, size_t> edge_to_number_of_visits;

    vector<BidirectionalPath*> tmp_path_storage;

    while (not paths_container.empty() and path_processing_iterations < max_paths_to_process) {
        DEBUG("Path container size: " << paths_container.size());
        DEBUG("Coverage map size: " << cover_map.size());
        auto current_path = paths_container.front();
        BidirectionalPath* path_copy = current_path;
        tmp_path_storage.push_back(path_copy);
        SubscribeCoverageMap(current_path, cover_map);
        paths_container.pop();
        EdgeId last_edge = current_path->Back();
        edge_to_number_of_visits[last_edge]++;
        size_t current_iterations = 0;
        if (edge_to_number_of_visits[last_edge] > max_edge_visits) {
            DEBUG("Edge was visited too often");
            continue;
        }
        PathContainer *empty_storage = nullptr;
        while (multi_extender.MakeGrowStep(*current_path, empty_storage)
            and current_iterations < max_path_growing_iterations) {
            ++current_iterations;
        }
        ++path_processing_iterations;
        if (current_iterations >= max_path_growing_iterations) {
            WARN("Stopped path growing because of too many steps");
            continue;
        }
        const auto &reached_vertices = multi_extender.GetReachedVertices();
        if (reached_vertices.find(target_vertex) != reached_vertices.end()) {
            DEBUG("Found target");
            return true;
        }
    }
    if (path_processing_iterations >= max_paths_to_process) {
        DEBUG("Had to process too many paths, returning");
        return true;
    }

    DEBUG("Could not find supported path");
    return false;
}

shared_ptr<path_extend::ExtensionChooser> CompositeConnectionPredicate::ConstructSimpleExtensionChooser() const {
    auto opts = params_.pe_params_.pset.extension_options;
    const size_t lib_index = params_.paired_lib_index_;
    const auto &lib = params_.dataset_info.reads[lib_index];
    shared_ptr<PairedInfoLibrary> paired_lib = MakeNewLib(gp_.g, lib, clustered_indices_[lib_index]);
    shared_ptr<CoverageAwareIdealInfoProvider> iip = nullptr;
    path_extend::PELaunchSupport support(params_.dataset_info, params_.pe_params_);
    if (opts.use_default_single_threshold) {
        if (params_.pe_params_.uneven_depth) {
            iip = make_shared<CoverageAwareIdealInfoProvider>(gp_.g, paired_lib, params_.dataset_info.RL());
        } else {
            double lib_cov = support.EstimateLibCoverage(lib_index);
            INFO("Estimated coverage of library #" << lib_index << " is " << lib_cov);
            iip = make_shared<GlobalCoverageAwareIdealInfoProvider>(gp_.g,
                                                                    paired_lib,
                                                                    params_.dataset_info.RL(),
                                                                    lib_cov);
        }
    }

    auto wc = make_shared<PathCoverWeightCounter>(gp_.g, paired_lib, params_.pe_params_.pset.normalize_weight,
                                                  support.SingleThresholdForLib(params_.pe_params_.pset,
                                                                                lib.data().pi_threshold),
                                                  iip);
    auto extension_chooser = make_shared<SimpleExtensionChooser>(gp_.g, wc, opts.weight_threshold, opts.priority_coeff);
    return extension_chooser;
}
CompositeConnectionPredicate::CompositeConnectionPredicate(
        const conj_graph_pack &gp_,
        shared_ptr<barcode_index::SimpleIntersectingScaffoldVertexExtractor> short_edge_extractor,
        shared_ptr<barcode_index::SimpleScaffoldVertexIndexInfoExtractor> barcode_extractor_,
        const ScaffoldingUniqueEdgeStorage &unique_storage_,
        const de::PairedInfoIndicesT<debruijn_graph::DeBruijnGraph> &clustered_indices_,
        const size_t length_bound_,
        const CompositeConnectionParams &params_,
        const LongEdgePairGapCloserParams &predicate_params_) :
    gp_(gp_),
    short_edge_extractor_(short_edge_extractor),
    long_edge_extractor_(barcode_extractor_),
    unique_storage_(unique_storage_),
    clustered_indices_(clustered_indices_),
    length_bound_(length_bound_),
    params_(params_),
    predicate_params_(predicate_params_) {}

CompositeConnectionParams::CompositeConnectionParams(const size_t paired_lib_index_,
                                                     const size_t prefix_length_,
                                                     const config::dataset &dataset_info,
                                                     const PathExtendParamsContainer &pe_params_) : paired_lib_index_(
    paired_lib_index_), prefix_length_(prefix_length_), dataset_info(dataset_info), pe_params_(pe_params_) {}

path_extend::ReadCloudMiddleDijkstraParams::ReadCloudMiddleDijkstraParams(const size_t count_threshold_,
                                                                          const size_t tail_threshold_,
                                                                          const size_t distance_,
                                                                          const path_extend::LongEdgePairGapCloserParams &edge_pair_gap_closer_params_)
    : count_threshold_(count_threshold_),
      tail_threshold_(tail_threshold_),
      distance_(distance_),
      edge_pair_gap_closer_params_(edge_pair_gap_closer_params_) {}
}

