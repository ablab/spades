#include <stack>
#include "read_cloud_path_extend/intermediate_scaffolding/scaffold_graph_polisher.hpp"
#include "common/assembly_graph/dijkstra/read_cloud_dijkstra/scaffold_graph_dijkstra.hpp"
#include "common/assembly_graph/dijkstra/read_cloud_dijkstra/path_extend_dijkstras.hpp"
#include "common/barcode_index/scaffold_vertex_index_builder.hpp"
#include "scaffold_graph_path_cleaner.hpp"
#include "modules/path_extend/read_cloud_path_extend/cluster_storage/path_cluster_storage_builder.hpp"

namespace path_extend {
namespace read_cloud {

GapCloserUtils::SimpleTransitionGraph GapCloserUtils::RemoveDisconnectedVertices(
        const ScaffoldGraphPolisher::SimpleTransitionGraph &graph, const ScaffoldVertex &source,
        const ScaffoldVertex &sink) const {
    SimpleTransitionGraph result;
    DEBUG("Removing disconnected vertices");
    ForwardReachabilityChecker forward_checker(graph);
    BackwardReachabilityChecker backward_checker(graph);
    forward_checker.Run(source, sink);
    auto passed_forward = forward_checker.GetPassedVertices();
    for (auto it = graph.begin(); it != graph.end(); ++it) {
        auto vertex = *it;
        DEBUG("Checking vertex: " << vertex.int_id());
        if (passed_forward.find(vertex) != passed_forward.end()) {
            DEBUG("Passed");
            result.AddVertex(vertex);
        }
    }
    for (auto it = result.begin(); it != result.end(); ++it) {
        auto vertex = *it;
        for (auto edge_it = graph.outcoming_begin(vertex); edge_it != graph.outcoming_end(vertex); ++edge_it) {
            auto next = *edge_it;
            if (passed_forward.find(next) != passed_forward.end()) {
                DEBUG("Adding edge: (" << vertex.int_id() << ", " << next.int_id() << ")");
                result.AddEdge(vertex, next);
            }
        }
    }
    return result;
}

CloudScaffoldSubgraphExtractor::SimpleGraphT CloudScaffoldSubgraphExtractor::ExtractSubgraphBetweenVertices(
    const CloudScaffoldSubgraphExtractor::ScaffoldGraph &scaffold_graph,
    const CloudScaffoldSubgraphExtractor::ScaffoldVertex &first,
    const CloudScaffoldSubgraphExtractor::ScaffoldVertex &second) const {
    DEBUG("Extracting scaffold subgraph");
    DEBUG("First: " << first.int_id());
    DEBUG("Second: " << second.int_id());
    DEBUG("First length: " << first.GetLengthFromGraph(g_));
    DEBUG("Second length: " << second.GetLengthFromGraph(g_));
    DEBUG("First coverage: " << first.GetCoverageFromGraph(g_));
    DEBUG("Second coverage: " << second.GetCoverageFromGraph(g_));
    SimpleGraphT result;
    unordered_set<ScaffoldVertex> forward_vertices;
    unordered_set<ScaffoldVertex> backward_vertices;
    unordered_set<ScaffoldVertex> subgraph_vertices;
    ScaffoldGraph::ScaffoldEdge edge(first, second);
    LongEdgePairGapCloserParams params(params_.count_threshold_, params_.large_length_threshold_,
                                       params_.share_threshold_, params.relative_coverage_threshold_,
                                       params_.small_length_threshold_, true);
    auto start = edge.getStart();
    auto end = edge.getEnd();

    DEBUG("Checking sizes");
    DEBUG("First size: " << scaff_vertex_extractor_->GetHeadSize(start));
    DEBUG("Second size: " << scaff_vertex_extractor_->GetHeadSize(end));
    DEBUG("Getting intersection");
    auto intersection_entry = scaff_vertex_extractor_->GetIntersection(start, end);
    DEBUG("Got intersection");
    auto pair_entry_extractor = make_shared<IntersectionBasedPairEntryProcessor>(intersection_entry,
                                                                                 scaff_vertex_extractor_);
    DEBUG("Got extractor");
    auto gap_closer_predicate = make_shared<LongEdgePairGapCloserPredicate>(g_,
                                                                            scaff_vertex_extractor_,
                                                                            params, start, end,
                                                                            pair_entry_extractor);
    DEBUG("Constructed predicates");
    omnigraph::ScaffoldDijkstraHelper helper;
    auto forward_dijkstra = helper.CreateForwardBoundedScaffoldDijkstra(scaffold_graph, first, second,
                                                                        params_.distance_threshold_,
                                                                        gap_closer_predicate);
    auto backward_dijkstra = helper.CreateBackwardBoundedScaffoldDijkstra(scaffold_graph, first, second,
                                                                          params_.distance_threshold_,
                                                                          gap_closer_predicate);
    DEBUG("Running dijkstra")
    forward_dijkstra.Run(first);
    //fixme avoid copying
    for (const auto &vertex: forward_dijkstra.ReachedVertices()) {
        TRACE("Adding forward vertex to subgraph: " << vertex.int_id());
        if (CheckSubGraphVertex(vertex, first, second)) {
            subgraph_vertices.insert(vertex);
        }
    }
    backward_dijkstra.Run(second);
    for (const auto &vertex: backward_dijkstra.ReachedVertices()) {
        TRACE("Adding backward vertex to subgraph: " << vertex.int_id());
        if (CheckSubGraphVertex(vertex, first, second)) {
            subgraph_vertices.insert(vertex);
        }
    }
    subgraph_vertices.insert(first);
    subgraph_vertices.insert(second);
    for (const auto &vertex: subgraph_vertices) {
        result.AddVertex(vertex);
    }
    unordered_set<ScaffoldVertex> intersection;
    for (const auto &vertex: forward_dijkstra.ReachedVertices()) {
        if (backward_dijkstra.DistanceCounted(vertex)) {
            intersection.insert(vertex);
        }
    }
    bool target_reached = intersection.size() > 0;
    DEBUG("Target reached: " << (target_reached ? "True" : "False"));
    DEBUG(subgraph_vertices.size() << " vertices in subgraph");
    for (const ScaffoldEdge &scaffold_edge: scaffold_graph.edges()) {
        if (CheckSubgraphEdge(scaffold_edge, first, second, subgraph_vertices)) {
            DEBUG("Adding edge: " << scaffold_edge.getStart().int_id() << ", " << scaffold_edge.getEnd().int_id());
            result.AddEdge(scaffold_edge.getStart(), scaffold_edge.getEnd());
        }
    }
//    GapCloserUtils utils;
//    auto cleaned_graph = utils.RemoveDisconnectedVertices(result, first, second);
    DEBUG(result.size() << " vertices in cleaned subgraph");
    DEBUG(result.GetEdgesCount() << " edges in cleaned subgraph");
    if (result.GetEdgesCount() + 1 > result.size()) {
        DEBUG("Complex subgraph");
    }
    if (result.GetEdgesCount() + 1 < result.size()) {
        DEBUG("Broken subgraph");
    }
    return result;
}
CloudScaffoldSubgraphExtractor::CloudScaffoldSubgraphExtractor(
    const Graph &g,
    shared_ptr<barcode_index::SimpleScaffoldVertexIndexInfoExtractor> extractor,
    const CloudSubgraphExtractorParams &params)
    : g_(g),
      scaff_vertex_extractor_(extractor),
      params_(params) {}
bool CloudScaffoldSubgraphExtractor::CheckSubgraphEdge(const ScaffoldEdge &edge,
                                                       const ScaffoldVertex &first,
                                                       const ScaffoldVertex &second,
                                                       const unordered_set<ScaffoldVertex> &subgraph_vertices) const {
    return subgraph_vertices.find(edge.getStart()) != subgraph_vertices.end() and
        subgraph_vertices.find(edge.getEnd()) != subgraph_vertices.end() and
        edge.getStart() != edge.getEnd() and edge.getStart() != second and edge.getEnd() != first;
}

bool CloudScaffoldSubgraphExtractor::CheckSubGraphVertex(const CloudScaffoldSubgraphExtractor::ScaffoldVertex &vertex,
                                                         const CloudScaffoldSubgraphExtractor::ScaffoldVertex &first,
                                                         const CloudScaffoldSubgraphExtractor::ScaffoldVertex &second) const {
    return vertex != first.GetConjugateFromGraph(g_) and vertex != second.GetConjugateFromGraph(g_);
}
ScaffoldGraph ScaffoldGraphPolisher::CleanSmallGraphUsingLargeGraph(
    const ScaffoldGraphPolisher::ScaffoldGraph &large_scaffold_graph,
    const ScaffoldGraphPolisher::ScaffoldGraph &small_scaffold_graph) const {
    ScaffoldGraphExtractor extractor;
    auto univocal_edges = extractor.ExtractReliableEdges(large_scaffold_graph);
    auto current_graph = small_scaffold_graph;
    DEBUG("Extracting paths");
    auto extracted_paths = ExtractPathsWithinUnivocal(current_graph, univocal_edges);
    INFO("Found " << extracted_paths.size() << " paths");
    DEBUG("Cleaning graph");
    ScaffoldGraphPathCleaner path_cleaner;
    path_cleaner.CleanScaffoldGraphUsingPaths(current_graph, extracted_paths);
    DEBUG("Cleaned graph");
    return current_graph;
}

ScaffoldGraphPolisher::InternalPaths ScaffoldGraphPolisher::ExtractPathsWithinUnivocal(
    const ScaffoldGraph &input_graph,
    const vector<ScaffoldGraphPolisher::ScaffoldEdge> &univocal_edges) const {
    DEBUG("Getting inserted path connections");
    InternalPaths result;
    for (const auto &edge: univocal_edges) {
        ScaffoldVertex source = edge.getStart();
        ScaffoldVertex sink = edge.getEnd();
        CloudScaffoldSubgraphExtractor subgraph_extractor(g_, scaff_vertex_extractor_, subgraph_extractor_params_);
        DEBUG("Extracting subgraph");
        const auto subgraph = subgraph_extractor.ExtractSubgraphBetweenVertices(input_graph, source, sink);
        DEBUG("Extracting paths");
        InternalPaths extracted_paths = path_extractor_->GetCorrectPaths(subgraph, source, sink);
        if (extracted_paths.size() == 1) {
            result.push_back(extracted_paths[0]);
        }
    }
    return result;
}

ScaffoldGraphPolisher::ScaffoldGraphPolisher(
    const Graph &g_,
    shared_ptr<barcode_index::SimpleScaffoldVertexIndexInfoExtractor> scaff_vertex_extractor,
    shared_ptr<CorrectPathExtractor> path_extractor,
    const CloudSubgraphExtractorParams &subgraph_extractor_params)
    : g_(g_),
      scaff_vertex_extractor_(scaff_vertex_extractor),
      path_extractor_(path_extractor),
      subgraph_extractor_params_(subgraph_extractor_params) {}

void ReachabilityChecker::Run(const VertexT &start, const VertexT &target) {
    std::unordered_set<VertexT> reached_vertices;
    DEBUG("Checking reachability for target: " << target.int_id());
    DEBUG("Starting processing from vertex " << start.int_id());
    ProcessVertex(start, target);
}

bool ReachabilityChecker::ProcessVertex(const ReachabilityChecker::VertexT &vertex,
                                        const ReachabilityChecker::VertexT &target) {
    DEBUG("Processing vertex: " << vertex.int_id());
    visited_.insert(vertex);
    bool result = false;
    if (vertex == target or passed_.find(vertex) != passed_.end()) {
        return true;
    }
    for (auto it = GetBeginIterator(vertex); it != GetEndIterator(vertex); ++it) {
        auto next = *it;
        DEBUG("Checking neighbour: " << next.int_id());
        if (next == target) {
            DEBUG("Found target");
            passed_.insert(vertex);
            DEBUG("Inserting " << vertex.int_id());
            passed_.insert(target);
            DEBUG("Inserting " << target.int_id());
            result = true;
        }
        if (visited_.find(next) == visited_.end()) {
            DEBUG("Not visited, processing");
            if (ProcessVertex(next, target)) {
                result = true;
            }
        } else {
            DEBUG("Visited");
            if (passed_.find(next) != passed_.end()) {
                passed_.insert(vertex);
                DEBUG("Inserting " << vertex.int_id());
                result = true;
            }
        }
    }
    if (result) {
        passed_.insert(vertex);
    }
    return result;
}
unordered_set<ReachabilityChecker::VertexT> ReachabilityChecker::GetPassedVertices() {
    return passed_;
}
ReachabilityChecker::~ReachabilityChecker() =
default;

ReachabilityChecker::ReachabilityChecker(const ReachabilityChecker::SimpleTransitionGraph &graph_)
    : visited_(), passed_(), graph_(graph_) {}
ReachabilityChecker::SimpleTransitionGraph::const_iterator ForwardReachabilityChecker::GetBeginIterator(
    const ReachabilityChecker::VertexT &vertex) const {
    return graph_.outcoming_begin(vertex);
}
ReachabilityChecker::SimpleTransitionGraph::const_iterator ForwardReachabilityChecker::GetEndIterator(
    const ReachabilityChecker::VertexT &vertex) const {
    return graph_.outcoming_end(vertex);
}

ForwardReachabilityChecker::ForwardReachabilityChecker(const ReachabilityChecker::SimpleTransitionGraph &graph_)
    : ReachabilityChecker(graph_) {}
ReachabilityChecker::SimpleTransitionGraph::const_iterator BackwardReachabilityChecker::GetBeginIterator(
    const ReachabilityChecker::VertexT &vertex) const {
    return graph_.incoming_begin(vertex);
}
ReachabilityChecker::SimpleTransitionGraph::const_iterator BackwardReachabilityChecker::GetEndIterator(
    const ReachabilityChecker::VertexT &vertex) const {
    return graph_.incoming_end(vertex);
}
BackwardReachabilityChecker::BackwardReachabilityChecker(const ReachabilityChecker::SimpleTransitionGraph &graph_)
    : ReachabilityChecker(graph_) {}

CloudSubgraphExtractorParams::CloudSubgraphExtractorParams(size_t distance_threshold_,
                                                           double share_threshold_,
                                                           size_t count_threshold_,
                                                           size_t small_length_threshold_,
                                                           size_t large_length_threshold_,
                                                           size_t min_length_for_barcode_collection) :
    distance_threshold_(distance_threshold_),
    share_threshold_(share_threshold_),
    count_threshold_(count_threshold_),
    small_length_threshold_(small_length_threshold_),
    large_length_threshold_(large_length_threshold_),
    min_length_for_barcode_collection_(min_length_for_barcode_collection) {}

ScaffoldGraph ScaffoldGraphPolisherLauncher::GetFinalScaffoldGraph(const conj_graph_pack &graph_pack,
                                                                   const ScaffoldGraphStorage &scaffold_graph_storage,
                                                                   bool path_scaffolding) {
    const auto &large_scaffold_graph = scaffold_graph_storage.GetLargeScaffoldGraph();
    const auto &small_scaffold_graph = scaffold_graph_storage.GetSmallScaffoldGraph();

    ScaffoldGraphGapCloserParamsConstructor params_constructor;
    auto subgraph_extractor_params =
        params_constructor.ConstructSubgraphExtractorParamsFromConfig(scaffold_graph_storage.GetLargeLengthThreshold(),
                                                                      configs_);
    auto path_extractor_params = params_constructor.ConstructPathExtractorParamsFromConfig(configs_);
    ScaffoldIndexInfoExtractorHelper scaffold_index_helper;
    auto scaffold_index_extractor = scaffold_index_helper.ConstructIndexExtractorFromParams(small_scaffold_graph,
                                                                                            graph_pack,
                                                                                            subgraph_extractor_params,
                                                                                            max_threads_);

    auto initial_cluster_storage = ConstructInitialStorage(graph_pack, small_scaffold_graph,
                                                           path_extractor_params, path_scaffolding);
    INFO("Initial cluster storage size: " << initial_cluster_storage->get_cluster_storage().Size());
    auto barcode_extractor_ptr = make_shared<barcode_index::FrameBarcodeIndexInfoExtractor>(graph_pack.barcode_mapper_ptr,
                                                                                            graph_pack.g);
    const size_t linkage_distance = path_extractor_params.linkage_distance_;
    const double relative_threshold = path_extractor_params.path_cluster_relative_threshold_;

    auto path_extractor = std::make_shared<CloudBasedPathExtractor>(graph_pack.g, initial_cluster_storage,
                                                                    barcode_extractor_ptr, linkage_distance,
                                                                    relative_threshold);

    ScaffoldGraphPolisher gap_closer(graph_pack.g, scaffold_index_extractor,
                                     path_extractor, subgraph_extractor_params);

    auto new_small_scaffold_graph =
        gap_closer.CleanSmallGraphUsingLargeGraph(large_scaffold_graph, small_scaffold_graph);
    return new_small_scaffold_graph;
}

shared_ptr<cluster_storage::InitialClusterStorage> ScaffoldGraphPolisherLauncher::ConstructInitialStorage(
        const conj_graph_pack &gp, const ScaffoldGraph &scaffold_graph,
        const PathExtractionParams &params, bool path_scaffolding) const {
    INFO("Constructing initial cluster storage");
    size_t cluster_storage_builder_threads = max_threads_;
    set<ScaffoldVertex> target_edges;
    std::copy(scaffold_graph.vbegin(), scaffold_graph.vend(), std::inserter(target_edges, target_edges.begin()));
    auto barcode_extractor_ptr = make_shared<barcode_index::FrameBarcodeIndexInfoExtractor>(gp.barcode_mapper_ptr, gp.g);
    size_t linkage_distance = params.linkage_distance_;
    size_t min_read_threshold = params.min_read_threshold_;
    if (not path_scaffolding) {
        auto edge_cluster_extractor =
            make_shared<cluster_storage::AccurateEdgeClusterExtractor>(gp.g, barcode_extractor_ptr,
                                                                       linkage_distance, min_read_threshold);
        auto storage_builder =
            std::make_shared<cluster_storage::EdgeInitialClusterStorageBuilder>(gp.g, edge_cluster_extractor,
                                                                                target_edges, linkage_distance,
                                                                                min_read_threshold,
                                                                                cluster_storage_builder_threads);
        auto result =
            std::make_shared<cluster_storage::InitialClusterStorage>(storage_builder->ConstructInitialClusterStorage());
        return result;
    }

    size_t edge_length_threshold = params.min_length_for_barcode_collection_;
    auto edge_cluster_extractor =
        make_shared<cluster_storage::AccurateEdgeClusterExtractor>(gp.g, barcode_extractor_ptr,
                                                                   linkage_distance, min_read_threshold);
    auto storage_builder =
        std::make_shared<cluster_storage::PathInitialClusterStorageBuilder>(gp.g, edge_cluster_extractor,
                                                                            target_edges, linkage_distance,
                                                                            min_read_threshold,
                                                                            cluster_storage_builder_threads,
                                                                            edge_length_threshold);
    auto result =
        std::make_shared<cluster_storage::InitialClusterStorage>(storage_builder->ConstructInitialClusterStorage());
    return result;
}
ScaffoldGraphPolisherLauncher::ScaffoldGraphPolisherLauncher(size_t max_threads, const ReadCloudConfigs &configs) :
    max_threads_(max_threads), configs_(configs) {}

CloudSubgraphExtractorParams ScaffoldGraphGapCloserParamsConstructor::ConstructSubgraphExtractorParamsFromConfig(
        size_t length_upper_bound,
        const ReadCloudConfigs &configs) {
    const size_t large_length_threshold = length_upper_bound;
    const size_t small_length_threshold = configs.long_edge_length_lower_bound;
    const size_t distance_threshold = configs.scaff_pol.max_scaffold_dijkstra_distance;
    const double share_threshold = configs.scaff_pol.share_threshold;
    const size_t count_threshold = configs.scaff_pol.read_count_threshold;
    const size_t min_length_for_barcode_collection = configs.scaff_con.min_edge_length_for_barcode_collection;
    CloudSubgraphExtractorParams subgraph_extractor_params(distance_threshold, share_threshold, count_threshold,
                                                           small_length_threshold, large_length_threshold,
                                                           min_length_for_barcode_collection);
    return subgraph_extractor_params;
}
PathExtractionParams ScaffoldGraphGapCloserParamsConstructor::ConstructPathExtractorParamsFromConfig(
        const ReadCloudConfigs &configs) {
    const size_t linkage_distance = configs.scaff_pol.path_cluster_linkage_distance;
    const double score_threshold = configs.scaff_pol.path_cluster_relative_threshold;
    const size_t min_read_threshold = configs.scaff_pol.path_cluster_min_reads;
    const size_t min_length_for_barcode_collection = configs.scaff_con.min_edge_length_for_barcode_collection;
    PathExtractionParams predicate_params(linkage_distance, score_threshold, min_read_threshold,
                                          min_length_for_barcode_collection);
    return predicate_params;
}

PathExtractionParams::PathExtractionParams(size_t linkage_distance,
                                           double path_cluster_relative_threshold,
                                           size_t min_read_threshold,
                                           size_t min_length_for_barcode_collection) :
    linkage_distance_(linkage_distance),
    path_cluster_relative_threshold_(path_cluster_relative_threshold),
    min_read_threshold_(min_read_threshold),
    min_length_for_barcode_collection_(min_length_for_barcode_collection) {}
bool GapCloserUtils::IsSimplePath(const GapCloserUtils::SimpleTransitionGraph &graph,
                                  const ScaffoldVertex &source,
                                  const ScaffoldVertex &sink) const {
    auto current_vertex = source;
    bool result = true;
    while (current_vertex != sink) {
        if (graph.GetOutdegree(current_vertex) != 1) {
            return false;
        }

        for (auto next_it = graph.outcoming_begin(current_vertex); next_it != graph.outcoming_end(current_vertex);
             ++next_it) {
            auto next = *next_it;
            current_vertex = next;
        }

    }
    return result;
}

bool CutVerticesExtractor::Check(const ScaffoldVertex &sink,
                                 const ScaffoldVertex &source,
                                 const ScaffoldVertex &candidate) {
    std::queue<ScaffoldVertex> vertex_queue;
    std::unordered_set<ScaffoldVertex> processed;
    vertex_queue.push(sink);
    while (not vertex_queue.empty()) {
        auto current_vertex = vertex_queue.front();
        vertex_queue.pop();
        if (current_vertex != candidate and processed.find(current_vertex) == processed.end()) {
            for (auto it = graph_.outcoming_begin(current_vertex); it != graph_.outcoming_end(current_vertex); ++it) {
                ScaffoldVertex next = *it;
                if (processed.find(next) == processed.end()) {
                    vertex_queue.push(next);
                }
            }
        }
        processed.insert(current_vertex);
    }
    return processed.find(source) == processed.end();
}
vector<CutVerticesExtractor::ScaffoldVertex> CutVerticesExtractor::GetCutVertices(const ScaffoldVertex &source,
                                                                                  const ScaffoldVertex &sink) {
    vector<CutVerticesExtractor::ScaffoldVertex> result;
    TRACE("Source: " << source.int_id());
    TRACE("Sink: " << sink.int_id());
    TRACE("Current graph: ");
    for (const auto &vertex: graph_) {
        for (auto it = graph_.outcoming_begin(vertex); it != graph_.outcoming_end(vertex); ++it) {
            TRACE(vertex.int_id() << " -> " << (*it).int_id());
        }
    }
    for (const auto &vertex: graph_) {
        if (vertex != source and vertex != sink and Check(source, sink, vertex)) {
            result.push_back(vertex);
        }
    }
    TRACE(result.size() << " cut vertices:");
    for (const auto &edge: result) {
        TRACE(edge.int_id())
    }
    return result;
}
CutVerticesExtractor::CutVerticesExtractor(const CutVerticesExtractor::SimpleTransitionGraph &graph_)
    : graph_(graph_) {}
shared_ptr<barcode_index::SimpleScaffoldVertexIndexInfoExtractor> ScaffoldIndexInfoExtractorHelper::ConstructIndexExtractorFromParams(
        const scaffold_graph::ScaffoldGraph &scaffold_graph,
        const conj_graph_pack &gp,
        const CloudSubgraphExtractorParams &subgraph_extractor_params,
        size_t max_threads) const {
    barcode_index::SimpleScaffoldVertexIndexBuilderHelper helper;
    const size_t tail_threshold = subgraph_extractor_params.large_length_threshold_;

    const size_t length_threshold = subgraph_extractor_params.min_length_for_barcode_collection_;
    const size_t count_threshold = subgraph_extractor_params.count_threshold_;

    auto barcode_extractor = make_shared<barcode_index::FrameBarcodeIndexInfoExtractor>(gp.barcode_mapper_ptr, gp.g);
    auto tail_threshold_getter = std::make_shared<barcode_index::ConstTailThresholdGetter>(tail_threshold);
    auto scaffold_vertex_index = helper.ConstructScaffoldVertexIndex(gp.g, *barcode_extractor, tail_threshold_getter,
                                                                     count_threshold, length_threshold, max_threads,
                                                                     scaffold_graph.vertices());
    auto scaffold_index_extractor =
        std::make_shared<barcode_index::SimpleScaffoldVertexIndexInfoExtractor>(scaffold_vertex_index);
    return scaffold_index_extractor;
}
}
}