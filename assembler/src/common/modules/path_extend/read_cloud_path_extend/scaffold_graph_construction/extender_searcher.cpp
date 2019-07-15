

#include "extender_searcher.hpp"
namespace path_extend {
namespace read_cloud {

bool ExtenderSearcher::IsReachable(VertexId target_vertex, path_extend::BidirectionalPath *start) const {
    auto paths_container = std::make_shared<path_extend::QueueContainer>();
    GraphCoverageMap cover_map(gp_.g);
    DEBUG("Initial path id: " << start->GetId());
    DEBUG("Initial path size: " << start->Size());
    auto start_copy = new BidirectionalPath(*start);
    paths_container->push(start_copy);
    DEBUG("Constructed path container");

    VertexId start_vertex = gp_.g.EdgeEnd(start->Back());
    DEBUG("Start: " << start_vertex.int_id());
    DEBUG("Target: " << target_vertex.int_id());
    if (start_vertex == target_vertex) {
        return true;
    }

    bool result = SearchTargetUsingExtenders(paths_container, cover_map, target_vertex);
    return result;
}

bool ExtenderSearcher::SearchTargetUsingExtenders(std::shared_ptr<path_extend::QueueContainer> paths_container,
                                                  GraphCoverageMap &cover_map,
                                                  const VertexId &target_vertex) const {
    size_t path_processing_iterations = 0;
    std::unordered_map<EdgeId, size_t> edge_to_number_of_visits;

    while (not paths_container->empty() and path_processing_iterations < search_params_.max_paths_to_process) {
        DEBUG("Path container size: " << paths_container->size());
        DEBUG("Coverage map size: " << cover_map.size());
        auto current_path = paths_container->front();
        SubscribeCoverageMap(current_path, cover_map);
        paths_container->pop();
        EdgeId last_edge = current_path->Back();
        edge_to_number_of_visits[last_edge]++;
        size_t current_iterations = 0;
        if (edge_to_number_of_visits[last_edge] > search_params_.max_edge_visits) {
            continue;
        }
        PathContainer *empty_storage = nullptr;
        auto multi_extender = ConstructSearchingExtender(paths_container, cover_map);
        while (multi_extender->MakeGrowStep(*current_path, empty_storage)
            and current_iterations < search_params_.max_path_growing_iterations) {
            DEBUG("Made grow step, iteration " << current_iterations);
            ++current_iterations;
        }
        ++path_processing_iterations;
        if (current_iterations >= search_params_.max_path_growing_iterations) {
            WARN("Stopped path growing because of too many iterations");
            continue;
        }
        const auto &reached_vertices = multi_extender->GetReachedVertices();
        if (reached_vertices.find(target_vertex) != reached_vertices.end()) {
            DEBUG("Found target");
            return true;
        }
    }
    if (path_processing_iterations >= search_params_.max_paths_to_process) {
        DEBUG("Had to process too many paths, returning");
        return true;
    }
    return false;
}
ExtenderSearcher::ExtenderSearcher(const conj_graph_pack &gp,
                                   std::shared_ptr<ExtensionChooser> extension_chooser,
                                   const SearchParams &search_params,
                                   const SearchingExtenderParams &extender_params, size_t length_bound) :
    gp_(gp), extension_chooser_(extension_chooser), search_params_(search_params),
    extender_params_(extender_params), length_bound_(length_bound) {}
std::shared_ptr<SearchingMultiExtender> ExtenderSearcher::ConstructSearchingExtender(
        std::shared_ptr<path_extend::QueueContainer> paths_container,
        GraphCoverageMap &cover_map) const {
    UsedUniqueStorage used_unique_storage(extender_params_.unique_storage, gp_.g);
    size_t insert_size = extender_params_.insert_size;
    bool investigate_loops = extender_params_.investigate_short_loops;
    bool short_loop_resolver = extender_params_.use_short_loop_cov_resolver;
    double weight_threshold = extender_params_.weight_threshold;

    auto multi_extender = make_shared<path_extend::SearchingMultiExtender>(gp_, cover_map, used_unique_storage,
                                                                           extension_chooser_, insert_size,
                                                                           investigate_loops, short_loop_resolver,
                                                                           weight_threshold, length_bound_,
                                                                           paths_container);
    return multi_extender;
}

SearchParams::SearchParams(size_t max_path_growing_iterations, size_t max_paths_to_process, size_t max_edge_visits) :
    max_path_growing_iterations(max_path_growing_iterations),
    max_paths_to_process(max_paths_to_process),
    max_edge_visits(max_edge_visits) {}
SearchParams::SearchParams() :
    max_path_growing_iterations(0),
    max_paths_to_process(0),
    max_edge_visits(0) {}
SearchingExtenderParams DefaultExtenderParamsConstructor::ConstructExtenderParams(size_t lib_index,
                                                                                  double weight_threshold) const {
    const auto &lib = dataset_info_.reads[lib_index];
    shared_ptr<PairedInfoLibrary> paired_lib = MakeNewLib(gp_.g, lib, gp_.clustered_indices[lib_index]);
    size_t insert_size = paired_lib->GetISMax();
    UsedUniqueStorage used_unique_storage(unique_storage_, gp_.g);
    const bool investigate_loops = true;
    const bool short_loop_resolver = false;
    const size_t length_bound = insert_size;

    SearchingExtenderParams params {unique_storage_, insert_size, investigate_loops,
                                    short_loop_resolver, weight_threshold, length_bound};
    return params;
}
DefaultExtenderParamsConstructor::DefaultExtenderParamsConstructor(const conj_graph_pack &gp,
                                                                   const config::dataset &dataset_info,
                                                                   const ScaffoldingUniqueEdgeStorage &unique_storage)
    : gp_(gp), dataset_info_(dataset_info), unique_storage_(unique_storage) {}
SearchingExtenderParams::SearchingExtenderParams(const ScaffoldingUniqueEdgeStorage &unique_storage) :
    unique_storage(unique_storage) {}
SearchingExtenderParams::SearchingExtenderParams(const ScaffoldingUniqueEdgeStorage &unique_storage,
                                                 size_t insert_size,
                                                 bool investigate_short_loops,
                                                 bool use_short_loop_cov_resolver,
                                                 double weight_threshold,
                                                 size_t length_bound) :
    unique_storage(unique_storage),
    insert_size(insert_size),
    investigate_short_loops(investigate_short_loops),
    use_short_loop_cov_resolver(use_short_loop_cov_resolver),
    weight_threshold(weight_threshold),
    length_bound(length_bound) {}
}
}
