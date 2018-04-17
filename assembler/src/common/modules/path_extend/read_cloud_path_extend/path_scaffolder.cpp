#include "common/pipeline/config_struct.hpp"
#include "read_cloud_path_extend/validation/scaffold_graph_validation.hpp"
#include "path_scaffolder.hpp"
#include "read_cloud_path_extend/scaffold_graph_construction/scaffold_graph_storage_constructor.hpp"
#include "scaffold_graph_extractor.hpp"

namespace path_extend {

void PathScaffolder::MergePaths(const PathContainer &old_paths) const {
    auto barcode_extractor = make_shared<barcode_index::FrameBarcodeIndexInfoExtractor>(gp_.barcode_mapper_ptr, gp_.g);
    ScaffoldGraphStorageConstructor storage_constructor(small_path_length_threshold_, large_path_length_threshold_, gp_);
    bool scaffolding_mode = true;

    INFO(small_path_length_threshold_);
    AnalyzePaths(old_paths, small_path_length_threshold_);

    size_t num_threads = cfg::get().max_threads;
    auto extractor = make_shared<barcode_index::FrameBarcodeIndexInfoExtractor>(gp_.barcode_mapper_ptr, gp_.g);
    CloudScaffoldGraphConstructor constructor(num_threads, gp_, extractor);
    auto path_scaffold_graph = constructor.ConstructScaffoldGraphFromPathContainer(old_paths, small_path_length_threshold_, scaffolding_mode);
//    ScaffoldGraphStorage storage(constructor.ConstructScaffoldGraphFromPathContainer(paths, large_length_threshold_,
//                                                                                     scaffolding_mode),
//
//    const ScaffoldGraphStorage storage = storage_constructor.ConstructStorageFromPaths(old_paths, scaffolding_mode);
//    ScaffoldGraphPolisher polisher(gp_);
//    bool path_polishing_mode = true;
//    auto path_scaffold_graph = polisher.GetScaffoldGraphFromStorage(storage, path_polishing_mode);

    INFO(path_scaffold_graph.VertexCount() << " vertices and " << path_scaffold_graph.EdgeCount()
                                           << " edges in path scaffold graph");

    //todo move validation somewhere else
    if (cfg::get().ts_res.debug_mode) {
        path_extend::validation::ScaffoldGraphValidator scaffold_graph_validator(gp_.g);
        const string path_to_reference = cfg::get().ts_res.statistics.genome_path;
        INFO("Path to reference: " << path_to_reference);
        INFO("Path exists: " << fs::check_existence(path_to_reference));
        const size_t small_length_threshold = small_path_length_threshold_;
        path_extend::validation::FilteredReferencePathHelper path_helper(gp_);
        auto reference_paths = path_helper.GetFilteredReferencePathsFromLength(path_to_reference, small_length_threshold);

        auto stats = scaffold_graph_validator.GetScaffoldGraphStats(path_scaffold_graph, reference_paths);
        INFO("False positive: " << stats.false_positive_);
        INFO("Single false transition: " << stats.single_false_transition_);
        INFO("False univocal edges: " << stats.false_univocal_edges_);
    }

    ScaffoldGraphExtractor graph_extractor;
    auto univocal_edges = graph_extractor.ExtractUnivocalEdges(path_scaffold_graph);
//    auto univocal_edges = graph_extractor.ExtractMaxScoreEdges(path_scaffold_graph);
    INFO("Found " << univocal_edges.size() << " univocal edges");
    MergeUnivocalEdges(univocal_edges);
}

PathScaffolder::PathScaffolder(const conj_graph_pack &gp_,
                               const ScaffoldingUniqueEdgeStorage &unique_storage_,
                               size_t small_path_length_threshold_, size_t large_path_length_threshold)
    : gp_(gp_), unique_storage_(unique_storage_),
      small_path_length_threshold_(small_path_length_threshold_),
      large_path_length_threshold_(large_path_length_threshold) {}
void PathScaffolder::ExtendPathAlongConnections(const PathScaffolder::ScaffoldVertex& start,
                                                const unordered_map<PathScaffolder::ScaffoldVertex,
                                                                    PathScaffolder::ScaffoldVertex> &merge_connections,
                                                const unordered_map<ScaffoldVertex, size_t> &start_to_distance) const {
    scaffold_graph::PathGetter path_getter;
    auto current = start;
    bool next_found = merge_connections.find(current) != merge_connections.end();
    auto start_path = path_getter.GetPathFromScaffoldVertex(start);
    while (next_found) {
        auto next = merge_connections.at(current);
        auto next_path = path_getter.GetPathFromScaffoldVertex(next);
        if (start_path->GetId() == next_path->GetId()) {
            break;
        }
        DEBUG("First path: " << start_path->GetId() << ", length : " << start_path->Length());
        DEBUG("Second path: " << next_path->GetId() << ", length: " << next_path->Length());
        DEBUG("First conj: " << start_path->GetConjPath()->GetId() << ", length : "
                             << start_path->GetConjPath()->Length());
        DEBUG("Second conj: " << next_path->GetConjPath()->GetId() << ", length: " << next_path->GetConjPath()->Length());
        DEBUG("Got paths")
        Gap path_distance_gap(static_cast<int>(start_to_distance.at(current)));
        DEBUG("Push back")
        start_path->PushBack(*next_path, path_distance_gap);
        DEBUG("Clear");
        next_path->Clear();
        DEBUG("Second path: " << next_path->GetId() << ", length: " << next_path->Length());
        DEBUG(next_path->Empty());
        DEBUG("Conjugate: " << next_path->GetConjPath()->GetId() << ", length: " << next_path->GetConjPath()->Length());
        DEBUG("Conjugate empty: " << next_path->GetConjPath()->Empty());
        current = next;
        next_found = merge_connections.find(current) != merge_connections.end();
    }
}

void PathScaffolder::MergeUnivocalEdges(const vector<PathScaffolder::ScaffoldEdge> &scaffold_edges) const {
    std::unordered_map<ScaffoldVertex, ScaffoldVertex> merge_connections;
    for (const auto &edge: scaffold_edges) {
        ScaffoldVertex start = edge.getStart();
        ScaffoldVertex end = edge.getEnd();
        DEBUG(start.int_id() << " -> " << end.int_id());
        DEBUG("Weight: " << edge.getWeight());
        VERIFY(merge_connections.find(start) == merge_connections.end());
        merge_connections.insert({start, end});
    }

    for (const auto &connection: merge_connections) {
        auto start = connection.first;
        auto end = connection.second;
        auto start_conjugate = start.getConjugateFromGraph(gp_.g);
        auto end_conjugate = end.getConjugateFromGraph(gp_.g);
        if (merge_connections.find(end_conjugate) == merge_connections.end() or
            merge_connections.at(end_conjugate) != start_conjugate) {
            WARN("Conjugate connection does not correspond to direct connection")
            merge_connections.at(end_conjugate) = start_conjugate;
        } else {
            merge_connections.insert({end_conjugate, start_conjugate});
        }
    }

    std::unordered_set<ScaffoldVertex> starts;
    std::unordered_set<ScaffoldVertex> used;
    for (const auto &connection: merge_connections) {
        auto current = connection.first;
        auto current_conjugate = current.getConjugateFromGraph(gp_.g);
        if (used.find(current) != used.end()) {
            continue;
        }
        bool prev_found = merge_connections.find(current_conjugate) != merge_connections.end();
        used.insert(current);
        used.insert(current_conjugate);
        bool prev_used = false;
        while (prev_found) {
            auto prev_conjugate = merge_connections.at(current_conjugate);
            if (used.find(prev_conjugate) != used.end()) {
                prev_used = true;
                break;
            }
            current = prev_conjugate.getConjugateFromGraph(gp_.g);
            current_conjugate = current.getConjugateFromGraph(gp_.g);
            prev_found = merge_connections.find(current_conjugate) != merge_connections.end();
            used.insert(current);
            used.insert(current_conjugate);
        }
        if (not prev_used) {
            starts.insert(current);
        }
    }
    std::unordered_map<ScaffoldVertex, size_t> start_to_distance;
    for (const auto& edge: scaffold_edges) {
        start_to_distance.insert({edge.getStart(), edge.getLength()});
        start_to_distance.insert({edge.getEnd().getConjugateFromGraph(gp_.g), edge.getLength()});
    }
    for (const auto& connection: merge_connections) {
        DEBUG(connection.first.int_id() << " -> " << connection.second.int_id());
    }
    scaffold_graph::PathGetter path_getter;
    INFO(starts.size() << " starts.");
    for (const auto &start: starts) {
        ScaffoldVertex current = start;
        bool next_found = merge_connections.find(current) != merge_connections.end();
        DEBUG("Start: " << current.int_id());
        while(next_found and merge_connections.at(current) != start) {
            current = merge_connections.at(current);
            next_found = merge_connections.find(current) != merge_connections.end();
            DEBUG(current.int_id());
        }
    }
    for (const auto &start: starts) {
        if (not path_getter.GetPathFromScaffoldVertex(start)->Empty()) {
            ExtendPathAlongConnections(start, merge_connections, start_to_distance);
        }
    }
}
void PathScaffolder::AnalyzePaths(const PathContainer &paths, size_t min_length) const {
    INFO("Analyzing long paths");
    auto dataset_info = cfg::get().ds;
    for (size_t lib_index = 0; lib_index < dataset_info.reads.lib_count(); ++lib_index) {
        const auto &lib = dataset_info.reads[lib_index];
        if (lib.type() != io::LibraryType::Clouds10x) {
            continue;
        }
        size_t no_candidates = 0;
        size_t single_candidate = 0;
        size_t simple_step_made = 0;
        size_t whole_step_made = 0;
        size_t rc_candidate = 0;
        size_t no_tip = 0;
        size_t total = 0;
        size_t multiple_candidates = 0;
        double total_no_candidate_coverage = 0;
        double total_single_coverage = 0;
        path_extend::PathExtendParamsContainer params(cfg::get().ds,
                                                      cfg::get().pe_params,
                                                      cfg::get().ss,
                                                      cfg::get().output_dir,
                                                      cfg::get().mode,
                                                      cfg::get().uneven_depth,
                                                      cfg::get().avoid_rc_connections,
                                                      cfg::get().use_scaffolder);
        const auto &pset = params.pset;
        shared_ptr<PairedInfoLibrary> paired_lib = MakeNewLib(gp_.g, lib, gp_.scaffolding_indices[lib_index]);

        shared_ptr<WeightCounter> counter = make_shared<ReadCountWeightCounter>(gp_.g, paired_lib);

        auto scaff_chooser = std::make_shared<ScaffoldingExtensionChooser>(gp_.g, counter,
                                                                           pset.scaffolder_options.cl_threshold,
                                                                           pset.scaffolder_options.var_coeff);
        GraphCoverageMap cover_map(gp_.g);
        ScaffoldingUniqueEdgeStorage unique_storage;
        UsedUniqueStorage used_unique_storage(unique_storage);
        auto scaff_extender = make_shared<ScaffoldingPathExtender>(gp_, cover_map,
                                                                   used_unique_storage, scaff_chooser,
                                                                   MakeGapAnalyzer(paired_lib->GetIsVar()),
                                                                   paired_lib->GetISMax(),
                                                                   false, /* investigate short loops */
                                                                   params.avoid_rc_connections);
        ExtensionChooser::EdgeContainer sources;
        PathContainer* empty;
        for (auto iter = gp_.g.ConstEdgeBegin(); !iter.IsEnd(); ++iter) {
            if (gp_.g.IncomingEdgeCount(gp_.g.EdgeStart(*iter)) == 0) {
                sources.push_back(EdgeWithDistance(*iter, 0));
            }
        }
        for (const auto& path: paths) {
            if (path.first->Length() <= min_length) {
                continue;
            }
            ++total;
            EdgeId back = path.first->Back();
            if (gp_.g.OutgoingEdgeCount(gp_.g.EdgeEnd(back)) != 0) {
                ++no_tip;
                continue;
            }
            auto candidates = scaff_chooser->Filter(*(path.first), sources);
            if (candidates.size() == 0) {
                ++no_candidates;
                total_no_candidate_coverage += gp_.g.coverage(back);

            } else if (candidates.size() == 1) {
                ++single_candidate;

                if (scaff_extender->MakeSimpleGrowStep(*(path.first), empty)) {
                    ++simple_step_made;
                }
                if (scaff_extender->MakeGrowStep((*(path.first)), empty)) {
                    ++whole_step_made;
                }
                total_single_coverage += gp_.g.coverage(back);
                if (gp_.g.conjugate(back) == candidates[0].e_ or back == candidates[0].e_) {
                    ++rc_candidate;
                }

            } else {
                ++multiple_candidates;
            }
        }
        INFO("Total paths: " << total);
        INFO("No tips: " << no_tip);
        INFO("Single candidate: " << single_candidate);
        if (single_candidate > 0) {
            double mean_single_candidate_coverage = total_single_coverage / static_cast<double>(single_candidate);
            INFO("Mean coverage of single candidate edge: " << mean_single_candidate_coverage);
        }
        INFO("Simple step made: " << simple_step_made);
        INFO("Whole step made: " << whole_step_made);
        INFO("RC candidate: " << rc_candidate);
        INFO("No candidates: " << no_candidates);
        if (no_candidates > 0) {
            double mean_no_candidate_coverage = total_no_candidate_coverage / static_cast<double>(no_candidates);
            INFO("Mean coverage of edge with no candidates: " << mean_no_candidate_coverage);
        }
        INFO("Multiple candidates: " << multiple_candidates);
    }

}
shared_ptr<GapAnalyzer> PathScaffolder::MakeGapAnalyzer(double is_variation) const {
    path_extend::PathExtendParamsContainer params(cfg::get().ds,
                                                  cfg::get().pe_params,
                                                  cfg::get().ss,
                                                  cfg::get().output_dir,
                                                  cfg::get().mode,
                                                  cfg::get().uneven_depth,
                                                  cfg::get().avoid_rc_connections,
                                                  cfg::get().use_scaffolder);
    const auto &pset = params.pset;

    vector<shared_ptr<GapAnalyzer>> joiners;
    if (pset.scaffolder_options.use_la_gap_joiner)
        joiners.push_back(std::make_shared<LAGapAnalyzer>(gp_.g, pset.scaffolder_options.min_overlap_length,
                                                          pset.scaffolder_options.flank_multiplication_coefficient,
                                                          pset.scaffolder_options.flank_addition_coefficient));


    joiners.push_back(std::make_shared<HammingGapAnalyzer>(gp_.g,
                                                           pset.scaffolder_options.min_gap_score,
                                                           pset.scaffolder_options.short_overlap,
                                                           (int) pset.scaffolder_options.basic_overlap_coeff
                                                               * cfg::get().ds.RL));

    //todo introduce explicit must_overlap_coeff and rename max_can_overlap -> can_overlap_coeff
    return std::make_shared<CompositeGapAnalyzer>(gp_.g,
                                                  joiners,
                                                  size_t(math::round(pset.scaffolder_options.max_can_overlap
                                                                         * is_variation)), /* may overlap threshold */
                                                  int(math::round(-pset.scaffolder_options.var_coeff * is_variation)), /* must overlap threshold */
                                                  pset.scaffolder_options.artificial_gap);
}
}