//#pragma once
//
//#include "statistics_processor.hpp"
//
//namespace transitions {
//
//struct ReferenceEdgeInfo {
//  size_t path_id_;
//  size_t edge_index_;
//  Range mapping_;
//};
//
//typedef std::unordered_map<EdgeId, ReferenceEdgeInfo> ReferenceEdgeIndex;
//
//struct LongEdgePathInfo {
//  size_t length_;
//  bool correct_;
//  vector<path_extend::validation::EdgeWithMapping> path_;
//};
//
//struct PathStats : public read_cloud_statistics::Statistic {
//    const size_t min_length_;
//    size_t incorrect_paths_;
//    size_t nonreference_paths_;
//    size_t number_of_paths_;
//    std::unordered_map<size_t, LongEdgePathInfo> path_to_info_;
//    ReferenceEdgeIndex reference_index_;
//
// public:
//    PathStats(size_t min_length) : Statistic("contig_paths_statistics"), min_length_(min_length), incorrect_paths_(0),
//                                   number_of_paths_(0), path_to_info_(), reference_index_() {}
//
//    void Serialize(const string& path) override {
//        ofstream fout(path);
//        const string sep = "\t";
//        fout << "Min length: " << min_length_ << std::endl;
//        fout << "Incorrect paths: " << incorrect_paths_ << std::endl;
//        fout << "Nonreference paths: " << nonreference_paths_ << std::endl;
//        fout << "Number of paths: " << number_of_paths_ << std::endl;
//        fout << "Length" << sep << "Correct" << std::endl;
//        for (const auto& entry: path_to_info_) {
//            fout << "Id: " << entry.first << std::endl;
//            fout << "Length: " << entry.second.length_ << std::endl;
//            string correct_output = entry.second.correct_ ? "True" : "False";
//            fout << "Correct: " << correct_output << std::endl;
//            for (const auto& ewm: entry.second.path_) {
//                EdgeId edge = ewm.edge_;
//                if (reference_index_.find(ewm.edge_) != reference_index_.end()) {
//                    fout << "Edge id: " << edge.int_id() << sep << "Path id: " << reference_index_.at(edge).path_id_
//                         << sep << "Edge index: " << reference_index_.at(edge).edge_index_ << std::endl;
//                }
//            }
//            fout << std::endl;
//            fout << std::endl;
//        }
//    }
//};
//
//
//class PathStatisticsExtractor : public read_cloud_statistics::StatisticProcessor {
// public:
//    typedef path_extend::validation::EdgeWithMapping EdgeWithMapping;
// private:
//
//    const size_t min_length_;
//    vector<vector<EdgeWithMapping>> reference_paths_;
//    vector<vector<EdgeWithMapping>> contig_paths_;
//    const Graph& g_;
//
// public:
//    PathStatisticsExtractor(size_t min_length,
//                            const vector<vector<EdgeWithMapping>>& reference_paths,
//                            const vector<vector<EdgeWithMapping>>& contig_paths, const Graph& g) :
//        StatisticProcessor("path_statistics_extractor"),
//        min_length_(min_length),
//        reference_paths_(reference_paths),
//        contig_paths_(contig_paths),
//        g_(g) {}
//
//    void FillStatistics() override {
//        auto path_stats_ptr = std::make_shared<PathStats>(GetPathStats(reference_paths_, contig_paths_));
//        AddStatistic(path_stats_ptr);
//    }
//
// private:
//
//    PathStats GetPathStats(const vector<vector<EdgeWithMapping>>& reference_paths,
//                           const vector<vector<EdgeWithMapping>>& contig_paths) const {
//        size_t path_counter = 0;
//        for (const auto& path: reference_paths) {
//            TRACE("Path " << path_counter);
//            for (const auto& ewm: path) {
//                TRACE(ewm.edge_.int_id() << " " << ewm.mapping_);
//            }
//        }
//
//        ReferenceEdgeIndex reference_edge_index = BuildReferenceEdgeIndex(reference_paths);
//        return BuildPathStats(contig_paths, reference_edge_index);
//    }
//
//    vector<vector<EdgeWithMapping>> FilterPaths(const vector<vector<EdgeWithMapping>>& contig_paths) const {
//        vector<vector<EdgeWithMapping>> result;
//        for (const auto& path: contig_paths) {
//            vector<EdgeWithMapping> current_path;
//            std::copy_if(path.begin(), path.end(), std::back_inserter(current_path), [this](const EdgeWithMapping& ewm) {
//              return g_.length(ewm.edge_) >= min_length_;
//            });
//            if (current_path.size() > 0) {
//                result.push_back(current_path);
//            }
//        }
//        return result;
//    }
//
//    ReferenceEdgeIndex BuildReferenceEdgeIndex(const vector<vector<EdgeWithMapping>>& reference_paths) const {
//        size_t path_id = 0;
//        ReferenceEdgeIndex result;
//        for (const auto& path: reference_paths) {
//            size_t edge_id = 0;
//            for (const auto& ewm: path) {
//                ReferenceEdgeInfo edge_info;
//                VERIFY(g_.length(ewm.edge_) >= min_length_);
//                edge_info.path_id_ = path_id;
//                edge_info.edge_index_ = edge_id;
//                edge_info.mapping_ = ewm.mapping_;
//                result.insert({ewm.edge_, edge_info});
//                ++edge_id;
//            }
//            ++path_id;
//        }
//        return result;
//    }
//
//    bool IsReferencePath(const vector<EdgeWithMapping>& contig_path, const ReferenceEdgeIndex& reference_index) const {
//        VERIFY(contig_path.size() > 0);
//        for (const auto& ewm: contig_path) {
//            if (reference_index.find(ewm.edge_) == reference_index.end()) {
//                return false;
//            }
//        }
//        return true;
//    }
//
//    bool IsCorrectPath(const vector<EdgeWithMapping>& contig_path, const ReferenceEdgeIndex& reference_index) const {
//        VERIFY(contig_path.size() > 0);
//        auto ewm = contig_path[0];
//        if (reference_index.find(ewm.edge_) == reference_index.end()) {
//            return true;
//        }
//        size_t current_path_id = reference_index.at(ewm.edge_).path_id_;
//        size_t current_edge_index = reference_index.at(ewm.edge_).edge_index_;
//        for (auto it = std::next(contig_path.begin()); it != contig_path.end(); ++it) {
//            EdgeWithMapping current_ewm = *it;
//            if (reference_index.find(current_ewm.edge_) == reference_index.end()) {
//                return true;
//            }
//            size_t path_id = reference_index.at(current_ewm.edge_).path_id_;
//            size_t edge_index = reference_index.at(current_ewm.edge_).edge_index_;
//            if (path_id != current_path_id) {
//                return false;
//            }
//            if (edge_index != current_edge_index + 1) {
//                return false;
//            }
//            current_edge_index = edge_index;
//        }
//        return true;
//    }
//
//    PathStats BuildPathStats(const vector<vector<EdgeWithMapping>>& contig_paths,
//                             const ReferenceEdgeIndex& reference_index) const {
//        PathStats result(min_length_);
//        size_t path_index = 0;
//        size_t incorrect_paths = 0;
//        size_t nonreference_paths = 0;
//        size_t number_of_paths = contig_paths.size();
//        unordered_map<size_t, LongEdgePathInfo> long_edge_path_info;
//        for (const auto& path: contig_paths) {
//            LongEdgePathInfo path_info;
//            bool is_correct = IsCorrectPath(path, reference_index);
//            bool is_reference_path = IsReferencePath(path, reference_index);
//            path_info.correct_ = is_correct;
//            if (not is_correct) {
//                ++incorrect_paths;
//            }
//            if (not is_reference_path) {
//                ++nonreference_paths;
//            }
//            path_info.length_ = path.size();
//            path_info.path_ = path;
//            long_edge_path_info.insert({path_index, path_info});
//            ++path_index;
//        }
//        result.number_of_paths_ = number_of_paths;
//        result.incorrect_paths_ = incorrect_paths;
//        result.nonreference_paths_ = nonreference_paths;
//        result.path_to_info_ = long_edge_path_info;
//        result.reference_index_ = reference_index;
//        return result;
//    }
//};
//
//struct CandidateStats {
//  ScaffoldVertex edge_;
//  size_t shared_barcodes_;
//  double cov_;
//  size_t barcodes_;
//  double score_;
//  size_t distance_;
//};
//
//struct PathEndInfo {
//  size_t candidates_;
//  double max_score_;
//  double second_score_;
//  EdgeId next_edge_;
//  EdgeId rc_edge_;
//  size_t barcodes_;
//  double coverage_;
//  vector<CandidateStats> candidate_stats_;
//};
//
//struct InitialFilterStats : public read_cloud_statistics::Statistic {
//  std::unordered_map<EdgeId, PathEndInfo> path_end_to_info_;
//
//  InitialFilterStats()
//      : Statistic("initial_filter_stats"), path_end_to_info_() {}
//
//  void Serialize(const string& path) override {
//      ofstream fout(path);
//      string sep = "\t";
//      fout << path_end_to_info_.size() << std::endl;
//      for (const auto& entry: path_end_to_info_) {
//          const auto& record = entry.second;
//          fout << entry.first.int_id() << sep << record.candidates_ << sep << record.max_score_
//               << sep << record.coverage_ << sep << record.barcodes_
//               << sep << record.next_edge_.int_id() << sep << record.rc_edge_.int_id() << std::endl;
//          for (const auto& stat: record.candidate_stats_) {
//              fout << stat.edge_.int_id() << sep << stat.shared_barcodes_ << sep << stat.cov_
//                   << sep << stat.barcodes_ << sep << stat.score_ << sep << stat.distance_ << std::endl;
//          }
//          fout << std::endl;
//      }
//  }
//};
//
//class InitialFilterStatisticsExtractor: public read_cloud_statistics::StatisticProcessor {
// public:
//    typedef std::unordered_map<EdgeId, EdgeId> NextInReferenceIndex;
//    typedef path_extend::scaffold_graph::ScaffoldGraph ScaffoldGraph;
//    typedef ScaffoldGraph::ScaffoldGraphVertex ScaffoldVertex;
//    typedef barcode_index::FrameBarcodeIndexInfoExtractor barcode_extractor_t;
//    typedef path_extend::validation::EdgeWithMapping EdgeWithMapping;
// private:
//    const ScaffoldGraph& scaffold_graph_;
//    const vector<vector<EdgeWithMapping>> reference_paths_;
//    shared_ptr<barcode_extractor_t> barcode_extractor_ptr_;
//    const path_extend::InitialTenXFilter initial_filter_;
//    const Graph& g_;
//    const path_extend::ScaffoldingUniqueEdgeStorage& unique_storage_;
//
// public:
//    InitialFilterStatisticsExtractor(const ScaffoldGraph& scaffold_graph_,
//                                     const vector<vector<EdgeWithMapping>>& reference_paths_,
//                                     const shared_ptr<barcode_extractor_t> barcode_extractor_ptr,
//                                     const path_extend::InitialTenXFilter initial_filter,
//                                     const Graph& g, const path_extend::ScaffoldingUniqueEdgeStorage& unique_storage)
//        : StatisticProcessor("initial_filter_stats"),
//          scaffold_graph_(scaffold_graph_), reference_paths_(reference_paths_),
//          barcode_extractor_ptr_(barcode_extractor_ptr), initial_filter_(initial_filter),
//          g_(g), unique_storage_(unique_storage)
//    {}
//
//    void FillStatistics() override {
//        path_extend::validation::ContigPathFilter contig_path_filter(unique_storage_);
//        auto filtered_paths = contig_path_filter.FilterPathsUsingUniqueStorage(reference_paths_);
//        auto initial_filter_stats = make_shared<InitialFilterStats>(GetInitialFilterStats(filtered_paths));
//        AddStatistic(initial_filter_stats);
//    }
//
// private:
//    InitialFilterStats GetInitialFilterStats(const vector<vector<EdgeWithMapping>>& reference_paths) const {
//        InitialFilterStats stats;
//        std::unordered_map<EdgeId, PathEndInfo> path_end_to_info;
//
//        path_extend::validation::ContigPathFilter contig_path_filter(unique_storage_);
//        auto filtered_reference_paths = contig_path_filter.FilterPathsUsingUniqueStorage(reference_paths);
//        auto next_in_reference_index = BuildNextInReferenceIndex(reference_paths);
//        size_t counter = 0;
//
//        vector<EdgeId> unique_edges;
//        std::copy(unique_storage_.begin(), unique_storage_.end(), std::back_inserter(unique_edges));
//        VERIFY(unique_edges.size() == unique_storage_.size());
//        const size_t threads = cfg::get().max_threads;
//#pragma omp parallel for num_threads(threads)
//        for (size_t i = 0; i < unique_edges.size(); ++i)
//        {
//            EdgeId edge = unique_edges[i];
//            auto out_edges = scaffold_graph_.OutgoingEdges(edge);
//            vector<CandidateStats> candidate_stats;
//            const size_t coverage_threshold = cfg::get().ts_res.tenx.initial_coverage_threshold;
//            const size_t tail_threshold = cfg::get().ts_res.tenx.tail_threshold;
//            for (const auto& scaffold_edge: out_edges) {
//                auto candidate = scaffold_edge.getEnd();
//                candidate_stats.push_back(GetCandidateStats(edge, candidate, scaffold_edge, coverage_threshold, tail_threshold));
//            }
//            vector<double> scores;
//            for (const auto& stats: candidate_stats) {
//                scores.push_back(stats.score_);
//            }
//            double max_score = std::accumulate(scores.begin(), scores.end(), 0, [](const size_t first, const size_t second) {
//              return std::max(first, second);
//            });
//            double second_score = 0;
//            if (scores.size() > 1) {
//                std::nth_element(scores.begin(), scores.begin() + 1, scores.end(), std::greater<double>());
//                second_score = scores[1];
//            }
//
//            PathEndInfo edge_info;
//            edge_info.max_score_ = max_score;
//            edge_info.second_score_ = second_score;
//            edge_info.candidates_ = out_edges.size();
//            edge_info.candidate_stats_ = candidate_stats;
//            vector<EdgeId> path_suffix({edge});
//            edge_info.barcodes_ = initial_filter_.ExtractBarcodesFromMultipleEdges(path_suffix, coverage_threshold, tail_threshold).size();
//            edge_info.coverage_ = g_.coverage(edge);
//            EdgeId next(nullptr);
//            if (next_in_reference_index.find(edge) != next_in_reference_index.end()) {
//                next = next_in_reference_index.at(edge);
//            }
//            edge_info.next_edge_ = next;
//            edge_info.rc_edge_ = g_.conjugate(edge);
//#pragma omp critical
//            {
//                path_end_to_info.insert({edge, edge_info});
//                ++counter;
//                if (counter % (unique_storage_.size() / 10) == 0) {
//                    INFO("Processed " << counter << " edges out of " << unique_storage_.size());
//                }
//            }
//        }
//        stats.path_end_to_info_ = path_end_to_info;
//        return stats;
//    }
//
//    NextInReferenceIndex BuildNextInReferenceIndex(const vector<vector<EdgeWithMapping>>& reference_paths) const {
//        NextInReferenceIndex next_in_reference_index;
//        std::unordered_set<EdgeId> reference_edges;
//        std::unordered_set<EdgeId> repeats;
//        for (const auto& path: reference_paths) {
//            if (path.size() != 0) {
//                for (auto first = path.begin(), second = std::next(path.begin()); second != path.end();
//                     ++first, ++second) {
//                    EdgeId e = (*first).edge_;
//                    EdgeId next = (*second).edge_;
//                    next_in_reference_index.insert({e, next});
//                }
//                EdgeId first = path[0].edge_;
//                EdgeId last = path.back().edge_;
//                next_in_reference_index.insert({last, first});
//            }
//        }
//        INFO("Next in reference index size: " << next_in_reference_index.size());
//        return next_in_reference_index;
//    }
//
//    CandidateStats GetCandidateStats(const ScaffoldGraph::ScaffoldGraphVertex& path_end,
//                                     const ScaffoldGraph::ScaffoldGraphVertex& candidate,
//                                     const ScaffoldGraph::ScaffoldEdge& scaffold_edge,
//                                     size_t coverage_threshold, size_t tail_threshold) const {
//        vector<ScaffoldVertex> path_edges({path_end});
//        auto barcodes_from_path = initial_filter_.ExtractBarcodesFromMultipleEdges(path_edges, coverage_threshold, tail_threshold);
//        auto barcodes_from_candidate = initial_filter_.ExtractBarcodesFromCandidate(candidate, coverage_threshold, tail_threshold);
//        double coverage_second = g_.coverage(candidate);
//        size_t distance = scaffold_edge.getLength();
//        vector<barcode_index::BarcodeId> shared_barcodes;
//        std::set_intersection(barcodes_from_candidate.begin(), barcodes_from_candidate.end(),
//                              barcodes_from_path.begin(), barcodes_from_path.end(),
//                              std::back_inserter(shared_barcodes));
//        double score = static_cast<double> (shared_barcodes.size()) / coverage_second;
//        CandidateStats stats;
//        stats.score_ = score;
//        stats.shared_barcodes_ = shared_barcodes.size();
//        stats.edge_ = candidate;
//        stats.barcodes_ = barcodes_from_candidate.size();
//        stats.cov_ = coverage_second;
//        stats.distance_ = distance;
//        return stats;
//    }
//
//    DECL_LOGGER("InitialFilterStatisticsExtractor");
//};
//}