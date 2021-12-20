#include "long_edge_statistics.hpp"

#include "utils/filesystem/path_helper.hpp"

namespace cont_index {
void LongEdgeStatisticsCounter::CountDoubleCoverageDistribution() const {
    using debruijn_graph::EdgeId;
    using debruijn_graph::Graph;
    using BarcodeExtractor = barcode_index::FrameBarcodeIndexInfoExtractor;
    auto barcode_extractor_ptr = std::make_shared<BarcodeExtractor>(barcode_index_, graph_);
    std::map<size_t, std::map<size_t, size_t>> double_cov_distribution;
    size_t lower_dc_bound = 10000;
    size_t upper_dc_bound = 50000;
    size_t dc_step = 10000;
    std::unordered_set<size_t> double_cov_gaps;
    for (size_t thr = lower_dc_bound; thr <= upper_dc_bound; thr+=dc_step) {
        double_cov_gaps.insert(thr);
    }
    INFO("Traversing");
    for (const EdgeId &edge: graph_.canonical_edges()) {
        if (graph_.length(edge) >= training_edge_length_) {
            typedef std::unordered_map<size_t, size_t> PerPosDistribution;
            std::unordered_map<size_t, PerPosDistribution> single_edge_double_cov;
            size_t frame_size = barcode_extractor_ptr->GetBinLength(edge);
            size_t num_of_frames = barcode_extractor_ptr->GetNumberOfBins(edge);
            for (const auto &gap: double_cov_gaps) {
                for (size_t i = 0; i < num_of_frames * frame_size - gap + frame_size; i += frame_size) {
                    single_edge_double_cov[gap][i] = 0;
                }
            }
            auto bar_begin = barcode_extractor_ptr->barcode_iterator_begin(edge);
            auto bar_end = barcode_extractor_ptr->barcode_iterator_end(edge);
            for (auto it = bar_begin; it != bar_end; ++it) {
                const auto &bit_positions = it->second.GetBitSet();
                std::vector<size_t> positions;
                size_t current_set_pos = bit_positions.find_first();
                while (current_set_pos != barcode_index::FrameBarcodeInfo::IsOnFrameT::npos) {
                    positions.push_back(current_set_pos * frame_size);
                    current_set_pos = bit_positions.find_next(current_set_pos);
                }
                if (positions.size() == 1) {
                    continue;
                }
                for (auto it1 = positions.begin(); it1 != positions.end(); ++it1) {
                    for (auto it2 = std::next(it1); it2 != positions.end(); ++it2) {
                        const size_t first_pos = *it1;
                        const size_t second_pos = *it2;
                        size_t diff = second_pos - first_pos;
                        if (diff > upper_dc_bound) {
                            break;
                        }
                        if (double_cov_gaps.count(diff) and first_pos < frame_size * num_of_frames - diff) {
                            single_edge_double_cov.at(diff).at(first_pos)++;
                        }
                    }
                }
            }
            for (const auto &entry: single_edge_double_cov) {
                const auto &diff = entry.first;
                const auto &pos_fragments = entry.second;
                for (const auto &entry: pos_fragments) {
                    size_t num_fragments = entry.second;
                    ++double_cov_distribution[diff][num_fragments];
                }
            }
//            INFO("Processed edge " << edge.int_id() << " , length " << graph.length(edge));
        }
    }
    std::string dc_output_path = fs::append_path(base_output_path_, "double_coverage_distribution.tsv");
    std::ofstream dc_stream(dc_output_path);
    dc_stream << "Gap\tFragments\tValue\n";
    for (const auto &gap_and_distr: double_cov_distribution) {
        const size_t &gap = gap_and_distr.first;
        const auto &distr = gap_and_distr.second;
        size_t total_fragments = 0;
        size_t total_pairs = 0;
        for (const auto &cov_and_number: distr) {
            dc_stream << gap << "\t" << cov_and_number.first << "\t" << cov_and_number.second << "\n";
            total_pairs += cov_and_number.first;
            total_fragments += cov_and_number.second;
        }
        double mean_double_cov = static_cast<double>(total_fragments) / static_cast<double>(total_pairs);
        INFO("Gap: " << gap << ", mean double coverage: " << mean_double_cov);
    }
}
LongEdgeStatisticsCounter::LongEdgeStatisticsCounter(const debruijn_graph::Graph &graph,
                                                     const barcode_index::FrameBarcodeIndex<debruijn_graph::Graph> &barcode_index,
                                                     const size_t &training_edge_length,
                                                     const size_t &training_edge_offset,
                                                     const size_t &read_count_threshold,
                                                     const size_t &read_linkage_distance,
                                                     const size_t &max_threads,
                                                     const std::string &base_output_path) :
    graph_(graph),
    barcode_index_(barcode_index),
    training_edge_length_(training_edge_length),
    training_edge_offset_(training_edge_offset),
    read_count_threshold_(read_count_threshold),
    read_linkage_distance_(read_linkage_distance),
    max_threads_(max_threads),
    base_output_path_(base_output_path) {}
void LongEdgeStatisticsCounter::CountClusterStatistics() const {
    using debruijn_graph::EdgeId;
    using debruijn_graph::Graph;
    using BarcodeExtractor = barcode_index::FrameBarcodeIndexInfoExtractor;
    using AccurateEdgeClusterExtractor = path_extend::read_cloud::cluster_storage::AccurateEdgeClusterExtractor;
    using EdgeInitialClusterStorageBuilder = path_extend::read_cloud::cluster_storage::EdgeInitialClusterStorageBuilder;
    auto barcode_extractor_ptr = std::make_shared<BarcodeExtractor>(barcode_index_, graph_);
    std::set<scaffold_graph::ScaffoldVertex> long_edges;
    for (const EdgeId &edge: graph_.canonical_edges()) {
        if (graph_.length(edge) >= training_edge_length_) {
            long_edges.insert(edge);
        }
    }
    INFO("Found " << long_edges.size() << " edges longer than " << training_edge_length_ << " in the graph");
    auto edge_cluster_extractor =
        std::make_shared<AccurateEdgeClusterExtractor>(graph_, barcode_extractor_ptr,
                                                       read_linkage_distance_, read_count_threshold_);
    EdgeInitialClusterStorageBuilder initial_builder(graph_, edge_cluster_extractor, long_edges,
                                                     read_linkage_distance_, read_count_threshold_, max_threads_);
    auto initial_cluster_storage = initial_builder.ConstructInitialClusterStorage();

    using DistributionExtractor = path_extend::read_cloud::fragment_statistics::ClusterDistributionExtractor;
    DistributionExtractor distribution_extractor(graph_, barcode_index_, read_count_threshold_, training_edge_length_,
                                                 training_edge_offset_, max_threads_);
    auto distribution_pack = distribution_extractor.GetDistributionsForStorage(initial_cluster_storage.get_cluster_storage());
    std::string length_output_path = fs::append_path(base_output_path_, "long_edge_fragment_lengths.tsv");
    std::ofstream length_stream(length_output_path);
    length_stream << "Length\tNumber\n";
    size_t total_long_edge_clusters = 0;
    for (const auto &entry: distribution_pack.length_distribution_) {
        length_stream << entry.first << "\t" << entry.second << std::endl;
        total_long_edge_clusters += entry.second;
    }
    INFO("Total long edge clusters: " << total_long_edge_clusters);
    std::string coverage_output_path = fs::append_path(base_output_path_, "long_edge_fragment_coverages.tsv");
    std::ofstream coverage_stream(coverage_output_path);
    coverage_stream << "Coverage\tNumber\n";
    for (const auto &entry: distribution_pack.coverage_distribution_) {
        coverage_stream << entry.first << "\t" << entry.second << std::endl;
    }
    std::string number_of_reads_output_path = fs::append_path(base_output_path_, "long_edge_fragment_reads.tsv");
    std::ofstream read_stream(number_of_reads_output_path);
    read_stream << "Number of reads\tNumber\n";
    for (const auto &entry: distribution_pack.num_reads_distribution_) {
        read_stream << entry.first << "\t" << entry.second << std::endl;
    }
}
}
