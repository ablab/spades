#pragma once

#include "distribution_extractor.hpp"

namespace path_extend {
namespace cluster_model {

template<class ClusterStatType>
struct DistributionStatistics {
  ClusterStatType min_;
  ClusterStatType max_;
  ClusterStatType mean_;
  ClusterStatType median_;

  DistributionStatistics(ClusterStatType min_, ClusterStatType max_, ClusterStatType mean_, ClusterStatType median_) :
      min_(min_), max_(max_), mean_(mean_), median_(median_) {}

  friend std::ostream& operator<< (std::ostream& stream, const DistributionStatistics<ClusterStatType> &statistics)
  {
      stream << "Min: " << statistics.min_ << std::endl
             << "Max: " << statistics.max_ << std::endl
             << "Mean: " << statistics.mean_ << std::endl
             << "Median: " << statistics.median_ << std::endl;
      return stream;
  }
};

struct StatisticsPack {
  DistributionStatistics<ClusterLength> length_statistics_;
  DistributionStatistics<ClusterCoverage> coverage_statistics_;

  StatisticsPack(const DistributionStatistics<ClusterLength> &length_statistics_,
                 const DistributionStatistics<ClusterCoverage> &coverage_statistics_) : length_statistics_(
      length_statistics_), coverage_statistics_(coverage_statistics_) {}
};

typedef std::map<size_t, StatisticsPack> StatisticsContainer;

class ClusterDistributionAnalyzer {
 private:
    const conj_graph_pack& gp_;
    size_t min_read_threshold_;
    size_t min_edge_length_;
    size_t min_cluster_offset_;
    size_t max_threads_;
 public:
    ClusterDistributionAnalyzer(const conj_graph_pack &gp_,
                                size_t min_read_threshold_,
                                size_t min_edge_length_,
                                size_t max_threads_)
        : gp_(gp_),
          min_read_threshold_(min_read_threshold_),
          min_edge_length_(min_edge_length_),
          max_threads_(max_threads_) {}

    DistributionPack GetDistributionsForDistance(size_t distance_threshold) {
        auto cluster_storage = GetInitialClusterStorage(distance_threshold);

//fixme remove code duplication by using boost::variant or smh else
//fixme lambdas everywhere
        auto cluster_predicate = [this](const cluster_storage::Cluster& cluster) {
          VERIFY_MSG(cluster.Size() == 1, "Cluster covers multiple edges!");
          path_extend::scaffold_graph::EdgeGetter edge_getter;
          auto map_info = cluster.GetMappings()[0];
          EdgeId cluster_edge = edge_getter.GetEdgeFromScaffoldVertex(map_info.GetEdge());
          size_t edge_length = this->gp_.g.length(cluster_edge);
          size_t left_pos = map_info.GetLeft();
          size_t right_pos = map_info.GetRight();
          return left_pos >= min_cluster_offset_ and right_pos + min_cluster_offset_ <= edge_length;
        };

        std::function<boost::optional<ClusterLength>(const cluster_storage::Cluster&)> length_extractor =
            [cluster_predicate](const cluster_storage::Cluster &cluster) {
                boost::optional<ClusterLength> result;
                if (cluster_predicate(cluster)) {
                    result = cluster.GetSpan();
                }
                return result;
            };
        std::function<boost::optional<ClusterCoverage>(const cluster_storage::Cluster&)> coverage_extractor =
            [cluster_predicate](const cluster_storage::Cluster &cluster) {
                boost::optional<ClusterCoverage> result;
                if (cluster_predicate(cluster)) {
                    result = cluster.GetCoverage();
                }
                return result;
            };
        SimpleDistributionExtractor distribution_extractor;

        DEBUG("Constructed extractors");

        auto length_distribution = distribution_extractor.ExtractDistribution(cluster_storage, length_extractor);
        auto coverage_distribution = distribution_extractor.ExtractDistribution(cluster_storage, coverage_extractor);

        DEBUG("Extracted distributions");

        DistributionPack result(length_distribution, coverage_distribution);
        cluster_storage.Clear();
        return result;
    }

    //fixme change method name or move serializing functionality elsewhere
    StatisticsContainer AnalyzeDistributions(const string& output_path, bool serialize_distributions=true) {
        fs::make_dir(output_path);
        const size_t min_distance = 1000;
        const size_t max_distance = 25000;
        const size_t distance_step = 4000;
        vector<size_t> distances;
        for (size_t distance = min_distance; distance <= max_distance; distance += distance_step) {
            distances.push_back(distance);
        }
        StatisticsContainer result;
        for (size_t distance: distances) {
            INFO("Getting distributions for distance " << distance);
            auto distributions = GetDistributionsForDistance(distance);
            if (serialize_distributions) {
                string distance_path = fs::append_path(output_path, std::to_string(distance));
                fs::make_dir(distance_path);
                PrintDistributions(distance_path, distributions);
            }
            auto length_statistics = GetDistributionStatistics(distributions.length_distribution_);
            auto coverage_statistics = GetDistributionStatistics(distributions.coverage_distribution_);
            StatisticsPack statistics(length_statistics, coverage_statistics);
            result.insert({distance, statistics});
        }
        return result;
    }

    void PrintStatistics(const StatisticsContainer& statistics, const string& output_path) {
        auto stats_stream = ofstream(fs::append_path(output_path, "statistics"));
        for (const auto& entry: statistics) {
            stats_stream << "Distance " << entry.first << ':' << std::endl;
            stats_stream << "Length" << std::endl;
            stats_stream << entry.second.length_statistics_ << std::endl;
            stats_stream << "Coverage" << std::endl;
            stats_stream << entry.second.coverage_statistics_;
            stats_stream << std::endl << "------------" << std::endl;
        }
    }

 private:
    template<class T>
    DistributionStatistics<T> GetDistributionStatistics(SimpleDistribution<T>& distribution) const {
        T min_element = std::numeric_limits<T>::max();
        T max_element = 0;
        T element_sum = 0;
        distribution.Sort();
        for (const T& element: distribution) {
            min_element = std::min(min_element, element);
            max_element = std::max(max_element, element);
            element_sum += element;
        }
        T mean = element_sum / static_cast<T>(distribution.Size());
        T median = distribution.At(distribution.Size() / 2);

        DistributionStatistics<T> result(min_element, max_element, mean, median);
        return result;
    }

    cluster_storage::ClusterStorage GetInitialClusterStorage(size_t distance_threshold) {
        auto barcode_extractor = make_shared<barcode_index::FrameBarcodeIndexInfoExtractor>(gp_.barcode_mapper_ptr, gp_.g);
        omnigraph::IterationHelper<Graph, EdgeId> edge_it_helper(gp_.g);
        std::set<scaffold_graph::ScaffoldVertex> long_edges;
        for (const auto& edge: edge_it_helper) {
            if (gp_.g.length(edge) >= min_edge_length_) {
                long_edges.insert(edge);
            }
        }
        INFO("Found " << long_edges.size() << " long edges in the graph");
        cluster_storage::EdgeInitialClusterStorageBuilder initial_builder(gp_.g, barcode_extractor, long_edges,
                                                                          distance_threshold, min_read_threshold_,
                                                                          max_threads_);
        auto initial_cluster_storage = initial_builder.ConstructInitialClusterStorage();
        return initial_cluster_storage.get_cluster_storage();
    }


    void PrintDistributions(const string& path, const DistributionPack& distributions) {
        auto length_stream = std::ofstream(fs::append_path(path, "length"));
        auto coverage_stream = std::ofstream(fs::append_path(path, "coverage"));
        length_stream << distributions.length_distribution_;
        coverage_stream << distributions.coverage_distribution_;
    }

    DECL_LOGGER("ClusterDistributionAnalyzer");
};
}
}