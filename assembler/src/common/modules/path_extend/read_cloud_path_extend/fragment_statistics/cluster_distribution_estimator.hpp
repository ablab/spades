//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "distribution_extractor.hpp"
#include "configs/pe_config_struct.hpp"
#include "modules/path_extend/read_cloud_path_extend/cluster_storage/initial_cluster_storage_builder.hpp"

namespace path_extend {
namespace read_cloud {
namespace fragment_statistics {

class MinTrainingLengthEstimator {
  public:
    using Graph = debruijn_graph::Graph;
    using EdgeId = debruijn_graph::EdgeId;

    MinTrainingLengthEstimator(const Graph &g, size_t min_total_length, size_t optimal_total_length, size_t min_edges)
        : g_(g),
          min_total_length_(min_total_length),
          optimal_total_length_(optimal_total_length),
          min_edges_(min_edges) {}

    boost::optional<size_t> EstimateTrainingLength() const;

  private:
    const Graph &g_;
    size_t min_total_length_;
    size_t optimal_total_length_;
    size_t min_edges_;
};

class MinTrainingLengthEstimatorHelper {
  public:
    using Graph = debruijn_graph::Graph;
    using ReadCloudConfigs = path_extend::pe_config::ReadCloud;

    MinTrainingLengthEstimatorHelper(const ReadCloudConfigs &configs) : configs_(configs) {}

    boost::optional<size_t> EstimateTrainingLength(const Graph &g) const;
  private:
    const ReadCloudConfigs &configs_;
};

struct StatisticsPack {
  typedef DistributionPack::ClusterLength ClusterLength;
  typedef DistributionPack::ClusterCoverage ClusterCoverage;

  DistributionStatistics<ClusterLength> length_statistics_;
  DistributionStatistics<ClusterCoverage> coverage_statistics_;

  StatisticsPack(const DistributionStatistics<ClusterLength> &length_statistics_,
                 const DistributionStatistics<ClusterCoverage> &coverage_statistics_) : length_statistics_(
      length_statistics_), coverage_statistics_(coverage_statistics_) {}
};

typedef std::map<size_t, StatisticsPack> StatisticsContainer;

class ClusterDistributionEstimator {
  public:
    typedef debruijn_graph::Graph Graph;
    typedef debruijn_graph::EdgeId EdgeId;
    typedef DistributionPack::ClusterLength ClusterLength;
    typedef DistributionPack::ClusterCoverage ClusterCoverage;
    typedef DistributionPack::ClusterLengthDistribution ClusterLengthDistribution;
    typedef DistributionPack::ClusterCoverageDistribution ClusterCoverageDistribution;
    typedef cluster_storage::ClusterStorage ClusterStorage;
    typedef cluster_storage::EdgeClusterExtractor EdgeClusterExtractor;
    typedef scaffold_graph::ScaffoldVertex ScaffoldVertex;
    typedef std::unordered_map<barcode_index::BarcodeId, std::vector<cluster_storage::Cluster>> ClusterCollection;

    ClusterDistributionEstimator(const Graph &g,
                                 const barcode_index::FrameBarcodeIndex<Graph> &barcode_mapper,
                                 size_t min_read_threshold,
                                 size_t min_edge_length,
                                 size_t min_cluster_offset,
                                 size_t max_threads)
        : g_(g),
          barcode_mapper_(barcode_mapper),
          min_read_threshold_(min_read_threshold),
          min_edge_length_(min_edge_length),
          min_cluster_offset_(min_cluster_offset),
          max_threads_(max_threads) {}

    DistributionPack GetClusterDistributions() const;

  private:
    DistributionPack EstimateFragmentDistributions(const std::unordered_set<ScaffoldVertex> &target_edges,
                                                   std::shared_ptr<EdgeClusterExtractor> edge_cluster_extractor) const;
    void UpdateDistributionFromEdge(DistributionPack& global_distribution, const DistributionPack& local_distribution) const;
    DistributionPack EstimateDistributionsFromClusters(const ClusterCollection &clusters) const;
    DistributionPack EstimateDistributionsForDistance(size_t distance_threshold) const;

    template<class DistributionT>
    DistributionStatistics<typename DistributionT::key_type> GetDistributionStatistics(
            const DistributionT &distribution) const;
    size_t EstimateDistance(const StatisticsContainer &statistics) const;

    const Graph &g_;
    const barcode_index::FrameBarcodeIndex<Graph> &barcode_mapper_;
    size_t min_read_threshold_;
    size_t min_edge_length_;
    size_t min_cluster_offset_;
    size_t max_threads_;

    DECL_LOGGER("ClusterDistributionExtractor");
};

class PercentileGetter {
  public:
    template<class DistributionT>
    typename DistributionT::key_type GetPercentile(const DistributionT &distribution, double percent) {
        size_t num_fragments = std::accumulate(distribution.begin(), distribution.end(), 0,
                                               [](const size_t prev, const auto &key_and_mult) {
                                                 return prev + key_and_mult.second;
                                               });

        auto index = static_cast<size_t>(static_cast<double>(num_fragments) * percent);
        size_t current_index = 0;
        typename DistributionT::key_type current_key = 0;
        for (const auto &key_and_mult: distribution) {
            current_index += key_and_mult.second;
            current_key = key_and_mult.first;
            if (current_index >= index) {
                return current_key;
            }
        }
        return current_key;
    }
};

class ClusterStatisticsExtractor {
  public:
    explicit ClusterStatisticsExtractor(const DistributionPack &cluster_distributions) :
        cluster_distributions_(cluster_distributions) {}

    DistributionPack GetDistributionPack() {
        return cluster_distributions_;
    }
    size_t GetLengthPercentile(double percent) {
        PercentileGetter percentile_getter;
        return percentile_getter.GetPercentile(cluster_distributions_.length_distribution_, percent);
    }
    double GetCoveragePercentile(double percent) {
        PercentileGetter percentile_getter;
        return percentile_getter.GetPercentile(cluster_distributions_.coverage_distribution_, percent);
    }

  private:
    DistributionPack cluster_distributions_;
};


class ClusterStatisticsExtractorHelper {
  public:
    typedef debruijn_graph::Graph Graph;
    typedef pe_config::ReadCloud ReadCloudConfigs;

    ClusterStatisticsExtractorHelper(const Graph &g,
                                     const barcode_index::FrameBarcodeIndex<Graph> &barcode_mapper,
                                     const ReadCloudConfigs &configs,
                                     size_t max_threads)
        : g_(g), barcode_mapper_(barcode_mapper), configs_(configs), max_threads_(max_threads) {}

    ClusterStatisticsExtractor GetStatisticsExtractor() const;

  private:
    const Graph &g_;
    const barcode_index::FrameBarcodeIndex<Graph> &barcode_mapper_;
    const ReadCloudConfigs &configs_;
    const size_t max_threads_;
};
}
}
}