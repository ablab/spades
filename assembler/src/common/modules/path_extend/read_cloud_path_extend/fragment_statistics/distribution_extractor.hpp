#pragma once

#include <map>

#include "common/utils/stl_utils.hpp"
#include "common/barcode_index/barcode_info_extractor.hpp"
#include "modules/path_extend/read_cloud_path_extend/cluster_storage/cluster_storage.hpp"

using namespace debruijn_graph;

namespace path_extend {
namespace read_cloud {
namespace fragment_statistics {

template<class T>
class SimplePrecisionComp {
  public:
    SimplePrecisionComp(const T epsilon = 1e-7) : epsilon_(epsilon) {}

    bool operator()(const T &first, const T &second) const {
        return math::ls(first + epsilon_, second);
    }
  private:
    T epsilon_;
};

struct DistributionPack {
  public:
    typedef uint64_t ClusterLength;
    typedef double ClusterCoverage;
    typedef std::map<ClusterLength, size_t> ClusterLengthDistribution;
    typedef std::map<ClusterCoverage, size_t, SimplePrecisionComp<ClusterCoverage>> ClusterCoverageDistribution;

    DistributionPack() :
        length_distribution_(), coverage_distribution_() {}

    explicit DistributionPack(const ClusterLengthDistribution &length_distribution) :
        length_distribution_(length_distribution), coverage_distribution_() {}

    DistributionPack(const ClusterLengthDistribution &length_distribution,
                     const ClusterCoverageDistribution &coverage_distribution) :
        length_distribution_(length_distribution),
        coverage_distribution_(coverage_distribution) {}

    ClusterLengthDistribution length_distribution_;
    ClusterCoverageDistribution coverage_distribution_;
};

class SimpleDistributionExtractor {
  public:
    typedef cluster_storage::ClusterStorage ClusterStorage;
    typedef cluster_storage::Cluster Cluster;

    template<class DistributionT>
    DistributionT ExtractDistribution(
        const ClusterStorage &cluster_storage,
        std::function<boost::optional<typename DistributionT::key_type>(const Cluster &)> cluster_stat_getter) {
        DistributionT result;
        for (const auto &entry: cluster_storage) {
            const cluster_storage::Cluster cluster = entry.second;
            auto statistic = cluster_stat_getter(cluster);
            if (statistic.is_initialized()) {
                ++result[statistic.get()];
            }
        }
        return result;
    }
};

}
}
}