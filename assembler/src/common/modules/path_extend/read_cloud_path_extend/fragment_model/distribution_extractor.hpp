#pragma once

#include "common/barcode_index/cluster_storage_extractor.hpp"
#include "common/barcode_index/barcode_info_extractor.hpp"

using namespace debruijn_graph;

namespace path_extend {
namespace cluster_model {

typedef size_t ClusterLength;
typedef double ClusterCoverage;

template<class T>
class SimpleDistribution {
    vector<T> distribution_;
    bool is_sorted_;

 public:
    SimpleDistribution()
        : distribution_(), is_sorted_(false) {}
    SimpleDistribution(const vector<T> &distribution_)
        : distribution_(distribution_), is_sorted_(false) {}

    void add(const boost::optional<T> &value) {
        if (value.is_initialized()) {
            distribution_.push_back(value.get());
        }
        is_sorted_ = false;
    }

    void sort() {
        std::sort(distribution_.begin(), distribution_.end());
        is_sorted_ = true;
    }

    bool is_sorted() const {
        return is_sorted_;
    };

    T at(size_t index) const {
        return distribution_.at(index);
    }

    size_t size() const {
        return distribution_.size();
    }

    auto begin() const {
        return distribution_.begin();
    }
    auto end() const {
        return distribution_.end();
    }

    friend std::ostream& operator<< (std::ostream& stream, const SimpleDistribution<T>& distribution) {
        for (const auto& element: distribution) {
            stream << element << '\t';
        }
        return stream;
    }
};

typedef SimpleDistribution<ClusterLength> ClusterLengthDistribution;
typedef SimpleDistribution<ClusterCoverage> ClusterCoverageDistribution;

struct DistributionPack {
  ClusterLengthDistribution length_distribution_;
  ClusterCoverageDistribution coverage_distribution_;

  DistributionPack(const ClusterLengthDistribution &length_distribution_,
                   const ClusterCoverageDistribution &coverage_distribution_) :
      length_distribution_(length_distribution_),
      coverage_distribution_(coverage_distribution_) {}
};

class SimpleDistributionExtractor {
 public:
    template<class ClusterStatType>
    SimpleDistribution<ClusterStatType> ExtractDistribution(
            const cluster_storage::ClusterStorage &cluster_storage,
            std::function<boost::optional<ClusterStatType>(const cluster_storage::Cluster &)> extractor_from_cluster) {
        SimpleDistribution<ClusterStatType> result;
        for (const auto &entry: cluster_storage) {
            const cluster_storage::Cluster cluster = entry.second;
            result.add(extractor_from_cluster(cluster));
        }
        return result;
    }
};

}
}