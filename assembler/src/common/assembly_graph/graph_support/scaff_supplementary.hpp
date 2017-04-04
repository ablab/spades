#pragma once

#include "assembly_graph/core/graph.hpp"
#include "pipeline/graph_pack.hpp"
#include "utils/logger/logger.hpp"
//FIXME
#include "modules/path_extend/pe_utils.hpp"
#include "modules/path_extend/pe_config_struct.hpp"
#include "modules/path_extend/paired_library.hpp"

namespace path_extend {
typedef debruijn_graph::EdgeId EdgeId;

/* Storage of presumably unique, relatively long edges. Filled by ScaffoldingUniqueEdgeAnalyzer
 *
 */
class ScaffoldingUniqueEdgeStorage {
    friend class ScaffoldingUniqueEdgeAnalyzer;
    set<EdgeId> unique_edges_;
    size_t min_unique_length_;

public:
    ScaffoldingUniqueEdgeStorage(): unique_edges_(), min_unique_length_(0) {
        DEBUG("storage created, empty");
    }

    ScaffoldingUniqueEdgeStorage(const ScaffoldingUniqueEdgeStorage&) = delete;
    ScaffoldingUniqueEdgeStorage& operator=(const ScaffoldingUniqueEdgeStorage&) = delete;

    ScaffoldingUniqueEdgeStorage(ScaffoldingUniqueEdgeStorage&&) = default;
    ScaffoldingUniqueEdgeStorage& operator=(ScaffoldingUniqueEdgeStorage&&) = default;

    bool IsUnique(EdgeId e) const {
        return (unique_edges_.find(e) != unique_edges_.end());
    }

    decltype(unique_edges_.begin()) begin() const {
        return unique_edges_.begin();
    }

    decltype(unique_edges_.end()) end() const {
        return unique_edges_.end();
    }

    decltype(unique_edges_.begin()) erase(decltype(unique_edges_.begin()) iter) {
        return unique_edges_.erase(iter);
    }

    size_t size() const {
        return unique_edges_.size();
    }
    size_t min_length() const {
        return min_unique_length_;
    }
    void set_min_length(size_t min_length) {
        min_unique_length_ = min_length;
    }

    const set<EdgeId>& unique_edges() const {
        return unique_edges_;
    }

protected:
    DECL_LOGGER("ScaffoldingUniqueEdgeStorage");
};

//Auxillary class required to fillin the unique edge storage.


class ScaffoldingUniqueEdgeAnalyzer {
    const debruijn_graph::conj_graph_pack &gp_;
    size_t length_cutoff_;
    double median_coverage_;
    double relative_coverage_variation_;
//for uniqueness detection
    static const size_t max_different_edges_ = 20;
    static const size_t max_dijkstra_depth_ = 1000;
    static const size_t max_dijkstra_vertices_ = 1000;
    static const size_t overwhelming_majority_ = 10;
    set<VertexId> GetChildren(VertexId v, map<VertexId, set<VertexId>> &dijkstra_cash) const;
    bool FindCommonChildren(EdgeId e1, EdgeId e2, map<VertexId, set<VertexId>> &dijkstra_cash) const;
    bool FindCommonChildren(const vector<pair<EdgeId, double>> &next_weights) const;
    bool FindCommonChildren(EdgeId from, size_t lib_index) const;
    map<EdgeId, size_t> FillNextEdgeVoting(BidirectionalPathMap<size_t>& active_paths, int direction) const;
    bool ConservativeByPaths(EdgeId e, const GraphCoverageMap &long_reads_cov_map,
                             const pe_config::LongReads &lr_config) const;
    bool ConservativeByPaths(EdgeId e, const GraphCoverageMap &long_reads_cov_map,
                             const pe_config::LongReads &lr_config, int direction) const;
    bool ConservativeByLength(EdgeId e);
    void CheckCorrectness(ScaffoldingUniqueEdgeStorage& unique_storage_pb);
protected:
    DECL_LOGGER("ScaffoldingUniqueEdgeAnalyzer")


    void SetCoverageBasedCutoff();
public:
    ScaffoldingUniqueEdgeAnalyzer(const debruijn_graph::conj_graph_pack &gp, size_t apriori_length_cutoff,
                                  double max_relative_coverage):
            gp_(gp),
            length_cutoff_(apriori_length_cutoff),
            relative_coverage_variation_(max_relative_coverage) {
        SetCoverageBasedCutoff();
    }
    void FillUniqueEdgeStorage(ScaffoldingUniqueEdgeStorage &storage);
    void ClearLongEdgesWithPairedLib(size_t lib_index, ScaffoldingUniqueEdgeStorage &storage) const;
    void FillUniqueEdgesWithLongReads(GraphCoverageMap &long_reads_cov_map,
                                      ScaffoldingUniqueEdgeStorage &unique_storage_pb,
                                      const pe_config::LongReads &lr_config);
};

class UsedUniqueStorage {
    set<EdgeId> used_;
    const ScaffoldingUniqueEdgeStorage& unique_;

public:
    UsedUniqueStorage(const UsedUniqueStorage&) = delete;
    UsedUniqueStorage& operator=(const UsedUniqueStorage&) = delete;

    UsedUniqueStorage(UsedUniqueStorage&&) = default;
    UsedUniqueStorage& operator=(UsedUniqueStorage&&) = default;

    explicit UsedUniqueStorage(const ScaffoldingUniqueEdgeStorage& unique):
            unique_(unique) {}

    void insert(EdgeId e) {
        if (unique_.IsUnique(e)) {
            used_.insert(e);
            used_.insert(e->conjugate());
        }
    }

//    const ScaffoldingUniqueEdgeStorage& unique_edge_storage() const {
//        return unique_;
//    }

    bool IsUsedAndUnique(EdgeId e) const {
        return (unique_.IsUnique(e) && used_.find(e) != used_.end());
    }

    bool UniqueCheckEnabled() const {
        return unique_.size() > 0;
    }

    bool TryUseEdge(BidirectionalPath &path, EdgeId e, const Gap &gap) {
        if (UniqueCheckEnabled()) {
            if (IsUsedAndUnique(e)) {
                return false;
            } else {
                insert(e);
            }
        }
        path.PushBack(e, gap);
        return true;
    }

};

//FIXME rename
struct UniqueData {
    size_t min_unique_length_;
    double unique_variation_;

    ScaffoldingUniqueEdgeStorage main_unique_storage_;
    vector<ScaffoldingUniqueEdgeStorage> unique_storages_;

    ScaffoldingUniqueEdgeStorage unique_pb_storage_;
    vector<PathContainer> long_reads_paths_;
    vector<GraphCoverageMap> long_reads_cov_map_;
};
}


