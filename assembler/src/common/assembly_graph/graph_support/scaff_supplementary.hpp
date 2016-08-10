#pragma once

#include "assembly_graph/core/graph.hpp"
#include "pipeline/graph_pack.hpp"
#include "utils/logger/logger.hpp"
//FIXME
#include "modules/path_extend/pe_utils.hpp"

namespace path_extend {
typedef debruijn_graph::EdgeId EdgeId;

/* Storage of presumably unique, relatively long edges. Filled by ScaffoldingUniqueEdgeAnalyzer
 *
 */
class ScaffoldingUniqueEdgeStorage {
    friend class ScaffoldingUniqueEdgeAnalyzer;
private:
    set <EdgeId> unique_edges_;
    size_t min_unique_length_;
public:
    ScaffoldingUniqueEdgeStorage(): unique_edges_(){
        DEBUG("storage created, empty");
    }

    bool IsUnique(EdgeId e) const {
        return (unique_edges_.find(e) != unique_edges_.end());
    }

    decltype(unique_edges_.begin()) begin() const {
        return unique_edges_.begin();
    }

    decltype(unique_edges_.end()) end() const {
        return unique_edges_.end();
    }

    size_t size() const {
        return unique_edges_.size();
    }
    size_t GetMinLength() const {
        return min_unique_length_;
    }
    void SetMinLength(size_t min_length)  {
        min_unique_length_ = min_length;
    }

    const set<EdgeId>& GetSet() const {
        return unique_edges_;
    }

protected:
    DECL_LOGGER("ScaffoldingUniqueEdgeStorage")

};

//Auxillary class required to fillin the unique edge storage.


class ScaffoldingUniqueEdgeAnalyzer {
    const debruijn_graph::conj_graph_pack &gp_;
    size_t length_cutoff_;
    double median_coverage_;
    double relative_coverage_variation_;

    bool CheckSuffixConservative(const BidirectionalPath& path1, size_t pos1, const BidirectionalPath& path2, size_t pos2) const;
    bool CheckPrefixConservative(const BidirectionalPath& path1, size_t pos1, const BidirectionalPath& path2, size_t pos2) const;
    bool ConsistentPath(const BidirectionalPath& path1, size_t pos1, const BidirectionalPath& path2, size_t pos2) const;
    bool ConservativeByPaths(EdgeId e, shared_ptr<GraphCoverageMap>& long_reads_cov_map);
    bool ConservativeByLength(EdgeId e);
    void CheckCorrectness(ScaffoldingUniqueEdgeStorage& unique_storage_pb);
protected:
    DECL_LOGGER("ScaffoldingUniqueEdgeAnalyzer")


    void SetCoverageBasedCutoff();
public:
    ScaffoldingUniqueEdgeAnalyzer(const debruijn_graph::conj_graph_pack &gp, size_t apriori_length_cutoff, double max_relative_coverage):gp_(gp), length_cutoff_(apriori_length_cutoff), relative_coverage_variation_(max_relative_coverage){
        SetCoverageBasedCutoff();
    }
    void FillUniqueEdgeStorage(ScaffoldingUniqueEdgeStorage &storage_);
    void FillUniqueEdgesWithLongReads(shared_ptr<GraphCoverageMap>& long_reads_cov_map, ScaffoldingUniqueEdgeStorage& unique_storage_pb);
};
}


