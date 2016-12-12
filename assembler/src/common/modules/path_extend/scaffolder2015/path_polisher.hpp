#pragma once

#include "assembly_graph/paths/path_processor.hpp"
#include "assembly_graph/paths/path_utils.hpp"
#include "assembly_graph/paths/bidirectional_path.hpp"
#include "assembly_graph/core/basic_graph_stats.hpp"
#include "modules/path_extend/paired_library.hpp"
#include "assembly_graph/graph_support/scaff_supplementary.hpp"
#include "common/pipeline/graph_pack.hpp"

namespace path_extend {

class PathGapCloser {
protected:
    const Graph& g_;
    size_t max_path_len_;
    int min_gap_;
public:
    virtual BidirectionalPath Polish(const BidirectionalPath& path) = 0;
//TODO:: config
    PathGapCloser(const Graph& g, size_t max_path_len): g_(g), max_path_len_(max_path_len), min_gap_(int(g.k() + 10)) {}

};

class MatePairGapCloser: public PathGapCloser {
    const shared_ptr<PairedInfoLibrary> lib_;
    const ScaffoldingUniqueEdgeStorage& storage_;

//TODO: config? somewhere else?
    static constexpr double weight_priority = 5;
public:
    EdgeId FindNext(const BidirectionalPath& path, size_t index,
                        const set<EdgeId>& present_in_paths, VertexId v) const;
    MatePairGapCloser(const Graph& g, size_t max_path_len, const shared_ptr<PairedInfoLibrary> lib, const ScaffoldingUniqueEdgeStorage& storage):
            PathGapCloser(g, max_path_len), lib_(lib), storage_(storage) {}
    BidirectionalPath Polish(const BidirectionalPath& path) override;
};

class DijkstraGapCloser: public PathGapCloser {

protected:

    BidirectionalPath Polish(const BidirectionalPath& path) override;

    size_t MinPathLength(const omnigraph::PathStorageCallback<Graph>& path_storage) const;

    bool FillWithMultiplePaths(const BidirectionalPath& path, size_t index,
                                       const omnigraph::PathStorageCallback<Graph>& path_storage,
                                       BidirectionalPath& result) const;

    bool FillWithBridge(const BidirectionalPath& path, size_t index,
                                                                  const omnigraph::PathStorageCallback<Graph>& path_storage,
                                                                  BidirectionalPath& result) const;

    size_t MinPathSize(const omnigraph::PathStorageCallback<Graph>& path_storage) const;

    vector<EdgeId> LCP(const omnigraph::PathStorageCallback<Graph>& path_storage) const;

    std::map<EdgeId, size_t> CountEdgesQuantity(const omnigraph::PathStorageCallback<Graph>& path_storage, size_t length_limit) const;

public:
    DijkstraGapCloser(const Graph& g, size_t max_path_len):
        PathGapCloser(g, max_path_len) {}


};

class PathPolisher {

private:
    const conj_graph_pack& gp_;
    vector<shared_ptr<PathGapCloser>> gap_closers;

private:
    void InfoAboutGaps(const PathContainer & result);
    BidirectionalPath Polish(const BidirectionalPath& path);

public:
    PathPolisher(const conj_graph_pack& gp, const config::dataset& dataset_info, const ScaffoldingUniqueEdgeStorage& storage, size_t max_resolvable_len);

    void PolishPaths(const PathContainer& paths, PathContainer& result);
};


}
