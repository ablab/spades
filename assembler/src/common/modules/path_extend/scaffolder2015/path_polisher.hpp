#pragma once

#include "assembly_graph/paths/path_processor.hpp"
#include "assembly_graph/paths/path_utils.hpp"
#include "assembly_graph/paths/bidirectional_path.hpp"
#include "assembly_graph/core/basic_graph_stats.hpp"
#include "modules/path_extend/paired_library.hpp"
#include "modules/path_extend/path_extender.hpp"
#include "assembly_graph/graph_support/scaff_supplementary.hpp"
#include "pipeline/graph_pack.hpp"

namespace path_extend {

class PathGapCloser {
protected:
    const Graph& g_;
    const size_t max_path_len_;
    const int min_gap_;

    virtual Gap CloseGap(const BidirectionalPath &original_path, size_t position,
                         BidirectionalPath &path) const = 0;
    DECL_LOGGER("PathGapCloser")
public:
    BidirectionalPath CloseGaps(const BidirectionalPath &path) const;

    PathGapCloser(const Graph& g, size_t max_path_len):
                  g_(g),
                  max_path_len_(max_path_len),
                  //TODO:: config
                  min_gap_(int(g.k() + 10)) {}

};

//Intermediate abstract class - majority of GapClosers needs only one next edge after gap, not all original path.
class TargetEdgeGapCloser : public PathGapCloser {
protected:
    //returns updated gap to target edge
    virtual Gap CloseGap(EdgeId target_edge, const Gap &gap, BidirectionalPath &path) const = 0;

    Gap CloseGap(const BidirectionalPath &original_path,
                 size_t position, BidirectionalPath &path) const final override {
        return CloseGap(original_path.At(position), original_path.GapAt(position), path);
    }

public:
    TargetEdgeGapCloser(const Graph& g, size_t max_path_len):
            PathGapCloser(g, max_path_len) {}

};

class PathExtenderGapCloser: public TargetEdgeGapCloser {
    shared_ptr<path_extend::PathExtender> extender_;

protected:
    Gap CloseGap(EdgeId target_edge, const Gap &gap, BidirectionalPath &path) const override;

public:
    PathExtenderGapCloser(const Graph& g, size_t max_path_len, shared_ptr<PathExtender> extender):
            TargetEdgeGapCloser(g, max_path_len), extender_(extender) {
        DEBUG("ext added");
    }
};

class MatePairGapCloser: public TargetEdgeGapCloser {
    const shared_ptr<PairedInfoLibrary> lib_;
    const ScaffoldingUniqueEdgeStorage& storage_;
//TODO: config? somewhere else?
    static constexpr double weight_priority = 5;

    EdgeId FindNext(const BidirectionalPath& path,
                    const set<EdgeId>& present_in_paths,
                    VertexId last_v, EdgeId target_edge) const;
protected:
    Gap CloseGap(EdgeId target_edge, const Gap &gap, BidirectionalPath &path) const override;

    DECL_LOGGER("MatePairGapCloser")

public:
    MatePairGapCloser(const Graph& g, size_t max_path_len,
                      const shared_ptr<PairedInfoLibrary> lib,
                      const ScaffoldingUniqueEdgeStorage& storage):
            TargetEdgeGapCloser(g, max_path_len), lib_(lib), storage_(storage) {}
};

//TODO switch to a different Callback, no need to store all paths
class DijkstraGapCloser: public TargetEdgeGapCloser {
    typedef vector<vector<EdgeId>> PathsT;

    Gap FillWithMultiplePaths(const PathsT& paths,
                              BidirectionalPath& result) const;

    Gap FillWithBridge(const Gap &orig_gap,
                       const PathsT& paths, BidirectionalPath& result) const;

    size_t MinPathLength(const PathsT& paths) const;

    size_t MinPathSize(const PathsT& paths) const;

    vector<EdgeId> LCP(const PathsT& paths) const;

    std::map<EdgeId, size_t> CountEdgesQuantity(const PathsT& paths, size_t length_limit) const;

protected:
    Gap CloseGap(EdgeId target_edge, const Gap &gap, BidirectionalPath &path) const override;

    DECL_LOGGER("DijkstraGapCloser")

public:
    DijkstraGapCloser(const Graph& g, size_t max_path_len):
        TargetEdgeGapCloser(g, max_path_len) {}

};

class PathPolisher {
    static const size_t MAX_POLISH_ATTEMPTS = 5;

    const conj_graph_pack &gp_;
    vector<shared_ptr<PathGapCloser>> gap_closers_;

    void InfoAboutGaps(const PathContainer& result);

    BidirectionalPath Polish(const BidirectionalPath& path);
    DECL_LOGGER("PathPolisher")

public:
    PathPolisher(const conj_graph_pack &gp,
                               const vector<shared_ptr<PathGapCloser>> &gap_closers):
            gp_(gp), gap_closers_(gap_closers) {
    }

    PathContainer PolishPaths(const PathContainer &paths);
};


}
