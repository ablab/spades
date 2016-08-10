//
// Created by andrey on 12.05.16.
//
#pragma once

#include "assembly_graph/paths/path_processor.hpp"
#include "assembly_graph/paths/path_utils.hpp"
#include "assembly_graph/paths/bidirectional_path.hpp"
#include "assembly_graph/core/basic_graph_stats.hpp"


namespace path_extend {

class PathGapCloser {

protected:

    const Graph& g_;

    virtual BidirectionalPath Polish(const BidirectionalPath& path) = 0;

public:
    PathGapCloser(const Graph& g): g_(g) {}

    virtual void PolishPaths(const PathContainer& paths, PathContainer& result);
};


class DijkstraGapCloser: public PathGapCloser {

protected:
    size_t max_path_len_;

    BidirectionalPath Polish(const BidirectionalPath& path) override;

    size_t MinPathLength(const omnigraph::PathStorageCallback<Graph>& path_storage) const;

    virtual void FillWithMultiplePaths(const BidirectionalPath& path, size_t index,
                                       const omnigraph::PathStorageCallback<Graph>& path_storage,
                                       BidirectionalPath& result);

public:
    DijkstraGapCloser(const Graph& g, size_t max_path_len):
        PathGapCloser(g), max_path_len_(max_path_len) {}


};


class CommonPrefixDijkstraGapCloser: public DijkstraGapCloser {

protected:

    size_t MinPathSize(const omnigraph::PathStorageCallback<Graph>& path_storage) const;

    vector<EdgeId> LCS(const omnigraph::PathStorageCallback<Graph>& path_storage) const;

    vector<EdgeId> LCP(const omnigraph::PathStorageCallback<Graph>& path_storage) const;

    void FillWithMultiplePaths(const BidirectionalPath& path, size_t index,
                               const omnigraph::PathStorageCallback<Graph>& path_storage,
                               BidirectionalPath& result) override ;

public:
    CommonPrefixDijkstraGapCloser(const Graph& g, size_t max_path_len):
        DijkstraGapCloser(g, max_path_len) {}


};

}
