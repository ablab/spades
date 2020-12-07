//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "assembly_graph/paths/bidirectional_path.hpp"
#include "overlap_remover.hpp"
#include "path_extender.hpp"

namespace path_extend {

inline void Deduplicate(const Graph &g, PathContainer &paths, GraphCoverageMap &coverage_map,
                        size_t min_edge_len, size_t max_path_diff,
                        bool equal_only = false) {
    //add sorting to guarantee survival of longest paths if max_path_diff used
    //paths.SortByLength(false);
    PathDeduplicator deduplicator(g, paths, coverage_map, min_edge_len, max_path_diff, equal_only);
    deduplicator.Deduplicate();
    paths.FilterEmptyPaths();
}

//Checks whether we are in a cycle of length 2, used only for seed selection.
inline bool InTwoEdgeCycle(EdgeId e, const Graph &g) {
    VertexId v = g.EdgeEnd(e);
    //allow to start from long edge with potential coverage one.
    if (g.OutgoingEdgeCount(v) >= 1 && (g.length(e) < 1000 || g.OutgoingEdgeCount(v) > 1)) {
        for (EdgeId e1 : g.OutgoingEdges(v)) {
            if (g.EdgeStart(e) == g.EdgeEnd(e1))
                return true;
        }
    }
    return false;
}

class PathExtendResolver {
    const Graph& g_;
    size_t k_;
public:
    PathExtendResolver(const Graph& g)
            : g_(g), k_(g.k()) {}
    
    PathContainer MakeSimpleSeeds() const {
        std::set<EdgeId> included;
        PathContainer edges;
        for (auto iter = g_.ConstEdgeBegin(/*canonical only*/true); !iter.IsEnd(); ++iter) {
            EdgeId e = *iter;
            if (g_.int_id(e) <= 0 || InTwoEdgeCycle(e, g_))
                continue;
            edges.CreatePair(g_, e);
        }
        return edges;
    }

    PathContainer ExtendSeeds(PathContainer &seeds, CompositeExtender &composite_extender) const {
        PathContainer paths;
        composite_extender.GrowAll(seeds, paths);
        return paths;
    }

    //Paths should be deduplicated first!
    void RemoveOverlaps(PathContainer &paths, GraphCoverageMap &coverage_map,
                        size_t min_edge_len, size_t max_path_diff,
                        bool end_start_only, bool cut_all) const {
        INFO("Removing overlaps");
        //VERIFY(min_edge_len == 0 && max_path_diff == 0);
        if (!cut_all) {
            INFO("Sorting paths");
            //sorting is currently needed to retain overlap instance in longest paths
            paths.SortByLength(false);
        }

        OverlapRemover overlap_remover(g_, paths, coverage_map,
                                       min_edge_len, max_path_diff);
        INFO("Marking overlaps");
        overlap_remover.MarkOverlaps(end_start_only, !cut_all);

        INFO("Splitting paths");
        PathSplitter splitter(overlap_remover.overlaps(), paths, coverage_map);
        splitter.Split();
        //splits are invalidated after this point

        INFO("Deduplicating paths");
        Deduplicate(g_, paths, coverage_map, min_edge_len, max_path_diff);
        INFO("Overlaps removed");
    }

    void AddUncoveredEdges(PathContainer &paths, GraphCoverageMap &coverageMap) const {
        for (auto iter = g_.ConstEdgeBegin(true); !iter.IsEnd(); ++iter) {
            EdgeId e = *iter;
            if (!coverageMap.IsCovered(e)) {
                AddPath(paths, BidirectionalPath::create(g_, e), coverageMap);
            }
        }
    }
protected:
    DECL_LOGGER("PEResolver")
};

}
