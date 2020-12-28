//***************************************************************************
//* Copyright (c) 2020 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "pe_resolver.hpp"
#include "overlap_remover.hpp"
#include "path_deduplicator.hpp"
#include "path_extender.hpp"

namespace path_extend {

using namespace debruijn_graph;

void Deduplicate(const Graph &g, PathContainer &paths, GraphCoverageMap &coverage_map,
                 size_t min_edge_len, size_t max_path_diff,
                 bool equal_only) {
    //add sorting to guarantee survival of longest paths if max_path_diff used
    //paths.SortByLength(false);
    PathDeduplicator deduplicator(g, paths, coverage_map, min_edge_len, max_path_diff, equal_only);
    deduplicator.Deduplicate();
    paths.FilterEmptyPaths();
}

//Checks whether we are in a cycle of length 2, used only for seed selection.
static bool InTwoEdgeCycle(EdgeId e, const Graph &g) {
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

PathContainer PathExtendResolver::MakeSimpleSeeds() const {
    std::set<EdgeId> included;
    PathContainer edges;
    for (EdgeId e : g_.canonical_edges()) {
        if (g_.int_id(e) <= 0 || InTwoEdgeCycle(e, g_))
            continue;
        edges.Create(g_, e);
    }
    return edges;
}

PathContainer PathExtendResolver::ExtendSeeds(PathContainer &seeds, CompositeExtender &composite_extender) const {
    PathContainer paths;
    composite_extender.GrowAll(seeds, paths);
    return paths;
}

//Paths should be deduplicated first!
void PathExtendResolver::RemoveOverlaps(PathContainer &paths, GraphCoverageMap &coverage_map,
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

void PathExtendResolver::AddUncoveredEdges(PathContainer &paths, GraphCoverageMap &coverageMap) const {
    for (EdgeId e : g_.canonical_edges()) {
        if (coverageMap.IsCovered(e))
            continue;

        CreatePath(paths, coverageMap,
                   g_, e);
    }
}

}
