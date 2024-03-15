//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "assembly_graph/paths/bidirectional_path.hpp"
#include "assembly_graph/paths/bidirectional_path_container.hpp"

namespace path_extend {

class CompositeExtender;
class GraphCoverageMap;

void Deduplicate(const debruijn_graph::Graph &g, PathContainer &paths, GraphCoverageMap &coverage_map,
                 size_t min_edge_len, size_t max_path_diff,
                 bool equal_only = false);

class PathExtendResolver {
    const debruijn_graph::Graph& g_;
    size_t k_;
public:
    PathExtendResolver(const debruijn_graph::Graph& g)
            : g_(g), k_(g.k()) {}
    
    PathContainer MakeSimpleSeeds() const;
    PathContainer ExtendSeeds(PathContainer &seeds, CompositeExtender &composite_extender) const;

    //Paths should be deduplicated first!
    void RemoveOverlaps(PathContainer &paths, GraphCoverageMap &coverage_map,
                        size_t min_edge_len, size_t max_path_diff,
                        bool end_start_only, bool cut_all) const;
    void AddUncoveredEdges(PathContainer &paths, GraphCoverageMap &coverageMap) const;
protected:
    DECL_LOGGER("PEResolver")
};

}
