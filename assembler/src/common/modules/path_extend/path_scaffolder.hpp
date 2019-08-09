//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "common/assembly_graph/paths/bidirectional_path_container.hpp"
#include "common/modules/path_extend/scaffolder2015/scaffold_graph.hpp"
#include "common/pipeline/graph_pack.hpp"
#include "common/pipeline/config_struct.hpp"

namespace path_extend {

class StartFinder {
    typedef path_extend::scaffold_graph::ScaffoldVertex ScaffoldVertex;
    typedef std::unordered_map<ScaffoldVertex, ScaffoldVertex> TransitionMap;

    const Graph &g_;
  public:
    StartFinder(const Graph &g);

    std::unordered_set<ScaffoldVertex> GetStarts(const TransitionMap &transition_map) const;
};

class PathScaffolder {
  public:
    typedef path_extend::scaffold_graph::ScaffoldGraph ScaffoldGraph;

    virtual ~PathScaffolder() = default;

    virtual void MergePaths(const ScaffoldGraph &scaffold_graph) const = 0;
};

class SimplePathScaffolder : public PathScaffolder {
  public:
    using PathScaffolder::ScaffoldGraph;
    typedef path_extend::scaffold_graph::ScaffoldVertex ScaffoldVertex;
    typedef ScaffoldGraph::ScaffoldEdge ScaffoldEdge;

    SimplePathScaffolder(const conj_graph_pack &gp, int default_gap);

    void MergePaths(const ScaffoldGraph &scaffold_graph) const override;

  private:

    void CondenseSimplePaths(const std::vector<ScaffoldEdge> &scaffold_edges) const;
    void ExtendPathAlongConnections(const ScaffoldVertex &start,
                                    const std::unordered_map<ScaffoldVertex, ScaffoldVertex> &merge_connections,
                                    const std::unordered_map<ScaffoldVertex, size_t> &start_to_length) const;

    const debruijn_graph::conj_graph_pack &gp_;
    const int default_gap_;

    DECL_LOGGER("SimplePathScaffolder");
};

}
