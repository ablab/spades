//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "assembly_graph/paths/bidirectional_path_container.hpp"
#include "auxiliary_graphs/scaffold_graph/scaffold_graph.hpp"
#include "configs/pe_config_struct.hpp"

namespace path_extend {

class StartFinder {
    typedef scaffold_graph::ScaffoldVertex ScaffoldVertex;
    typedef std::unordered_map<ScaffoldVertex, ScaffoldVertex> TransitionMap;
    typedef debruijn_graph::Graph Graph;

    const debruijn_graph::Graph &g_;
  public:
    StartFinder(const debruijn_graph::Graph &g);

    std::unordered_set<ScaffoldVertex> GetStarts(const TransitionMap &transition_map) const;
};

class PathScaffolder {
  public:
    typedef scaffold_graph::ScaffoldGraph ScaffoldGraph;

    virtual ~PathScaffolder() = default;

    virtual void MergePaths(const ScaffoldGraph &scaffold_graph) const = 0;
};

class SimplePathScaffolder : public PathScaffolder {
  public:
    using PathScaffolder::ScaffoldGraph;
    typedef scaffold_graph::ScaffoldVertex ScaffoldVertex;
    typedef ScaffoldGraph::ScaffoldEdge ScaffoldEdge;
    typedef debruijn_graph::Graph Graph;

    SimplePathScaffolder(const debruijn_graph::Graph &g, int default_gap);

    void MergePaths(const ScaffoldGraph &scaffold_graph) const override;

  private:

    void CondenseSimplePaths(const std::vector<ScaffoldEdge> &scaffold_edges) const;
    void ExtendPathAlongConnections(const ScaffoldVertex &start,
                                    const std::unordered_map<ScaffoldVertex, ScaffoldVertex> &merge_connections,
                                    const std::unordered_map<ScaffoldVertex, size_t> &start_to_length) const;

    const debruijn_graph::Graph &g_;
    const int default_gap_;

    DECL_LOGGER("SimplePathScaffolder");
};

}
