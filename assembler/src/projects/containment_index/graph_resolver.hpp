//***************************************************************************
//* Copyright (c) 2021-2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "vertex_resolver.hpp"

#include "assembly_graph/core/observable_graph.hpp"
#include "assembly_graph/paths/bidirectional_path_container.hpp"

#pragma once

namespace cont_index {
class GraphResolver {
  public:
    using EdgeId = debruijn_graph::EdgeId;
    using VertexId = debruijn_graph::VertexId;
    using LinkId = debruijn_graph::LinkId;
    using LinkMap = std::unordered_map<EdgeId, LinkId>;

    struct GraphResolverInfo {
      public:
        using VertexMap = std::unordered_map<VertexId, VertexId>;
        using EdgeMap = std::unordered_map<EdgeId, EdgeId>;

        GraphResolverInfo(const VertexMap &transformed_vertex_to_original, const EdgeMap &original_edge_to_transformed)
            : transformed_vertex_to_original(transformed_vertex_to_original),
              original_edge_to_transformed(original_edge_to_transformed) {}
        VertexMap transformed_vertex_to_original;
        EdgeMap original_edge_to_transformed;
    };

    GraphResolverInfo TransformGraph(debruijn_graph::Graph &graph,
                                     const path_extend::PathContainer &paths,
                                     const VertexResults &vertex_results) const;
  private:
    GraphResolverInfo::VertexMap SplitVertices(debruijn_graph::Graph &graph,
                                               const VertexResults &vertex_results) const;
    GraphResolverInfo::EdgeMap MergePaths(debruijn_graph::Graph &graph, const path_extend::PathContainer &paths) const;
    LinkMap GetLinkMap(const debruijn_graph::Graph &graph,
                       const VertexId &vertex,
                       const VertexResult &vertex_result) const;
};
}
