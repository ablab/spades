//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2018-2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "assembly_graph/components/graph_component.hpp"
#include "assembly_graph/core/graph.hpp"
#include "io/utils/edge_namer.hpp"

#include <ostream>

namespace omnigraph {

template <class Graph>
class GraphComponent;
}

namespace gfa {

class GFAWriter {
  protected:
    typedef debruijn_graph::DeBruijnGraph Graph;
    typedef debruijn_graph::VertexId VertexId;
    typedef omnigraph::GraphComponent<Graph> Component;

public:
    GFAWriter(const Graph &graph, std::ostream &os,
              io::EdgeNamingF<Graph> naming_f = io::IdNamingF<Graph>())
            : graph_(graph),
              edge_namer_(graph_, naming_f),
              os_(os) {}

    virtual ~GFAWriter() = default;

    void WriteSegmentsAndLinks() {
        WriteHeader();
        WriteSegments();
        WriteLinks();
    }
    void WriteSegmentsAndLinks(const Component &gc);

  protected:
    virtual void WriteHeader();

  private:
    void WriteSegments();
    void WriteLinks();

    void WriteSegments(const Component &gc);
    void WriteLinks(const Component &gc);

    void WriteVertexLinks(const VertexId &vertex);
    void WriteVertexLinks(const VertexId &vertex, const Component &gc);

  protected:
    const Graph &graph_;
    io::CanonicalEdgeHelper<Graph> edge_namer_;
    std::ostream &os_;
};

}
