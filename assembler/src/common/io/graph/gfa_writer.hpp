//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "assembly_graph/components/graph_component.hpp"
#include "assembly_graph/core/graph.hpp"
#include "io/utils/edge_namer.hpp"
#include "io/utils/id_mapper.hpp"

#include <memory>
#include <string>
#include <ostream>

namespace omnigraph {

template <class Graph>
class GraphComponent;

}

namespace gfa {

class GFAWriter {
  protected:
    typedef debruijn_graph::DeBruijnGraph Graph;
    typedef omnigraph::GraphComponent<Graph> Component;

public:
    GFAWriter(const Graph &graph, std::ostream &os,
              io::EdgeNamingF<Graph> naming_f = io::IdNamingF<Graph>())
            : graph_(graph),
              edge_namer_(graph_, naming_f),
              os_(os) {
    }

    void WriteSegmentsAndLinks() {
        WriteSegments();
        WriteLinks();
    }
    void WriteSegmentsAndLinks(const Component &gc);

  private:
    void WriteSegments();
    void WriteLinks();

    void WriteSegments(const Component &gc);
    void WriteLinks(const Component &gc);

  protected:
    const Graph &graph_;
    io::CanonicalEdgeHelper<Graph> edge_namer_;
    std::ostream &os_;
};

class GFAComponentWriter {
protected:
    typedef debruijn_graph::DeBruijnGraph Graph;
public:
    GFAComponentWriter(const omnigraph::GraphComponent<Graph> &component, std::ostream &os,
              io::EdgeNamingF<Graph> naming_f = io::IdNamingF<Graph>())
            : component_(component),
              edge_namer_(component_.g(), naming_f),
              os_(os) {
    }

    void WriteSegmentsAndLinks() {
        WriteSegments();
        WriteLinks();
    }

private:
    void WriteSegments();
    void WriteLinks();

protected:
    const omnigraph::GraphComponent<Graph> &component_;
    io::CanonicalEdgeHelper<Graph> edge_namer_;
    std::ostream &os_;

};

}

