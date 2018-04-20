//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "assembly_graph/core/graph.hpp"
#include "common/io/utils/edge_namer.hpp"

#include <memory>
#include <string>
#include <ostream>

namespace gfa {

class GFAWriter {
  protected:
    typedef debruijn_graph::DeBruijnGraph Graph;

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

  private:
    void WriteSegments();
    void WriteLinks();

  protected:
    const Graph &graph_;
    io::CanonicalEdgeHelper<Graph> edge_namer_;
    std::ostream &os_;

};

}

