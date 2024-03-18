//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2015-2022 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "assembly_graph/core/action_handlers.hpp"
#include "assembly_graph/core/debruijn_data.hpp"
#include "utils/verify.hpp"

#include <unordered_set>

namespace omnigraph {

template<class Graph>
class OverlapHandler : public omnigraph::GraphActionHandler<Graph> {

 public:
  typedef typename Graph::VertexId VertexId;
  typedef typename Graph::EdgeId EdgeId;
  typedef typename Graph::LinkId LinkId;

  OverlapHandler(Graph &g)
      : base(g, "OverlapHandler"), g_(g) {}

  virtual void HandleDelete(VertexId v) {
      if (g_.is_complex(v)) {
          g_.clear_links(v);
      }
  }
  virtual void HandleDelete(EdgeId e) {
      if (g_.is_complex(g_.EdgeEnd(e))) {
          g_.erase_links_with_inedge(g_.EdgeEnd(e), e);
      }
      if (g_.is_complex(g_.EdgeStart(e))) {
          g_.erase_links_with_outedge(g_.EdgeStart(e), e);
      }
  }

 private:
  typedef omnigraph::GraphActionHandler<Graph> base;

  Graph &g_;
};
}