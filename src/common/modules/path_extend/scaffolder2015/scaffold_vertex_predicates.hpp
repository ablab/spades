//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "auxiliary_graphs/scaffold_graph/scaffold_vertex.hpp"
#include "assembly_graph/core/graph.hpp"

#include <memory>

namespace path_extend {
namespace scaffolder {

class ScaffoldVertexPredicate
  : public func::AbstractPredicate<const scaffold_graph::ScaffoldVertex &> {
  public:
    virtual ~ScaffoldVertexPredicate() = default;
};

class LengthChecker : public ScaffoldVertexPredicate {
  public:
    typedef scaffold_graph::ScaffoldVertex ScaffoldVertex;

    LengthChecker(size_t length_threshold, const Graph &g);

    bool Check(const ScaffoldVertex &vertex) const override;

  private:
    const size_t length_threshold_;
    const Graph &g_;
};

class AndPredicate : public ScaffoldVertexPredicate {
  public:
    typedef scaffold_graph::ScaffoldVertex ScaffoldVertex;

    AndPredicate(std::shared_ptr<ScaffoldVertexPredicate> first, std::shared_ptr<ScaffoldVertexPredicate> second);

    bool Check(const ScaffoldVertex &scaffold_vertex) const override;

  private:
    std::shared_ptr<ScaffoldVertexPredicate> first_;
    std::shared_ptr<ScaffoldVertexPredicate> second_;
};

}
}