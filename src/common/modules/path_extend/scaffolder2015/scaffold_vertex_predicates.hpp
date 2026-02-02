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

}
}
