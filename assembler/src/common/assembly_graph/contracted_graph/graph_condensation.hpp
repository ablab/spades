//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "common/assembly_graph/contracted_graph/contracted_graph.hpp"

namespace contracted_graph {

class UnbranchingPathExtractor {
  public:
    typedef path_extend::scaffold_graph::ScaffoldVertex ScaffoldVertex;
    typedef std::vector<ScaffoldVertex> SimplePath;
    std::vector<SimplePath> ExtractUnbranchingPaths(const ContractedGraph &graph) const;
};

}