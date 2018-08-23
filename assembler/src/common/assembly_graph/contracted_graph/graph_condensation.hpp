#pragma once
#include "contracted_graph.hpp"

namespace contracted_graph {

class UnbranchingPathExtractor {
    typedef path_extend::scaffold_graph::ScaffoldVertex ScaffoldVertex;
    typedef vector<ScaffoldVertex> SimplePath;
 public:
    vector<SimplePath> ExtractUnbranchingPaths(const ContractedGraph &graph) const;
};

}