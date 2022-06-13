//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "transition_extractor.hpp"

namespace path_extend {
namespace read_cloud {
namespace validation {

class SimpleTransitionGraphValidator {
  public:
    typedef scaffold_graph::ScaffoldVertex ScaffoldVertex;
    typedef SimpleGraph<ScaffoldVertex> SimpleTransitionGraph;

    SimpleTransitionGraphValidator(const ContigTransitionStorage &reference_transition_storage,
                                   size_t length_threshold);
    boost::optional<std::vector<ScaffoldVertex>> GetCorrectPath(const SimpleTransitionGraph &graph,
                                                                ScaffoldVertex source, ScaffoldVertex sink) const;
  private:
    const ContigTransitionStorage reference_transition_storage_;
    const size_t length_threshold_;

    DECL_LOGGER("SimpleTransitionGraphValidator");
};

class SimpleTransitionGraphValidatorConstructor {
  public:
    SimpleTransitionGraphValidatorConstructor(const graph_pack::GraphPack &gp,
                                              size_t length_threshold);
    SimpleTransitionGraphValidator GetValidator(const string &path_to_reference) const;

  private:
    const graph_pack::GraphPack &gp_;
    const Graph &g_;
    const debruijn_graph::EdgeIndex<Graph> &index_;
    const debruijn_graph::KmerMapper<Graph> &kmer_mapper_;
    const size_t length_threshold_;
};
}
}
}