//***************************************************************************
//* Copyright (c) 2021-2023 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "graph_resolver.hpp"

namespace cont_index {
class TransformedGraphIO {
  public:
    explicit TransformedGraphIO(io::IdMapper<string> *id_mapper) : id_mapper_(id_mapper) {}
    void PrintGraph(const debruijn_graph::Graph &graph,
                    const GraphResolver::GraphResolverInfo &resolver_info,
                    const std::filesystem::path &output_base) const;

  private:
    io::IdMapper<std::string> *id_mapper_;
};
}
