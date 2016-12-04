//***************************************************************************
//* Copyright (c) 2016 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "assembly_graph/core/graph.hpp"
#include "assembly_graph/paths/mapping_path.hpp"

extern "C" {
struct bwaidx_s;
typedef struct bwaidx_s bwaidx_t;

struct mem_opt_s;
typedef struct mem_opt_s mem_opt_t;
};

namespace alignment {

class BWAIndex {
  public:
    // bwaidx / memopt are incomplete below, therefore we need to outline ctor
    // and dtor.
    BWAIndex(const debruijn_graph::Graph& g);
    ~BWAIndex();

    omnigraph::MappingPath<debruijn_graph::EdgeId> AlignSequence(const Sequence &sequence) const;
  private:
    void Init();

    const debruijn_graph::Graph& g_;

    // Store the options in memory
    std::unique_ptr<mem_opt_t, void(*)(void*)> memopt_;

    // hold the full index structure
    std::unique_ptr<bwaidx_t, void(*)(bwaidx_t*)> idx_;

    std::vector<debruijn_graph::EdgeId> ids_;
};
    
}
