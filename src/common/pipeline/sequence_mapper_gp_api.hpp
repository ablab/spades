//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2015-2022 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "graph_pack.hpp"

#include "alignment/sequence_mapper.hpp"

namespace debruijn_graph {
std::shared_ptr<BasicSequenceMapper<Graph, EdgeIndex<Graph>>> MapperInstance(const graph_pack::GraphPack &gp);

std::shared_ptr<BasicSequenceMapper<Graph, EdgeIndex<Graph>>> MapperInstance(const graph_pack::GraphPack &gp,
                                                                             const EdgeIndex<Graph> &index);
}
