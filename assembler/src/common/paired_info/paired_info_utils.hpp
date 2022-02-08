//***************************************************************************
//* Copyright (c) 2021 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "paired_info.hpp"

#include "alignment/sequence_mapper_fwd.hpp"
#include "assembly_graph/core/graph.hpp"
#include "library/library_fwd.hpp"
#include "library/library_data.hpp"
#include "adt/bf.hpp"

namespace paired_info {

using SequencingLib = io::SequencingLibrary<debruijn_graph::config::LibraryData>;
using PairedInfoFilter = bf::counting_bloom_filter<std::pair<debruijn_graph::Graph::EdgeId,
                                                             debruijn_graph::Graph::EdgeId>, 2>;
using PairedIndex = omnigraph::de::UnclusteredPairedInfoIndexT<debruijn_graph::Graph>;

bool CollectLibInformation(const debruijn_graph::Graph &gp,
                           const debruijn_graph::SequenceMapper<debruijn_graph::Graph> &mapper,
                           size_t &edgepairs, SequencingLib &reads,
                           size_t edge_length_threshold);

void FillPairedIndex(const debruijn_graph::Graph &gp,
                     const debruijn_graph::SequenceMapper<debruijn_graph::Graph> &mapper,
                     SequencingLib &reads,
                     PairedIndex &index,
                     std::unique_ptr<PairedInfoFilter> filter, unsigned filter_threshold,
                     unsigned round_thr = 0, bool use_binary = true);

std::unique_ptr<PairedInfoFilter> FillEdgePairFilter(const debruijn_graph::Graph &gp,
                                                     const debruijn_graph::SequenceMapper<debruijn_graph::Graph> &mapper,
                                                     SequencingLib &reads,
                                                     size_t edgepairs);
}

