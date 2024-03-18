//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2021-2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "kmer_sequence_mapper.hpp"

#include "assembly_graph/core/graph.hpp"
#include "edge_index.hpp"

#include <parallel_hashmap/phmap.h>
#include <string>
#include <utility>
#include <vector>

namespace alignment {

void ShortKMerReadMapper::Init(const std::string &workdir) {
    index_.reset(new debruijn_graph::EdgeIndex<debruijn_graph::Graph>(g_, workdir));
    index_->Refill(k_);
}

omnigraph::MappingPath<debruijn_graph::EdgeId> ShortKMerReadMapper::MapSequence(const Sequence &sequence,
                                                                                bool) const {
    if (sequence.size() < k_)
        return { };

    // Extract all k-mers from Sequence and count how many times each edge appears.
    // The top one (but with at least min_occ) hits wins
    phmap::flat_hash_map<EdgeId, unsigned> occ;
    EdgeIndex::KMer kmer = sequence.start<RtSeq>(k_) >> 'A';
    for (size_t j = k_ - 1; j < sequence.size(); ++j) {
        uint8_t inchar = sequence[j];
        kmer <<= inchar;

        auto pos = index_->get(kmer);
        if (pos.second == EdgeIndex::NOT_FOUND)
            continue;

        occ[pos.first] += 1;
    }

    if (occ.empty())
        return { };
        
    std::vector<std::pair<EdgeId, unsigned>> occv(occ.begin(), occ.end());
    std::sort(occv.begin(), occv.end(),
              [](const auto &lhs, const auto &rhs) {
                  if (lhs.second == rhs.second)
                      return lhs.first < rhs.first;

                  return lhs.second > rhs.second;
              });

    if (occv.front().second < min_occ_)
        return { };

    return { occv.front().first, {} };
}


}
