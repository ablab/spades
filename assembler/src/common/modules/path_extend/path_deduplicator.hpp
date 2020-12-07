//***************************************************************************
//* Copyright (c) 2020 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "overlap_remover.hpp"
#include "pe_utils.hpp"
#include "assembly_graph/paths/bidirectional_path.hpp"

namespace path_extend {

class PathDeduplicator {
    PathContainer &paths_;
    const bool equal_only_;
    const OverlapFindingHelper helper_;

    bool IsRedundant(const BidirectionalPath &path) const {
        TRACE("Checking if path redundant " << path.GetId());
        for (const BidirectionalPath *candidate : helper_.FindCandidatePaths(path)) {
            TRACE("Considering candidate " << candidate->GetId());
//                VERIFY(candidate != path && candidate != path->GetConjPath());
            if (candidate->GetId() == path.GetId() ||
                candidate->GetId() == path.GetConjPath()->GetId())
                continue;

            if (equal_only_ ? helper_.IsEqual(path, *candidate) : helper_.IsSubpath(path, *candidate))
                return true;
        }
        return false;
    }
public:
    PathDeduplicator(const Graph &g,
                     PathContainer &paths,
                     GraphCoverageMap &coverage_map,
                     size_t min_edge_len,
                     size_t max_diff,
                     bool equal_only) :
            paths_(paths),
            equal_only_(equal_only),
            helper_(g, coverage_map, min_edge_len, max_diff) {}

    //TODO use path container filtering?
    void Deduplicate() {
        for (auto & path_pair : paths_) {
            auto &path = path_pair.first;
            if (IsRedundant(*path)) {
                TRACE("Clearing path " << path->str());
                path->Clear();
            }
        }
    }

private:
    DECL_LOGGER("PathDeduplicator");
};

}
