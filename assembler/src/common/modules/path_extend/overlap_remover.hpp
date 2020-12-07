//***************************************************************************
//* Copyright (c) 2020 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "assembly_graph/core/graph.hpp"
#include "assembly_graph/paths/bidirectional_path.hpp"
#include "sequence/range.hpp"

namespace path_extend {

class GraphCoverageMap;
class PathContainer;
typedef std::unordered_map<uint64_t, std::set<size_t>> SplitsStorage;

//TODO think about symmetry and what if it breaks?
class OverlapFindingHelper {
    const debruijn_graph::Graph &g_;
    const GraphCoverageMap &coverage_map_;
    const size_t min_edge_len_;
    const size_t max_diff_;
    const bool try_extend_;

    //TODO think of the cases when (gap + length) < 0
    //Changes second argument on success
    void TryExtendToEnd(const BidirectionalPath &path, size_t &pos) const {
        if (pos < path.Size() &&
            path.GapAt(pos).gap + path.LengthAt(pos) <= max_diff_)
            pos = path.Size();
    }

    //Changes second argument on success
    void TryExtendToStart(const BidirectionalPath &path, size_t &pos) const {
        if (pos > 0 && path.Length() - path.LengthAt(pos) <= max_diff_)
            pos = 0;
    }

    std::pair<Range, Range> ComparePaths(const BidirectionalPath &path1,
                                         const BidirectionalPath &path2,
                                         size_t start2) const;
public:
    OverlapFindingHelper(const debruijn_graph::Graph &g,
                         const GraphCoverageMap &coverage_map,
                         size_t min_edge_len,
                         size_t max_diff) :
            g_(g),
            coverage_map_(coverage_map),
            min_edge_len_(min_edge_len),
            max_diff_(max_diff),
            //had to enable try_extend, otherwise equality lost symmetry
            try_extend_(max_diff_ > 0) {
    }

    bool IsSubpath(const BidirectionalPath &path,
                   const BidirectionalPath &other) const;

    //NB! Equality is not transitive if max_diff is > 0
    bool IsEqual(const BidirectionalPath &path,
                 const BidirectionalPath &other) const;

    std::pair<size_t, size_t> CommonPrefix(const BidirectionalPath &path1,
                                           const BidirectionalPath &path2) const;

    //overlap is forced to start from the beginning of path1
    std::pair<Range, Range> FindOverlap(const BidirectionalPath &path1,
                                        const BidirectionalPath &path2,
                                        bool end_start_only) const;

    std::vector<const BidirectionalPath*> FindCandidatePaths(const BidirectionalPath &path) const;
private:
    DECL_LOGGER("OverlapFindingHelper");
};

class OverlapRemover {
    const PathContainer &paths_;
    const OverlapFindingHelper helper_;
    SplitsStorage splits_;

    bool AlreadyAdded(const BidirectionalPath &p, size_t pos) const {
        auto it = splits_.find(p.GetId());
        return it != splits_.end() && it->second.count(pos);
    }

    //TODO if situation start ==0 && end==p.Size is not interesting then code can be simplified
    bool AlreadyAdded(const BidirectionalPath &p, size_t start, size_t end) const {
        if (start == 0 && AlreadyAdded(p, end))
            return true;
        if (end == p.Size() && AlreadyAdded(*p.GetConjPath(), p.Size() - start))
            return true;
        return false;
    }

    //NB! This can only be launched over paths taken from path container!
    size_t AnalyzeOverlaps(const BidirectionalPath &path, const BidirectionalPath &other,
                           bool end_start_only, bool retain_one_copy) const;
    void MarkStartOverlaps(const BidirectionalPath &path, bool end_start_only, bool retain_one_copy);
    void InnerMarkOverlaps(bool end_start_only, bool retain_one_copy);

public:
    OverlapRemover(const debruijn_graph::Graph &g,
                   const PathContainer &paths,
                   GraphCoverageMap &coverage_map,
                   size_t min_edge_len,// = 0,
                   size_t max_diff) // = 0)
            :  paths_(paths),
               helper_(g, coverage_map,
                       min_edge_len, max_diff) {
    }

    //Note that during start/end removal all repeat instance have to be cut
    void MarkOverlaps(bool end_start_only, bool retain_one_copy) {
        VERIFY(!end_start_only || !retain_one_copy);
        INFO("Marking start/end overlaps");
        InnerMarkOverlaps(/*end/start overlaps only*/ true, /*retain one copy*/ false);
        if (!end_start_only) {
            INFO("Marking remaining overlaps");
            InnerMarkOverlaps(/*end/start overlaps only*/ false, retain_one_copy);
        }
    }

    const SplitsStorage& overlaps() const { return splits_; }

private:
    DECL_LOGGER("OverlapRemover");
};

class PathSplitter {
    const SplitsStorage splits_;
    PathContainer &paths_;
    GraphCoverageMap &coverage_map_;

    std::set<size_t> TransformConjSplits(const BidirectionalPath &p) const;
    std::set<size_t> GatherAllSplits(const BidirectionalPath &p,
                                     const BidirectionalPath &cp) const;
    void SplitPath(BidirectionalPath * const p, const std::set<size_t> &path_splits);
public:
    PathSplitter(const SplitsStorage &splits,
                 PathContainer &paths,
                 GraphCoverageMap &coverage_map) :
            splits_(splits),
            paths_(paths),
            coverage_map_(coverage_map) {}

    void Split();

private:
    DECL_LOGGER("PathSplitter");
};


}
