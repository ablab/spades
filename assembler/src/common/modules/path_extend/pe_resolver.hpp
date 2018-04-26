//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

/*
 * pe_resolver.hpp
 *
 *  Created on: Mar 12, 2012
 *      Author: andrey
 */

#pragma once

#include "path_extender.hpp"

namespace path_extend {

typedef const BidirectionalPath * PathPtr;
typedef unordered_map<PathPtr, set<size_t>> SplitsStorage;

inline void PopFront(BidirectionalPath * const path, size_t cnt) {
    path->GetConjPath()->PopBack(cnt);
}

class OverlapRemover {
    const PathContainer &paths_;
    const OverlapFindingHelper helper_;
    SplitsStorage splits_;

    bool AlreadyAdded(PathPtr ptr, size_t pos) const {
        auto it = splits_.find(ptr);
        return it != splits_.end() && it->second.count(pos);
    }

    //TODO if situation start ==0 && end==p.Size is not interesting then code can be simplified
    bool AlreadyAdded(const BidirectionalPath &p, size_t start, size_t end) const {
        if (start == 0 && AlreadyAdded(&p, end))
            return true;
        if (end == p.Size() && AlreadyAdded(p.GetConjPath(), p.Size() - start))
            return true;
        return false;
    }

    //NB! This can only be launched over paths taken from path container!
    size_t AnalyzeOverlaps(const BidirectionalPath &path, const BidirectionalPath &other,
                           bool end_start_only, bool retain_one_copy) const {
        VERIFY(!retain_one_copy || !end_start_only);
        auto range_pair = helper_.FindOverlap(path, other, end_start_only);
        size_t overlap = range_pair.first.size();
        auto other_range = range_pair.second;

        if (overlap == 0) {
            return 0;
        }

        //checking if region on the other path has not been already added
        //TODO discuss if the logic is needed/correct. It complicates the procedure and prevents trivial parallelism.
        if (retain_one_copy &&
                AlreadyAdded(other,
                             other_range.start_pos,
                             other_range.end_pos) &&
                /*forcing "cut_all" behavior on conjugate paths*/
                &other != path.GetConjPath() &&
                /*certain overkill*/
                &other != &path) {
            return 0;
        }

        if (&other == &path) {
            if (overlap == path.Size())
                return 0;
            overlap = std::min(overlap, other_range.start_pos);
        }

        if (&other == path.GetConjPath()) {
            overlap = std::min(overlap, other.Size() - other_range.end_pos);
        }

        DEBUG("First " << overlap << " edges of the path will be removed");
        DEBUG(path.str());
        DEBUG("Due to overlap with path");
        DEBUG(other.str());
        DEBUG("Range " << other_range);

        return overlap;
    }

    void MarkStartOverlaps(const BidirectionalPath &path, bool end_start_only, bool retain_one_copy) {
        set<size_t> overlap_poss;
        for (PathPtr candidate : helper_.FindCandidatePaths(path)) {
            size_t overlap = AnalyzeOverlaps(path, *candidate,
                                             end_start_only, retain_one_copy);
            if (overlap > 0) {
                overlap_poss.insert(overlap);
            }
        }
        if (!overlap_poss.empty()) {
            utils::insert_all(splits_[&path], overlap_poss);
        }
    }

    void InnerMarkOverlaps(bool end_start_only, bool retain_one_copy) {
        for (auto path_pair: paths_) {
            //TODO think if this "optimization" is necessary
            if (path_pair.first->Size() == 0)
                continue;
            MarkStartOverlaps(*path_pair.first, end_start_only, retain_one_copy);
            MarkStartOverlaps(*path_pair.second, end_start_only, retain_one_copy);
        }
    }

public:
    OverlapRemover(const Graph &g,
                         const PathContainer &paths,
                         GraphCoverageMap &coverage_map,
                         size_t min_edge_len,// = 0,
                         size_t max_diff) :// = 0) :
            paths_(paths),
            helper_(g, coverage_map,
                    min_edge_len, max_diff) {
    }

    //Note that during start/end removal all repeat instance have to be cut
//    void MarkOverlaps(bool end_start_only = false, bool retain_one_copy = true) {
    void MarkOverlaps(bool end_start_only, bool retain_one_copy) {
        VERIFY(!end_start_only || !retain_one_copy);
        INFO("Marking start/end overlaps");
        InnerMarkOverlaps(/*end/start overlaps only*/ true, /*retain one copy*/ false);
        if (!end_start_only) {
            INFO("Marking remaining overlaps");
            InnerMarkOverlaps(/*end/start overlaps only*/ false, retain_one_copy);
        }
    }

    const SplitsStorage& overlaps() const {
        return splits_;
    }

private:
    DECL_LOGGER("OverlapRemover");
};

class PathSplitter {
    const SplitsStorage splits_;
    PathContainer &paths_;
    GraphCoverageMap &coverage_map_;

    set<size_t> TransformConjSplits(PathPtr p) const {
        set<size_t> path_splits;
        size_t path_len = p->Size();
        auto it = splits_.find(p);
        if (it != splits_.end()) {
//                std::transform(it->second.begin(), it->second.end(),
//                               std::inserter(path_splits, path_splits.end()),
//                               [=] (size_t pos) {return path_len - pos;});
            for (size_t pos : it->second) {
                path_splits.insert(path_len - pos);
            }
        }
        return path_splits;
    }

    set<size_t> GatherAllSplits(const PathPair &pp) const {
        VERIFY(pp.first->Size() == pp.second->Size());
        set<size_t> path_splits = TransformConjSplits(pp.second);
        auto it = splits_.find(pp.first);
        if (it != splits_.end()) {
            utils::insert_all(path_splits, it->second);
        }
        return path_splits;
    }

    void SplitPath(BidirectionalPath * const p, const set<size_t> &path_splits) {
        size_t start_pos = 0;
        for (size_t split_pos : path_splits) {
            if (split_pos == 0)
                continue;
            if (split_pos == p->Size())
                break;
            AddPath(paths_, p->SubPath(start_pos, split_pos), coverage_map_);
            start_pos = split_pos;
        }
        PopFront(p, start_pos);
    }

public:
    PathSplitter(const SplitsStorage &splits,
                 PathContainer &paths,
                 GraphCoverageMap &coverage_map) :
            splits_(splits),
            paths_(paths),
            coverage_map_(coverage_map) {}

     void Split() {
         vector<PathPair> tmp_paths(paths_.begin(), paths_.end());
         for (auto path_pair: tmp_paths) {
             SplitPath(path_pair.first, GatherAllSplits(path_pair));
         }
     }

private:
    DECL_LOGGER("PathSplitter");
};

class PathDeduplicator {
    PathContainer &paths_;
    const bool equal_only_;
    const OverlapFindingHelper helper_;

    bool IsRedundant(PathPtr path) const {
        TRACE("Checking if path redundant " << path->GetId());
        for (auto candidate : helper_.FindCandidatePaths(*path)) {
            TRACE("Considering candidate " << candidate->GetId());
//                VERIFY(candidate != path && candidate != path->GetConjPath());
            if (candidate == path || candidate == path->GetConjPath())
                continue;
            if (equal_only_ ? helper_.IsEqual(*path, *candidate) : helper_.IsSubpath(*path, *candidate)) {
                return true;
            }
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
        for (auto path_pair : paths_) {
            auto path = path_pair.first;
            if (IsRedundant(path)) {
                TRACE("Clearing path " << path->str());
                path->Clear();
            }
        }
    }

private:
    DECL_LOGGER("PathDeduplicator");
};

inline void Deduplicate(const Graph &g, PathContainer &paths, GraphCoverageMap &coverage_map,
                 size_t min_edge_len, size_t max_path_diff,
                 bool equal_only = false) {
    //add sorting to guarantee survival of longest paths if max_path_diff used
    //paths.SortByLength(false);
    PathDeduplicator deduplicator(g, paths, coverage_map, min_edge_len, max_path_diff, equal_only);
    deduplicator.Deduplicate();
    paths.FilterEmptyPaths();
}

class PathExtendResolver {

    const Graph& g_;
    size_t k_;

public:
    PathExtendResolver(const Graph& g): g_(g), k_(g.k()) {
    }

    PathContainer MakeSimpleSeeds() const {
        std::set<EdgeId> included;
        PathContainer edges;
        for (auto iter = g_.ConstEdgeBegin(/*canonical only*/true); !iter.IsEnd(); ++iter) {
            EdgeId e = *iter;
            if (g_.int_id(e) <= 0 || InTwoEdgeCycle(e, g_))
                continue;
            edges.AddPair(new BidirectionalPath(g_, e), new BidirectionalPath(g_, g_.conjugate(e)));
        }
        return edges;
    }

    PathContainer ExtendSeeds(PathContainer &seeds, CompositeExtender &composite_extender) const {
        PathContainer paths;
        composite_extender.GrowAll(seeds, paths);
        return paths;
    }

    //Paths should be deduplicated first!
    void RemoveOverlaps(PathContainer &paths, GraphCoverageMap &coverage_map,
                        size_t min_edge_len, size_t max_path_diff,
                        bool end_start_only, bool cut_all) const {
        INFO("Removing overlaps");
        //VERIFY(min_edge_len == 0 && max_path_diff == 0);
        if (!cut_all) {
            INFO("Sorting paths");
            //sorting is currently needed to retain overlap instance in longest paths
            paths.SortByLength(false);
        }

        OverlapRemover overlap_remover(g_, paths, coverage_map,
                                             min_edge_len, max_path_diff);
        INFO("Marking overlaps");
        overlap_remover.MarkOverlaps(end_start_only, !cut_all);

        INFO("Splitting paths");
        PathSplitter splitter(overlap_remover.overlaps(), paths, coverage_map);
        splitter.Split();
        //splits are invalidated after this point

        INFO("Deduplicating paths");
        Deduplicate(g_, paths, coverage_map, min_edge_len, max_path_diff);
        INFO("Overlaps removed");
    }

    void AddUncoveredEdges(PathContainer &paths, GraphCoverageMap &coverageMap) const {
        for (auto iter = g_.ConstEdgeBegin(true); !iter.IsEnd(); ++iter) {
            EdgeId e = *iter;
            if (!coverageMap.IsCovered(e)) {
                AddPath(paths, BidirectionalPath(g_, e), coverageMap);
            }
        }
    }

protected:
    DECL_LOGGER("PEResolver")
};

}
