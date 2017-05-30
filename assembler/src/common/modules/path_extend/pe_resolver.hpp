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

//FIXME think about overlaps of path with itself and with its conjugate

//FIXME think about symmetry and what if it breaks?
class OverlapFindingHelper {
    const Graph &g_;
    const GraphCoverageMap &coverage_map_;
    const size_t min_edge_len_;
    const size_t max_diff_;
    const bool try_extend_;

    //Changes second argument on success
    bool TryExtendToEnd(const BidirectionalPath &path, size_t &pos) const {
        //FIXME rewrite via LengthAt
        size_t cum = 0;
        for (size_t j = pos; j < path.Size(); ++j) {
            cum += path.ShiftLength(j);
            if (cum > max_diff_) {
                return false;
            }
        }
        pos = path.Size();
        return true;
    }

    //Changes second argument on success
    bool TryExtendToStart(const BidirectionalPath &path, size_t &pos) const {
        size_t tmp_pos = path.Size() - pos;
        if (TryExtendToEnd(*path.GetConjPath(), tmp_pos)) {
            pos = path.Size() - tmp_pos;
            return true;
        }
        return false;
    }


    pair<Range, Range> ComparePaths(const BidirectionalPath &path1,
                                      const BidirectionalPath &path2,
                                      size_t start2) const {
        //TODO change to edit distance?
        //FIXME what to do with gaps here? Currently process stops if gap > max_diff_
        size_t shift1 = 0;
        //path1 is always matched from the start
        const size_t start1 = 0;
        size_t end1 = start1;
        size_t end2 = start2;

        for (size_t i = start1; i < path1.Size(); ++i) {
            shift1 += path1.GapAt(i).gap; //First gap is always zero
            if (shift1 > max_diff_)
                break;

            bool match = false;
            size_t j = end2;
            size_t shift2 = 0;
            for (; j < path2.Size(); ++j) {
                if (end1 == 0) {
                    //Force first match to start with pos2
                    if (j > start2) {
                        break;
                    }
                } else {
                    shift2 += path2.GapAt(j).gap;
                }

                if (shift2 > max_diff_) //shift1 + max_diff)
                    break;
                if (path1.At(i) == path2.At(j) &&
                    (&path1 != &path2 || i != j)) {
                    match = true;
                    break;
                } else {
//                    shift2 += path2.ShiftLength(j);
                    shift2 += g_.length(path2.At(j));
                }
            }
            if (match) {
                end1 = i+1;
                end2 = j+1;
                shift1 = 0;
            } else {
                shift1 += g_.length(path1.At(i));
            }
        }

        //FIXME!!! starts should be extended too!!!
        //Might want to ignore if equal/conjugate paths
        //Extending the ends of the paths if possible
        if (try_extend_ && end1 > 0) {
            TryExtendToEnd(path1, end1);
            TryExtendToEnd(path2, end2);
            TryExtendToStart(path2, start2);
        }
        return make_pair(Range(start1, end1), Range(start2, end2));
    }

public:
    OverlapFindingHelper(const Graph &g,
                         const GraphCoverageMap &coverage_map,
                         size_t min_edge_len,
                         size_t max_diff,
                         bool try_extend = true) :
            g_(g),
            coverage_map_(coverage_map),
            min_edge_len_(min_edge_len),
            max_diff_(max_diff),
            try_extend_(try_extend) {

    }

    bool IsSubpath(const BidirectionalPath &path,
                   const BidirectionalPath &other) const {
        for (size_t j = 0; j < other.Size(); ++j) {
            auto range_pair = ComparePaths(path, other, j);
            if (range_pair.first.end_pos == path.Size()) {
                return true;
            }
        }
        return false;
    }

    bool IsEqual(const BidirectionalPath &path,
                 const BidirectionalPath &other) const {
        auto range_pair = ComparePaths(path, other, 0);
        return range_pair.first.end_pos == path.Size()
               && range_pair.second.end_pos == other.Size();
    }

    //overlap is forced to start from the beginning of path1
    pair<Range, Range> FindOverlap(const BidirectionalPath &path1,
                                   const BidirectionalPath &path2,
                                   bool end_start_only) const {
        size_t max_overlap = 0;
        pair<Range, Range> matching_ranges;
        for (size_t j = 0; j < path2.Size(); ++j) {
            auto range_pair = ComparePaths(path1, path2, j);
            VERIFY(range_pair.first.start_pos == 0);
            //checking if overlap is valid
            if (end_start_only && range_pair.second.end_pos != path2.Size())
                continue;

            size_t overlap_size = range_pair.first.size();
            if (overlap_size > max_overlap ||
                //prefer overlaps with end of path2
                (overlap_size == max_overlap &&
                        range_pair.second.end_pos == path2.Size())) {
                max_overlap = overlap_size;
                matching_ranges = range_pair;
            }
        }
        return matching_ranges;
    }

    vector<PathPtr> FindCandidatePaths(const BidirectionalPath &path) const {
        set<PathPtr> candidates;
        //FIXME needs discussion
        if (min_edge_len_ == 0) {
            size_t cum_len = 0;
            for (size_t i = 0; i < path.Size(); ++i) {
                EdgeId e = path.At(i);
                cum_len += path.GapAt(i).gap;
                if (cum_len > max_diff_)
                    break;
                utils::insert_all(candidates, coverage_map_.GetCoveringPaths(e));
                cum_len += g_.length(e);
            }
            VERIFY(path.Size() == 0 || candidates.count(&path));
        } else {
            for (size_t i = 0; i < path.Size(); ++i) {
                EdgeId e = path.At(i);
                if (g_.length(e) >= min_edge_len_) {
                    utils::insert_all(candidates, coverage_map_.GetCoveringPaths(e));
                }
            }
        }
        return vector<PathPtr>(candidates.begin(), candidates.end());
    }

};

class DecentOverlapRemover {
    const Graph &g_;
    const PathContainer &paths_;
    const GraphCoverageMap &coverage_map_;
    SplitsStorage &splits_;
    const OverlapFindingHelper helper_;

    bool AlreadyAdded(PathPtr ptr, size_t pos) const {
        auto it = splits_.find(ptr);
        return it != splits_.end() && it->second.count(pos);
    }

    //FIXME if situation start ==0 && end==p.Size is not interesting then code can be simplified
    bool AlreadyAdded(const BidirectionalPath &p, size_t start, size_t end) const {
        if (start == 0 && AlreadyAdded(&p, end))
            return true;
        if (end == p.Size() && AlreadyAdded(p.GetConjPath(), p.Size() - start))
            return true;
        return false;
    }

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
                /*forcing cut all behavior on conjugate paths*/
                &other != path.GetConjPath() &&
                /*certain overkill*/
                &other != &path) {
            return 0;
        }

        if (&other == &path) {
            overlap = std::min(overlap, other_range.start_pos);
        }

        if (&other == path.GetConjPath()) {
            overlap = std::min(overlap, other.Size() - other_range.end_pos);
        }

        DEBUG("First " << overlap << " edges of the path will be removed");
        path.PrintDEBUG();
        DEBUG("Due to overlap with path");
        other.PrintDEBUG();

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
        if (!overlap_poss.empty())
            utils::insert_all(splits_[&path], overlap_poss);
    }

public:
    //FIXME First need to remove subpaths
    //FIXME Order path container in order of increasing length for retain_one to work properly!!!
    DecentOverlapRemover(const Graph &g,
                         const PathContainer &paths,
                         GraphCoverageMap &coverage_map,
                         SplitsStorage &splits,
                         size_t min_edge_len,// = 0,
                         size_t max_diff) :// = 0) :
            g_(g),
            paths_(paths),
            coverage_map_(coverage_map),
            splits_(splits),
            helper_(g, coverage_map,
                    min_edge_len, max_diff) {
    }

    void MarkOverlaps(bool end_start_only, bool retain_one_copy) {
        for (auto path_pair: paths_) {
            //FIXME think if this "optimization" is necessary
            if (path_pair.first->Size() == 0)
                continue;
            MarkStartOverlaps(*path_pair.first, end_start_only, retain_one_copy);
            MarkStartOverlaps(*path_pair.second, end_start_only, retain_one_copy);
        }
    }

//    const SplitsStorage& overlaps() const {
//        return splits_;
//    }
private:
    DECL_LOGGER("DecentOverlapRemover");
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

};

class PathDeduplicator {
    const Graph& g_;
    PathContainer &paths_;
    const bool equal_only_;
    const OverlapFindingHelper helper_;

    bool IsRedundant(PathPtr path) const {
        for (auto candidate : helper_.FindCandidatePaths(*path)) {
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
                     size_t max_diff = 0,
                     bool equal_only = false) :
            g_(g),
            paths_(paths),
            equal_only_(equal_only),
            //FIXME discuss min_edge_len
            helper_(g, coverage_map, 0, max_diff) {}

    //FIXME use path container filtering?
    void Deduplicate() {
        for (auto path_pair : paths_) {
            auto path = path_pair.first;
            if (IsRedundant(path)) {
                DEBUG("Clearing path");
                path->PrintDEBUG();
                path->Clear();
            }
        }
    }
};

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

    PathContainer ExtendSeeds(PathContainer &seeds, ContigsMaker &pathExtender) const {
        PathContainer paths;
        pathExtender.GrowAll(seeds, paths);
        return paths;
    }

    void Deduplicate(PathContainer &paths, GraphCoverageMap &coverage_map,
                     size_t /*min_edge_len*/, size_t max_path_diff,
                     bool equal_only = false) const {
        PathDeduplicator deduplicator(g_, paths, coverage_map, max_path_diff, equal_only);
        deduplicator.Deduplicate();
    }

    void RemoveRNAOverlaps(PathContainer& /*paths*/, GraphCoverageMap& /*coverage_map*/,
                          size_t /*min_edge_len*/, size_t /*max_path_diff*/) const {

        VERIFY(false);
    }

    void RemoveOverlaps(PathContainer &paths, GraphCoverageMap &coverage_map,
                        size_t min_edge_len, size_t max_path_diff,
                        bool cut_all) const {
        INFO("Removing overlaps");
        VERIFY(min_edge_len == 0 && max_path_diff == 0);

        INFO("Deduplicating paths");
        Deduplicate(paths, coverage_map, min_edge_len, max_path_diff);

        INFO("Deduplicated");

        SplitsStorage splits;

        DecentOverlapRemover overlap_remover(g_, paths, coverage_map, splits,
                                             min_edge_len, max_path_diff);
        //FIXME can't we simplify logic?
        //DO NOT CHANGE ORDER!
        INFO("Marking start/end overlaps");
        overlap_remover.MarkOverlaps(/*retain one copy*/ false,/*end/start overlaps only*/ true);
        INFO("Marking remaining overlaps");
        overlap_remover.MarkOverlaps(/*retain one copy*/ !cut_all,/*end/start overlaps only*/ false);

        INFO("Splitting paths");
        PathSplitter splitter(splits, paths, coverage_map);
        splitter.Split();
        //splits are invalidated after this point

        INFO("Deduplicating paths");
        Deduplicate(paths, coverage_map, min_edge_len, max_path_diff);
        INFO("Overlaps removed");
        //FIXME Add removal of empty paths here
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
