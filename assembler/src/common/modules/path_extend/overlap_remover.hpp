//***************************************************************************
//* Copyright (c) 2020 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "pe_utils.hpp"
#include "path_extender.hpp" // FIXME: Temporary
#include "assembly_graph/core/graph.hpp"
#include "assembly_graph/paths/bidirectional_path.hpp"
#include "sequence/range.hpp"

namespace path_extend {

typedef const BidirectionalPath * PathPtr;
typedef std::unordered_map<PathPtr, std::set<size_t>> SplitsStorage;

inline void PopFront(BidirectionalPath * const path, size_t cnt) {
    path->GetConjPath()->PopBack(cnt);
}

//TODO think about symmetry and what if it breaks?
class OverlapFindingHelper {
    const Graph &g_;
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
                                         size_t start2) const {
        TRACE("Comparing paths " << path1.GetId() << " and " << path2.GetId());
        //TODO change to edit distance?
        int shift1 = 0;
        //path1 is always matched from the start
        const size_t start1 = 0;
        size_t end1 = start1;
        size_t end2 = start2;

        for (size_t i = start1; i < path1.Size(); ++i) {
            if (abs(shift1) > int(max_diff_))
                break;

            bool match = false;
            size_t j = end2;
            int shift2 = 0;
            for (; j < path2.Size(); ++j) {
                if (end1 == 0) {
                    //Force first match to start with pos2
                    if (j > start2) {
                        break;
                    }
                }

                if (abs(shift2) > int(max_diff_))
                    break;
                if (path1.At(i) == path2.At(j) &&
                        (end1 == 0 ||
                            abs(shift1 + path1.GapAt(i).gap - shift2 - path2.GapAt(j).gap) <= int(max_diff_))) {
                    match = true;
                    break;
                } else {
                    shift2 += path2.ShiftLength(j);
                }
            }
            if (match) {
                end1 = i+1;
                end2 = j+1;
                shift1 = 0;
            } else {
                shift1 += path1.ShiftLength(i);
            }
        }

        //Extending the ends of the paths if possible
        if (try_extend_ && end1 > 0) {
            TryExtendToEnd(path1, end1);
            TryExtendToEnd(path2, end2);
            //no need to extend path1 left
            VERIFY(start1 == 0);
            TryExtendToStart(path2, start2);
        }

        return {Range(start1, end1), Range(start2, end2)};
    }

public:
    OverlapFindingHelper(const Graph &g,
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
                   const BidirectionalPath &other) const {
        for (size_t j = 0; j < other.Size(); ++j) {
            auto range_pair = ComparePaths(path, other, j);
            if (range_pair.first.end_pos == path.Size()) {
                return true;
            }
        }
        return false;
    }

    //NB! Equality is not transitive if max_diff is > 0
    bool IsEqual(const BidirectionalPath &path,
                 const BidirectionalPath &other) const {
        auto ends_pair = CommonPrefix(path, other);
        return ends_pair.first == path.Size()
               && ends_pair.second == other.Size();
    }


    std::pair<size_t, size_t> CommonPrefix(const BidirectionalPath &path1,
                                           const BidirectionalPath &path2) const {
        std::pair<size_t, size_t> answer(0, 0);
        size_t cum = 0;
        size_t max_overlap = 0;
        for (size_t j = 0; j < path2.Size(); ++j) {
            auto range_pair = ComparePaths(path1, path2, j);
            if (range_pair.second.start_pos == 0 && range_pair.first.size() > max_overlap) {
                answer = {range_pair.first.end_pos, range_pair.second.end_pos};
                max_overlap = range_pair.first.size();
            }

            if (!try_extend_)
                break;

            cum += path2.ShiftLength(j);
            if (cum > max_diff_)
                break;
        }
        return answer;
    };

    //overlap is forced to start from the beginning of path1
    std::pair<Range, Range> FindOverlap(const BidirectionalPath &path1,
                                        const BidirectionalPath &path2,
                                        bool end_start_only) const {
        size_t max_overlap = 0;
        std::pair<Range, Range> matching_ranges;
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

    std::vector<const BidirectionalPath*> FindCandidatePaths(const BidirectionalPath &path) const {
        std::set<const BidirectionalPath*> candidates;
        size_t cum_len = 0;
        for (size_t i = 0; i < path.Size(); ++i) {
            if (cum_len > max_diff_)
                break;
            EdgeId e = path.At(i);
            if (g_.length(e) >= min_edge_len_) {
                utils::insert_all(candidates, coverage_map_.GetCoveringPaths(e));
                cum_len += path.ShiftLength(i);
            }
        }
        return std::vector<const BidirectionalPath*>(candidates.begin(), candidates.end());
    }

private:
    DECL_LOGGER("OverlapFindingHelper");
};

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
        std::set<size_t> overlap_poss;
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
        for (auto &path_pair : paths_) {
            //TODO think if this "optimization" is necessary
            if (path_pair.first->Size() == 0)
                continue;
            
            if (path_pair.first->IsCycle()) {
                VERIFY(path_pair.first->GetCycleOverlapping() == path_pair.second->GetCycleOverlapping());
                auto overlapping = path_pair.first->GetCycleOverlapping();
                if (overlapping > 0)
                    splits_[path_pair.first.get()].insert(overlapping);
            } else {
                MarkStartOverlaps(*path_pair.first, end_start_only, retain_one_copy);
                MarkStartOverlaps(*path_pair.second, end_start_only, retain_one_copy);
            }
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

    std::set<size_t> TransformConjSplits(const BidirectionalPath *p) const {
        std::set<size_t> path_splits;
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

    std::set<size_t> GatherAllSplits(const BidirectionalPath *p,
                                     const BidirectionalPath *cp) const {
        VERIFY(p->Size() == cp->Size());
        std::set<size_t> path_splits = TransformConjSplits(cp);
        auto it = splits_.find(p);
        if (it != splits_.end()) {
            utils::insert_all(path_splits, it->second);
        }
        return path_splits;
    }

    void SplitPath(BidirectionalPath * const p, const std::set<size_t> &path_splits) {
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
         std::vector<std::pair<BidirectionalPath*, BidirectionalPath*>> tmp_paths;
         for (const auto &entry : paths_)
             tmp_paths.emplace_back(entry.first.get(), entry.second.get());
         for (auto & path_pair : tmp_paths) {
             SplitPath(path_pair.first, GatherAllSplits(path_pair.first,
                                                        path_pair.second));
         }
     }

private:
    DECL_LOGGER("PathSplitter");
};

class PathDeduplicator {
    PathContainer &paths_;
    const bool equal_only_;
    const OverlapFindingHelper helper_;

    bool IsRedundant(BidirectionalPath *path) const {
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
        for (auto & path_pair : paths_) {
            auto &path = path_pair.first;
            if (IsRedundant(path.get())) {
                TRACE("Clearing path " << path->str());
                path->Clear();
            }
        }
    }

private:
    DECL_LOGGER("PathDeduplicator");
};

}
