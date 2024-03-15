//***************************************************************************
//* Copyright (c) 2020 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "overlap_remover.hpp"
#include "path_extender.hpp" // FIXME: Temporary

namespace path_extend {

static void PopFront(BidirectionalPath &path, size_t cnt) {
    path.GetConjPath()->PopBack(cnt);
}

std::pair<Range, Range> OverlapFindingHelper::ComparePaths(const BidirectionalPath &path1,
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

bool OverlapFindingHelper::IsSubpath(const BidirectionalPath &path,
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
bool OverlapFindingHelper::IsEqual(const BidirectionalPath &path,
                                   const BidirectionalPath &other) const {
    auto ends_pair = CommonPrefix(path, other);
    return ends_pair.first == path.Size() && ends_pair.second == other.Size();
}


std::pair<size_t, size_t> OverlapFindingHelper::CommonPrefix(const BidirectionalPath &path1,
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
std::pair<Range, Range> OverlapFindingHelper::FindOverlap(const BidirectionalPath &path1,
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

std::vector<const BidirectionalPath*> OverlapFindingHelper::FindCandidatePaths(const BidirectionalPath &path) const {
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

size_t OverlapRemover::AnalyzeOverlaps(const BidirectionalPath &path, const BidirectionalPath &other,
                                       bool end_start_only, bool retain_one_copy) const {
    VERIFY(!retain_one_copy || !end_start_only);
    auto range_pair = helper_.FindOverlap(path, other, end_start_only);
    size_t overlap = range_pair.first.size();
    auto other_range = range_pair.second;

    if (overlap == 0)
        return 0;

    //checking if region on the other path has not been already added
    //TODO discuss if the logic is needed/correct. It complicates the procedure and prevents trivial parallelism.
    if (retain_one_copy &&
        AlreadyAdded(other, other_range.start_pos, other_range.end_pos) &&
        /*forcing "cut_all" behavior on conjugate paths*/
        other.GetId() != path.GetConjPath()->GetId() &&
        /*certain overkill*/
        other.GetId() != path.GetId()) {
        return 0;
    }

    if (other.GetId() == path.GetId()) {
        if (overlap == path.Size())
            return 0;
        overlap = std::min(overlap, other_range.start_pos);
    }

    if (other.GetId() == path.GetConjPath()->GetId()) {
        overlap = std::min(overlap, other.Size() - other_range.end_pos);
    }

    DEBUG("First " << overlap << " edges of the path will be removed");
    DEBUG(path.str());
    DEBUG("Due to overlap with path");
    DEBUG(other.str());
    DEBUG("Range " << other_range);

    return overlap;
}

void OverlapRemover::MarkStartOverlaps(const BidirectionalPath &path, bool end_start_only, bool retain_one_copy) {
    std::set<size_t> overlap_poss;
    for (const BidirectionalPath *candidate : helper_.FindCandidatePaths(path)) {
        size_t overlap = AnalyzeOverlaps(path, *candidate,
                                         end_start_only, retain_one_copy);
        if (overlap > 0)
            overlap_poss.insert(overlap);
    }

    if (!overlap_poss.empty()) {
        utils::insert_all(splits_[path.GetId()], overlap_poss);
    }
}

void OverlapRemover::InnerMarkOverlaps(bool end_start_only, bool retain_one_copy) {
    for (auto &path_pair : paths_) {
        //TODO think if this "optimization" is necessary
        if (path_pair.first->Size() == 0)
            continue;

        if (path_pair.first->IsCycle()) {
            VERIFY(path_pair.first->GetCycleOverlapping() == path_pair.second->GetCycleOverlapping());
            auto overlapping = path_pair.first->GetCycleOverlapping();
            if (overlapping > 0)
                splits_[path_pair.first->GetId()].insert(overlapping);
        } else {
            MarkStartOverlaps(*path_pair.first, end_start_only, retain_one_copy);
            MarkStartOverlaps(*path_pair.second, end_start_only, retain_one_copy);
        }
    }
}

std::set<size_t> PathSplitter::TransformConjSplits(const BidirectionalPath &p) const {
    std::set<size_t> path_splits;
    size_t path_len = p.Size();
    auto it = splits_.find(p.GetId());
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

std::set<size_t> PathSplitter::GatherAllSplits(const BidirectionalPath &p,
                                               const BidirectionalPath &cp) const {
    VERIFY(p.Size() == cp.Size());
    std::set<size_t> path_splits = TransformConjSplits(cp);
    auto it = splits_.find(p.GetId());
    if (it != splits_.end()) {
        utils::insert_all(path_splits, it->second);
    }
    return path_splits;
}

void PathSplitter::SplitPath(BidirectionalPath * const p, const std::set<size_t> &path_splits) {
    size_t start_pos = 0;
    for (size_t split_pos : path_splits) {
        if (split_pos == 0)
            continue;
        if (split_pos == p->Size())
            break;
        CreatePath(paths_, coverage_map_,
                   p->SubPath(start_pos, split_pos));
        start_pos = split_pos;
    }
    PopFront(*p, start_pos);
}

void PathSplitter::Split() {
    std::vector<std::pair<BidirectionalPath*, BidirectionalPath*>> tmp_paths;
    for (const auto &entry : paths_)
        tmp_paths.emplace_back(entry.first.get(), entry.second.get());
    for (auto & path_pair : tmp_paths) {
        SplitPath(path_pair.first,
                  GatherAllSplits(*path_pair.first, *path_pair.second));
    }
}

}
