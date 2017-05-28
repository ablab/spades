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

#ifndef PE_RESOLVER_HPP_
#define PE_RESOLVER_HPP_

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
        for (auto path_pair: paths_) {
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

//class OverlapRemover {
//    const Graph& g_;
//    GraphCoverageMap& coverage_map_;
//
//    bool IsSamePath(BidirectionalPath * path1,
//                    BidirectionalPath * path2) const {
//        return *path2 == *path1 || *path2 == *path1->GetConjPath();
//    }
//
//    void RemoveOverlap(PathContainer& paths, BidirectionalPath* path1,
//                       BidirectionalPath* path2, size_t overlap_size) const {
//        DEBUG("Path 1");
//        path1->PrintDEBUG();
//        DEBUG("Path 2");
//        path2->PrintDEBUG();
//        BidirectionalPath* conj2 = path2->GetConjPath();
//        if (path1->IsOverlap() && overlap_size == path1->Size()) {
//            DEBUG("Detaching overlap from path 2 " << path2->GetConjPath()->GetId()
//                                                   << " because of path 1 " << path1->GetId());
//            conj2->PopBack(overlap_size);
//        } else if (path2->IsOverlap() && path2->Size() == overlap_size) {
//            DEBUG("Detaching overlap from path 1 " << path1->GetId() << " because of path 2" << path2->GetId());
//            path1->PopBack(overlap_size);
//        } else if (overlap_size < path2->Size()
//                   && overlap_size < path1->Size()) {
//            AddPath(paths, path1->SubPath(path1->Size() - overlap_size), coverage_map_)->SetOverlap(true);
//            DEBUG("Detaching overlap " << path1->GetId() << " and " << conj2->GetId());
//            path1->PopBack(overlap_size);
//            conj2->PopBack(overlap_size);
//        }
//    }
//
//    void FindAndRemovePathOverlap(PathContainer& all_paths,
//                                  BidirectionalPath* path) const {
//        //FIXME why case of single edge is skipped?
//        if (path->Size() <= 1 ||
//            coverage_map_.GetCoverage(path->Back()) <= 1) {
//            return;
//        }
//        size_t max_overlap_size = 0;
//        BidirectionalPath* overlap_path = nullptr;
//        for (auto p : coverage_map_.GetCoveringPaths(path->Back())) {
//            if (IsSamePath(p, path)) {
//                continue;
//            }
//            size_t overlap_size = OverlapSize(*path, *p);
//            if (overlap_size > max_overlap_size) {
//                max_overlap_size = overlap_size;
//                overlap_path = p;
//            } else if (overlap_size > 0
//                       && overlap_size == max_overlap_size
//                       && p->GetId() < overlap_path->GetId()) {
//                //FIXME check if this logic is needed
//                overlap_path = p;
//            }
//        }
//        if (max_overlap_size > 0)
//            RemoveOverlap(all_paths, path, overlap_path, max_overlap_size);
//    }
//
//public:
//    OverlapRemover(const Graph& g, GraphCoverageMap& cm)
//            : g_(g), coverage_map_(cm) {
//    }
//
//    void RemoveOverlaps(PathContainer& paths) const {
//        for (size_t i = 0; i < paths.size(); i++) {
//            FindAndRemovePathOverlap(paths, paths.Get(i));
//            FindAndRemovePathOverlap(paths, paths.GetConjugate(i));
//        }
//    }
//
//};
//

//inline size_t CommonForward(const BidirectionalPath &path, size_t pos1, size_t pos2) {
//    VERIFY(pos1 <= pos2);
//    if (pos1 == pos2)
//        return 0;
//    size_t i = 0;
//    while (pos2 + i < path.Size() && path.At(pos1 + i) == path.At(pos2 + i)) {
//        ++i;
//    }
//    return i;
//}

//inline void MarkMaximumNonUniquePrefix(const BidirectionalPath &path, SplitsStorage &storage) {
//    if (path.Size() == 0) {
//        return;
//    }
//
//    size_t answer = 0;
//    for (size_t pos : path.FindAll(path.Front())) {
//        answer = std::max(answer, CommonForward(path, 0, pos));
//    }
//    storage[&path].insert(answer);
//}
//
////FIXME why this has been done via the tmp vector?
//inline void MarkNonUniquePrefixes(const PathContainer &paths, SplitsStorage &storage) {
//    for (const auto &path_pair : paths) {
//        MarkMaximumNonUniquePrefix(*path_pair.first, storage);
//        MarkMaximumNonUniquePrefix(*path_pair.first, storage);
//    }
//}

////FIXME seems to work out of the box now, but need some thought
//inline void MarkPseudoSelfConjugatePaths(PathContainer &paths, SplitsStorage &storage) {
//    for (const auto &path_pair : paths) {
//        BidirectionalPath * path1 = path_pair.first;
//        BidirectionalPath * path2 = path_pair.second;
//        if (path1 != path2) {
//            size_t i = 0;
//            while (i < path1->Size() && path1->At(i) == path2->At(i)) {
//                i++;
//            }
//            if (i > 0) {
//                storage[path1].insert(i);
//                storage[path2].insert(i);
//            }
//        }
//    }
//}

//FIXME remove on cleanup
//class SimpleOverlapRemover {
//
//public:
//    SimpleOverlapRemover(const Graph& g,
//                         GraphCoverageMap& cm,
//                         size_t min_edge_len,
//                         size_t max_path_diff)
//            : g_(g), coverage_map_(cm),
//              min_edge_len_(min_edge_len),
//              max_path_diff_(max_path_diff) {
//    }
//
//    void RemoveSimilarPaths(PathContainer& paths,
//                            bool del_only_equal, bool del_subpaths, bool del_begins,
//                            bool del_all, bool add_overlap_begins) const {
//        DEBUG("== Removing similar paths ==");
//        DEBUG("Min edge len " << min_edge_len_ << ", max path diff " << max_path_diff_);
//        DEBUG("Only equal " << del_only_equal << ", subpaths " << del_subpaths
//                            << ", starts " << del_begins << ", all " << del_all
//                            << ", add starts " << add_overlap_begins);
//        //FIXME loop over edges looks a really bad idea
//        for (EdgeId edge : GetSortedEdges()) {
//            BidirectionalPathSet cov_paths = coverage_map_.GetCoveringPaths(edge);
//            std::vector<BidirectionalPath*> cov_vect(cov_paths.begin(), cov_paths.end());
//            std::sort(cov_vect.begin(), cov_vect.end(), [] (BidirectionalPath* p1, BidirectionalPath* p2) {
//                return p1->GetId() < p2->GetId();
//            });
//            for (size_t vect_i = 0; vect_i < cov_vect.size(); ++vect_i) {
//                BidirectionalPath* path1 = cov_vect.at(vect_i);
//                if (cov_paths.find(path1) == cov_paths.end()) {
//                    continue;
//                }
//                VERIFY(cov_paths.find(path1) != cov_paths.end());
//                for (size_t vect_i1 = vect_i + 1; vect_i1 < cov_vect.size(); ++vect_i1) {
//                    BidirectionalPath* path2 = cov_vect.at(vect_i1);
//                    if (path1 == path2 || path1 == path2->GetConjPath()) {
//                        continue;
//                    }
//                    if (cov_paths.find(path2) == cov_paths.end())
//                        continue;
//                    VERIFY(cov_paths.find(path2) != cov_paths.end());
//                    if ((*path1) == (*path2)) {
//                        if (path2->IsOverlap()) {
//                            path1->SetOverlap(true);
//                        }
//                        DEBUG("Removing path " << path2->GetId() << " because of path " << path1->GetId());
//                        path2->PrintDEBUG();
//                        path1->PrintDEBUG();
//                        path2->Clear();
//                        cov_paths = coverage_map_.GetCoveringPaths(edge);
//                        continue;
//                    }
//                    if (g_.length(edge) <= min_edge_len_ || path1->IsOverlap() || path2->IsOverlap() || del_only_equal) {
//                        continue;
//                    }
//                    CompareAndCut(paths, edge, path1, path2,
//                                  del_subpaths, del_begins, del_all, add_overlap_begins);
//                    cov_paths = coverage_map_.GetCoveringPaths(edge);
//                }
//            }
//        }
//        DEBUG("== Emd removing similar paths ==");
//    }
//
//private:
//
//    void CompareAndCut(PathContainer& paths, EdgeId edge, BidirectionalPath* path1, BidirectionalPath* path2,
//                       bool del_subpaths, bool del_begins,
//                       bool del_all, bool add_overlap_begins) const {
//        vector<size_t> positions1 = path1->FindAll(edge);
//        vector<size_t> positions2 = path2->FindAll(edge);
//        size_t i1 = 0;
//        size_t i2 = 0;
//        bool renewed = false;
//        while (i1 < positions1.size()) {
//            while (i2 < positions2.size()) {
//                DEBUG("CompareAndCutFromPos paths " << g_.int_id(edge));
//                CompareAndCutFromPos(paths, path1, (int) positions1[i1], path2,
//                                     (int) positions2[i2],
//                                     del_subpaths, del_begins, del_all, add_overlap_begins);
//
//                if (positions1[i1] >= path1->Size() || path1->At(positions1[i1]) != edge || positions2[i2] >= path2->Size() || path2->At(positions2[i2]) != edge) {
//                    vector<size_t> new_positions1 = path1->FindAll(edge);
//                    vector<size_t> new_positions2 = path2->FindAll(edge);
//
//                    if (new_positions1.size() == positions1.size() && new_positions2.size() == positions2.size()) {
//                        return;
//                    }
//                    else {
//                        positions1 = new_positions1;
//                        positions2 = new_positions2;
//                        i1 = 0;
//                        i2 = 0;
//                        renewed = true;
//                        break;
//                    }
//                    ++i2;
//                }
//                ++i2;
//            }
//
//            if (renewed) {
//                renewed = false;
//                continue;
//            }
//            ++i1;
//        }
//    }
//
//    void CompareAndCutFromPos(PathContainer& paths, BidirectionalPath* path1, int pos1,
//                              BidirectionalPath* path2, int pos2,
//                              bool delete_subpaths, bool delete_begins,
//                              bool delete_all, bool add_overlap_begins) const {
//        int last2 = pos2;
//        int last1 = pos1;
//        if (last1 >= (int) path1->Size() || last2 >= (int) path2->Size()) {
//            return;
//        }
//        vector<int> other_path_end;
//        pair<int, int> posRes = ComparePaths(last1, last2, *path1, *path2, max_path_diff_);
//        last1 = posRes.first;
//        last2 = posRes.second;
//        BidirectionalPath* conj1 = path1->GetConjPath();
//        BidirectionalPath* conj2 = path2->GetConjPath();
//        size_t first1 = conj1->Size() - pos1 - 1;
//        size_t first2 = conj2->Size() - pos2 - 1;
//        posRes = ComparePaths(first1, first2, *conj1, *conj2, max_path_diff_);
//        first2 = conj2->Size() - posRes.second - 1;
//        first1 = conj1->Size() - posRes.first - 1;
//        if ((int)path2->LengthAt(last2) - (int)g_.length(path2->At(last2)) < (int) max_path_diff_) {
//            last2 = (int)path2->Size() - 1;
//        }
//        if ((int)path2->Length() - (int)path2->LengthAt(first2) < (int) max_path_diff_) {
//            first2 = 0;
//        }
//        if ((int)path1->LengthAt(last1) - (int)g_.length(path1->At(last1)) < (int) max_path_diff_) {
//            last1 = (int)path1->Size() - 1;
//        }
//        if ((int)path1->Length() - (int)path1->LengthAt(first1) < (int) max_path_diff_) {
//            first1 = 0;
//        }
//
//        CutOverlaps(paths, path1, first1, last1, path1->Size(), path2,
//                         first2, last2, path2->Size(), delete_subpaths,
//                         delete_begins, delete_all, add_overlap_begins);
//    }
//
//    void AddOverlap(PathContainer& paths, BidirectionalPath* path, size_t first, size_t last) const {
//        AddPath(paths, path->SubPath(first, last + 1), coverage_map_);
//    }
//
//    bool CutOverlaps(PathContainer& paths, BidirectionalPath* path1,
//                     size_t first1, size_t last1, size_t size1,
//                     BidirectionalPath* path2, size_t first2,
//                     size_t last2, size_t size2,
//                     bool del_subpaths, bool del_begins,
//                     bool del_all, bool add_overlap_begins) const {
//        DEBUG("Path1");
//        path1->PrintDEBUG();
//        DEBUG("Path2");
//        path2->PrintDEBUG();
//        if (first1 == 0 && last1 == size1 - 1 && del_subpaths) {
//            DEBUG("Removing path1 " << path1->GetId() << " because of path2 " << path2->GetId());
//            path1->Clear();
//        } else if (first2 == 0 && last2 == size2 - 1 && del_subpaths) {
//            DEBUG("Removing path2 " << path2->GetId() << " because of path1 " << path1->GetId());
//            path2->Clear();
//        } else if (first2 == 0 && first1 == 0 && del_begins) {
//            DEBUG("Path " << path1->GetId() << ", len " << path1->Length()
//                          << " and path " << path2->GetId() << ", len "
//                          << path2->Length() <<  " have similar starts");
//            DEBUG("Path 1: " << last1 << " edges of length "
//                             << path1->Length() - path1->LengthAt(min(last1 + 1, path1->Size() - 1)));
//            DEBUG("Path 2: " << last2 << " edges of length "
//                             << path2->Length() - path2->LengthAt(min(last2 + 1, path2->Size() - 1)));
//
//            if (add_overlap_begins) {
//                AddOverlap(paths, path1, first1, last1);
//                DEBUG("Detaching overlap " << path2->GetId() << " and " << path1->GetId());
//                PopFront(path1, last1 + 1);
//                PopFront(path2, last2 + 1);
//            } else if (path1->Length() < path2->Length()) {
//                DEBUG("Detaching overlap from " << path1->GetId() << " because of "<< path2->GetId());
//                PopFront(path1, last1 + 1);
//            } else {
//                DEBUG("Detaching overlap from " << path2->GetId() << " because of "<< path1->GetId());
//                PopFront(path2, last2 + 1);
//            }
//        } else if ((last1 == size1 - 1 && last2 == size2 - 1) && del_begins) {
//            DEBUG("Path " << path1->GetId() << ", len " << path1->Length() << " and path " << path2->GetId() << ", len " << path2->Length() << " have similar ends");
//            DEBUG("Path 1: " << path1->Size() - first1 << " edges of length " << path1->LengthAt(first1));
//            DEBUG("Path 2: " << path2->Size() - first2 << " edges of length " << path2->LengthAt(first2));
//
//            if (add_overlap_begins){
//                AddOverlap(paths, path1, first1, last1);
//                DEBUG("Detaching overlap " << path2->GetId() << " and " << path1->GetId());
//                path1->PopBack(last1 + 1 - first1);
//                path2->PopBack(last2 + 1 - first2);
//            }
//            if (path1->Length() < path2->Length()) {
//                DEBUG("Detaching overlap from " << path1->GetId() << " because of "<< path2->GetId());
//                path1->PopBack(last1 + 1 - first1);
//            } else {
//                DEBUG("Detaching overlap from " << path2->GetId() << " because of "<< path1->GetId());
//                path2->PopBack(last2 + 1 - first2);
//            }
//        } else if (first2 == 0 && del_all) {
//            DEBUG("Detaching overlap from " << path2->GetConjPath()->GetId() << " because of "<< path1->GetId());
//            PopFront(path2, last2 + 1);
//        } else if (last2 == size2 - 1 && del_all) {
//            DEBUG("Detaching overlap from " << path2->GetId() << " because of "<< path1->GetId());
//            path2->PopBack(last1 + 1 - first1);
//        } else if (first1 == 0 && del_all) {
//            DEBUG("Detaching overlap from " << path1->GetConjPath()->GetId() << " because of "<< path2->GetId());
//            path1->GetConjPath()->PopBack(last1 + 1);
//        } else if (last1 == size1 - 1 && del_all) {
//            DEBUG("Detaching overlap from " << path1->GetId() << " because of "<< path2->GetId());
//            path1->PopBack(last1 + 1 - first1);
//        } else {
//            return false;
//        }
//        return true;
//    }
//
//    std::vector<EdgeId> GetSortedEdges() const {
//        std::vector<EdgeId> edges(GraphEdgeIterator<Graph>(g_, g_.begin()),
//                                  GraphEdgeIterator<Graph>(g_, g_.end()));
//        std::sort(edges.begin(), edges.end(), [&] (const EdgeId& e1, const EdgeId& e2) {
//            if (g_.length(e1) == g_.length(e2)) {
//                return e1.int_id() < e2.int_id();
//            }
//            return g_.length(e1) > g_.length(e2);
//        });
//        return edges;
//    }
//
//    const Graph& g_;
//    GraphCoverageMap& coverage_map_;
//    size_t min_edge_len_;
//    size_t max_path_diff_;
//protected:
//    DECL_LOGGER("PEResolver")
//};

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

//    void RemoveEqualPaths(PathContainer &paths, GraphCoverageMap &coverage_map,
//                          size_t /*min_edge_len*/, size_t max_path_diff) const {
//        PathDeduplicator deduplicator(g_, paths, coverage_map, max_path_diff, /*equal only*/true);
//        deduplicator.Deduplicate();
//    }

    void RemoveRNAOverlaps(PathContainer& /*paths*/, GraphCoverageMap& /*coverage_map*/,
                          size_t /*min_edge_len*/, size_t /*max_path_diff*/) const {

        VERIFY(false);
        //FIXME remove on cleanup
//        SimpleOverlapRemover remover(g_, coverage_map, min_edge_len, max_path_diff);
//        remover.RemoveSimilarPaths(paths, true, false, false, false, false);
//
//        remover.RemoveSimilarPaths(paths, false, true, false, false, false);
//
//        OverlapRemover(g_, coverage_map).RemoveOverlaps(paths);
//
//        remover.RemoveSimilarPaths(paths, true, false, false, false, false);
    }

    void RemoveOverlaps(PathContainer &paths, GraphCoverageMap &coverage_map,
                        size_t min_edge_len, size_t max_path_diff,
                        bool cut_all) const {
        INFO("Removing overlaps");
        VERIFY(min_edge_len == 0 && max_path_diff == 0);

        INFO("Deduplicating paths");
        Deduplicate(paths, coverage_map, min_edge_len, max_path_diff);

//        DEBUG("Initial paths");
//        for (auto path_pair: paths) {
//            path_pair.first->PrintDEBUG();
//            path_pair.second->PrintDEBUG();
//        }

        //FIXME remove on cleanup
//        //writer.WritePathsToFASTA(paths, output_dir + "/before.fasta");
//        //DEBUG("Removing subpaths");
//        //delete not only eq,
//        remover.RemoveSimilarPaths(paths, false, true, false, false, add_overlaps_begin);
//        //writer.WritePathsToFASTA(paths, output_dir + "/remove_similar.fasta");
//        //DEBUG("Remove overlaps")
//        OverlapRemover(g_, coverage_map).RemoveOverlaps(paths);
//        //writer.WritePathsToFASTA(paths, output_dir + "/after_remove_overlaps.fasta");
//        remover.RemoveSimilarPaths(paths, true, false, false, false, add_overlaps_begin);
//        //writer.WritePathsToFASTA(paths, output_dir + "/remove_equal.fasta");
//        //DEBUG("remove similar path. Max difference " << max_overlap);
//        remover.RemoveSimilarPaths(paths, false, true, true, true, add_overlaps_begin);

        SplitsStorage splits;

//        if (cut_pseudo_self_conjugate)
//            MarkPseudoSelfConjugatePaths(paths, splits);

//        MarkNonUniquePrefixes(paths, splits);
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

        INFO("Deduplicating paths");
        Deduplicate(paths, coverage_map, min_edge_len, max_path_diff);
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

} /* PE_RESOLVER_HPP_ */

#endif
