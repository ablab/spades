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

class OverlapFindingHelper {
    const Graph &g_;
    const GraphCoverageMap &coverage_map_;
    const size_t min_edge_len_;
    const size_t max_diff_;
public:
    OverlapFindingHelper(const Graph &g,
                         const GraphCoverageMap &coverage_map,
                         size_t min_edge_len,
                         size_t max_diff) :
            g_(g),
            coverage_map_(coverage_map),
            min_edge_len_(min_edge_len),
            max_diff_(max_diff) {

    }

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

    pair<size_t, size_t> ComparePaths(const BidirectionalPath &path1,
                                      const BidirectionalPath &path2,
                                      size_t pos2,
                                      bool try_extend = true) const {
        //FIXME change to edit distance
        size_t shift1 = 0;
        size_t i = 0;
        size_t end1 = 0;
        size_t end2 = pos2;

        for (; i < path1.Size(); ++i) {
            shift1 += path1.GapAt(i).gap;
            if (shift1 > max_diff_)
                break;

            bool match = false;
            size_t j = end2;
            size_t shift2 = 0;
            for (; j < path2.Size(); ++j) {
                shift2 += path2.GapAt(j).gap;
                if (shift2 > max_diff_) //shift1 + max_diff)
                    break;
                if (path1.At(i) == path2.At(j)) {
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

        //Shifting end2 to the end of the path if close
        if (try_extend && end1 > 0) {
            TryExtendToEnd(path1, end1);
            TryExtendToEnd(path2, end2);
        }
        return make_pair(end1, end2);
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
        } else {
            for (size_t i = 0; i < path.Size(); ++i) {
                EdgeId e = path.At(i);
                if (g_.length(e) >= min_edge_len_) {
                    utils::insert_all(candidates, coverage_map_.GetCoveringPaths(e));
                }
            }
        }
        VERIFY(candidates.count(&path));
        candidates.erase(&path);
        //FIXME think if this is necessary; related to Anton's code on removal of pseudo-self-conj paths
        candidates.erase(path.GetConjPath());
        return vector<PathPtr>(candidates.begin(), candidates.end());
    }

};

class DecentOverlapRemover {
    const Graph &g_;
    const PathContainer &paths_;
    const GraphCoverageMap &coverage_map_;
    const bool retain_one_copy_;
    const bool end_start_only_;
    const OverlapFindingHelper helper_;
    SplitsStorage splits_;

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

    size_t AnalyzeOverlaps(const BidirectionalPath &path, const BidirectionalPath &other) {
        size_t max_overlap = 0;
        size_t other_start = size_t(-1);
        size_t other_end = size_t(-1);
        for (size_t j = 0; j < other.Size(); ++j) {
            auto ends_pair = helper_.ComparePaths(path, other, j);
            size_t end = ends_pair.first;
            if (end > max_overlap ||
                    (end == max_overlap && ends_pair.second == other.Size())) {
                max_overlap = end;
                other_start = j;
                other_end = ends_pair.second;
            }
        }

        VERIFY(!retain_one_copy_ || !end_start_only_);

        if (end_start_only_ && other_end != other.Size()) {
            return 0;
        }

        //checking if region on the other path has not been already added
        //TODO discuss if the logic is needed/correct. It complicates the procedure and prevents trivial parallelism.
        if (max_overlap == 0 || (retain_one_copy_ && AlreadyAdded(other, other_start, other_end))) {
            return 0;
        }

        return max_overlap;
    }

    void MarkStartOverlaps(const BidirectionalPath &path) {
        set<size_t> overlap_poss;
        for (PathPtr candidate : helper_.FindCandidatePaths(path)) {
            size_t overlap = AnalyzeOverlaps(path, *candidate);
            if (overlap > 0) {
                overlap_poss.insert(overlap);
            }
        }
        splits_[&path] = overlap_poss;
    }

public:
    //FIXME First need to remove subpaths
    //FIXME Order path container in order of increasing length for retain_one to work properly!!!
    DecentOverlapRemover(const Graph &g,
                         const PathContainer &paths,
                         GraphCoverageMap &coverage_map,
                         bool retain_one,
                         bool end_start_only,
                         size_t min_edge_len,// = 0,
                         size_t max_diff) :// = 0) :
            g_(g),
            paths_(paths),
            coverage_map_(coverage_map),
            retain_one_copy_(retain_one),
            end_start_only_(end_start_only),
            helper_(g, coverage_map,
                    min_edge_len, max_diff) {
    }

    void MarkOverlaps() {
        for (auto path_pair: paths_) {
            MarkStartOverlaps(*path_pair.first);
            MarkStartOverlaps(*path_pair.second);
        }
    }

    const SplitsStorage& overlaps() const {
        return splits_;
    }
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

    bool IsSubpath(const BidirectionalPath &path,
                   const BidirectionalPath &other) const {
        for (size_t j = 0; j < other.Size(); ++j) {
            auto ends_pair = helper_.ComparePaths(path, other, j);
            if (ends_pair.first == path.Size()) {
                return true;
            }
        }
        return false;
    }

    bool IsEqual(const BidirectionalPath &path,
                 const BidirectionalPath &other) const {
        auto ends_pair = helper_.ComparePaths(path, other, 0);
        return ends_pair.first == path.Size()
               && ends_pair.second == other.Size();
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

    void Deduplicate() {
        for (auto path_pair : paths_) {
            auto path = path_pair.first;
            for (auto candidate : helper_.FindCandidatePaths(*path)) {
                VERIFY(candidate != path && candidate != path->GetConjPath());
                if (equal_only_ ? IsEqual(*path, *candidate) : IsSubpath(*path, *candidate)) {
                    path->Clear();
                }
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
inline size_t NonUniqueCommon(BidirectionalPath * path, int pos1, int pos2) {
    if (pos1 == pos2)
        return 0;
    size_t answer = 0;
    while (pos1 >= 0) {
        if (path->At(pos1) == path->At(pos2)) {
            pos1--;
            pos2--;
            answer++;
        } else {
            break;
        }
    }
    return answer;
}

inline size_t MaximumNonUniqueSuffix(BidirectionalPath * path) {
    if (path->Size() == 0) {
        return 0;
    }

    size_t answer = 0;
    EdgeId back = path->Back();
    for (size_t pos : path->FindAll(back)) {
        answer = std::max(answer, NonUniqueCommon(path, (int) pos, (int) path->Size() - 1));
    }
    return answer;
}

//FIXME why this has been done via the tmp vector?
inline void CutNonUniqueSuffix(const PathContainer &paths) {
    for (const auto &path_pair : paths) {
        BidirectionalPath * path1 = path_pair.first;
        BidirectionalPath * path2 = path_pair.second;
        path1->PopBack(MaximumNonUniqueSuffix(path1));
        path2->PopBack(MaximumNonUniqueSuffix(path2));
    }
}

//FIXME seems to work out of the box now, but need some thought
inline void CutPseudoSelfConjugatePaths(PathContainer &paths, GraphCoverageMap &cov_map) {
    vector<pair<BidirectionalPath *, BidirectionalPath *>> tmp_paths(paths.begin(), paths.end());
    for (auto it = tmp_paths.begin(); it != tmp_paths.end(); ++it) {
        BidirectionalPath * path1 = it->first;
        BidirectionalPath * path2 = it->second;
        if(path1 != path2) {
            size_t last = 0;
            while (last < path1->Size() && path1->At(last) == path2->At(last)) {
                last++;
            }
            if(last > 0) {
                AddPath(paths, path1->SubPath(0, last), cov_map);
                path1->PopBack(last);
                path2->PopBack(last);
            }
        }
    }
}

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

protected:
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

    void RemoveEqualPaths(PathContainer &paths, GraphCoverageMap &coverage_map,
                          size_t /*min_edge_len*/, size_t max_path_diff) const {
        PathDeduplicator deduplicator(g_, paths, coverage_map, max_path_diff, /*equal only*/true);
        deduplicator.Deduplicate();
    }

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
                        bool cut_all,
                        bool cut_pseudo_self_conjugate) const {
        if (cut_pseudo_self_conjugate)
            CutPseudoSelfConjugatePaths(paths, coverage_map);

        CutNonUniqueSuffix(paths);
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

        VERIFY(min_edge_len == 0 && max_path_diff == 0 && cut_all);
        DecentOverlapRemover overlap_remover(g_, paths, coverage_map,
                                             /*retain one copy*/!cut_all,
                                             /*end/start overlaps only*/false,
                                             min_edge_len, max_path_diff);
        overlap_remover.MarkOverlaps();

        PathSplitter splitter(overlap_remover.overlaps(), paths, coverage_map);
        splitter.Split();

        PathDeduplicator deduplicator(g_, paths, coverage_map, max_path_diff, /*equal only*/false);
        deduplicator.Deduplicate();

        DEBUG("end removing");
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
