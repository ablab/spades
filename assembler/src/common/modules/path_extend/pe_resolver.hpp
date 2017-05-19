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

class OverlapRemover {
    const Graph& g_;
    GraphCoverageMap& coverage_map_;

    bool IsSamePath(BidirectionalPath * path1,
                    BidirectionalPath * path2) const {
        return *path2 == *path1 || *path2 == *path1->GetConjPath();
    }

    void RemoveOverlap(PathContainer& paths, BidirectionalPath* path1,
                       BidirectionalPath* path2, size_t overlap_size) const {
        DEBUG("Path 1");
        path1->PrintDEBUG();
        DEBUG("Path 2");
        path2->PrintDEBUG();
        BidirectionalPath* conj2 = path2->GetConjPath();
        if (path1->IsOverlap() && overlap_size == path1->Size()) {
            DEBUG("Detaching overlap from path 2 " << path2->GetConjPath()->GetId()
                                                   << " because of path 1 " << path1->GetId());
            conj2->PopBack(overlap_size);
        } else if (path2->IsOverlap() && path2->Size() == overlap_size) {
            DEBUG("Detaching overlap from path 1 " << path1->GetId() << " because of path 2" << path2->GetId());
            path1->PopBack(overlap_size);
        } else if (overlap_size < path2->Size()
                   && overlap_size < path1->Size()) {
            AddPath(paths, path1->SubPath(path1->Size() - overlap_size), coverage_map_)->SetOverlap(true);
            DEBUG("Detaching overlap " << path1->GetId() << " and " << conj2->GetId());
            path1->PopBack(overlap_size);
            conj2->PopBack(overlap_size);
        }
    }

    void FindAndRemovePathOverlap(PathContainer& all_paths,
                                  BidirectionalPath* path) const {
        //FIXME why case of single edge is skipped?
        if (path->Size() <= 1 ||
            coverage_map_.GetCoverage(path->Back()) <= 1) {
            return;
        }
        size_t max_overlap_size = 0;
        BidirectionalPath* overlap_path = nullptr;
        for (auto p : coverage_map_.GetCoveringPaths(path->Back())) {
            if (IsSamePath(p, path)) {
                continue;
            }
            size_t overlap_size = OverlapSize(*path, *p);
            if (overlap_size > max_overlap_size) {
                max_overlap_size = overlap_size;
                overlap_path = p;
            } else if (overlap_size > 0
                       && overlap_size == max_overlap_size
                       && p->GetId() < overlap_path->GetId()) {
                //FIXME check if this logic is needed
                overlap_path = p;
            }
        }
        if (max_overlap_size > 0)
            RemoveOverlap(all_paths, path, overlap_path, max_overlap_size);
    }

public:
    OverlapRemover(const Graph& g, GraphCoverageMap& cm)
            : g_(g), coverage_map_(cm) {
    }

    void RemoveOverlaps(PathContainer& paths) const {
        for (size_t i = 0; i < paths.size(); i++) {
            FindAndRemovePathOverlap(paths, paths.Get(i));
            FindAndRemovePathOverlap(paths, paths.GetConjugate(i));
        }
    }

};

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

class SimpleOverlapRemover {

public:
    SimpleOverlapRemover(const Graph& g,
                         GraphCoverageMap& cm,
                         size_t min_edge_len,
                         size_t max_path_diff)
            : g_(g), coverage_map_(cm),
              min_edge_len_(min_edge_len),
              max_path_diff_(max_path_diff) {
    }

    void CutPseudoSelfConjugatePaths(PathContainer& paths) {
        vector<pair<BidirectionalPath *, BidirectionalPath *>> tmp_paths(paths.begin(), paths.end());
        for (auto it = tmp_paths.begin(); it != tmp_paths.end(); ++it) {
            BidirectionalPath * path1 = it->first;
            BidirectionalPath * path2 = it->second;
            if(path1 != path2) {
                size_t last = 0;
                while(last < path1->Size() && path1->operator [](last) == path2->operator [](last)) {
                    last++;
                }
                if(last > 0) {
                    AddOverlap(paths, path1, 0, last - 1);
                    //FIXME is it even correct?!
                    path1->PopBack(last);
                    path2->PopBack(last);
                }
            }
        }
    }

    void RemoveSimilarPaths(PathContainer& paths,
                            bool del_only_equal, bool del_subpaths, bool del_begins,
                            bool del_all, bool add_overlap_begins) const {
        DEBUG("== Removing similar paths ==");
        DEBUG("Min edge len " << min_edge_len_ << ", max path diff " << max_path_diff_);
        DEBUG("Only equal " << del_only_equal << ", subpaths " << del_subpaths
                            << ", starts " << del_begins << ", all " << del_all
                            << ", add starts " << add_overlap_begins);
        //FIXME loop over edges looks a really bad idea
        for (EdgeId edge : GetSortedEdges()) {
            BidirectionalPathSet cov_paths = coverage_map_.GetCoveringPaths(edge);
            std::vector<BidirectionalPath*> cov_vect(cov_paths.begin(), cov_paths.end());
            std::sort(cov_vect.begin(), cov_vect.end(), [] (BidirectionalPath* p1, BidirectionalPath* p2) {
                return p1->GetId() < p2->GetId();
            });
            for (size_t vect_i = 0; vect_i < cov_vect.size(); ++vect_i) {
                BidirectionalPath* path1 = cov_vect.at(vect_i);
                if (cov_paths.find(path1) == cov_paths.end()) {
                    continue;
                }
                VERIFY(cov_paths.find(path1) != cov_paths.end());
                for (size_t vect_i1 = vect_i + 1; vect_i1 < cov_vect.size(); ++vect_i1) {
                    BidirectionalPath* path2 = cov_vect.at(vect_i1);
                    if (path1 == path2 || path1 == path2->GetConjPath()) {
                        continue;
                    }
                    if (cov_paths.find(path2) == cov_paths.end())
                        continue;
                    VERIFY(cov_paths.find(path2) != cov_paths.end());
                    if ((*path1) == (*path2)) {
                        if (path2->IsOverlap()) {
                            path1->SetOverlap(true);
                        }
                        DEBUG("Removing path " << path2->GetId() << " because of path " << path1->GetId());
                        path2->PrintDEBUG();
                        path1->PrintDEBUG();
                        path2->Clear();
                        cov_paths = coverage_map_.GetCoveringPaths(edge);
                        continue;
                    }
                    if (g_.length(edge) <= min_edge_len_ || path1->IsOverlap() || path2->IsOverlap() || del_only_equal) {
                        continue;
                    }
                    CompareAndCut(paths, edge, path1, path2,
                                  del_subpaths, del_begins, del_all, add_overlap_begins);
                    cov_paths = coverage_map_.GetCoveringPaths(edge);
                }
            }
        }
        DEBUG("== Emd removing similar paths ==");
    }

private:
    
    void CompareAndCut(PathContainer& paths, EdgeId edge, BidirectionalPath* path1, BidirectionalPath* path2,
                       bool del_subpaths, bool del_begins,
                       bool del_all, bool add_overlap_begins) const {
        vector<size_t> positions1 = path1->FindAll(edge);
        vector<size_t> positions2 = path2->FindAll(edge);
        size_t i1 = 0;
        size_t i2 = 0;
        bool renewed = false;
        while (i1 < positions1.size()) {
            while (i2 < positions2.size()) {
                DEBUG("CompareAndCutFromPos paths " << g_.int_id(edge));
                CompareAndCutFromPos(paths, path1, (int) positions1[i1], path2,
                                     (int) positions2[i2],
                                     del_subpaths, del_begins, del_all, add_overlap_begins);

                if (positions1[i1] >= path1->Size() || path1->At(positions1[i1]) != edge || positions2[i2] >= path2->Size() || path2->At(positions2[i2]) != edge) {
                    vector<size_t> new_positions1 = path1->FindAll(edge);
                    vector<size_t> new_positions2 = path2->FindAll(edge);

                    if (new_positions1.size() == positions1.size() && new_positions2.size() == positions2.size()) {
                        return;
                    }
                    else {
                        positions1 = new_positions1;
                        positions2 = new_positions2;
                        i1 = 0;
                        i2 = 0;
                        renewed = true;
                        break;
                    }
                    ++i2;
                }
                ++i2;
            }

            if (renewed) {
                renewed = false;
                continue;
            }
            ++i1;
        }
    }

    void CompareAndCutFromPos(PathContainer& paths, BidirectionalPath* path1, int pos1,
                              BidirectionalPath* path2, int pos2,
                              bool delete_subpaths, bool delete_begins,
                              bool delete_all, bool add_overlap_begins) const {
        int last2 = pos2;
        int last1 = pos1;
        if (last1 >= (int) path1->Size() || last2 >= (int) path2->Size()) {
            return;
        }
        vector<int> other_path_end;
        pair<int, int> posRes = ComparePaths(last1, last2, *path1, *path2, max_path_diff_);
        last1 = posRes.first;
        last2 = posRes.second;
        BidirectionalPath* conj1 = path1->GetConjPath();
        BidirectionalPath* conj2 = path2->GetConjPath();
        size_t first1 = conj1->Size() - pos1 - 1;
        size_t first2 = conj2->Size() - pos2 - 1;
        posRes = ComparePaths(first1, first2, *conj1, *conj2, max_path_diff_);
        first2 = conj2->Size() - posRes.second - 1;
        first1 = conj1->Size() - posRes.first - 1;
        if ((int)path2->LengthAt(last2) - (int)g_.length(path2->At(last2)) < (int) max_path_diff_) {
            last2 = (int)path2->Size() - 1;
        }
        if ((int)path2->Length() - (int)path2->LengthAt(first2) < (int) max_path_diff_) {
            first2 = 0;
        }
        if ((int)path1->LengthAt(last1) - (int)g_.length(path1->At(last1)) < (int) max_path_diff_) {
            last1 = (int)path1->Size() - 1;
        }
        if ((int)path1->Length() - (int)path1->LengthAt(first1) < (int) max_path_diff_) {
            first1 = 0;
        }

        CutOverlaps(paths, path1, first1, last1, path1->Size(), path2,
                         first2, last2, path2->Size(), delete_subpaths,
                         delete_begins, delete_all, add_overlap_begins);
    }

    void AddOverlap(PathContainer& paths, BidirectionalPath* path, size_t first, size_t last) const {
        AddPath(paths, path->SubPath(first, last + 1), coverage_map_);
    }

    void TrimFront(BidirectionalPath *path, size_t cnt) const {
        path->GetConjPath()->PopBack(cnt);
    }

    bool CutOverlaps(PathContainer& paths, BidirectionalPath* path1,
                     size_t first1, size_t last1, size_t size1,
                     BidirectionalPath* path2, size_t first2,
                     size_t last2, size_t size2,
                     bool del_subpaths, bool del_begins,
                     bool del_all, bool add_overlap_begins) const {
        DEBUG("Path1");
        path1->PrintDEBUG();
        DEBUG("Path2");
        path2->PrintDEBUG();
        if (first1 == 0 && last1 == size1 - 1 && del_subpaths) {
            DEBUG("Removing path1 " << path1->GetId() << " because of path2 " << path2->GetId());
            path1->Clear();
        } else if (first2 == 0 && last2 == size2 - 1 && del_subpaths) {
            DEBUG("Removing path2 " << path2->GetId() << " because of path1 " << path1->GetId());
            path2->Clear();
        } else if (first2 == 0 && first1 == 0 && del_begins) {
            DEBUG("Path " << path1->GetId() << ", len " << path1->Length()
                          << " and path " << path2->GetId() << ", len "
                          << path2->Length() <<  " have similar starts");
            DEBUG("Path 1: " << last1 << " edges of length "
                             << path1->Length() - path1->LengthAt(min(last1 + 1, path1->Size() - 1)));
            DEBUG("Path 2: " << last2 << " edges of length "
                             << path2->Length() - path2->LengthAt(min(last2 + 1, path2->Size() - 1)));

            if (add_overlap_begins) {
                AddOverlap(paths, path1, first1, last1);
                DEBUG("Detaching overlap " << path2->GetId() << " and " << path1->GetId());
                TrimFront(path1, last1 + 1);
                TrimFront(path2, last2 + 1);
            } else if (path1->Length() < path2->Length()) {
                DEBUG("Detaching overlap from " << path1->GetId() << " because of "<< path2->GetId());
                TrimFront(path1, last1 + 1);
            } else {
                DEBUG("Detaching overlap from " << path2->GetId() << " because of "<< path1->GetId());
                TrimFront(path2, last2 + 1);
            }
        } else if ((last1 == size1 - 1 && last2 == size2 - 1) && del_begins) {
            DEBUG("Path " << path1->GetId() << ", len " << path1->Length() << " and path " << path2->GetId() << ", len " << path2->Length() << " have similar ends");
            DEBUG("Path 1: " << path1->Size() - first1 << " edges of length " << path1->LengthAt(first1));
            DEBUG("Path 2: " << path2->Size() - first2 << " edges of length " << path2->LengthAt(first2));

            if (add_overlap_begins){
                AddOverlap(paths, path1, first1, last1);
                DEBUG("Detaching overlap " << path2->GetId() << " and " << path1->GetId());
                path1->PopBack(last1 + 1 - first1);
                path2->PopBack(last2 + 1 - first2);
            }
            if (path1->Length() < path2->Length()) {
                DEBUG("Detaching overlap from " << path1->GetId() << " because of "<< path2->GetId());
                path1->PopBack(last1 + 1 - first1);
            } else {
                DEBUG("Detaching overlap from " << path2->GetId() << " because of "<< path1->GetId());
                path2->PopBack(last2 + 1 - first2);
            }
        } else if (first2 == 0 && del_all) {
            DEBUG("Detaching overlap from " << path2->GetConjPath()->GetId() << " because of "<< path1->GetId());
            TrimFront(path2, last2 + 1);
        } else if (last2 == size2 - 1 && del_all) {
            DEBUG("Detaching overlap from " << path2->GetId() << " because of "<< path1->GetId());
            path2->PopBack(last1 + 1 - first1);
        } else if (first1 == 0 && del_all) {
            DEBUG("Detaching overlap from " << path1->GetConjPath()->GetId() << " because of "<< path2->GetId());
            path1->GetConjPath()->PopBack(last1 + 1);
        } else if (last1 == size1 - 1 && del_all) {
            DEBUG("Detaching overlap from " << path1->GetId() << " because of "<< path2->GetId());
            path1->PopBack(last1 + 1 - first1);
        } else {
            return false;
        }
        return true;
    }

    std::vector<EdgeId> GetSortedEdges() const {
        std::vector<EdgeId> edges(GraphEdgeIterator<Graph>(g_, g_.begin()),
                                  GraphEdgeIterator<Graph>(g_, g_.end()));
        std::sort(edges.begin(), edges.end(), [&] (const EdgeId& e1, const EdgeId& e2) {
            if (g_.length(e1) == g_.length(e2)) {
                return e1.int_id() < e2.int_id();
            }
            return g_.length(e1) > g_.length(e2);
        });
        return edges;
    }

    const Graph& g_;
    GraphCoverageMap& coverage_map_;
    size_t min_edge_len_;
    size_t max_path_diff_;
protected:
    DECL_LOGGER("PEResolver")
};

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
                          size_t min_edge_len) const {

        SimpleOverlapRemover remover(g_, coverage_map, min_edge_len, min_edge_len);
        remover.RemoveSimilarPaths(paths, true, false, false, false, false);
    }

    void RemoveRNAOverlaps(PathContainer& paths, GraphCoverageMap& coverage_map,
                          size_t min_edge_len, size_t max_path_diff) const {

        SimpleOverlapRemover remover(g_, coverage_map, min_edge_len, max_path_diff);
        remover.RemoveSimilarPaths(paths, true, false, false, false, false);

        remover.RemoveSimilarPaths(paths, false, true, false, false, false);

        OverlapRemover(g_, coverage_map).RemoveOverlaps(paths);

        remover.RemoveSimilarPaths(paths, true, false, false, false, false);
    }

    void RemoveOverlaps(PathContainer &paths, GraphCoverageMap &coverage_map,
                        size_t min_edge_len, size_t max_path_diff,
                        bool add_overlaps_begin,
                        bool cut_preudo_self_conjugate) const {
        SimpleOverlapRemover remover(g_, coverage_map, min_edge_len, max_path_diff);
        if (cut_preudo_self_conjugate)
            remover.CutPseudoSelfConjugatePaths(paths);

        CutNonUniqueSuffix(paths);
        //writer.WritePathsToFASTA(paths, output_dir + "/before.fasta");
        //DEBUG("Removing subpaths");
        //delete not only eq,
        remover.RemoveSimilarPaths(paths, false, true, false, false, add_overlaps_begin);
        //writer.WritePathsToFASTA(paths, output_dir + "/remove_similar.fasta");
        //DEBUG("Remove overlaps")
        OverlapRemover(g_, coverage_map).RemoveOverlaps(paths);
        //writer.WritePathsToFASTA(paths, output_dir + "/after_remove_overlaps.fasta");
        remover.RemoveSimilarPaths(paths, true, false, false, false, add_overlaps_begin);
        //writer.WritePathsToFASTA(paths, output_dir + "/remove_equal.fasta");
        //DEBUG("remove similar path. Max difference " << max_overlap);
        remover.RemoveSimilarPaths(paths, false, true, true, true, add_overlaps_begin);
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
