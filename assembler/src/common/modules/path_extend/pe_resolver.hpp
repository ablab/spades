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


class SimpleOverlapRemover {

public:
    SimpleOverlapRemover(const Graph& g, GraphCoverageMap& cm)
            : g_(g), coverage_map_(cm) {
    }

    void RemoveOverlaps(PathContainer& paths) const {
        for (size_t i = 0; i < paths.size(); i++) {
            FindAndRemovePathOverlap(paths, paths.Get(i));
            FindAndRemovePathOverlap(paths, paths.GetConjugate(i));
        }
    }

    size_t NonUniqueCommon(BidirectionalPath * path, int pos1, int pos2) {
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

    size_t MaximumNonUniqueSuffix(BidirectionalPath * path) {
        if (path->Size() == 0) {
            return 0;
        }

        size_t answer = 0;
        EdgeId back = path->Back();
        vector<size_t> all_pos = path->FindAll(back);
        for (size_t i = 0; i < all_pos.size() - 1; ++i) {
            answer = std::max(answer, NonUniqueCommon(path, (int) all_pos[i], (int) path->Size() - 1));
        }
        return answer;
    }

    void CutNonUniqueSuffix(PathContainer& paths) {
        vector<pair<BidirectionalPath *, BidirectionalPath *>> tmp_paths(paths.begin(), paths.end());
        for (auto it = tmp_paths.begin(); it != tmp_paths.end(); ++it) {
            BidirectionalPath * path1 = it->first;
            BidirectionalPath * path2 = it->second;
            size_t longest_suffix1 = MaximumNonUniqueSuffix(path1);
            path1->PopBack(longest_suffix1);
            size_t longest_suffix2 = MaximumNonUniqueSuffix(path2);
            path2->PopBack(longest_suffix2);
        }
    }

    void CutPseudoSelfConjugatePaths(PathContainer& paths) {
        vector<pair<BidirectionalPath *, BidirectionalPath *>> tmp_paths(paths.begin(), paths.end());
        for (auto it = tmp_paths.begin(); it != tmp_paths.end(); ++it) {
            BidirectionalPath * path1 = it->first;
            BidirectionalPath * path2 = it->second;
            bool ups = false;
            if(path1 != path2) {
                size_t last = 0;
                while(last < path1->Size() && path1->operator [](last) == path2->operator [](last)) {
                    last++;
                }
                if(last > 0) {
                    AddOverlap(paths, path1, 0, last - 1);
                    path1->PopBack(last);
                    path2->PopBack(last);
                }
            }
            if(ups) path1->Print();
        }
    }

    void RemoveSimilarPaths(PathContainer& paths, size_t min_edge_len, size_t max_path_diff, bool del_only_equal, bool del_subpaths, bool del_begins, bool del_all, bool add_overlap_begins) const {
        DEBUG("== Removing similar paths ==");
        DEBUG("Min edge len " << min_edge_len << ", max path diff " << max_path_diff)
        DEBUG("Only equal " << del_only_equal << ", subpaths " << del_subpaths << ", starts " << del_begins << ", all " << del_all << ", add starts " << add_overlap_begins);
        std::vector<EdgeId> edges = GetSortedEdges();
        for (size_t edgeIndex = 0; edgeIndex < edges.size(); ++edgeIndex) {
            EdgeId edge = edges.at(edgeIndex);
            BidirectionalPathSet cov_paths = coverage_map_.GetCoveringPaths(edge);
            std::vector<BidirectionalPath*> cov_vect(cov_paths.begin(), cov_paths.end());
            std::sort(cov_vect.begin(), cov_vect.end(), PathIdCompare);
            for (size_t vect_i = 0; vect_i < cov_vect.size(); ++vect_i) {
                BidirectionalPath* path1 = cov_vect.at(vect_i);
                if (cov_paths.find(path1) == cov_paths.end()) {
                    continue;
                }
                for (size_t vect_i1 = vect_i + 1; vect_i1 < cov_vect.size(); ++vect_i1) {
                    BidirectionalPath* path2 = cov_vect.at(vect_i1);
                    if (path1 == path2 || path1 == path2->GetConjPath()) {
                        continue;
                    }
                    if (cov_paths.find(path2) == cov_paths.end())
                        continue;
                    if ((*path1) == (*path2)) {
                        if (path2->IsOverlap()) {
                            path1->SetOverlap(true);
                        }
                        DEBUG("Removing path " << path2->GetId() << " because of path " << path1->GetId());
                        path2->Print();
                        path1->Print();
                        path2->Clear();
                        cov_paths = coverage_map_.GetCoveringPaths(edge);
                        continue;
                    }
                    if (g_.length(edge) <= min_edge_len || path1->IsOverlap() || path2->IsOverlap() || del_only_equal) {
                        continue;
                    }
                    CompareAndCut(paths, edge, path1, path2, max_path_diff,
                                  del_subpaths, del_begins, del_all, add_overlap_begins);
                    cov_paths = coverage_map_.GetCoveringPaths(edge);
                }
            }
        }
        DEBUG("== Emd removing similar paths ==");
    }

private:
    
    void SubscribeCoverageMap(BidirectionalPath* path) const {
        path->Subscribe(&coverage_map_);
        for (size_t i = 0; i < path->Size(); ++i) {
            coverage_map_.BackEdgeAdded(path->At(i), path, path->GapAt(i));
        }
    }

    void CompareAndCut(PathContainer& paths, EdgeId edge, BidirectionalPath* path1, BidirectionalPath* path2,
                       size_t max_path_diff,
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
                                     (int) positions2[i2], max_path_diff,
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
                              size_t max_path_diff,
                              bool delete_subpaths, bool delete_begins,
                              bool delete_all, bool add_overlap_begins) const {
        int last2 = pos2;
        int last1 = pos1;
        if (last1 >= (int) path1->Size() || last2 >= (int) path2->Size()) {
            return;
        }
        vector<int> other_path_end;
        pair<int, int> posRes = ComparePaths(last1, last2, *path1, *path2, max_path_diff);
        last1 = posRes.first;
        last2 = posRes.second;
        BidirectionalPath* conj1 = path1->GetConjPath();
        BidirectionalPath* conj2 = path2->GetConjPath();
        size_t first1 = conj1->Size() - pos1 - 1;
        size_t first2 = conj2->Size() - pos2 - 1;
        posRes = ComparePaths(first1, first2, *conj1, *conj2, max_path_diff);
        first2 = conj2->Size() - posRes.second - 1;
        first1 = conj1->Size() - posRes.first - 1;
        if ((int)path2->LengthAt(last2) - (int)g_.length(path2->At(last2)) < (int) max_path_diff) {
            last2 = (int)path2->Size() - 1;
        }
        if ((int)path2->Length() - (int)path2->LengthAt(first2) < (int) max_path_diff) {
            first2 = 0;
        }
        if ((int)path1->LengthAt(last1) - (int)g_.length(path1->At(last1)) < (int) max_path_diff) {
            last1 = (int)path1->Size() - 1;
        }
        if ((int)path1->Length() - (int)path1->LengthAt(first1) < (int) max_path_diff) {
            first1 = 0;
        }

        CutOverlaps(paths, path1, first1, last1, path1->Size(), path2,
                         first2, last2, path2->Size(), delete_subpaths,
                         delete_begins, delete_all, add_overlap_begins);
    }

    void AddOverlap(PathContainer& paths, BidirectionalPath* path1, size_t first1, size_t last1) const {
        BidirectionalPath* overlap = new BidirectionalPath(path1->SubPath(first1, last1 + 1));
        BidirectionalPath* conj_overlap = new BidirectionalPath(overlap->Conjugate());
        SubscribeCoverageMap(overlap);
        SubscribeCoverageMap(conj_overlap);
        paths.AddPair(overlap, conj_overlap);
    }

    bool CutOverlaps(PathContainer& paths, BidirectionalPath* path1, size_t first1, size_t last1, size_t size1, BidirectionalPath* path2, size_t first2,
                     size_t last2, size_t size2, bool del_subpaths, bool del_begins, bool del_all, bool add_overlap_begins) const {
        if (first1 == 0 && last1 == size1 - 1 && del_subpaths) {
            DEBUG("Removing path " << path1->GetId() << " because of path " << path2->GetId());
            path1->Print();
            path2->Print();
            path1->Clear();
        } else if (first2 == 0 && last2 == size2 - 1 && del_subpaths) {
            DEBUG("Removing path " << path2->GetId() << " because of path " << path1->GetId());
            path2->Print();
            path1->Print();
            path2->Clear();
        } else if (first2 == 0 && first1 == 0 && del_begins) {
            DEBUG("Path " << path1->GetId() << ", len " << path1->Length() << " and path " << path2->GetId() << ", len " << path2->Length() <<  " have similar starts");
            DEBUG("Path 1: " << last1 << " edges of length " << path1->Length() - path1->LengthAt(min(last1 + 1, path1->Size() - 1)));
            DEBUG("Path 2: " << last2 << " edges of length " << path2->Length() - path2->LengthAt(min(last2 + 1, path2->Size() - 1)));
            DEBUG("Path 1 has overlap start " << path1->HasOverlapedBegin() << ", path 2 has overlap start " <<  path2->HasOverlapedBegin());

            if (add_overlap_begins) {
                AddOverlap(paths, path1, first1, last1);
                DEBUG("Detaching overlap " << path2->GetId() << " and " << path1->GetId());
                path2->Print();
                path1->Print();
                path1->GetConjPath()->PopBack(last1 + 1);
                path2->GetConjPath()->PopBack(last2 + 1);
            } else if (path1->Length() < path2->Length()) {
                DEBUG("Detaching overlap from " << path1->GetId() << " because of "<< path2->GetId());
                path1->Print();
                path2->Print();
                path1->GetConjPath()->PopBack(last1 + 1);
            } else {
                DEBUG("Detaching overlap from " << path2->GetId() << " because of "<< path1->GetId());
                path2->Print();
                path1->Print();
                path2->GetConjPath()->PopBack(last2 + 1);
            }
        } else if ((last1 == size1 - 1 && last2 == size2 - 1) && del_begins) {
            DEBUG("Path " << path1->GetId() << ", len " << path1->Length() << " and path " << path2->GetId() << ", len " << path2->Length() << " have similar ends");
            DEBUG("Path 1: " << path1->Size() - first1 << " edges of length " << path1->LengthAt(first1));
            DEBUG("Path 2: " << path2->Size() - first2 << " edges of length " << path2->LengthAt(first2));
            DEBUG("Path 1 has overlap end " << path1->HasOverlapedEnd() << ", path 2 has overlap end " <<  path2->HasOverlapedEnd());

            if (add_overlap_begins){
                AddOverlap(paths, path1, first1, last1);
                DEBUG("Detaching overlap " << path2->GetId() << " and " << path1->GetId());
                path2->Print();
                path1->Print();
                path1->PopBack(last1 + 1 - first1);
                path2->PopBack(last2 + 1 - first2);
            }
            if (path1->Length() < path2->Length()) {
                DEBUG("Detaching overlap from " << path1->GetId() << " because of "<< path2->GetId());
                path1->Print();
                path2->Print();
                path1->PopBack(last1 + 1 - first1);
            } else {
                DEBUG("Detaching overlap from " << path2->GetId() << " because of "<< path1->GetId());
                path2->Print();
                path1->Print();
                path2->PopBack(last2 + 1 - first2);
            }
        } else if (first2 == 0 && del_all) {
            DEBUG("Detaching overlap from " << path2->GetConjPath()->GetId() << " because of "<< path1->GetId());
            DEBUG("Does it have overlap in the beginning: " << path2->HasOverlapedBegin());
            path2->Print();
            DEBUG(" >>>> ")
            path1->Print();
            DEBUG(" ==== ");
            path2->GetConjPath()->PopBack(last2 + 1);
        } else if (last2 == size2 - 1 && del_all) {
            DEBUG("Detaching overlap from " << path2->GetId() << " because of "<< path1->GetId());
            DEBUG("Does it have overlap in the end: " << path2->HasOverlapedEnd());
            path2->Print();
            DEBUG(" >>>> ")
            path1->Print();
            DEBUG(" ==== ");
            path2->PopBack(last1 + 1 - first1);
        } else if (first1 == 0 && del_all) {
            DEBUG("Detaching overlap from " << path1->GetConjPath()->GetId() << " because of "<< path2->GetId());
            DEBUG("Does it have overlap in the end: " << path1->HasOverlapedBegin());
            path1->Print();
            DEBUG(" >>>> ")
            path2->Print();
            DEBUG(" ==== ");
            path1->GetConjPath()->PopBack(last1 + 1);
        } else if (last1 == size1 - 1 && del_all) {
            DEBUG("Detaching overlap from " << path1->GetId() << " because of "<< path2->GetId());
            DEBUG("Does it have overlap in the end: " << path1->HasOverlapedBegin());
            path1->Print();
            DEBUG(" >>>> ")
            path2->Print();
            DEBUG(" ==== ");
            path1->PopBack(last1 + 1 - first1);
        } else {
            return false;
        }
        return true;
    }

    std::vector<EdgeId> GetSortedEdges() const {
        std::set<EdgeId> edges_set;
        for (auto iter = g_.ConstEdgeBegin(); !iter.IsEnd(); ++iter) {
            edges_set.insert(*iter);
            edges_set.insert(g_.conjugate(*iter));
        }
        std::vector<EdgeId> edges(edges_set.begin(), edges_set.end());
        std::sort(edges.begin(), edges.end(), EdgeLengthAndIdComparator(g_));
        return edges;
    }

    bool HasAlreadyOverlapedEnd(BidirectionalPath * path) const {
        return !path->IsOverlap() and path->HasOverlapedEnd();
    }

    bool HasAlreadyOverlapedBegin(BidirectionalPath * path) const {
        return !path->IsOverlap() and path->HasOverlapedBegin();
    }

    bool IsSamePath(BidirectionalPath * path1,
                    BidirectionalPath * path2) const {
        return *path2 == *path1 or *path2 == *path1->GetConjPath();
    }

    void RemoveOverlap(PathContainer& paths, BidirectionalPath* path1,
                       BidirectionalPath* path2, size_t overlap_size) const {
        BidirectionalPath* conj2 = path2->GetConjPath();
        if (path1->IsOverlap() && overlap_size == path1->Size()) {
            DEBUG("Detaching overlap from " << path2->GetConjPath()->GetId() << " because of "<< path1->GetId());
            path2->Print();
            path1->Print();
            conj2->PopBack(overlap_size);
            path2->SetOverlapedBeginTo(path1);
        } else if (path2->IsOverlap() && path2->Size() == overlap_size) {
            DEBUG("Detaching overlap from " << path1->GetId() << " because of "<< path2->GetId());
            path1->Print();
            path2->Print();
            path1->PopBack(overlap_size);
            path1->SetOverlapedEndTo(path2);
        } else if (overlap_size < path2->Size()
                && overlap_size < path1->Size()) {
            BidirectionalPath *overlap = new BidirectionalPath(g_, path1->Back());
            BidirectionalPath *conj_overlap = new BidirectionalPath(g_, g_.conjugate(path1->Back()));
            SubscribeCoverageMap(overlap);
            SubscribeCoverageMap(conj_overlap);
            paths.AddPair(overlap, conj_overlap);
            DEBUG("Detaching overlap " << path1->GetId() << " and " << conj2->GetId());
            path1->Print();
            conj2->Print();
            path1->PopBack();
            conj2->PopBack();

            for (size_t i = 1; i < overlap_size; ++i) {
                conj_overlap->PushBack(g_.conjugate(path1->Back()));
                path1->PopBack();
                conj2->PopBack();
            }
            overlap->SetOverlap(true);
            path1->SetOverlapedEndTo(overlap);
            path2->SetOverlapedBeginTo(overlap);
        }
    }

    void FindAndRemovePathOverlap(PathContainer& all_paths,
                                  BidirectionalPath* path1) const {
        int last = (int) path1->Size() - 1;
        if (last <= 0 or coverage_map_.GetCoverage(path1->At(last)) <= 1) {
            return;
        }
        BidirectionalPathSet paths =
                coverage_map_.GetCoveringPaths(path1->At(last));
        BidirectionalPath* overlap_path = NULL;
        size_t overlap_size = 0;
        for (auto path_iter = paths.begin(); path_iter != paths.end();
                ++path_iter) {
            if (IsSamePath(*path_iter, path1)) {
                continue;
            }
            size_t over_size = path1->OverlapEndSize(*path_iter);
            if (over_size > overlap_size) {
                overlap_size = over_size;
                overlap_path = *path_iter;
            } else if (over_size == overlap_size &&
                    (overlap_path == NULL || (*path_iter)->GetId() < overlap_path->GetId())) {
                overlap_path = *path_iter;
            }
        }
        if (overlap_path == NULL) {
            return;
        }
        if (overlap_size > 0) {
            RemoveOverlap(all_paths, path1, overlap_path, overlap_size);
        }
    }

    class EdgeLengthAndIdComparator {
    public:
        EdgeLengthAndIdComparator(const Graph& g)
                : g_(g) {
        }
        bool operator()(const EdgeId& e1, const EdgeId& e2) const {
            if (g_.length(e1) > g_.length(e2)) {
                return true;
            }
            if (g_.length(e2) > g_.length(e1)) {
                return false;
            }
            return e1.int_id() < e2.int_id();
        }
    private:
        const Graph& g_;
    };

    const Graph& g_;
    GraphCoverageMap& coverage_map_;
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
        for (auto iter = g_.ConstEdgeBegin(); !iter.IsEnd(); ++iter) {
            if (g_.int_id(*iter) <= 0 or InTwoEdgeCycle(*iter, g_))
                continue;
            if (included.count(*iter) == 0) {
                BidirectionalPath * first = new BidirectionalPath(g_, *iter);
                BidirectionalPath * second = new BidirectionalPath(g_, g_.conjugate(*iter));
                edges.AddPair(first,second);
                included.insert(*iter);
                included.insert(g_.conjugate(*iter));
            }
        }
        return edges;
    }

    PathContainer ExtendSeeds(PathContainer &seeds, ContigsMaker &pathExtender) const {
        PathContainer paths;
        pathExtender.GrowAll(seeds, paths);
        return paths;
    }

    void RemoveEqualPaths(PathContainer &paths, GraphCoverageMap &coverage_map,
                          size_t min_edge_len) const  {

        SimpleOverlapRemover remover(g_, coverage_map);
        remover.RemoveSimilarPaths(paths, min_edge_len, min_edge_len, true, false, false, false, false);
    }

    void RemoveRNAOverlaps(PathContainer& paths, GraphCoverageMap& coverage_map,
                          size_t min_edge_len, size_t max_path_diff) const  {

        SimpleOverlapRemover remover(g_, coverage_map);
        remover.RemoveSimilarPaths(paths, min_edge_len, max_path_diff, true, false, false, false, false);

        remover.RemoveSimilarPaths(paths, min_edge_len, max_path_diff, false, true, false, false, false);

        remover.RemoveOverlaps(paths);

        remover.RemoveSimilarPaths(paths, min_edge_len, max_path_diff, true, false, false, false, false);
    }

    void RemoveOverlaps(PathContainer &paths, GraphCoverageMap &coverage_map,
                        size_t min_edge_len, size_t max_path_diff,
                        bool add_overlaps_begin,
                        bool cut_preudo_self_conjugate) const {
        SimpleOverlapRemover remover(g_, coverage_map);
        if (cut_preudo_self_conjugate)
            remover.CutPseudoSelfConjugatePaths(paths);

        remover.CutNonUniqueSuffix(paths);
        //writer.WritePathsToFASTA(paths, output_dir + "/before.fasta");
        //DEBUG("Removing subpaths");
        //delete not only eq,
        remover.RemoveSimilarPaths(paths, min_edge_len, max_path_diff, false, true, false, false, add_overlaps_begin);
        //writer.WritePathsToFASTA(paths, output_dir + "/remove_similar.fasta");
        //DEBUG("Remove overlaps")
        remover.RemoveOverlaps(paths);
        //writer.WritePathsToFASTA(paths, output_dir + "/after_remove_overlaps.fasta");
        remover.RemoveSimilarPaths(paths, min_edge_len, max_path_diff, true, false, false, false, add_overlaps_begin);
        //writer.WritePathsToFASTA(paths, output_dir + "/remove_equal.fasta");
        //DEBUG("remove similar path. Max difference " << max_overlap);
        remover.RemoveSimilarPaths(paths, min_edge_len, max_path_diff, false, true, true, true, add_overlaps_begin);
        DEBUG("end removing");
    }

    void RemoveMatePairEnds(PathContainer& paths, size_t min_edge_len) const {
        DEBUG("remove mp ends");
        for (size_t i = 0; i < paths.size(); ++i) {
            RemoveMatePairEnd(*paths.Get(i), min_edge_len);
            RemoveMatePairEnd(*paths.GetConjugate(i), min_edge_len);
        }
    }

    void AddUncoveredEdges(PathContainer &paths, GraphCoverageMap &coverageMap) const {
        std::set<EdgeId> included;
        for (auto iter = g_.ConstEdgeBegin(); !iter.IsEnd(); ++iter) {
            if (included.count(*iter) == 0 && !coverageMap.IsCovered(*iter)) {
                BidirectionalPath* path = new BidirectionalPath(g_, *iter);
                BidirectionalPath* conj = new BidirectionalPath(g_, g_.conjugate(*iter));
                path->Subscribe(&coverageMap);
                conj->Subscribe(&coverageMap);
                coverageMap.BackEdgeAdded(path->At(0), path, path->GapAt(0));
                coverageMap.BackEdgeAdded(conj->At(0), conj, conj->GapAt(0));
                paths.AddPair(path, conj);
                included.insert(*iter);
                included.insert(g_.conjugate(*iter));
            }
        }
    }

private:
    void RemoveMatePairEnd(BidirectionalPath& path, size_t min_edge_len) const {
        int pos = int(path.Size()) - 1;
        while (pos > 0 and g_.length(path.At(pos)) < min_edge_len) {
            path.PopBack();
            pos--;
        }
    }
protected:
    DECL_LOGGER("PEResolver")
};

} /* PE_RESOLVER_HPP_ */

#endif
