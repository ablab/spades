/*
 * mate_pair_scaffolding.hpp
 *
 *  Created on: Jul 31, 2013
 *      Author: ira
 */

#ifndef MATE_PAIR_SCAFFOLDING_HPP_
#define MATE_PAIR_SCAFFOLDING_HPP_

#include <map>
#include <vector>
#include "bidirectional_path.hpp"
#include "paired_library.hpp"
#include "pe_utils.hpp"
#include "../graph_pack.hpp"
#include "../../include/de/paired_info.hpp"
#include "../../include/omni/path_processor.hpp"

namespace path_extend {

/*using debruijn_graph::EdgeId;

struct IndexPair {
    size_t index1_;
    size_t index2_;
    IndexPair(size_t index1, size_t index2)
            : index1_(index1),
              index2_(index2) {
    }
    bool operator<(const IndexPair& pair) const {
        if (index1_ < pair.index1_) {
            return true;
        } else if (index1_ == pair.index1_) {
            return index2_ < pair.index2_;
        }
        return false;
    }
};

double CountNewDistance(double old_d, double d, double old_w, double w) {
    if (math::gr(w + old_w, 0.0)) {
        return d = (old_d * old_w + d * w) / (w + old_w);
    }
    return 0.0;
}

class PathsEdgesPairInfo {
public:
    PathsEdgesPairInfo()
            : indexes_(0, 0),
              w_(0.0),
              d_(0.0) {

    }
    PathsEdgesPairInfo(size_t index1, size_t index2, double w, double d)
            : indexes_(index1, index2),
              w_(w),
              d_(d) {
    }
    PathsEdgesPairInfo(const PathsEdgesPairInfo& pi)
            : indexes_(pi.index1(), pi.index2()),
              w_(pi.w()),
              d_(pi.d()),
              distances_(pi.distances_) {

    }
    size_t index1() const {
        return indexes_.index1_;
    }
    size_t index2() const {
        return indexes_.index2_;
    }

    double w() const {
        return w_;
    }

    double d() const {
        return d_;
    }

    void AddPairInfo(double w, double d) {
        d_ = CountNewDistance(d_, d, w_, w);
        w_ += w;
        distances_[d] = w;
    }

    bool operator<(const PathsEdgesPairInfo& pair) const {
        return indexes_ < pair.indexes_;
    }
private:
    IndexPair indexes_;
    double w_;
    double d_;
    std::map<double, double> distances_;
};

class PairPathsInfo {
public:
    PairPathsInfo(const BidirectionalPath* path1,
                  const BidirectionalPath* path2)
            : path1_(path1),
              path2_(path2),
              weight_(0.0),
              dist_(0.0) {
    }
    void AddPairInfo(size_t index1, size_t index2, double w, double d) {
        double dist_between_paths = d - path1_->LengthAt(index1)
                - (path2_->Length() - path2_->LengthAt(index2));
        dist_ = CountNewDistance(dist_, dist_between_paths, weight_, w);
        weight_ += w;
        IndexPair pair(index1, index2);
        if (pair_edges_.find(pair) == pair_edges_.end()) {
            PathsEdgesPairInfo pair_info(index1, index2, w, d);
            pair_edges_[pair] = pair_info;
        } else {
            pair_edges_[pair].AddPairInfo(w, d);
        }
    }
    bool operator<(const PairPathsInfo& pair) const {
        if (path1_->GetId() < pair.path1_->GetId()) {
            return true;
        } else if (path1_->GetId() == pair.path1_->GetId()) {
            return path2_->GetId() < pair.path2_->GetId();
        }
        return false;
    }

    double d() const {
        return dist_;
    }

    double weight() const {
        return weight_;
    }

    bool CanBeValid(size_t max_dist, size_t min_dist) const {
        const Graph& g = path1_->graph();
        PathContainer paths = FindPathBetween(g, max_dist);
        if (paths.size() == 0) {
            return false;
        }
        for (size_t i = 0; i < paths.size(); ++i) {
            size_t curr_dist = paths[i].Length();
            std::set<PathsEdgesPairInfo> ideal_pi = GenerateIdealPairs(
                    g, curr_dist, max_dist, min_dist);
            for (auto pi = ideal_pi.begin(); pi != ideal_pi.end(); ++pi) {
                IndexPair index(pi->index1(), pi->index2());
                if (pair_edges_.find(index) != pair_edges_.end()
                        && pair_edges_.find(index)->second.w() > 2) {
                    return true;
                }
            }
        }
        return false;
    }

    void AnalyzeAndPrint(size_t max_dist, size_t min_dist) const {
        const Graph& g = path1_->graph();
        PathContainer paths = FindPathBetween(g, max_dist);
        std::stringstream str1;
        str1 << "\t\tPaths size " << paths.size() << ", their lengths : ";
        std::set<size_t> lengths;
        for (size_t i = 0; i < paths.size(); ++i) {
            if (lengths.find(paths[i].Length()) == lengths.end()) {
                //str << paths[i].Length() << " ";
                lengths.insert(paths[i].Length());
            }
        }DEBUG(str1.str());
        set<size_t> different_length;
        for (auto iter = lengths.begin(); iter != lengths.end(); ++iter) {
            std::stringstream str;
            size_t curr_dist = *iter;
            str << "\t\t\tdist " << curr_dist << " : ";
            std::set<PathsEdgesPairInfo> ideal_pi = GenerateIdealPairs(
                    g, curr_dist, max_dist, min_dist);
            size_t ideal_supported_pairs = ideal_pi.size();
            size_t supported = 0;
            std::set<size_t> ideal_supported_edges;
            std::set<size_t> edges_supported;
            for (auto pi = ideal_pi.begin(); pi != ideal_pi.end(); ++pi) {
                IndexPair index(pi->index1(), pi->index2());
                double weight = 0.;
                ideal_supported_edges.insert(pi->index1());
                if (pair_edges_.find(index) != pair_edges_.end()) {
                    weight = pair_edges_.find(index)->second.w();
                    if (weight > 2) {
                        supported += 1;
                        edges_supported.insert(pi->index1());
                    }
                }
                //str << " (" << pi->index1() << ", " << pi->index2() << " , "
                //      << pi->w() << ", " << weight << "), ";
            }

            if (supported > 0
                    && different_length.find(ideal_supported_pairs)
                            == different_length.end()) {
                different_length.insert(ideal_supported_pairs);
                str << "paires should be supported " << ideal_supported_pairs
                        << " supported " << supported
                        << " edges should be supported "
                        << ideal_supported_edges.size() << " supported "
                        << edges_supported.size();
                DEBUG(str.str());
            }

        }
    }

    void Print() const {
        std::stringstream str;
        for (auto iter = pair_edges_.begin(); iter != pair_edges_.end();
                ++iter) {
            str << "( " << iter->first.index1_ << "("
                    << path1_->graph().length(path1_->At(iter->first.index1_))
                    << ")" << ", " << iter->first.index2_ << "("
                    << path1_->graph().length(path2_->At(iter->first.index2_))
                    << ")" << "), ";
        }DEBUG("\t\tSuppoting edges : " << str.str());
    }
private:
    void ConvertPaths(const Graph& g,
                      const PathStorageCallback<Graph>& callback,
                      PathContainer& result) const {
        for (size_t i = 0; i < callback.size(); ++i) {
            const vector<EdgeId>& path = callback.paths()[i];
            BidirectionalPath* our_path = new BidirectionalPath(g, path);
            BidirectionalPath* conj_our_path = new BidirectionalPath(
                    our_path->Conjugate());
            result.AddPair(our_path, conj_our_path);
        }
    }

    std::set<PathsEdgesPairInfo> GenerateIdealPairs(const Graph& g, size_t dist,
                                                    size_t max_dist, size_t min_dist) const {
        std::set<PathsEdgesPairInfo> result;
        for (size_t i1 = 0; i1 < path1_->Size(); ++i1) {
            for (size_t i2 = 0; i2 < path2_->Size(); ++i2) {
                double dist2 = 0;
                if (i2 > 0){
                    dist2 = path2_->Length() -  path2_->LengthAt(i2) ;
                }
                double curr_dist =
                        dist + path1_->LengthAt(i1) + dist2;
                if (max_dist + g.length(path2_->At(i2)) > curr_dist && curr_dist + g.length(path1_->At(i1)) > min_dist) {
                    double weight = std::min(g.length(path1_->At(i1)),
                                             g.length(path2_->At(i2)));
                    weight = std::min(weight, (double) max_dist - curr_dist);
                    PathsEdgesPairInfo pi(i1, i2, weight, curr_dist);
                    result.insert(pi);
                }
            }
        }
        return result;
    }
    PathContainer FindPathBetween(const Graph& g, size_t max_dist) const {
        PathContainer result;
        using namespace omnigraph;
        PathStorageCallback<Graph> path_store(g);
        VertexId begin = g.EdgeEnd(path1_->At(path1_->Size() - 1));
        VertexId end = g.EdgeStart(path2_->At(0));
        PathProcessor<Graph> path_processor(g, 0, max_dist, begin, end,
                                            path_store);
        path_processor.Process();
        ConvertPaths(g, path_store, result);
        return result;
    }

    const BidirectionalPath* path1_;
    const BidirectionalPath* path2_;
    double weight_;
    double dist_;  //distance between the end of the first path and the begin of the second path
    std::map<IndexPair, PathsEdgesPairInfo> pair_edges_;
};

class PathAllPairInfo {
public:
    PathAllPairInfo(const BidirectionalPath* path, size_t max_dist, size_t min_dist)
            : path_(path),
              max_dist_(max_dist),
              min_dist_(min_dist){
    }

    ~PathAllPairInfo() {
        for (auto iter = neighbors_.begin(); iter != neighbors_.end(); ++iter) {
            delete iter->second;
        }
        neighbors_.clear();
    }
    void AddPairInfo(const BidirectionalPath* neighbor, size_t index1,
                     size_t index2, double w, double d) {
        if (neighbors_.find(neighbor) == neighbors_.end()) {
            neighbors_[neighbor] = new PairPathsInfo(path_, neighbor);
        }
        neighbors_[neighbor]->AddPairInfo(index1, index2, w, d);
    }
    void Print() {
        DEBUG("count of neighbors " << neighbors_.size());
        size_t good_size = 0;
        for (auto iter = neighbors_.begin(); iter != neighbors_.end(); ++iter) {
            if (iter->second->weight() > 3
                    && iter->second->CanBeValid(max_dist_, min_dist_)
                    && path_->Length() > max_dist_
                     && iter->first->Length() > max_dist_) {
                good_size++;
                DEBUG("\t index " << iter->first->GetId() << " its length "<< iter->first->Length()<<" distance " << iter->second->d() << " weight " << iter->second->weight());
                iter->second->AnalyzeAndPrint(max_dist_, min_dist_);
                //iter->second->Print();
            }
        }

        if (good_size > 0) {
            DEBUG("Extension supported " << good_size);
        }
    }
private:
    const BidirectionalPath* path_;
    size_t max_dist_;
    size_t min_dist_;
    std::map<const BidirectionalPath*, PairPathsInfo*> neighbors_;
};

class PathsPairInfoContainer {
public:
    PathsPairInfoContainer(const debruijn_graph::conj_graph_pack& gp,
                           PathContainer& paths, const PairedInfoLibrary& lib)
            : gp_(gp),
              paths_(paths),
              lib_(lib),
              max_dist_(lib.GetISMax()),
              min_dist_(lib.GetISMin()),
              coverage_map_(gp.g, paths) {

    }
    ~PathsPairInfoContainer() {
        for (auto iter = pair_info_.begin(); iter != pair_info_.end(); ++iter) {
            delete iter->second;
        }
        pair_info_.clear();
    }

    void FillPairInfo() {
        DEBUG("Fill pair info");
        for (auto iter = paths_.begin(); iter != paths_.end(); ++iter) {
            BidirectionalPath* path = iter.get();
            size_t end_length = 0;
            for (int i = (int) path->Size() - 1;
                    i >= 0 && end_length < max_dist_; i--) {
                EdgeId edge1 = path->At(i);
                std::vector<omnigraph::de::PairInfo<EdgeId> > edge_pair_infos = lib_
                        .index_.GetEdgeInfo(edge1);
                for (size_t pi_i = 0; pi_i < edge_pair_infos.size(); ++pi_i) {
                    AnalyzeTwoEdges(path, edge1, i,
                                    edge_pair_infos[pi_i].second,
                                    edge_pair_infos[pi_i]);
                }
                end_length = path->LengthAt(i);
            }
        }
        Print();
    }

private:
    void AnalyzeTwoEdges(BidirectionalPath* path, EdgeId e1, size_t index1,
                         EdgeId e2, const omnigraph::de::PairInfo<EdgeId>& pi) {
        std::set<BidirectionalPath*> paths2 = coverage_map_.GetCoveringPaths(
                e2);
        for (auto iter = paths2.begin(); iter != paths2.end(); ++iter) {
            vector<size_t> positions2 = (*iter)->FindAll(e2);
            for (size_t i = 0; i < positions2.size(); ++i) {
                AnalyzeTwoPaths(path, e1, index1, *iter, e2, positions2[i], pi);
            }
        }

    }

    void AnalyzeTwoPaths(BidirectionalPath* path1, EdgeId e1, size_t index1,
                         BidirectionalPath* path2, EdgeId e2, size_t index2,
                         const omnigraph::de::PairInfo<EdgeId>& pi) {
        if (path1 == path2 || path1->GetConjPath() == path2) {
            return;
        }
        if (pi.d()
                >= path1->LengthAt(index1) + path2->Length()
                        - path2->LengthAt(index2)) {
            if (pair_info_.find(path1) == pair_info_.end()) {
                pair_info_[path1] = new PathAllPairInfo(path1, max_dist_, min_dist_);
            }
            pair_info_[path1]->AddPairInfo(path2, index1, index2, pi.weight(),
                                           pi.d());
        }
    }

    void Print() {
        for (auto iter = pair_info_.begin(); iter != pair_info_.end(); ++iter) {
            DEBUG("path " << iter->first->GetId());
            iter->second->Print();
        }
    }

    const debruijn_graph::conj_graph_pack& gp_;
    PathContainer& paths_;
    const PairedInfoLibrary& lib_;
    size_t max_dist_;
    size_t min_dist_;
    std::map<const BidirectionalPath*, PathAllPairInfo*> pair_info_;
    GraphCoverageMap coverage_map_;
};*/

}  // namespace path_extend

#endif /* MATE_PAIR_SCAFFOLDING_HPP_ */
