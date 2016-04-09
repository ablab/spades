//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "../utils/bulge_utils.hpp"
#include "../utils/element_printers.hpp"

#include "gluing_vertices_definer.hpp"

using namespace debruijn_graph;
using namespace io;

namespace dipspades {

// both values inclusive
// first index points on index of the start edge in subpath
// last index points on index of the last edge in subpath
typedef pair<size_t, size_t> subpath_range;

class SplitResult {
    Graph &graph_;
    vector<EdgeId> path1_;
    vector<EdgeId> path2_;

    bool CheckExtention(vector<EdgeId> &old_path, vector<EdgeId> &new_path){
        if(old_path.size() != 0 && new_path.size() != 0)
            return graph_.EdgeEnd(old_path[old_path.size() - 1]) ==
                    graph_.EdgeStart(new_path[0]);
        return true;
    }

    void ExtendPath(vector<EdgeId> &old_path, vector<EdgeId> new_path){
        VERIFY(CheckExtention(old_path, new_path));
        old_path.insert(old_path.end(), new_path.begin(), new_path.end());
    }

public:
    SplitResult(Graph &graph, vector<EdgeId> path1, vector<EdgeId> path2) :
        graph_(graph), path1_(path1), path2_(path2) { }

    SplitResult(Graph &graph) : graph_(graph) { }

    void ExtendPaths(SplitResult new_results) {
        ExtendPath(path1_, new_results.path1());
        ExtendPath(path2_, new_results.path2());
    }

    vector<EdgeId> path1() { return path1_; }

    vector<EdgeId> path2() { return path2_; }

    bool IsEmpty() { return path1_.size() == 0 || path2_.size() == 0; }
};

class SubpathsSplitter {

    typedef vector<EdgeId> edge_path;

    Graph &graph_;
    pair<edge_path, edge_path> old_paths_;
    pair<vector<size_t>, vector<size_t> > part_lens_;
    pair<size_t, size_t> num_splits_;
    pair<edge_path, edge_path > split_paths_;
    pair<vector<size_t>, vector<size_t> > partlen_split_paths_;
    pair<size_t, size_t> spath_lens_;
    pair<subpath_range, subpath_range> ranges_;

    enum owner { no_owner, first_path, second_path };
    struct vertex_to_split {
        double rel_dist_;     // relative distance of vertex from the start of subpath
        owner owner_path_;     // path-owner of this vertex
        size_t edge_ind_;     // index of edge end of which matches with vertex

        vertex_to_split(double rel_dist, owner owner_path, size_t edge_ind) :
            rel_dist_(rel_dist), owner_path_(owner_path), edge_ind_(edge_ind) { }
    };

    vector<vertex_to_split> vertices_to_split_;

    struct split_find_result {
        EdgeId edge_;
        size_t pos_;
        bool correct_;

        split_find_result() : edge_(), pos_(), correct_(false) { }

        split_find_result(EdgeId edge, size_t pos, bool correct) :
            edge_(edge),
            pos_(pos),
            correct_(correct) { }

        void initialize(EdgeId edge, size_t pos, bool correct){
            edge_ = edge;
            pos_ = pos;
            correct_ = correct;
        }
    };

    void clear(){
        split_paths_.first.clear();
        split_paths_.second.clear();
        partlen_split_paths_.first.clear();
        partlen_split_paths_.second.clear();
        vertices_to_split_.clear();
    }

    edge_path CutSubpath(edge_path path, subpath_range range){
        return edge_path(path.begin() + range.first, path.begin() + range.second + 1);
    }

    size_t DefineStartEdgePosition(size_t num_splits, subpath_range range){
        if(num_splits == 0)
            return range.second;
        return range.first;
    }

    size_t DefineSpathLenBefore(vector<size_t> &part_len, subpath_range range){
        if(range.first == 0)
            return 0;
        return part_len[range.first - 1];
    }

    size_t DefineSubpathLen(vector<size_t> &part_len, subpath_range range){
        return part_len[range.second] - DefineSpathLenBefore(part_len, range);
    }

    void InitializeSplitVectors(size_t num_splits, edge_path &split_path, vector<size_t> &split_lens,
            edge_path old_path, subpath_range range){
        if(num_splits == 0){
            split_path = CutSubpath(old_path, range);
            split_lens = CalculatePathPartLens(graph_, split_path);
            return;
        }
        split_path.push_back(old_path[range.first]);
        split_lens.push_back(graph_.length(old_path[range.first]));
    }

    size_t IndexFirstOppositeEdge(size_t split_index, set<size_t> &processed_indices){
        for(size_t index = split_index + 1; index < vertices_to_split_.size(); index++)
            if(vertices_to_split_[index].owner_path_ != vertices_to_split_[split_index].owner_path_)
                if(processed_indices.find(index) == processed_indices.end())
                    return index;
        return vertices_to_split_.size() - 1;
    }

    bool VerticesMergePossible(size_t ind1, size_t ind2){
        VERIFY(vertices_to_split_[ind1].owner_path_ != vertices_to_split_[ind2].owner_path_);
        // todo replace magic const to config file
        return fabs(vertices_to_split_[ind1].rel_dist_ - vertices_to_split_[ind2].rel_dist_) < .01;
    }

    bool OwnersMatch(size_t ind, owner owner_path){
        return vertices_to_split_[ind].owner_path_ == owner_path;
    }

    pair<size_t, size_t> OrderByPaths(size_t ind1, size_t ind2){
        VERIFY(vertices_to_split_[ind1].owner_path_ != vertices_to_split_[ind2].owner_path_);
        if(OwnersMatch(ind1, first_path))
            return pair<size_t, size_t>(ind1, ind2);
        return pair<size_t, size_t>(ind2, ind1);
    }

    void PerformMerge(pair<size_t, size_t> indexes){
        size_t edge_ind1 = vertices_to_split_[indexes.first].edge_ind_ + 1;
        split_paths_.first.push_back(old_paths_.first[edge_ind1]);
        partlen_split_paths_.first.push_back(partlen_split_paths_.first[partlen_split_paths_.first.size() - 1] +
                graph_.length(old_paths_.first[edge_ind1]));

        size_t edge_ind2 = vertices_to_split_[indexes.second].edge_ind_ + 1;
        split_paths_.second.push_back(old_paths_.second[edge_ind2]);
        partlen_split_paths_.second.push_back(partlen_split_paths_.second[partlen_split_paths_.second.size() - 1] +
                graph_.length(old_paths_.second[edge_ind2]));
    }

    EdgeId get_last_edge_by_owner(owner path_owner){
        if(path_owner == first_path)
            return split_paths_.first[split_paths_.first.size() - 1];
        return split_paths_.second[split_paths_.second.size() - 1];
//        size_t path_index = vertices_to_split_[index].edge_ind_;
//        owner path_owner = vertices_to_split_[index].owner_path_;
//        if(path_owner == first_path)
//            return old_paths_.first[path_index];
//        return old_paths_.second[path_index];
    }

    split_find_result FindSplitPosition(pair<size_t, size_t> indices, size_t oppos_spaths_len,
            vector<size_t> &oppos_pathlen){
        EdgeId edge_to_split = get_last_edge_by_owner(
                vertices_to_split_[indices.second].owner_path_);
        size_t split_pos = size_t(vertices_to_split_[indices.first].rel_dist_ * double(oppos_spaths_len));
        TRACE("Absolute split position " << split_pos);
        TRACE("oppos_pathlen[oppos_pathlen.size() - 2] - " << oppos_pathlen[oppos_pathlen.size() - 2]);
        if(oppos_pathlen.size() != 1 && split_pos >= oppos_pathlen[oppos_pathlen.size() - 2])
            split_pos -= oppos_pathlen[oppos_pathlen.size() - 2];

        if(split_pos == 0) split_pos++;

        TRACE("Edge before split - " << graph_.str(edge_to_split) <<
                ", split pos - " << split_pos);

        return split_find_result(edge_to_split, split_pos, split_pos < graph_.length(edge_to_split));
    }

    void UpdateSplittedPath(edge_path &path, pair<EdgeId, EdgeId> splitted_edges){
        if(path.size() == 0)
            path.push_back(splitted_edges.first);
        else
            path[path.size() - 1] = splitted_edges.first;
        path.push_back(splitted_edges.second);
    }

    void UpdatesplittedPartLens(vector<size_t> &part_lens, pair<EdgeId, EdgeId> splitted_edges){
        if(part_lens.size() == 0)
            part_lens.push_back(graph_.length(splitted_edges.first));
        else if(part_lens.size() == 1)
            part_lens[0] = graph_.length(splitted_edges.first);
        else
            part_lens[part_lens.size() - 1] = part_lens[part_lens.size() - 2] +
            graph_.length(splitted_edges.first);
        part_lens.push_back(part_lens[part_lens.size() - 1] + graph_.length(splitted_edges.second));
    }

    void SplitOppositeEdge(split_find_result split_res, edge_path &oppos_path,
            vector<size_t> &oppos_partlen){
        if(!split_res.correct_ || graph_.length(split_res.edge_) < split_res.pos_)
            return;
        pair<EdgeId, EdgeId> splitted_edges = graph_.SplitEdge(split_res.edge_, split_res.pos_);
        TRACE("Edges after split - " << graph_.str(splitted_edges.first) << " " <<
                graph_.str(splitted_edges.second));
        UpdateSplittedPath(oppos_path, splitted_edges);
        UpdatesplittedPartLens(oppos_partlen, splitted_edges);
    }

    // first from pairs - splitted, second - opposite
    bool PerformSplit(pair<size_t, size_t> indexes,
            pair<edge_path&, edge_path& > split_paths,
            pair<vector<size_t>&, vector<size_t>& > split_partlens,
            pair<edge_path&, edge_path& > default_paths,
            pair<size_t, size_t> spaths_len,
            pair<size_t, size_t> num_splits){

        TRACE("New path1 before: " << SimplePathWithVerticesToString(graph_, split_paths.first));
        TRACE("New path2 before: " << SimplePathWithVerticesToString(graph_, split_paths.second));

        TRACE("FindEdgeAndSplitPosition");
        split_find_result split_res = FindSplitPosition(indexes, spaths_len.second,
                split_partlens.second);

        if(!split_res.correct_){
            TRACE("Split was not performed");
            return false;
        }

        TRACE("SplitOppositeEdge");
        SplitOppositeEdge(split_res, split_paths.second, split_partlens.second);

        // update non splitted path
        TRACE("Update non splitted path");
        if(num_splits.second != 0){
            size_t edge_ind = vertices_to_split_[indexes.first].edge_ind_ + 1;
            split_paths.first.push_back(default_paths.first[edge_ind]);
            split_partlens.first.push_back(split_partlens.first[split_partlens.first.size() - 1] +
                    graph_.length(default_paths.first[edge_ind]));
        }

        TRACE("New path1 after: " << SimplePathWithVerticesToString(graph_, split_paths.first));
        TRACE("New path2 after: " << SimplePathWithVerticesToString(graph_, split_paths.second));

        return true;
    }

    // function expect that order in pair_to_order matches with (first_path, second_path)
    template<typename T>
    pair<T&, T&> OrderBySplitAndOpposite(size_t split_ind, size_t oppos_ind, pair<T, T> &pair_to_order){
        VERIFY(vertices_to_split_[split_ind].owner_path_ !=
                vertices_to_split_[oppos_ind].owner_path_);
        if(vertices_to_split_[split_ind].owner_path_ == first_path)
            return pair<T&, T&>(pair_to_order.first, pair_to_order.second);
        return pair<T&, T&>(pair_to_order.second, pair_to_order.first);
    }

    bool PerformSplitting(subpath_range range1, subpath_range range2){
        TRACE("Vector initialization");
        InitializeSplitVectors(num_splits_.second, split_paths_.first,
                partlen_split_paths_.first,  old_paths_.first, range1);
        InitializeSplitVectors(num_splits_.first, split_paths_.second,
                partlen_split_paths_.second, old_paths_.second, range2);

        size_t num_done_splits = 0;
        size_t split_index = 1;

        TRACE("Splitting cycle starts");

        set<size_t> processed_indices;
        while(num_done_splits < num_splits_.first + num_splits_.second){
            TRACE("Splitted index - " << split_index << " , owner - " <<
                    vertices_to_split_[split_index].owner_path_);

            size_t opposite_index = IndexFirstOppositeEdge(split_index, processed_indices);
            TRACE("Opposite index - " << opposite_index << ", owner - " <<
                    vertices_to_split_[opposite_index].owner_path_);

            if(processed_indices.find(split_index) == processed_indices.end()){
                if(VerticesMergePossible(split_index, opposite_index) &&
                        (opposite_index != vertices_to_split_.size() - 2) &&
                                (opposite_index != vertices_to_split_.size() - 1)){

                    TRACE("Merge starts");
                    PerformMerge(OrderByPaths(split_index, opposite_index));
                    num_done_splits += 2;
                    processed_indices.insert(opposite_index);

                    TRACE("Merge was performed");
                }
                else{
                    TRACE("Split starts");

                    bool split_res = PerformSplit(pair<size_t, size_t>(split_index, opposite_index),
                            OrderBySplitAndOpposite<edge_path>(split_index, opposite_index, split_paths_),
                            OrderBySplitAndOpposite<vector<size_t> >(split_index, opposite_index, partlen_split_paths_),
                            OrderBySplitAndOpposite<edge_path>(split_index, opposite_index, old_paths_),
                            OrderBySplitAndOpposite<size_t>(split_index, opposite_index, spath_lens_),
                            OrderBySplitAndOpposite<size_t>(split_index, opposite_index, num_splits_));

                    if(!split_res)
                        return false;

                    num_done_splits++;
                    TRACE("Split was performed");
                }

                processed_indices.insert(split_index);
            }
            TRACE("Number done splittings - " << num_done_splits);
            split_index ++;
            TRACE("-------------------------");
        }
        TRACE("Splitting cycle ends");
        TRACE("-------------------------");
        return true;
    }

    void CreateVectorSplitVertices(subpath_range range1, subpath_range range2){
        pair<size_t, size_t> lens_before_spath(DefineSpathLenBefore(part_lens_.first, range1),
                DefineSpathLenBefore(part_lens_.second, range2));
        spath_lens_ = pair<size_t, size_t>(DefineSubpathLen(part_lens_.first, range1),
                DefineSubpathLen(part_lens_.second, range2));
        vertices_to_split_.push_back(vertex_to_split(0, no_owner, 0));
        size_t iter1 = range1.first;
        size_t iter2 = range2.first;

        TRACE("Partlens for 1st vector: " << VectorToString<size_t>(part_lens_.first));
        TRACE("Partlens for 2nd vector: " << VectorToString<size_t>(part_lens_.second));
        TRACE("Slens before - " << lens_before_spath.first << " " << lens_before_spath.second);

        for(size_t i = 0; i < num_splits_.first + num_splits_.second; i++){
            double rel_dist1 = double(part_lens_.first[iter1] - lens_before_spath.first) /
                    double(spath_lens_.first);
            double rel_dist2 = double(part_lens_.second[iter2] - lens_before_spath.second) /
                    double(spath_lens_.second);
            if(rel_dist1 < rel_dist2){
                vertices_to_split_.push_back(vertex_to_split(rel_dist1, first_path, iter1));
                iter1++;
            }
            else{
                vertices_to_split_.push_back(vertex_to_split(rel_dist2, second_path, iter2));
                iter2++;
            }
        }
        vertices_to_split_.push_back(vertex_to_split(1.0, second_path, iter2));
        vertices_to_split_.push_back(vertex_to_split(1.0, first_path, iter1));
    }

public:
    SubpathsSplitter(Graph &graph, shared_ptr<BaseBulge> bulge) :
        graph_(graph),
        old_paths_(pair<edge_path, edge_path>(bulge->path1(), bulge->path2())),
        part_lens_(make_pair(CalculatePathPartLens(graph_, old_paths_.first),
                CalculatePathPartLens(graph_, old_paths_.second))),
        num_splits_(),
        split_paths_(),
        partlen_split_paths_() { }

    SplitResult SplitSubpaths(subpath_range range1, subpath_range range2) {
        clear();

        // number of splits on the 1st and the 2nd subpaths
        ranges_.first = range1;
        ranges_.second = range2;

        num_splits_.first = range1.second - range1.first;
        num_splits_.second = range2.second - range2.first;

        TRACE("Range 1: " << range1.first << " - " << range1.second);
        TRACE("Range 2: " << range2.first << " - " << range2.second);
        TRACE("Num splits 1 - " << num_splits_.first << ", num splits 2 - " << num_splits_.second);

        TRACE("Subpath to split1 - " << SimplePathWithVerticesToString(graph_, CutSubpath(old_paths_.first, range1)));
        TRACE("Subpath to split2 - " << SimplePathWithVerticesToString(graph_, CutSubpath(old_paths_.second, range2)));

        if(num_splits_.first + num_splits_.second == 0)
            return SplitResult(graph_, CutSubpath(old_paths_.first, range1),
                    CutSubpath(old_paths_.second, range2));

        CreateVectorSplitVertices(range1, range2);

        TRACE("Vertices to split:");
        for(auto it = vertices_to_split_.begin(); it != vertices_to_split_.end(); it++)
            TRACE(it->rel_dist_ << " " << it->owner_path_ << " " << it->edge_ind_ );

        TRACE("Auxiliary vectors were created");

        if(!PerformSplitting(range1, range2))
            return SplitResult(graph_);

        TRACE("Splitted spath1 - " << SimplePathWithVerticesToString(graph_, split_paths_.first));
        TRACE("Splitted spath2 - " << SimplePathWithVerticesToString(graph_, split_paths_.second));
        return SplitResult(graph_, split_paths_.first, split_paths_.second);
    }

private:
    DECL_LOGGER("SubpathSplitter");
};

class BulgeSplitter {
    Graph &graph_;
public:
    BulgeSplitter(Graph &graph) : graph_(graph) { }

    shared_ptr<BaseBulge> SplitBulge(shared_ptr<BaseBulge> bulge, GluingVericesDefinerResults gluing_def_results) {
        if(bulge->IsSimple()){
            TRACE("Bulge is simple. Splitting was not performed");
            return shared_ptr<BaseBulge>(new Bulge(graph_, graph_.k(), bulge->path1(), bulge->path2()));
        }

        SubpathsSplitter spaths_splitter(graph_, bulge);
        if(gluing_def_results.size() == 0){
            // one big split
            TRACE("No gluing vertices. Split will perform between start and end vertices");
            auto split_res = spaths_splitter.SplitSubpaths(
                    subpath_range(0, bulge->path1().size() - 1),
                    subpath_range(0, bulge->path2().size() - 1));
            TRACE("bulge was splitted");
            TRACE("1st new bulge side - " << SimplePathWithVerticesToString(graph_, split_res.path1()));
            TRACE("2nd new bulge side - " << SimplePathWithVerticesToString(graph_, split_res.path2()));
            return shared_ptr<BaseBulge>(new Bulge(graph_, graph_.k(), split_res.path1(), split_res.path2()));
        }
        TRACE(gluing_def_results.size() << " - number of gluing pairs");
        // splitting before first gluing pair
        TRACE("Splitting before first gluing pair");
        auto split_result = spaths_splitter.SplitSubpaths(
                subpath_range(0, gluing_def_results.begin()->first),
                subpath_range(0, gluing_def_results.begin()->second));

        if(split_result.IsEmpty())
            return shared_ptr<BaseBulge>(new Bulge(graph_));

        // perform all intermediate splittings
        TRACE("All intermediate splittings");
        for(auto iter1 = gluing_def_results.begin(), iter2 = ++gluing_def_results.begin();
                iter2 != gluing_def_results.end(); iter1++, iter2++){
            TRACE("Gluing pairs - (" << iter1->first << " " << iter1->second << ") (" <<
                    iter2->first << " " << iter2->second << ")");
            auto new_split_res = spaths_splitter.SplitSubpaths(
                    subpath_range(iter1->first + 1, iter2->first),
                    subpath_range(iter1->second + 1, iter2->second));
            if(new_split_res.IsEmpty())
                return shared_ptr<BaseBulge>(new Bulge(graph_));
            split_result.ExtendPaths(new_split_res);
        }

        // splitting after last gluing last pair
        TRACE("Splitting after last gluing last pair");
        auto last_split_res = spaths_splitter.SplitSubpaths(
                subpath_range((--gluing_def_results.end())->first + 1, bulge->path1().size() - 1),
                subpath_range((--gluing_def_results.end())->second + 1, bulge->path2().size() - 1));
        if(last_split_res.IsEmpty())
            return shared_ptr<BaseBulge>(new Bulge(graph_));
        split_result.ExtendPaths(last_split_res);

        TRACE("New bulge path1 - " << SimplePathWithVerticesToString(graph_, split_result.path1()));
        TRACE("New bulge path2 - " << SimplePathWithVerticesToString(graph_, split_result.path2()));
        TRACE("Splitting completed");

        return shared_ptr<BaseBulge>(new Bulge(graph_, graph_.k(), split_result.path1(), split_result.path2()));
    }

private:
    DECL_LOGGER("BulgeSplitter");
};

}
