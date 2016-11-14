//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "coloring.hpp"

#include <vector>
#include <map>
#include <common/visualization/graph_labeler.hpp>
#include "common/adt/bag.hpp"

namespace cap {

template<class Graph>
class GenomePath: public GraphActionHandler<Graph> {
    typedef GraphActionHandler<Graph> base;
public:
    typedef typename Graph::EdgeId EdgeId;
    typedef typename vector<EdgeId>::const_iterator iterator;
private:
    vector<EdgeId> path_;
    bool cyclic_;
    bag<EdgeId> edge_mult_;

    int Idx(size_t idx) {
        if (cyclic_) {
            return idx % path_.size();
        }
        if (idx >= path_.size())
            return -1;
        return idx;
    }

    bool CheckAllMatch(const vector<EdgeId>& edges, size_t pos) {
        for (size_t j = 0; j < edges.size(); ++j) {
            int idx = Idx(pos + j);
            if (idx > 0) {
                if (path_[idx] != edges[j])
                    return false;
            } else {
                return false;
            }
        }
        return true;
    }

    bool CheckNoneMatch(const vector<EdgeId>& edges, size_t pos) {
        for (size_t j = 0; j < edges.size(); ++j) {
            int idx = Idx(pos + j);
            if (idx > 0) {
                if (path_[idx] == edges[j])
                    return false;
            } else {
                return true;
            }
        }
        return true;
    }

    bool CheckConsistency(const vector<EdgeId>& edges, size_t pos) {
        return CheckAllMatch(edges, pos) || CheckNoneMatch(edges, pos);
//        if (Idx(pos) < 0)
//            return true;
//        bool first_matched = (path_[idx] == edges[0]);
//        for (size_t j = 0; j < edges.size(); ++j) {
//            int idx = Idx(pos + j);
//            if (idx > 0) {
//                if (first_matched ^ (path_[idx] == edges[j]) != 0)
//                    return false;
//            } else {
//                return !first_matched;
//            }
//        }
//        return true;
    }

    bool CheckConsistency(const vector<EdgeId>& edges) {
        VERIFY(!edges.empty());
        size_t mult = edge_mult_.mult(edges[0]);
        DEBUG("Mult of " << this->g().str(edges[0]) << " is " << mult);
        for (size_t i = 1; i < edges.size(); ++i) {
            DEBUG(
                    "Mult of " << this->g().str(edges[i]) << " is " << edge_mult_.mult(edges[i]));
            if (!CheckConsistency(edges, i)
                    || edge_mult_.mult(edges[i]) != mult) {
                return false;
            }
        }
        return true;
    }

    void MovePrefixBack(size_t prefix_length) {
        VERIFY(cyclic_);
        vector<EdgeId> tmp(path_.begin(), path_.begin() + prefix_length);
        path_.erase(path_.begin(), path_.begin() + prefix_length);
        path_.insert(path_.end(), tmp.begin(), tmp.end());
    }

    void SubstituteNonBorderFragment(size_t start_pos, size_t end_pos,
            const vector<EdgeId>& subst) {
        VERIFY(start_pos < end_pos && end_pos <= path_.size());
        ChangeMult(start_pos, end_pos, subst);
        path_.insert(
                path_.erase(path_.begin() + start_pos, path_.begin() + end_pos),
                subst.begin(), subst.end());
    }

    void FillEdgeMult() {
        for (auto it = path_.begin(); it != path_.end(); ++it) {
            DEBUG("Edge " << this->g().str(*it) << " is genomic")
            edge_mult_.put(*it);
//            edge_mult_.put(this->g().conjugate(*it));
        }
    }

    void ChangeMult(size_t start_pos, size_t end_pos,
            const vector<EdgeId>& subst) {
        for (size_t i = start_pos; i < end_pos; ++i) {
            bool could_take = edge_mult_.take(path_[i]);
            VERIFY(could_take);
        }
        for (auto it = subst.begin(); it != subst.end(); ++it) {
            edge_mult_.put(*it);
        }
    }

public:
    GenomePath(const Graph& g, const vector<EdgeId>& path, bool cyclic = false) :
            base(g, "GenomePath"), path_(path), cyclic_(cyclic) {
        FillEdgeMult();
    }

    /*virtual*/
    void HandleMerge(const vector<EdgeId>& old_edges, EdgeId new_edge) {
//        DEBUG(
//                "Handling merge of edges " << this->g().str(old_edges) << " into edge " << this->g().str(new_edge));
        VERIFY(CheckConsistency(old_edges));
//        DEBUG("Path before: " << this->g().str(path_));
        auto it = find(path_.begin(), path_.end(), old_edges.front());
        while (it != path_.end()) {
            size_t start = it - path_.begin();
            size_t end = start + old_edges.size();
            Substitute(start, end, vector<EdgeId> { new_edge });
//            DEBUG("Path after find: " << this->g().str(path_));
            it = find(path_.begin(), path_.end(), old_edges.front());
        }
        //debug
        for (auto it2 = old_edges.begin(); it2 != old_edges.end(); ++it2) {
//            DEBUG("Checking " << this->g().str(*it2))
            VERIFY(find(path_.begin(), path_.end(), *it2) == path_.end());
            VERIFY(edge_mult_.mult(*it2) == 0);
        }
        //debug
//        DEBUG("Path final: " << this->g().str(path_));
        DEBUG("Merge handled");
    }

    /*virtual*/
    void HandleDelete(EdgeId e) {
        DEBUG(
                "Multiplicity of edge " << this->g().str(e) << " in delete " << edge_mult_.mult(e));
        VERIFY(edge_mult_.mult(e) == 0);
    }

    //for cyclic paths, end_pos might be > path.size()
    //might change indices unexpectedly
    void Substitute(size_t start_pos, size_t end_pos,
            const vector<EdgeId>& subst) {
        DEBUG("Substitute called");
        VERIFY(cyclic_ || end_pos <= path_.size());
        if (end_pos <= path_.size()) {
            SubstituteNonBorderFragment(start_pos, end_pos, subst);
        } else {
            size_t prefix_length = end_pos - path_.size();
            VERIFY(start_pos >= prefix_length);
            MovePrefixBack(prefix_length);
            SubstituteNonBorderFragment(start_pos - prefix_length, path_.size(),
                    subst);
        }
    }

    size_t size() const {
        return path_.size();
    }

    iterator begin() const {
        return path_.begin();
    }

    iterator end() const {
        return path_.end();
    }

    size_t mult(EdgeId e) {
        return edge_mult_.mult(e);
    }

private:
    DECL_LOGGER("GenomePath")
    ;
};

template<class Graph>
class AssemblyPathCallback: public PathProcessor<Graph>::Callback {
    typedef typename Graph::EdgeId EdgeId;
    typedef typename std::vector<EdgeId> Path;

private:
    const Graph& g_;
    const ColorHandler<Graph>& coloring_;
    const TColorSet assembly_color_;
    size_t edge_count_;

    std::vector<Path> paths_;

    bool CheckPath(const vector<EdgeId>& path) const {
        DEBUG("Checking path " << g_.str(path));
        if (path.size() > edge_count_) {
            DEBUG("false");
            return false;
        }
        for (auto it = path.begin(); it != path.end(); ++it) {
            if ((coloring_.Color(*it) & assembly_color_) == kEmptyColorSet) {
                DEBUG("false");
                return false;
            }
        }
        DEBUG("true");
        return true;
    }

public:
    AssemblyPathCallback(const Graph& g, const ColorHandler<Graph>& coloring,
            TColorSet assembly_color, size_t edge_count) :
            g_(g), coloring_(coloring), assembly_color_(assembly_color), edge_count_(
                    edge_count) {
    }

    virtual void HandleReversedPath(const Path& rev_path) {
        Path path = this->ReversePath(rev_path);
        if (CheckPath(path)) {
            paths_.push_back(path);
        }
    }

    size_t size() const {
        return paths_.size();
    }

    vector<Path> paths() const {
        return paths_;
    }
};

template<class Graph>
class SimpleInDelCorrector {
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;
    Graph& g_;
    ColorHandler<Graph>& coloring_;

    //become invalidated during process
//    const EdgesPositionHandler<Graph>& edge_pos_;
    GenomePath<Graph> genome_path_;
    const TColorSet genome_color_;
    const TColorSet assembly_color_;

    vector<EdgeId> FindAssemblyPath(VertexId start, VertexId end,
            size_t edge_count_bound, size_t min_length, size_t max_length) {
        AssemblyPathCallback<Graph> assembly_callback(g_, coloring_,
                assembly_color_, edge_count_bound);
        PathProcessor<Graph> path_finder(g_, min_length, max_length, start, end,
                assembly_callback);
        path_finder.Process();
        if (assembly_callback.size() == 0) {
            DEBUG("Couldn't find assembly path");
        } else if (assembly_callback.size() > 1) {
            DEBUG("Found several assembly paths");
            DEBUG("Taking first");
            return assembly_callback.paths().front();
        } else {
            DEBUG("Found unique assembly path");
            return assembly_callback.paths().front();
        }
        return {};
    }

    int TryFindGenomePath(size_t pos, VertexId end, size_t edge_count_bound) {
        for (size_t i = 0;
                i + pos < genome_path_.size() && i < edge_count_bound; ++i) {
            if (g_.EdgeEnd(*(genome_path_.begin() + pos + i)) == end) {
                return pos + i + 1;
            }
        }
        return -1;
    }

//    bag<TColorSet> ColorLengths(const vector<EdgeId>& edges) {
//        bag<TColorSet> answer;
//        for (size_t i = 0; i < edges.size(); ++i) {
//            answer.put(coloring_.Color(edges[i]), g_.length(edges[i]));
//        }
//        return answer;
//    }

    size_t VioletLengthOfGenomeUnique(const vector<EdgeId>& edges) {
        size_t answer = 0;
        for (size_t i = 0; i < edges.size(); ++i) {
            // TODO make reference, do not copy!!!
            TColorSet edge_color_set = coloring_.Color(edges[i]);
            if (edge_color_set[0] && edge_color_set[1]
                    && (genome_path_.mult(edges[i]) == 1)) {
                answer += g_.length(edges[i]);
            }
        }
        return answer;
    }

    //genome pos exclusive
//    size_t CumulativeGenomeLengthToPos(size_t pos) {
//        size_t answer = 0;
//        for (size_t i = 0; i < pos; ++i) {
//            answer += g_.length(genome_path_[i]);
//        }
//        return answer;
//    }

    bool CheckGenomePath(size_t genome_start, size_t genome_end) {
        return VioletLengthOfGenomeUnique(
                vector<EdgeId>(genome_path_.begin() + genome_start,
                        genome_path_.begin() + genome_end)) < 25;
    }

    optional<pair<size_t, size_t>> FindGenomePath(VertexId start, VertexId end,
            size_t edge_count_bound) {
        for (size_t i = 0; i < genome_path_.size(); ++i) {
            if (g_.EdgeStart(*(genome_path_.begin() + i)) == start) {
                int path_end = TryFindGenomePath(i, end, edge_count_bound);
                if (path_end > 0 && CheckGenomePath(i, path_end))
                    return make_optional(make_pair(size_t(i), size_t(path_end)));
            }
        }
        return boost::none;
    }

    void RemoveObsoleteEdges(const vector<EdgeId>& edges) {
        for (auto it = SmartSetIterator<Graph, EdgeId>(g_, edges.begin(),
                edges.end()); !it.IsEnd(); ++it) {
            if (coloring_.Color(*it) == genome_color_
                    && genome_path_.mult(*it) == 0
                    && genome_path_.mult(g_.conjugate(*it)) == 0) {
                DEBUG("Removing edge " << g_.str(*it) << " as obsolete");
                VertexId start = g_.EdgeStart(*it);
                VertexId end = g_.EdgeEnd(*it);
                g_.DeleteEdge(*it);
                DEBUG("Comressing start");
                g_.CompressVertex(start);
                if (!g_.RelatedVertices(start, end)) {
                    DEBUG("Comressing end");
                    g_.CompressVertex(end);
                }
                DEBUG("Edge removed");
            }
        }
    }

    string GenomePathStr(size_t genome_start, size_t genome_end) const {
        return g_.str(
                vector<EdgeId>(genome_path_.begin() + genome_start,
                        genome_path_.begin() + genome_end));
    }

    void GenPicAlongPath(const vector<EdgeId> path, size_t cnt) {
    utils::MakeDirPath("ref_correction");
        WriteComponentsAlongPath(g_, visualization::graph_labeler::StrGraphLabeler<Graph>(g_),
                "ref_correction/" + ToString(cnt) + ".dot", 100000, 10,
                TrivialMappingPath(g_, path), *ConstructColorer(coloring_));
    }

    void GenPicAroundEdge(EdgeId e, size_t cnt) {
        utils::MakeDirPath("ref_correction");
        GraphComponent<Graph> component = omnigraph::EdgeNeighborhood(g_, e, 10, 100000);
        visualization::visualization_utils::WriteComponent(g_, "ref_correction/" + ToString(cnt) + ".dot", component, coloring_.GetInstance(),
                                      visualization::graph_labeler::StrGraphLabeler<Graph>(g_));
    }

    void CorrectGenomePath(size_t genome_start, size_t genome_end,
            const vector<EdgeId>& assembly_path) {
        static size_t cnt = 0;
        DEBUG(
                "Case " << ++cnt << " Substituting genome path " << GenomePathStr(genome_start, genome_end) << " with assembly path " << g_.str(assembly_path));
        vector<EdgeId> genomic_edges;
        for (size_t i = genome_start; i < genome_end; ++i) {
            genomic_edges.push_back(*(genome_path_.begin() + i));
        }
        GenPicAlongPath(genomic_edges, cnt * 100);
        GenPicAlongPath(assembly_path, cnt * 100 + 1);
        for (size_t i = 0; i < assembly_path.size(); ++i) {
            coloring_.PaintEdge(assembly_path[i], genome_color_);
        }
        genome_path_.Substitute(genome_start, genome_end, assembly_path);
        RemoveObsoleteEdges(genomic_edges);
        GenPicAroundEdge(
                *((genome_start < genome_path_.size()) ?
                        (genome_path_.begin() + genome_start) :
                        genome_path_.end() - 1), cnt * 100 + 2);
    }

//    pair<string, pair<size_t, size_t>> ContigIdAndPositions(EdgeId e) {
//        vector<EdgePosition> poss = edge_pos_.GetEdgePositions(e);
//        VERIFY(!poss.empty());
//        if (poss.size() > 1) {
//            WARN("Something strange with assembly positions");
//            return make_pair("", make_pair(0, 0));
//        }
//        EdgePosition pos = poss.front();
//        return make_pair(pos.contigId_, make_pair(pos.start(), pos.end()));
//    }

//    void WriteAltPath(EdgeId e, const vector<EdgeId>& genome_path) {
//        LengthIdGraphLabeler<Graph> basic_labeler(g_);
//        EdgePosGraphLabeler<Graph> pos_labeler(g_, edge_pos_);
//
//        CompositeLabeler<Graph> labeler(basic_labeler, pos_labeler);
//
//        string alt_path_folder = folder_ + ToString(g_.int_id(e)) + "/";
//        make_dir(alt_path_folder);
//        WriteComponentsAlongPath(g_, labeler, alt_path_folder + "path.dot", /*split_length*/
//        1000, /*vertex_number*/15, TrivialMappingPath(g_, genome_path),
//                *ConstructBorderColorer(g_, coloring_));
//    }

//todo use contig constraints here!!!
    void AnalyzeGenomeEdge(EdgeId e) {
        DEBUG("Analysing shortcut genome edge " << g_.str(e));
//        VERIFY(genome_path_.mult(e) > 0);
        DEBUG("Multiplicity " << genome_path_.mult(e));
        if (genome_path_.mult(e) == 1) {
            vector<EdgeId> assembly_path = FindAssemblyPath(g_.EdgeStart(e),
                    g_.EdgeEnd(e), 100, 0, g_.length(e) + 1000);
            if (!assembly_path.empty()) {
                DEBUG("Assembly path " << g_.str(assembly_path));
                auto it = std::find(genome_path_.begin(), genome_path_.end(),
                        e);
                VERIFY(it != genome_path_.end());
                size_t pos = it - genome_path_.begin();
                CorrectGenomePath(pos, pos + 1, assembly_path);
            } else {
                DEBUG("Couldn't find assembly path");
            }
        }
    }

    void AnalyzeAssemblyEdge(EdgeId e) {
        DEBUG("Analysing shortcut assembly edge " << g_.str(e));
        optional < pair < size_t, size_t >> genome_path = FindGenomePath(
                g_.EdgeStart(e), g_.EdgeEnd(e), /*edge count bound*//*100*/
                300);
        if (genome_path) {
            CorrectGenomePath(genome_path->first, genome_path->second,
                    vector<EdgeId> { e });
        } else {
            DEBUG("Empty genome path");
        }
    }

public:
    SimpleInDelCorrector(Graph& g, ColorHandler<Graph>& coloring,
            const vector<EdgeId>& genome_path, TColorSet genome_color,
            TColorSet assembly_color) :
            g_(g), coloring_(coloring), genome_path_(g_, genome_path), genome_color_(
                    genome_color), assembly_color_(assembly_color) {
    }

    void Analyze() {
        //remove_dir("ref_correction");
        for (auto it = g_.SmartEdgeBegin(); !it.IsEnd(); ++it) {
            if (coloring_.Color(*it) == genome_color_
                    && genome_path_.mult(*it) > 0) {
                AnalyzeGenomeEdge(*it);
            }
            if (coloring_.Color(*it) == assembly_color_) {
                AnalyzeAssemblyEdge(*it);
            }
        }
    }

private:
    DECL_LOGGER("SimpleInDelCorrector")
    ;
};

}
