//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* Copyright (c) 2019 University of Warwick
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "stop_condon_finder.hpp"
#include "assembly_graph/core/graph.hpp"
#include "assembly_graph/paths/path_utils.hpp"
#include "assembly_graph/components/graph_component.hpp"
#include "modules/alignment/sequence_mapper.hpp"
#include "io/utils/edge_namer.hpp"

namespace cds_subgraphs {

using namespace debruijn_graph;

static size_t RoundedProduct(size_t l, double coeff) {
    return size_t(math::round(coeff * double(l)));
}

//template<class Graph>
//class EdgeTrackingCallback : public omnigraph::PathProcessor<Graph>::Callback {
//    typedef typename Graph::EdgeId EdgeId;
//    typedef typename Graph::VertexId VertexId;
//
//public:
//    EdgeTrackingCallback(const Graph &g) :
//            g_(g) {
//    }
//
//    void HandleReversedPath(const std::vector<EdgeId> &path) override {
//        utils::insert_all(edges_, path);
//    }
//
//    const std::set<EdgeId> &edges() const {
//        return edges_;
//    }
//
//private:
//    const Graph &g_;
//    std::set<EdgeId> edges_;
//};

//class PathFindingRelevantComponentFinder {
//    const Graph &g_;
//    const double min_len_frac_;
//    const double max_len_frac_;
//
//public:
//    PathFindingRelevantComponentFinder(const Graph &g, double min_len_frac = 0.5,
//                                       double max_len_frac = 1.5) :
//            g_(g), min_len_frac_(min_len_frac), max_len_frac_(max_len_frac) {}
//
//    GraphComponent<Graph> RelevantComponent(size_t base_len,
//                                            const std::unordered_map<VertexId, size_t> &starts,
//                                            const std::unordered_map<VertexId, size_t> &ends) const {
//        EdgeTrackingCallback<Graph> callback(g_);
//        for (auto s_v_d : starts) {
//            VertexId start = s_v_d.first;
//            DEBUG("Searching paths for source " << g_.str(start));
//            size_t max_length = 0;
//            for (auto e_v_d : ends) {
//                size_t total_length = s_v_d.second + e_v_d.second + base_len;
//                if (total_length > max_length)
//                    max_length = total_length;
//            }
//            max_length = size_t(math::round(max_len_frac_ * double(max_length)));
//            DEBUG("Creating path processor for max path length " << max_length);
//            //todo increase dijkstra limits
//            PathProcessor<Graph> processor(g_, start, max_length);
//            for (auto e_v_d : ends) {
//                VertexId end = e_v_d.first;
//                DEBUG("Finding paths to sink " << g_.str(end));
//                size_t total_length = s_v_d.second + e_v_d.second + base_len;
//                processor.Process(end,
//                                  size_t(math::round(min_len_frac_ * double(total_length))),
//                                  size_t(math::round(max_len_frac_ * double(total_length))),
//                                  callback
//                );
//            }
//        }
//
//        const auto &edges = callback.edges();
//        return GraphComponent<Graph>::FromEdges(g_, edges.begin(), edges.end());
//    }
//
//};

//FIXME consider using edges rather than vertices
class MinDistRelevantComponentFinder {
    //TODO use throughout file
    using DistInfo = std::unordered_map<GraphPos, size_t>;
    using BaseDistF = std::function<size_t(GraphPos, GraphPos)>;
    const Graph &g_;
    const io::EdgeNamingF<Graph> edge_naming_f_;
    const double max_len_coeff_;

    size_t MaxDist(const DistInfo &v_ds) const {
        size_t max = 0;
        for (const auto &v_d : v_ds) {
            if (v_d.second > max)
                max = v_d.second;
        }
        return max;
    }

    //returns min distance to exit among the appropriate paths or max value if none
    size_t MinExitDist(const DistInfo &s_dist, const DistInfo &e_dist, const BaseDistF &base_dist_f) const {
        size_t min_e_dist = std::numeric_limits<size_t>::max();
        //TODO seems like can be simplified
        for (const auto &s_d : s_dist) {
            for (const auto &e_d : e_dist) {
                if (base_dist_f(s_d.first, e_d.first) >= s_d.second + e_d.second) {
                    if (e_d.second < min_e_dist) {
                        min_e_dist = e_d.second;
                    }
                }
            }
        }
        return min_e_dist;
    }

    bool CheckConnectedTo(VertexId v, const std::set<VertexId> &vertices) const {
        for (EdgeId e : g_.IncomingEdges(v))
            if (vertices.count(g_.EdgeStart(e)))
                return true;
        return false;
    }

public:
    //TODO NB: max_len_frac used twice with slightly different meaning
    MinDistRelevantComponentFinder(const Graph &g,
                                   io::EdgeNamingF<Graph> edge_naming_f,
                                   double max_len_frac = 1.5) :
            g_(g), edge_naming_f_(edge_naming_f), max_len_coeff_(max_len_frac) {}

    //TODO check how base_len is defined
    omnigraph::GraphComponent<Graph> RelevantComponent(size_t cds_len_est,
                                                       size_t base_len,
                                                       const DistInfo &starts,
                                                       const DistInfo &ends) const;
};

class PartialGenePathProcessor {
    const Graph &g_;
    const MinDistRelevantComponentFinder rel_comp_finder_;
    const io::EdgeNamingF<Graph> edge_naming_f_;

    std::string PrintEdgePath(const EdgePath &path) const {
        std::stringstream ss;
        std::string delim = "";
        ss << "[";
        for (EdgeId e : path.sequence()) {
            ss << delim << edge_naming_f_(g_, e);
            delim = ", ";
        }
        ss << "], ";
        ss << "start: " << path.start_pos() << ", end: " << path.end_pos();
        return ss.str();
    }

    Sequence PathSeq(const EdgePath &p) const {
        return PathSequence(g_, p);
    }

    EdgePath RCEdgePath(const EdgePath &p) const {
        const auto &edges = p.sequence();
        std::vector<EdgeId> rc_edges(edges.size());
        std::transform(edges.rbegin(), edges.rend(), rc_edges.begin(), [this](EdgeId e) { return g_.conjugate(e); });
        return EdgePath(rc_edges, g_.length(edges.back()) - p.end_pos(), g_.length(edges.front()) - p.start_pos());
    }

    std::vector<EdgePath> RCEdgePaths(const std::vector<EdgePath> &ps) const {
        std::vector<EdgePath> paths(ps.size());
        std::transform(ps.begin(), ps.end(), paths.begin(), [this](const EdgePath &p) { return RCEdgePath(p); });
        return paths;
    }

    GraphPos RCPos(GraphPos p) const {
        return std::make_pair(g_.conjugate(p.first), g_.length(p.first) - p.second);
    }

    std::unordered_map<GraphPos, size_t>
    RCTerminates(const std::unordered_map<GraphPos, size_t> &terminate_dists) const {
        std::unordered_map<GraphPos, size_t> ans;
        for (auto v_d : terminate_dists) {
            ans[RCPos(v_d.first)] = v_d.second;
        }
        return ans;
    };

//    EdgePath CropRight(const EdgePath &path, size_t min_len) const {
//        if (path.size() == 0)
//            return EdgePath();
//
//        const size_t start_pos = path.start_pos();
//        std::vector<EdgeId> edges = {path.sequence().front()};
//        size_t total_len = g_.length(path.sequence().front()) - start_pos;
//        size_t i = 1;
//
//        while (i < path.size() && total_len < min_len) {
//            EdgeId e = path[i];
//            edges.push_back(e);
//            total_len += g_.length(e);
//            ++i;
//        }
//
//        size_t end_pos;
//        if (i == path.size())
//            end_pos = path.end_pos();
//        else
//            end_pos = g_.length(path[i - 1]);
//
//        return EdgePath(edges, start_pos, end_pos);
//    }
//
//    EdgePath CropLeft(const EdgePath &path, size_t min_len) const {
//        return RCEdgePath(CropRight(RCEdgePath(path), min_len));
//    }

public:
    PartialGenePathProcessor(const Graph &g,
                             io::EdgeNamingF<Graph> edge_naming_f,
                             double max_len_coeff = 1.5) :
            g_(g),
            rel_comp_finder_(g_, edge_naming_f, max_len_coeff),
            edge_naming_f_(edge_naming_f) {}

    //potential stop codons of the gene will be added to stop_codon_poss
    omnigraph::GraphComponent<Graph> ProcessPartialGenePath(const EdgePath &gene_path,
                                                            size_t cds_len_est,
                                                            std::set<GraphPos> *stop_codon_poss) const;

};

class CDSSubgraphExtractor {
    const Graph &g_;
    const SequenceMapper<Graph> &mapper_;
    const PartialGenePathProcessor &cds_path_processor_;
    const ReadPathFinder<Graph> path_finder_;

    EdgePath ExtractPath(const std::string &partial_cds, const SequenceMapper<Graph> &mapper) const;

public:
    CDSSubgraphExtractor(const Graph &g,
                         const SequenceMapper<Graph> &mapper,
                         const PartialGenePathProcessor &cds_path_processor):
        g_(g),
        mapper_(mapper),
        cds_path_processor_(cds_path_processor),
        path_finder_(g, false, /*max gap length*/100) {

    }

    //TODO configure + might need to use query length instead of path length
    //TODO think about interface
    omnigraph::GraphComponent<Graph> ProcessPartialCDS(const std::string &partial_cds,
                                                       size_t cds_length_estimate,
                                                       std::set<GraphPos> *stop_codon_poss,
                                                       size_t min_len_to_explore = 200,
                                                       double frac_to_explore = 0.5) const;
};

}
