//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "visualization/visualization.hpp"
#include "assembly_graph/graph_support/basic_edge_conditions.hpp"
#include "modules/alignment/sequence_mapper.hpp"
#include "assembly_graph/core/action_handlers.hpp"

namespace debruijn_graph {

template<class Graph>
class EdgeQuality: public visualization::graph_labeler::GraphLabeler<Graph>, public omnigraph::GraphActionHandler<Graph> {
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    map<EdgeId, size_t> quality_;
    size_t k_;

    template<class Index>
    void FillQuality(const Index &index
            , const KmerMapper<Graph>& kmer_mapper, const Sequence &genome) {
        if (genome.size() < k_)
            return;
        RtSeq cur = genome.start<RtSeq>(k_);
        cur >>= 0;
        for (size_t i = 0; i + k_ - 1 < genome.size(); i++) {
            cur <<= genome[i + k_ - 1];
            auto corr_cur = kmer_mapper.Substitute(cur);
            if (index.contains(corr_cur)) {
                quality_[index.get(corr_cur).first]++;
            }
        }
    }

public:

    template<class Index>
    void Fill(const Index &index
            , const KmerMapper<Graph>& kmer_mapper
            , const Sequence &genome) {
        DEBUG("Filling quality values");
        FillQuality(index, kmer_mapper, genome);
        FillQuality(index, kmer_mapper, !genome);
        DEBUG(quality_.size() << " edges have non-zero quality");
    }

    EdgeQuality(const Graph &graph) :
            omnigraph::GraphActionHandler<Graph>(graph, "EdgeQuality"),
            k_(graph.k() + 1) {
    }

    virtual void HandleAdd(EdgeId /*e*/) {
    }

    virtual void HandleDelete(EdgeId e) {
        quality_.erase(e);
    }

    virtual void HandleMerge(const vector<EdgeId>& old_edges, EdgeId new_edge) {
        size_t res = 0;
        for (size_t i = 0; i < old_edges.size(); i++) {
            res += quality_[old_edges[i]];
        }
        quality_[new_edge] += res;
    }

    virtual void HandleGlue(EdgeId new_edge, EdgeId edge1, EdgeId edge2) {
        quality_[new_edge] += quality_[edge2];
        quality_[new_edge] += quality_[edge1];
    }

    virtual void HandleSplit(EdgeId old_edge, EdgeId new_edge1,
            EdgeId new_edge2) {
        if (old_edge == this->g().conjugate(old_edge)) {
            WARN("EdgeQuality does not support self-conjugate splits");
            return;
        }
        VERIFY(old_edge != this->g().conjugate(old_edge));
        quality_[new_edge1] = quality_[old_edge] * this->g().length(new_edge1)
                / (this->g().length(new_edge1) + this->g().length(new_edge2));
        quality_[new_edge2] = quality_[old_edge] * this->g().length(new_edge2)
                / (this->g().length(new_edge1) + this->g().length(new_edge2));
    }

    double quality(EdgeId edge) const {
        auto it = quality_.find(edge);
        if (it == quality_.end())
            return 0.;
        else
            return 1. * (double) it->second / (double) this->g().length(edge);
    }

    bool IsPositiveQuality(EdgeId edge) const {
        return math::gr(quality(edge), 0.);
    }

    bool IsZeroQuality(EdgeId edge) const {
        return math::eq(quality(edge), 0.);
    }

    virtual std::string label(VertexId /*vertexId*/) const {
        return "";
    }

    virtual std::string label(EdgeId edge) const {
        double q = quality(edge);
        return (q == 0) ? "" : "quality: " + std::to_string(q);
    }

    void clear() {
        quality_.clear();
    }

private:
    DECL_LOGGER("EdgeQuality");
};

template<class Graph>
class QualityLoggingRemovalHandler {
    typedef typename Graph::EdgeId EdgeId;
    const Graph& g_;
    const EdgeQuality<Graph>& quality_handler_;
    size_t black_removed_;
    size_t total_;
    bool handle_all_;

    virtual void HandlePositiveQuality(EdgeId /*e*/) {

    }

public:
    QualityLoggingRemovalHandler(const Graph& g, const EdgeQuality<Graph>& quality_handler,
                                 bool handle_all = false) :
            g_(g), quality_handler_(quality_handler), black_removed_(0), total_(0), handle_all_(handle_all) {
    }

    void HandleDelete(EdgeId e) {
        total_++;
        if (handle_all_ || math::gr(quality_handler_.quality(e), 0.)) {
            TRACE("Deleting good edge id = " << g_.int_id(e)
                  << "; length = " << g_.length(e)
                  << "; quality = " << quality_handler_.quality(e)
                  << "; cov = " << g_.coverage(e));
            HandlePositiveQuality(e);
        } else {
            black_removed_++;
        }
    }

    const Graph& g() const {
        return g_;
    }

    const EdgeQuality<Graph>& quality_handler() const {
        return quality_handler_;
    }

    virtual ~QualityLoggingRemovalHandler() {
        TRACE("Overall stats: total removed = " << total_
              << "; bad removed = " << black_removed_
              << "; good removed = " << total_ - black_removed_);
    }

private:
    DECL_LOGGER("QualityLoggingRemovalHandler");
};

template<class Graph>
class QualityEdgeLocalityPrintingRH : public QualityLoggingRemovalHandler<Graph> {
    typedef QualityLoggingRemovalHandler<Graph> base;
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    visualization::visualization_utils::LocalityPrintingRH<Graph> printing_rh_;
public:
    QualityEdgeLocalityPrintingRH(const Graph& g
            , const EdgeQuality<Graph>& quality_handler
            , const visualization::graph_labeler::GraphLabeler<Graph>& labeler
            , std::shared_ptr<visualization::graph_colorer::GraphColorer<Graph>> colorer
            , const string& output_folder, bool handle_all = false) :
            base(g, quality_handler, handle_all),
            printing_rh_(g, labeler, colorer, output_folder)
    {}

    virtual void HandlePositiveQuality(EdgeId e) {
        printing_rh_.HandleDelete(e, "_" + std::to_string(this->quality_handler().quality(e)));
    }

private:
    DECL_LOGGER("QualityEdgeLocalityPrintingRH");
};

//earlier version from rel_cov branch
//template<class Graph>
//class EdgeNeighborhoodFinder: public omnigraph::GraphSplitter<Graph> {
//private:
//  typedef typename Graph::EdgeId EdgeId;
//  typedef typename Graph::VertexId VertexId;
//  EdgeId edge_;
//  size_t max_size_;
//  size_t edge_length_bound_;
//  bool finished_;
//public:
//  EdgeNeighborhoodFinder(const Graph &graph, EdgeId edge, size_t max_size
//          , size_t edge_length_bound) :
//          GraphSplitter<Graph>(graph), edge_(edge), max_size_(
//                  max_size), edge_length_bound_(edge_length_bound), finished_(
//                  false) {
//  }
//
//  GraphComponent<Graph> NextComponent() {
//      CountingDijkstra<Graph> cf(this->graph(), max_size_,
//              edge_length_bound_);
//      set<VertexId> result_set;
//      cf.run(this->graph().EdgeStart(edge_));
//      vector<VertexId> result_start = cf.ReachedVertices();
//      result_set.insert(result_start.begin(), result_start.end());
//      cf.run(this->graph().EdgeEnd(edge_));
//      vector<VertexId> result_end = cf.ReachedVertices();
//      result_set.insert(result_end.begin(), result_end.end());
//
//      ComponentCloser<Graph> cc(this->graph(), edge_length_bound_);
//      cc.CloseComponent(result_set);
//
//      finished_ = true;
//      return GraphComponent<Graph>(this->graph(), result_set.begin(), result_set.end());
//  }
//
//  /*virtual*/ bool Finished() {
//      return finished_;
//  }
//};
//
//template<class Graph>
//class EdgeLocalityPrintingRH {
//  typedef typename Graph::EdgeId EdgeId;
//  typedef typename Graph::VertexId VertexId;
//  const Graph& g_;
//  const GraphLabeler<Graph>& labeler_;
//  const string& output_folder_;
//    std::function<double (EdgeId)>& quality_f_;
////    size_t black_removed_;
////    size_t colored_removed_;
//public:
//  EdgeLocalityPrintingRH(const Graph& g
//          , const GraphLabeler<Graph>& labeler
//          , const string& output_folder
//            , std::function<double (EdgeId)> quality_f = 0) :
//          g_(g),
//          labeler_(labeler), output_folder_(output_folder),
//            quality_f_(quality_f){
//  }
//
//  void HandleDelete(EdgeId edge) {
//            TRACE("Deleting edge " << g_.str(edge));
//            if (quality_f_ && math::gr(quality_f_(edge), 0.))
//                INFO("EdgeLocalityPrintRH handling the edge with positive quality : " << quality_f_(edge) << " " << g_.str(edge));
//
//            string folder = output_folder_ + "edges_deleted/";
//            path::make_dir(folder);
//            //todo magic constant
//            map<EdgeId, string> empty_coloring;
//            visualization::visualization_utils::WriteComponent(g_, EdgeNeighborhood<Graph>(g_, edge, 50, 250),
//                  folder + "edge_" +  std::to_string(g_.int_id(edge)) + ".dot", empty_coloring, labeler_);
//  }
//
//private:
//  DECL_LOGGER("QualityEdgeLocalityPrintingRH")
//  ;
//};

//template<class Graph, class Index>
//class EdgeQuality: public GraphLabeler<Graph>, public GraphActionHandler<Graph> {
//    typedef typename Graph::EdgeId EdgeId;
//    typedef typename Graph::VertexId VertexId;
//    map<EdgeId, size_t> quality_;
//    size_t k_;
//
//public:
//
//    void FillQuality(const Index &index
//            , const KmerMapper<Graph>& kmer_mapper, const Sequence &genome) {
//        if (genome.size() < k_)
//            return;
//        RtSeq cur = genome.start<RtSeq>(k_);
//        cur >>= 0;
//        for (size_t i = 0; i + k_ - 1 < genome.size(); i++) {
//            cur <<= genome[i + k_ - 1];
//            auto corr_cur = kmer_mapper.Substitute(cur);
//            if (index.contains(corr_cur)) {
//                quality_[index.get(corr_cur).first]++;
//            }
//        }
//    }
//
//    EdgeQuality(const Graph &graph, const Index &index,
//    const KmerMapper<Graph>& kmer_mapper,
//    const Sequence &genome) :
//
//            GraphActionHandler<Graph>(graph, "EdgeQualityLabeler"),
//            k_(kmer_mapper.get_k()) {
//        FillQuality(index, kmer_mapper, genome);
//        FillQuality(index, kmer_mapper, !genome);
//    }
//
//    virtual ~EdgeQuality() {
//    }
//
//    virtual void HandleAdd(EdgeId /*e*/) {
//    }
//
//    virtual void HandleDelete(EdgeId e) {
//        quality_.erase(e);
//    }
//
//    virtual void HandleMerge(const vector<EdgeId>& old_edges, EdgeId new_edge) {
//        size_t res = 0;
//        for (size_t i = 0; i < old_edges.size(); i++) {
//            res += quality_[old_edges[i]];
//        }
//        quality_[new_edge] += res;
//    }
//
//    virtual void HandleGlue(EdgeId new_edge, EdgeId edge1, EdgeId edge2) {
//        quality_[new_edge] += quality_[edge2];
//        quality_[new_edge] += quality_[edge1];
//    }
//
//    virtual void HandleSplit(EdgeId old_edge, EdgeId new_edge1,
//            EdgeId new_edge2) {
//        quality_[new_edge1] = quality_[old_edge] * this->g().length(new_edge1)
//                / (this->g().length(new_edge1) + this->g().length(new_edge2));
//        quality_[new_edge2] = quality_[old_edge] * this->g().length(new_edge2)
//                / (this->g().length(new_edge1) + this->g().length(new_edge2));
//    }
//
//    double quality(EdgeId edge) const {
//        auto it = quality_.find(edge);
//        if (it == quality_.end())
//            return 0.;
//        else
//            return 1. * (double) it->second / (double) this->g().length(edge);
//    }
//
//    bool IsPositiveQuality(EdgeId edge) const {
//        return math::gr(quality(edge), 0.);
//    }
//
//    virtual std::string label(VertexId /*vertexId*/) const {
//        return "";
//    }
//
//    virtual std::string label(EdgeId edge) const {
//        double q = quality(edge);
//        return (q == 0) ? "" : "quality: " + std::to_string(q);
//    }
//
//};
//
//template<class Graph, class Index>
//class QualityLoggingRemovalHandler {
//    typedef typename Graph::EdgeId EdgeId;
//    const Graph& g_;
//    const EdgeQuality<Graph, Index>& quality_handler_;
////  size_t black_removed_;
////  size_t colored_removed_;
//public:
//    QualityLoggingRemovalHandler(const Graph& g, const EdgeQuality<Graph, Index>& quality_handler) :
//            g_(g), quality_handler_(quality_handler)/*, black_removed_(0), colored_removed_(
//     0)*/{
//
//    }
//
//    void HandleDelete(EdgeId edge) {
//        if (math::gr(quality_handler_.quality(edge), 0.)) {
//            TRACE("Deleting edge " << g_.str(edge) << " with quality " << quality_handler_.quality(edge));
//        } else {
////          TRACE("Deleting edge " << g_.int_id(edge) << " with zero quality");
//        }
////      if (math::gr(quality_handler_.quality(edge), 0.))
////          colored_removed_++;
////      else
////          black_removed_++;
//    }
//
//private:
//    DECL_LOGGER("QualityLoggingRemovalHandler")
//    ;
//};
//
//template<class Graph, class Index>
//class QualityLoggingRemovalCountHandler {
//    typedef typename Graph::EdgeId EdgeId;
//    const Graph& g_;
//    const EdgeQuality<Graph, Index>& quality_handler_;
//    size_t black_removed_;
//    size_t total;
//
//public:
//    QualityLoggingRemovalCountHandler(const Graph& g, const EdgeQuality<Graph, Index>& quality_handler) :
//            g_(g), quality_handler_(quality_handler)/*, black_removed_(0), colored_removed_(
//     0)*/{
//        black_removed_ = 0;
//        total = 0;
//    }
//
//    void HandleDelete(EdgeId edge) {
//        total++;
//        if (math::gr(quality_handler_.quality(edge), 0.)) {
//            TRACE("Deleting good edge " << g_.int_id(edge) << " with quality " << quality_handler_.quality(edge) << " cov " << g_.coverage(edge) << " length " << g_.length(edge));
//        }else{
//            black_removed_++;
//        }
//        if ((total % (1<<10)) != 0)
//            TRACE("Removed still " << black_removed_ << " " << total);
//    }
//
//private:
//};
//
//template<class Graph, class Index>
//class QualityEdgeLocalityPrintingRH {
//    typedef typename Graph::EdgeId EdgeId;
//    typedef typename Graph::VertexId VertexId;
//    const Graph& g_;
//    const EdgeQuality<Graph, Index>& quality_handler_;
//    const omnigraph::GraphLabeler<Graph>& labeler_;
//    const visualization::graph_colorer::GraphColorer<Graph>& colorer_;
//    const string& output_folder_;
////  size_t black_removed_;
////  size_t colored_removed_;
//public:
//    QualityEdgeLocalityPrintingRH(const Graph& g
//            , const EdgeQuality<Graph, Index>& quality_handler
//            , const visualization::graph_labeler::GraphLabeler<Graph>& labeler
//            , const visualization::graph_colorer::GraphColorer<Graph>& colorer
//            , const string& output_folder) :
//            g_(g), quality_handler_(quality_handler),
//            labeler_(labeler), colorer_(colorer), output_folder_(output_folder){
//    }
//
//    void HandleDelete(EdgeId edge) {
//        if (quality_handler_.IsPositiveQuality(edge)) {
//            DEBUG("Deleting edge " << g_.str(edge) << " with quality " << quality_handler_.quality(edge));
//            string folder = output_folder_ + "colored_edges_deleted/";
//            path::make_dir(folder);
//            //todo magic constant
////          map<EdgeId, string> empty_coloring;
//            shared_ptr<GraphSplitter<Graph>> splitter = EdgeNeighborhoodFinder<Graph>(g_, edge, 50, 250);
//            visualization::visualization_utils::WriteComponents(g_, *splitter/*, "locality_of_edge_" + std::to_string(g_.int_id(edge))*/
//                    , folder + "edge_" +  std::to_string(g_.int_id(edge)) + "_" + std::to_string(quality_handler_.quality(edge)) + ".dot"
//                    , colorer_, labeler_);
//        } else {
//            TRACE("Deleting edge " << g_.str(edge) << " with zero quality");
//        }
//    }
//
//private:
//    DECL_LOGGER("QualityEdgeLocalityPrintingRH")
//    ;
//};
//
//template<class Graph, class Index>
//class QualityPairInfoHandler {
//    typedef typename Graph::EdgeId EdgeId;
//    typedef typename Graph::VertexId VertexId;
//    typedef omnigraph::PairInfo<EdgeId> PairInfo;
//    typedef vector<PairInfo> PairInfos;
//    const Graph& g_;
//    const EdgeQuality<Graph, Index>& quality_handler_;
//    const GraphLabeler<Graph>& labeler_;
//    const string& output_folder_;
//    const PairedInfoIndex<ConjugateDeBruijnGraph>& index_;
////  size_t black_removed_;
////  size_t colored_removed_;
//public:
//    QualityPairInfoHandler(const Graph& g
//            , const EdgeQuality<Graph, Index>& quality_handler
//            , const GraphLabeler<Graph>& labeler
//            , const string& output_folder
//            , const PairedInfoIndex<ConjugateDeBruijnGraph>& index) :
//            g_(g), quality_handler_(quality_handler),
//            labeler_(labeler), output_folder_(output_folder), index_(index) {
//    }
//
//    void HandleDelete(EdgeId edge) {
//        if (quality_handler_.IsPositiveQuality(edge)) {
//            cout << "Deleting edge " << g_.str(edge) << " with quality " << quality_handler_.quality(edge) << endl;
//            string folder = output_folder_ + "colored_edges_deleted/";
//            path::make_dir(folder);
//            //todo magic constant
//            PairInfos infos = index_.GetEdgeInfo(edge);
//            if (infos.size() > 0){
//                for (size_t i = 0; i<infos.size(); i++){
//                    cout << "Tip Info " << g_.int_id(infos[i].first) << " " << g_.int_id(infos[i].second) << " " << infos[i].d << " " << infos[i].weight << " " << infos[i].variance << endl;
//                }
//            }
//            map<EdgeId, string> empty_coloring;
//            shared_ptr<GraphSplitter<Graph>> splitter = EdgeNeighborhoodFinder<Graph>(g_, edge, 50,
//                    250);
//
//            visualization::visualization_utils::WriteComponents(g_, *splitter, TrueFilter<vector<VertexId>>(), "locality_of_edge_" + std::to_string(g_.int_id(edge))
//                    , folder + "edge_" +  std::to_string(g_.int_id(edge)) + "_" + std::to_string(quality_handler_.quality(edge)) + ".dot"
//                    , empty_coloring, labeler_);
//        }
//    }
//
//private:
//};
//
////todo what is the difference with QELPRH?!
//template<class Graph>
//class EdgeLocalityPrintingRH {
//    typedef typename Graph::EdgeId EdgeId;
//    typedef typename Graph::VertexId VertexId;
//    const Graph& g_;
//    const GraphLabeler<Graph>& labeler_;
//    const string& output_folder_;
//    std::function<double (EdgeId)>& quality_f_;
////  size_t black_removed_;
////  size_t colored_removed_;
//public:
//    EdgeLocalityPrintingRH(const Graph& g
//            , const GraphLabeler<Graph>& labeler
//            , const string& output_folder
//            , std::function<double (EdgeId)> quality_f = 0) :
//            g_(g),
//            labeler_(labeler), output_folder_(output_folder),
//            quality_f_(quality_f){
//    }
//
//    void HandleDelete(EdgeId edge) {
//            TRACE("Deleting edge " << g_.str(edge));
//            if (quality_f_ && math::gr(quality_f_(edge), 0.))
//                INFO("Handling the edge with positive quality : " << quality_f_(edge) << " " << g_.str(edge));
//
//            string folder = output_folder_ + "edges_deleted/";
//            path::make_dir(folder);
//            //todo magic constant
//            map<EdgeId, string> empty_coloring;
//            shared_ptr<GraphSplitter<Graph>> splitter = EdgeNeighborhoodFinder<Graph>(g_, edge, 50, 250);
//            visualization::visualization_utils::WriteComponents(g_, *splitter, TrueFilter<vector<VertexId>>(), "locality_of_edge_" + std::to_string(g_.int_id(edge))
//                    , folder + "edge_" +  std::to_string(g_.int_id(edge)) + ".dot", empty_coloring, labeler_);
//    }
//
//private:
//    DECL_LOGGER("EdgeLocalityPrintingRH")
//    ;
//};

}
