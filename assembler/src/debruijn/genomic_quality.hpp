#pragma once

#include "omni/visualization/visualization.hpp"
#include "omni/basic_edge_conditions.hpp"

namespace debruijn_graph {

template<class Graph, class Index>
class EdgeQuality: public GraphLabeler<Graph>, public GraphActionHandler<Graph> {
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    map<EdgeId, size_t> quality_;
    size_t k_;

public:

    void FillQuality(const Index &index
            , const KmerMapper<Graph>& kmer_mapper, const Sequence &genome) {
        if (genome.size() < k_)
            return;
        runtime_k::RtSeq cur = genome.start<runtime_k::RtSeq>(k_);
        cur >>= 0;
        for (size_t i = 0; i + k_ - 1 < genome.size(); i++) {
            cur <<= genome[i + k_ - 1];
            auto corr_cur = kmer_mapper.Substitute(cur);
            if (index.contains(corr_cur)) {
                quality_[index.get(corr_cur).first]++;
            }
        }
    }

    EdgeQuality(const Graph &graph, const Index &index,
    const KmerMapper<Graph>& kmer_mapper,
    const Sequence &genome) :
            GraphActionHandler<Graph>(graph, "EdgeQualityLabeler"),
            k_(kmer_mapper.get_k()) {
        FillQuality(index, kmer_mapper, genome);
        FillQuality(index, kmer_mapper, !genome);
    }

    virtual ~EdgeQuality() {
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
        return (q == 0) ? "" : "quality: " + ToString(q);
    }

};

//template<class Graph, class Index>
//class ZeroQualityCondition : public EdgeCondition<Graph> {
//    typedef EdgeCondition<Graph> base;
//    const EdgeQuality<Graph, Index>& edge_qual_;
//
//public:
//    ZeroQualityCondition(const Graph& g, const EdgeQuality<Graph, Index>& edge_qual)
//            : base(g),
//              edge_qual_(edge_qual) {
//
//    }
//
//    bool Check(EdgeId e) const {
//        return edge_qual_.IsZeroQuality(e);
//    }
//
//};

template<class Graph, class Index>
class QualityLoggingRemovalHandler {
    typedef typename Graph::EdgeId EdgeId;
    const Graph& g_;
    const EdgeQuality<Graph, Index>& quality_handler_;
//  size_t black_removed_;
//  size_t colored_removed_;
public:
    QualityLoggingRemovalHandler(const Graph& g, const EdgeQuality<Graph, Index>& quality_handler) :
            g_(g), quality_handler_(quality_handler)/*, black_removed_(0), colored_removed_(
     0)*/{

    }

    void HandleDelete(EdgeId edge) {
        if (math::gr(quality_handler_.quality(edge), 0.)) {
            TRACE("Deleting edge " << g_.int_id(edge) << " with quality " << quality_handler_.quality(edge));
        } else {
//          TRACE("Deleting edge " << g_.int_id(edge) << " with zero quality");
        }
//      if (math::gr(quality_handler_.quality(edge), 0.))
//          colored_removed_++;
//      else
//          black_removed_++;
    }

private:
    DECL_LOGGER("QualityLoggingRemovalHandler")
    ;
};

template<class Graph, class Index>
class QualityLoggingRemovalCountHandler {
    typedef typename Graph::EdgeId EdgeId;
    const Graph& g_;
    const EdgeQuality<Graph, Index>& quality_handler_;
    size_t black_removed_;
    size_t total;

public:
    QualityLoggingRemovalCountHandler(const Graph& g, const EdgeQuality<Graph, Index>& quality_handler) :
            g_(g), quality_handler_(quality_handler)/*, black_removed_(0), colored_removed_(
     0)*/{
        black_removed_ = 0;
        total = 0;
    }

    void HandleDelete(EdgeId edge) {
        total++;
        if (math::gr(quality_handler_.quality(edge), 0.)) {
            TRACE("Deleting good edge " << g_.int_id(edge) << " with quality " << quality_handler_.quality(edge) << " cov " << g_.coverage(edge) << " length " << g_.length(edge));
        }else{
            black_removed_++;
        }
        if ((total % (1<<10)) != 0)
            TRACE("Removed still " << black_removed_ << " " << total);
    }

private:
};

template<class Graph, class Index>
class QualityEdgeLocalityPrintingRH {
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    const Graph& g_;
    const EdgeQuality<Graph, Index>& quality_handler_;
    const GraphLabeler<Graph>& labeler_;
    const string& output_folder_;
//  size_t black_removed_;
//  size_t colored_removed_;
public:
    QualityEdgeLocalityPrintingRH(const Graph& g
            , const EdgeQuality<Graph, Index>& quality_handler
            , const GraphLabeler<Graph>& labeler
            , const string& output_folder) :
            g_(g), quality_handler_(quality_handler),
            labeler_(labeler), output_folder_(output_folder){
    }

    void HandleDelete(EdgeId edge) {
        if (quality_handler_.IsPositiveQuality(edge)) {
            DEBUG("Deleting edge " << g_.str(edge) << " with quality " << quality_handler_.quality(edge));
            string folder = output_folder_ + "colored_edges_deleted/";
            path::make_dir(folder);
            //todo magic constant
//          map<EdgeId, string> empty_coloring;
            omnigraph::visualization::WriteComponent(g_, omnigraph::EdgeNeighborhood<Graph>(g_, edge, 50, 250)
                    , folder + "edge_" +  ToString(g_.int_id(edge)) + "_" + ToString(quality_handler_.quality(edge)) + ".dot"
                    , omnigraph::visualization::DefaultColorer(g_), labeler_);
        } else {
            TRACE("Deleting edge " << g_.str(edge) << " with quality " << quality_handler_.quality(edge));
        }
    }

private:
    DECL_LOGGER("QualityEdgeLocalityPrintingRH")
    ;
};

}
