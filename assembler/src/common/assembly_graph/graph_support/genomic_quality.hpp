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

#include <algorithm>

namespace debruijn_graph {

template<class Graph>
class EdgeQuality: public visualization::graph_labeler::GraphLabeler<Graph>,
                   public omnigraph::GraphActionHandler<Graph> {
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    std::map<EdgeId, size_t> quality_;
    size_t k_;

    template<class Index>
    void FillQuality(const Index &index, const KmerMapper<Graph> &kmer_mapper, const Sequence &genome) {
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
    void Fill(const Index &index, const KmerMapper<Graph> &kmer_mapper, const Sequence &genome) {
        DEBUG("Filling quality values");
        FillQuality(index, kmer_mapper, genome);
        FillQuality(index, kmer_mapper, !genome);
        DEBUG(quality_.size() << " edges have non-zero quality");
    }

    template<class SequenceMapper, class SingleReadStream>
    void Fill(const SequenceMapper &mapper, SingleReadStream &stream) {
        io::SingleRead r;
        while (!stream.eof()) {
            stream >> r;
            for (auto e_mr : mapper.MapRead(r)) {
                quality_[e_mr.first] += e_mr.second.mapped_range.size();
            }
        }
    }

    void Fill(const EdgeQuality &edge_qual) {
        VERIFY(&(this->g()) == &edge_qual.g());
        for (const auto &e_q : edge_qual.quality_) {
            quality_[e_q.first] += e_q.second;
        }
    }

    EdgeQuality(const Graph &graph) :
            omnigraph::GraphActionHandler<Graph>(graph, "EdgeQuality"),
            k_(graph.k() + 1) {
    }

    void HandleAdd(EdgeId /*e*/) override {
    }

    void HandleDelete(EdgeId e) override {
        quality_.erase(e);
    }

    void HandleMerge(const std::vector<EdgeId> &old_edges, EdgeId new_edge) override {
        size_t res = 0;
        for (size_t i = 0; i < old_edges.size(); i++) {
            res += quality_[old_edges[i]];
        }
        quality_[new_edge] += res;
    }

    void HandleGlue(EdgeId new_edge, EdgeId edge1, EdgeId edge2) override {
        quality_[new_edge] += quality_[edge2];
        quality_[new_edge] += quality_[edge1];
    }

    void HandleSplit(EdgeId old_edge, EdgeId new_edge1, EdgeId new_edge2) override {
        if (!quality_.count(old_edge))
            return;

        const auto &g = this->g();
        auto qual_frac_f = [&] (EdgeId e) {
            return size_t(std::max(1., math::round(
                   double(quality_[old_edge]) * double(g.length(e)) / double(g.length(old_edge)))));
        };
        if (old_edge == g.conjugate(old_edge)) {
            quality_[new_edge1] = qual_frac_f(new_edge1);
            quality_[g.conjugate(new_edge1)] = qual_frac_f(new_edge1);
            quality_[g.conjugate(new_edge2)] = qual_frac_f(new_edge2);
        } else {
            quality_[new_edge1] = qual_frac_f(new_edge1);
            quality_[new_edge2] = qual_frac_f(new_edge2);
        }
    }

    double quality(EdgeId edge) const {
        auto it = quality_.find(edge);
        if (it == quality_.end())
            return 0.;
        else
            return (double) it->second / (double) this->g().length(edge);
    }

    void AddQuality(EdgeId edge, double quality) {
        //VERIFY(!IsPositiveQuality(e));
        quality_[edge] += (size_t) math::round(quality * (double) this->g().length(edge));
    }

    bool IsPositiveQuality(EdgeId edge) const {
        return math::gr(quality(edge), 0.);
    }

    bool IsZeroQuality(EdgeId edge) const {
        return math::eq(quality(edge), 0.);
    }

    std::string label(VertexId /*vertexId*/) const override {
        return "";
    }

    std::string label(EdgeId edge) const override {
        double q = quality(edge);
        return (q == 0) ? "" : "quality: " + std::to_string(q);
    }

    void clear() {
        //VERIFY_MSG(false, "Quality was cleared unexpectedly");
        quality_.clear();
    }

    std::set<EdgeId> PositiveQualEdges() const {
        std::set<EdgeId> answer;
        for (const auto & e_q : quality_) {
            if (e_q.second > 0) {
                answer.insert(e_q.first);
            }
        }
        return answer;
    }

private:
    DECL_LOGGER("EdgeQuality");
};

template<class Graph>
class QualityLoggingRemovalHandler {
    typedef typename Graph::EdgeId EdgeId;
    const Graph &g_;
    const EdgeQuality<Graph> &quality_handler_;
    size_t black_removed_;
    size_t total_;
    bool handle_all_;

    virtual void HandlePositiveQuality(EdgeId /*e*/) {

    }

public:
    QualityLoggingRemovalHandler(const Graph &g, const EdgeQuality<Graph> &quality_handler,
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

    const Graph &g() const {
        return g_;
    }

    const EdgeQuality<Graph> &quality_handler() const {
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
    QualityEdgeLocalityPrintingRH(const Graph &g
            , const EdgeQuality<Graph> &quality_handler
            , const visualization::graph_labeler::GraphLabeler<Graph> &labeler
            , std::shared_ptr<visualization::graph_colorer::GraphColorer<Graph>> colorer
            , const std::string &output_folder, bool handle_all = false) :
            base(g, quality_handler, handle_all),
            printing_rh_(g, labeler, colorer, output_folder)
    {}

    void HandlePositiveQuality(EdgeId e) override {
        printing_rh_.HandleDelete(e, "_" + std::to_string(this->quality_handler().quality(e)));
    }

private:
    DECL_LOGGER("QualityEdgeLocalityPrintingRH");
};

}
