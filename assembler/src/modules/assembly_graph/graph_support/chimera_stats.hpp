//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "assembly_graph/stats/statistics.hpp"
#include "assembly_graph/graph_support/genomic_quality.hpp"

namespace debruijn_graph {

namespace stats {

template<class Graph>
class ChimericEdgeClassifier {
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;

    const Graph& g_;
    size_t length_bound_;
    const EdgeQuality<Graph>& edge_qual_;
    bool real_edges_mode_;

    template<class EdgeContainer>
    vector<EdgeId> FilterNotEqual(const EdgeContainer& edges,
            EdgeId edge) const {
        vector<EdgeId> answer;
        for (EdgeId e : edges) {
            if (e != edge) {
                answer.push_back(e);
            }
        }
        return answer;
    }

    bool TopologyAndQualCheck(const vector<EdgeId>& edges) const {
        return edges.size() == 1 && edge_qual_.IsPositiveQuality(edges.front());
    }

    bool TopologyAndQualCheck(VertexId v, EdgeId e) const {
        return TopologyAndQualCheck(
                FilterNotEqual(g_.OutgoingEdges(v), e))
                && TopologyAndQualCheck(
                        FilterNotEqual(g_.IncomingEdges(v), e));
    }

    bool TopologyAndQualCheck(EdgeId e) const {
        return TopologyAndQualCheck(g_.EdgeStart(e), e)
                && TopologyAndQualCheck(g_.EdgeEnd(e), e);
    }

public:
    ChimericEdgeClassifier(const Graph& g, size_t length_bound, const EdgeQuality<Graph>& edge_qual, bool real_edges_mode = false)
    : g_(g),
      length_bound_(length_bound),
      edge_qual_(edge_qual),
      real_edges_mode_(real_edges_mode) {
    }

    bool IsTrivialChimeric(EdgeId e) const {
        bool correct_qual = real_edges_mode_ ? edge_qual_.IsPositiveQuality(e) : edge_qual_.IsZeroQuality(e);
        return correct_qual && g_.length(e) <= length_bound_
                && TopologyAndQualCheck(e);
    }

private:
    DECL_LOGGER("ChimericEdgeClassifier");
};

template<class Graph>
class InterstrandAnalyzer {
    const static size_t infinity = -1u;
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;

    const Graph& g_;
    size_t dist_bound_;
    const MappingPath<EdgeId> genome_path_;

    bool Relax(size_t& a, size_t b) const {
        if (b < a) {
            a = b;
            return true;
        }
        return false;
    }

    size_t GenomicDistance(size_t genome_path_pos, EdgeId e2,
            size_t distance_bound) const {
        for (size_t i = genome_path_pos + 1; i < genome_path_.size(); ++i) {
            int gap =
                    (int)(genome_path_[i].second.initial_range.start_pos
                            - genome_path_[genome_path_pos].second.initial_range.end_pos);
            VERIFY(gap >= 0);
            if (size_t(gap) > distance_bound)
                return infinity;
            if (genome_path_[i].first == e2)
                return gap;
        }
        return infinity;
    }

    size_t ShortestGenomicDistance(EdgeId e1, EdgeId e2,
            size_t distance_bound) const {
        size_t best = infinity;
        for (size_t i = 0; i < genome_path_.size(); ++i) {
            if (genome_path_[i].first == e1) {
                Relax(best, GenomicDistance(i, e2, distance_bound));
            }
        }
        return best;
    }

    size_t InnerInterstrandDistance(EdgeId e) const {
        size_t answer = infinity;
        EdgeId e1 = g_.GetUniqueIncomingEdge(g_.EdgeStart(e));
        EdgeId e2 = g_.GetUniqueOutgoingEdge(g_.EdgeEnd(e));
        if (g_.length(e2) > dist_bound_)
            return -1;
        Relax(answer,
                ShortestGenomicDistance(e1, g_.conjugate(e2),
                        dist_bound_ - g_.length(e2)));
        Relax(answer,
                ShortestGenomicDistance(e2, g_.conjugate(e1),
                        dist_bound_ - g_.length(e2)));
        return answer + g_.length(e2);
    }


public:
    InterstrandAnalyzer(const Graph& g, size_t dist_bound, const MappingPath<EdgeId> genome_path)
            : g_(g),
              dist_bound_(dist_bound),
              genome_path_(genome_path) {
    }

    //todo rewrite and think of additionally detecting thorns with no path
    //returns -1u if no interstrand path or interstrand distance > dist_bound
    size_t InterstrandDistance(EdgeId e) const {
        size_t answer = infinity;
        Relax(answer, InnerInterstrandDistance(e));
        Relax(answer, InnerInterstrandDistance(g_.conjugate(e)));
        //todo maybe unnecessary check
        return answer <= dist_bound_ ? answer : -1u;
    }

private:
    DECL_LOGGER("InterstrandAnalyzer");
};

template<class Graph>
class ChimericEdgeStats {
    const static size_t infinity = -1u;
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;

    const Graph& g_;
    const ChimericEdgeClassifier<Graph>& chimeric_edge_classifier_;
    const InterstrandAnalyzer<Graph>& interstrand_analyzer_;
    ostream& out_;

protected:
    virtual string Head() {
        std::stringstream ss;
        ss << "int_id\t"
                << "length\t"
                << "coverage\t"
                << "interstrand_dist"
                << endl;
        return ss.str();
    }

    virtual string ReportChimera(EdgeId e, size_t interstrand_dist) {
        std::stringstream ss;
        ss << g_.int_id(e) << "\t"
                << g_.length(e) << "\t"
                << g_.coverage(e) << "\t";
        if (interstrand_dist < infinity) {
            ss << interstrand_dist;
        } else {
            ss << -1;
        }
        ss << endl;
        return ss.str();
    }

    const Graph& g() const {
        return g_;
    }

public:
    ChimericEdgeStats(const Graph& g,
                      const ChimericEdgeClassifier<Graph>& chimeric_edge_classifier,
                      const InterstrandAnalyzer<Graph>& interstrand_analyzer,
                      ostream& out)
            : g_(g),
              chimeric_edge_classifier_(chimeric_edge_classifier),
              interstrand_analyzer_(interstrand_analyzer),
              out_(out) {
    }

    virtual ~ChimericEdgeStats() {
    }

    void operator()() {
        out_ << Head() << endl;
        set<EdgeId> visited;
        for (auto it = g_.SmartEdgeBegin(); !it.IsEnd(); ++it) {
            if (visited.count(*it) > 0)
                continue;
            visited.insert(*it);
            visited.insert(g_.conjugate(*it));
            if (chimeric_edge_classifier_.IsTrivialChimeric(*it)) {
                out_ << ReportChimera(*it, interstrand_analyzer_.InterstrandDistance(*it)) << endl;
            }
        }
    }
};

template<class Graph>
class ChimeraRelativeCoverageStats : public ChimericEdgeStats<Graph> {
    typedef ChimericEdgeStats<Graph> base;
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    typedef std::function<double(EdgeId, VertexId)> LocalCoverageFT;

    simplification::relative_coverage::RelativeCoverageHelper<Graph> rel_helper_;

    double RelativeCoverage(VertexId v, EdgeId base_edge) {
        return rel_helper_.RelativeCoverageToReport(v, rel_helper_.LocalCoverage(base_edge, v));
    }

public:
    ChimeraRelativeCoverageStats(const Graph& g,
                                 const ChimericEdgeClassifier<Graph>& edge_classifier,
                                 const InterstrandAnalyzer<Graph>& interstrand_analyzer,
                                 LocalCoverageFT local_coverage_f,
                                 ostream& out)
            : base(g, edge_classifier, interstrand_analyzer, out),
              rel_helper_(g, local_coverage_f, 2.0/*any value works here*/) {
    }

protected:
    virtual string Head() {
        return base::Head() + "\tmin_rel_cov\tmax_rel_cov";
    }

    virtual string ReportChimera(EdgeId e, size_t interstrand_dist) {
        double start_cov = RelativeCoverage(this->g().EdgeStart(e), e);
        double end_cov = RelativeCoverage(this->g().EdgeEnd(e), e);
        stringstream ss;
        ss << base::ReportChimera(e, interstrand_dist) << "\t"
                << std::min(start_cov, end_cov) << "\t"
                << std::max(start_cov, end_cov);
        return ss.str();
    }

private:
    DECL_LOGGER("ChimeraRelativeCoverageStats");
};

}
}
