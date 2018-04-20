//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

/*
 * erroneous_connection_remover.hpp
 *
 *  Created on: May 31, 2011
 *      Author: sergey
 */

#pragma once

#include "assembly_graph/graph_support/graph_processing_algorithm.hpp"
#include "assembly_graph/graph_support/basic_edge_conditions.hpp"
#include "math/xmath.h"
#include "assembly_graph/dijkstra/dijkstra_helper.hpp"
#include "assembly_graph/core/coverage.hpp"
#include "assembly_graph/graph_support/detail_coverage.hpp"
#include "modules/simplification/topological_edge_conditions.hpp"

namespace omnigraph {

//todo move to rnaSPAdes project
template<class Graph>
class RelativeCoverageECCondition: public EdgeCondition<Graph> {
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    typedef EdgeCondition<Graph> base;

    const double rcec_ratio_;

    template<class ContainerType>
    double SumCompetitorCoverage(EdgeId ec_edge, const ContainerType& edges) const {
        const Graph &g = this->g();
        double sum = 0;
        for (EdgeId e : edges) {
            //update if competitor edge is not loop
            if (e != ec_edge && g.EdgeStart(e) != g.EdgeEnd(e))
                sum += g.coverage(e);
        }
        return sum;
    }

    double AvgLocalityCoverage(EdgeId ec_edge) const {
        const Graph &g = this->g();
        VertexId start = g.EdgeStart(ec_edge), end = g.EdgeEnd(ec_edge);
        auto in_start = g.IncomingEdges(start);
        auto out_start = g.OutgoingEdges(start);
        auto in_end = g.IncomingEdges(end);
        auto out_end = g.OutgoingEdges(end);
        double total_edges = double(g.IncomingEdgeCount(start) + g.OutgoingEdgeCount(start) +
            g.IncomingEdgeCount(end) + g.OutgoingEdgeCount(end) - 2);
        return (SumCompetitorCoverage(ec_edge, in_start) +
                SumCompetitorCoverage(ec_edge, out_start) +
                SumCompetitorCoverage(ec_edge, in_end) +
                SumCompetitorCoverage(ec_edge, out_end)) / total_edges;
    }

    template<class ContainerType>
    double MaxCompetitorCoverage(EdgeId ec_edge, const ContainerType& edges) const {
        const Graph &g = this->g();
        double result = 0;
        for (EdgeId e : edges) {
            //update if competitor edge is not loop
            if (e != ec_edge && g.EdgeStart(e) != g.EdgeEnd(e))
                result = std::max(result, g.coverage(e));
        }
        return result;
    }

    double MaxCompetitorCoverage(EdgeId ec_edge) const {
        const Graph &g = this->g();
        VertexId start = g.EdgeStart(ec_edge), end = g.EdgeEnd(ec_edge);
        auto in_start = g.IncomingEdges(start);
        auto out_start = g.OutgoingEdges(start);
        auto in_end = g.IncomingEdges(end);
        auto out_end = g.OutgoingEdges(end);
        return std::max(
                std::max(MaxCompetitorCoverage(ec_edge, in_start),
                         MaxCompetitorCoverage(ec_edge, out_start)),
                std::max(MaxCompetitorCoverage(ec_edge, in_end),
                         MaxCompetitorCoverage(ec_edge, out_end)));
    }

public:

    RelativeCoverageECCondition(const Graph& g, double rcec_ratio) :
            base(g), rcec_ratio_(rcec_ratio) {
    }

    bool Check(EdgeId e) const override {
        //+1 is a trick to deal with edges of 0 coverage from iterative run
        double locality_coverage = AvgLocalityCoverage(e) + 1;
        return math::le(this->g().coverage(e), rcec_ratio_ * locality_coverage);
    }

};

//todo move to rnaSPAdes project
template<class Graph>
func::TypedPredicate<typename Graph::EdgeId> AddRelativeCoverageECCondition(const Graph &g, double rcec_ratio,
                                                                            func::TypedPredicate<typename Graph::EdgeId> condition) {
    return func::And(RelativeCoverageECCondition<Graph>(g, rcec_ratio), condition);
}

//todo move to rnaSPAdes project
template<class Graph>
inline bool IsSimpleBulge(const Graph &g, typename Graph::EdgeId e){
    size_t edge_count = g.GetEdgesBetween(g.EdgeStart(e), g.EdgeEnd(e)).size();

    return edge_count == g.OutgoingEdgeCount(g.EdgeStart(e)) &&
           edge_count == g.IncomingEdgeCount(g.EdgeEnd(e)) &&
           edge_count >= 2;
}

template<class Graph>
inline bool IsAlternativePathExist(const Graph &g, typename Graph::EdgeId e){
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;

    MostCoveredSimpleAlternativePathChooser<Graph> path_chooser(g, e);

    VertexId start = g.EdgeStart(e);
    TRACE("Start " << g.str(start));
    VertexId end = g.EdgeEnd(e);
    TRACE("End " << g.str(end));

    ProcessPaths(g, 0, std::numeric_limits<std::size_t>::max(), start, end, path_chooser);

    const vector<EdgeId>& path = path_chooser.most_covered_path();
    double path_coverage = path_chooser.max_coverage();
    if (!path.empty() && math::gr(path_coverage, 0.)) {
        VERIFY(g.EdgeStart(path[0]) == start);
        VERIFY(g.EdgeEnd(path.back()) == end);

        return true;
    }
    else
        return false;
}

template<class Graph>
inline bool IsAlternativeInclusivePathExist(const Graph &g, typename Graph::EdgeId forbidden_edge, typename Graph::EdgeId compulsory_edge){
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;

    MostCoveredSimpleAlternativePathChooser<Graph> path_chooser(g, forbidden_edge);

    VertexId start = g.EdgeStart(forbidden_edge);
    TRACE("Start " << g.str(start));
    VertexId end = g.EdgeEnd(forbidden_edge);
    TRACE("End " << g.str(end));

    ProcessPaths(g, 0, std::numeric_limits<std::size_t>::max(), start, end, path_chooser);

    const vector<EdgeId>& path = path_chooser.most_covered_path();
    double path_coverage = path_chooser.max_coverage();
    if (!path.empty() && math::gr(path_coverage, 0.)) {
        VERIFY(g.EdgeStart(path[0]) == start);
        VERIFY(g.EdgeEnd(path.back()) == end);

        if(std::find(path.begin(), path.end(), compulsory_edge) != path.end()){
            return true;
        }
    }
    return false;
}

template<class Graph>
inline bool IsReachableBulge(const Graph &g, typename Graph::EdgeId e){
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;

    bool res = IsAlternativePathExist(g, e);
    if(res)
        return res;
    else{
        VertexId start = g.EdgeStart(e), end = g.EdgeEnd(e);
        vector<EdgeId> incident;
        utils::push_back_all(incident, g.IncomingEdges(end));
        utils::push_back_all(incident, g.OutgoingEdges(start));
        for (auto it = incident.begin(); it != incident.end(); ++it){
            res = IsAlternativeInclusivePathExist(g, *it, e);
            if(res){
                return res;
            }
        }
    }
    return false;
}

//todo move to rnaSPAdes project
template<class Graph>
class NotBulgeECCondition : public EdgeCondition<Graph> {
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    typedef EdgeCondition<Graph> base;

public:

    NotBulgeECCondition(const Graph &g)
            : base(g) {

    }

    bool Check(EdgeId e) const {
        if (HasAlternatives(this->g(), e) && !IsSimpleBulge(this->g(), e)){
            DEBUG("edge id = " << this->g().int_id(e)
                 << " between = " << this->g().GetEdgesBetween(this->g().EdgeStart(e), this->g().EdgeEnd(e)).size()
                 << " between ids: " << this->g().GetEdgesBetween(this->g().EdgeStart(e), this->g().EdgeEnd(e))
                 << " outgoing s = " << this->g().OutgoingEdgeCount(this->g().EdgeStart(e))
                 << " incoming e = " << this->g().IncomingEdgeCount(this->g().EdgeEnd(e)));
        }
//        return !IsSimpleBulge(this->g(), e);
        return !IsReachableBulge(this->g(), e);
    }

private:
    DECL_LOGGER("NotBulgeECCondition");

};

//todo move to rnaSPAdes project
template<class Graph>
func::TypedPredicate<typename Graph::EdgeId> AddNotBulgeECCondition(const Graph &g,
                                                                    func::TypedPredicate<typename Graph::EdgeId> condition) {
    return func::And(NotBulgeECCondition<Graph>(g), condition);
}

template<class Graph>
class TopologicalThornCondition : public EdgeCondition<Graph> {
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    typedef EdgeCondition<Graph> base;
    typedef std::vector<EdgeId> Path;

    size_t max_jump_distance_;
    size_t max_edge_cnt_;

    bool CheckEdgeCounts(EdgeId e) const {
        if (this->g().EdgeStart(e) == this->g().EdgeEnd(e))
            return false;
        if (this->g().OutgoingEdgeCount(this->g().EdgeStart(e)) != 2)
            return false;
        if (this->g().IncomingEdgeCount(this->g().EdgeStart(e)) != 1)
            return false;
        if (this->g().OutgoingEdgeCount(this->g().EdgeEnd(e)) != 1)
            return false;
        if (this->g().IncomingEdgeCount(this->g().EdgeEnd(e)) != 2)
            return false;
        return true;
    }

public:

    TopologicalThornCondition(Graph& g,
                              size_t max_jump_dist,
                              size_t max_edge_cnt = std::numeric_limits<size_t>::max())
            : base(g),
              max_jump_distance_(max_jump_dist),
              max_edge_cnt_(max_edge_cnt) {
    }

    bool Check(EdgeId e) const override {
        const Graph& g = this->g();
        if (!CheckEdgeCounts(e))
            return false;

        //fixme micro-optimization to be removed
        if (g.conjugate(g.EdgeStart(e)) == g.EdgeEnd(e)) {
            return true;
        }

        auto comparator = [](const Path& a, const Path& b) {return a.size() >= b.size();};

        BestPathStorage<Graph, decltype(comparator)> callback(g, comparator);
        ProcessPaths(g, 0, max_jump_distance_, g.EdgeStart(e), g.conjugate(g.EdgeEnd(e)), callback, max_edge_cnt_);
        return (bool) callback.best_path();
    }
};

template<class Graph>
class AdditionalMDAThornCondition : public EdgeCondition<Graph> {
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    typedef EdgeCondition<Graph> base;
    typedef std::vector<EdgeId> Path;

    size_t uniqueness_length_;

    bool CheckUniqueness(EdgeId e) const {
        return this->g().length(e) >= uniqueness_length_;
    }

    bool CheckUnique(VertexId v) const {
        return this->g().CheckUniqueIncomingEdge(v) &&
                CheckUniqueness(this->g().GetUniqueIncomingEdge(v));
    }

    bool CheckUniqueCondition(EdgeId e) const {
        TRACE("Checking conditions for edge start");
        return CheckUnique(this->g().EdgeStart(e)) ||
                CheckUnique(this->g().conjugate(this->g().EdgeEnd(e)));
    }

    template<class EdgeContainer>
    bool CheckAlternativesForEC(const EdgeContainer& edges, EdgeId base) const {
        for (EdgeId e: edges) {
            if (e != base && this->g().length(e) < 400
                    && math::ls(this->g().coverage(e) / this->g().coverage(base), 15.)) {
                return false;
            }
        }
        return true;
    }

    bool CheckForECAround(EdgeId e) const {
        return CheckAlternativesForEC(
                this->g().IncidentEdges(this->g().EdgeStart(e)), e)
                && CheckAlternativesForEC(
                this->g().IncidentEdges(this->g().EdgeEnd(e)), e);
    }

 public:

    AdditionalMDAThornCondition(Graph& g, size_t uniqueness_length)
            : base(g),
              uniqueness_length_(uniqueness_length) {
    }

    bool Check(EdgeId e) const override {
        return CheckUniqueCondition(e) || CheckForECAround(e);
    }

 private:
    DECL_LOGGER("AdditionalMDAThornCondition");
};

//todo move to rnaSPAdes simplification
template<class Graph>
class ECLoopRemover : public EdgeProcessingAlgorithm<Graph> {
    typedef std::less<typename Graph::EdgeId> Comparator;
    typedef EdgeProcessingAlgorithm<Graph> base;
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;

    double ec_threshold_;
    double relative_threshold_;
    const FlankingCoverage<Graph> &flanking_coverage_;
    EdgeRemover<Graph> edge_remover_;
    size_t coverage_loops_removed = 0;
    size_t dead_loops_removed = 0;
    size_t not_dead_loops_removed = 0;
    size_t coverage_rc_loops_removed = 0;
    size_t dead_rc_loops_removed = 0;
    size_t not_dead_rc_loops_removed = 0;

    bool IsLoop(EdgeId e) {
        return this->g().EdgeStart(e) == this->g().EdgeEnd(e);
    }

    bool IsRCLoop(EdgeId e) {
        return this->g().EdgeStart(e) == this->g().conjugate(this->g().EdgeEnd(e));
    }

    bool IsAnyLoop(EdgeId e) {
        return IsRCLoop(e) || IsLoop(e);
    }

    void RemoveHiddenLoopEC(EdgeId e, bool break_on_end) {
        if (IsLoop(e))
            coverage_loops_removed++;
        else
            coverage_rc_loops_removed++;
        if (this->g().length(e) <= this->g().k())
            edge_remover_.DeleteEdge(e);
        else {
            if (break_on_end) {
                auto split_result = this->g().SplitEdge(e, this->g().length(e) - this->g().k());
                edge_remover_.DeleteEdge(split_result.second);
            } else {
                auto split_result = this->g().SplitEdge(e, this->g().k());
                edge_remover_.DeleteEdge(split_result.first);
            }
        }

    }
    void RemoveLoopWithNoCheck(EdgeId e) {
        if (IsLoop(e)) {
            if (this->g().IncomingEdgeCount(this->g().EdgeStart(e)) == 1 || this->g().OutgoingEdgeCount(this->g().EdgeStart(e)) == 1)
                dead_loops_removed++;
            else
                not_dead_loops_removed++;
        } else {
            if (this->g().IncomingEdgeCount(this->g().EdgeStart(e)) == 2)
                dead_rc_loops_removed++;
            else
                not_dead_rc_loops_removed++;

        }
        edge_remover_.DeleteEdge(e);
    }

    bool FindHiddenLoopEC(EdgeId e) {
        if(flanking_coverage_.GetInCov(e) * relative_threshold_ < flanking_coverage_.GetOutCov(e) && flanking_coverage_.GetInCov(e) < ec_threshold_) {
            //start is bad, end is OK.
            RemoveHiddenLoopEC(e, false);
            return true;
        } else if(flanking_coverage_.GetOutCov(e) * relative_threshold_ < flanking_coverage_.GetInCov(e) && flanking_coverage_.GetOutCov(e) < ec_threshold_) {
            //end is bad, start is OK.
            RemoveHiddenLoopEC(e, true);
            return true;
        }
        RemoveLoopWithNoCheck(e);
        return false;
    }

    bool ProcessEdge(EdgeId e) {
        if (IsAnyLoop(e)) {
            DEBUG("Susp loop: " << this->g().int_id(e) << endl);
            bool res = FindHiddenLoopEC(e);
            if (res) {DEBUG ("was removed");} else {DEBUG("was not removed"); }
            return res;
        }
        return false;
    }


public:
    ECLoopRemover(Graph &g, const FlankingCoverage<Graph> &flanking_coverage, double ec_threshold, double relative_threshold,
                  EdgeRemovalHandlerF<Graph> removal_handler = 0): base(g),ec_threshold_(ec_threshold),
                                                                            relative_threshold_(relative_threshold), flanking_coverage_(flanking_coverage),
                                                                            edge_remover_(g, removal_handler){
    }
    void PrintLoopStats(){
        INFO("Loops: accurately removed/deadend removed/other: "<< coverage_loops_removed <<"/" << dead_loops_removed << "/" <<not_dead_loops_removed);
        INFO("RC loops: accurately removed/deadend removed/other: "<< coverage_rc_loops_removed <<"/" << dead_rc_loops_removed << "/" <<not_dead_rc_loops_removed);
    }
private:
    DECL_LOGGER("ECLoopRemover");
};

template<class Graph>
class MetaHiddenECRemover: public PersistentProcessingAlgorithm<Graph, typename Graph::VertexId> {
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    typedef PersistentProcessingAlgorithm<Graph, VertexId> base;
    const FlankingCoverage<Graph>& flanking_coverage_;
    size_t uniqueness_length_;
    double relative_threshold_;

    EdgeDisconnector<Graph> disconnector_;

    void DisconnectEdges(VertexId v) {
        while (!this->g().IsDeadEnd(v)) {
            disconnector_(*(this->g().out_begin(v)), /*compress*/false);
        }
    }

    bool CheckUniqueness(EdgeId e) {
        return UniquePathLengthLowerBound(this->g(), uniqueness_length_)(e);
    }

    void ProcessHiddenEC(VertexId v) {
        VERIFY(this->g().OutgoingEdgeCount(v) == 2);
        vector<EdgeId> edges(this->g().out_begin(v), this->g().out_end(v));
        if (math::gr(flanking_coverage_.CoverageOfStart(edges.front()),
                    flanking_coverage_.CoverageOfStart(edges.back()))) {
            std::swap(edges.front(), edges.back());
        }
        double c1 = flanking_coverage_.CoverageOfStart(edges.front());
        double c2 = flanking_coverage_.CoverageOfStart(edges.back());
        TRACE("c1 " << c1 << "; c2 " << c2);
        if (math::ls(c1 * relative_threshold_, c2)) {
            TRACE("Disconnecting " << this->g().str(edges.front()));
            disconnector_(edges.front());
        } else {
            TRACE("Disconnecting " << this->g().str(edges.front()) << " and " << this->g().str(edges.back()));
            DisconnectEdges(v);
        }
    }

    bool CheckSuspicious(VertexId v) {
        if (this->g().IncomingEdgeCount(v) != 1 || this->g().OutgoingEdgeCount(v) != 2) {
            return false;
        }
        vector<EdgeId> edges;
        utils::push_back_all(edges, this->g().OutgoingEdges(v));
        VERIFY(edges.size() == 2);
        if (this->g().conjugate(edges[0]) != edges[1]) {
            return false;
        }
        return CheckUniqueness(this->g().GetUniqueIncomingEdge(v));
    }

protected:

    bool Process(VertexId v) override {
        if (CheckSuspicious(v)) {
            ProcessHiddenEC(v);
            return true;
        }
        return false;
    }

public:
    MetaHiddenECRemover(Graph& g, size_t chunk_cnt,
                    const FlankingCoverage<Graph> &flanking_coverage,
                    size_t uniqueness_length,
                    double relative_threshold,
                    EdgeRemovalHandlerF<Graph> removal_handler = 0)
            : base(g, nullptr, /*canonical only*/ false, std::less<VertexId>(), /*track changes*/false), 
              flanking_coverage_(flanking_coverage),
              uniqueness_length_(uniqueness_length),
              relative_threshold_(relative_threshold),
              disconnector_(g, removal_handler, g.k() + 1) {
        this->interest_el_finder_ = std::make_shared<ParallelInterestingElementFinder<Graph, VertexId>>(
                [&](VertexId v) {
                    return CheckSuspicious(v);
                }, chunk_cnt);
    }

private:
    DECL_LOGGER("MetaHiddenECRemover");
};

//be careful unreliability_threshold_ is dependent on ec_threshold_!
template<class Graph>
class HiddenECRemover: public PersistentProcessingAlgorithm<Graph, typename Graph::VertexId> {
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    typedef PersistentProcessingAlgorithm<Graph, VertexId> base;
    const FlankingCoverage<Graph>& flanking_coverage_;
    size_t uniqueness_length_;
    double unreliability_threshold_;
    double ec_threshold_;
    double relative_threshold_;

    EdgeDisconnector<Graph> disconnector_;

    void DisconnectEdges(VertexId v) {
        while (!this->g().IsDeadEnd(v)) {
            disconnector_(*(this->g().out_begin(v)), /*compress*/false);
        }
    }

    bool CheckUniqueness(EdgeId e) {
        //todo why 8???
        omnigraph::MultiplicityCounter<Graph> mult_counter(this->g(), uniqueness_length_, 8);

        vector<EdgeId> edges;
        utils::push_back_all(edges, this->g().OutgoingEdges(this->g().EdgeEnd(e)));
        VERIFY(edges.size() == 2);
        return (this->g().conjugate(edges[0]) == edges[1] && mult_counter.count(e, this->g().EdgeStart(e)) <= 1) ||
                this->g().length(e) >= uniqueness_length_;
    }

    bool ProcessHiddenEC(VertexId v) {
        TRACE("Processing outgoing edges for vertex " << this->g().str(v));
        VERIFY(this->g().OutgoingEdgeCount(v) == 2)
        vector<EdgeId> edges(this->g().out_begin(v), this->g().out_end(v));
        if (math::gr(flanking_coverage_.CoverageOfStart(edges.front()),
                    flanking_coverage_.CoverageOfStart(edges.back()))) {
            std::swap(edges.front(), edges.back());
        }
        double c1 = flanking_coverage_.CoverageOfStart(edges.front());
        TRACE("Flank start of e1 " << this->g().str(edges.front()) << ": " << c1);
        double c2 = flanking_coverage_.CoverageOfStart(edges.back());
        TRACE("Flank start of e1 " << this->g().str(edges.back()) << ": " << c2);
        if (math::ls(c2, unreliability_threshold_)) {
            TRACE("Disconnecting both edges from vertex " << this->g().str(v));
            DisconnectEdges(v);
            return true;
        }
        if (math::ls(c1 * relative_threshold_, c2) && math::ls(c1, ec_threshold_)) {
            TRACE("Disconnecting edge " << this->g().str(edges.front()) << " from vertex " << this->g().str(v));
            disconnector_(edges.front());
            return true;
        }
        return false;
    }

    bool CheckSuspicious(VertexId v) {
        if (this->g().IncomingEdgeCount(v) != 1 || this->g().OutgoingEdgeCount(v) != 2) {
            return false;
        }
        return CheckUniqueness(this->g().GetUniqueIncomingEdge(v));
    }

protected:

    bool Process(VertexId v) override {
        if (CheckSuspicious(v)) {
            return ProcessHiddenEC(v);
        }
        return false;
    }

public:
    HiddenECRemover(Graph& g, size_t chunk_cnt,
                    const FlankingCoverage<Graph> &flanking_coverage,
                    size_t uniqueness_length,
                    double unreliability_coeff,
                    double ec_threshold, double relative_threshold,
                    EdgeRemovalHandlerF<Graph> removal_handler = 0)
            : base(g, nullptr, /*canonical only*/ false, std::less<VertexId>(), /*track changes*/false), 
              flanking_coverage_(flanking_coverage),
              uniqueness_length_(uniqueness_length),
              unreliability_threshold_(unreliability_coeff * ec_threshold), ec_threshold_(ec_threshold),
              relative_threshold_(relative_threshold),
              disconnector_(g, removal_handler, g.k() + 1) {
        VERIFY(math::gr(unreliability_coeff, 0.));
        this->interest_el_finder_ = std::make_shared<ParallelInterestingElementFinder<Graph, VertexId>>(
                [&](VertexId v) {
                    return CheckSuspicious(v);
                }, chunk_cnt);
    }

private:
    DECL_LOGGER("HiddenECRemover");
};

template<class Graph>
class SelfConjugateDisruptor: public EdgeProcessingAlgorithm<Graph> {
    typedef EdgeProcessingAlgorithm<Graph> base;
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    size_t max_repeat_len_;
    EdgeRemover<Graph> edge_remover_;

protected:

    bool ProcessEdge(EdgeId e) override {
        Graph& g = this->g();
        if (e == g.conjugate(e)) {
            TRACE("Disrupting self-conjugate edge " << g.str(e));
            UniquePathFinder<Graph> unique_path_finder(g);
            size_t e_len = g.length(e);
            size_t induced_repeat_len =
                    CumulativeLength(g, unique_path_finder.UniquePathBackward(e)) - (e_len / 2);
            if (induced_repeat_len > max_repeat_len_ || g.OutgoingEdgeCount(g.EdgeEnd(e)) == 0) {
                EdgeId to_del = e;
                if (e_len > 1) {
                    to_del = g.SplitEdge(e, e_len / 2).second;
                }
                edge_remover_.DeleteEdge(to_del);
                return true;
            }
        }
        return false;
    }

public:
    SelfConjugateDisruptor(Graph& g, size_t max_repeat_len,
                           std::function<void(EdgeId)> removal_handler = 0)
            : base(g, true),
              max_repeat_len_(max_repeat_len),
              edge_remover_(g, removal_handler) {
    }

private:
    DECL_LOGGER("SelfConjugateDisruptor");
};
}
