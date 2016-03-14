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

#include "graph_processing_algorithm.hpp"
#include "basic_edge_conditions.hpp"
#include "omni_tools.hpp"
#include "omni_utils.hpp"
#include "func.hpp"
#include "xmath.h"
#include "dijkstra_tools/dijkstra_helper.hpp"

namespace omnigraph {

template<class Graph>
pred::TypedPredicate<typename Graph::EdgeId>
NecessaryECCondition(const Graph& g, size_t max_length, double max_coverage) {
    return AddAlternativesPresenceCondition(g, pred::And(LengthUpperBound<Graph>(g, max_length),
                                                        CoverageUpperBound<Graph>(g, max_coverage)));
}

template<class Graph>
bool RemoveErroneousEdgesInCoverageOrder(Graph &g,
                                         pred::TypedPredicate<typename Graph::EdgeId> removal_condition,
                                         double max_coverage,
                                         std::function<void(typename Graph::EdgeId)> removal_handler) {
    omnigraph::EdgeRemovingAlgorithm<Graph> erroneous_edge_remover(g,
                                                                   AddAlternativesPresenceCondition(g, removal_condition),
                                                                   removal_handler);

    return erroneous_edge_remover.Run(CoverageComparator<Graph>(g),
                                      CoverageUpperBound<Graph>(g, max_coverage));
}

template<class Graph>
bool RemoveErroneousEdgesInLengthOrder(Graph &g,
                                       pred::TypedPredicate<typename Graph::EdgeId> removal_condition,
                                       size_t max_length,
                                       std::function<void(typename Graph::EdgeId)> removal_handler) {
    omnigraph::EdgeRemovingAlgorithm<Graph> erroneous_edge_remover(g,
                                                                   AddAlternativesPresenceCondition(g, removal_condition),
                                                                   removal_handler);

    return erroneous_edge_remover.Run(LengthComparator<Graph>(g),
                                      LengthUpperBound<Graph>(g, max_length));
}

template<class Graph>
class SelfConjugateCondition : public EdgeCondition<Graph> {
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    typedef EdgeCondition<Graph> base;

 public:

    SelfConjugateCondition(const Graph& g)
            : base(g) {
    }

    bool Check(EdgeId e) const {
        return e == this->g().conjugate(e);
    }

 private:
    DECL_LOGGER("SelfConjugateCondition");
};

//coverage comparator
//template<class Graph>
//class RelativeCoverageCondition : public EdgeCondition<Graph> {
//    typedef typename Graph::EdgeId EdgeId;
//    typedef typename Graph::VertexId VertexId;
//    typedef EdgeCondition<Graph> base;
//
//    double min_coverage_gap_;
//
//    bool StrongNeighbourCondition(EdgeId neighbour_edge,
//                                  EdgeId possible_ec) const {
//        return neighbour_edge == possible_ec
//                || math::gr(this->g().coverage(neighbour_edge),
//                            this->g().coverage(possible_ec) * min_coverage_gap_);
////	              || this->g().length(neighbour_edge)
////	                      >= neighbour_length_threshold_;
//    }
//
//    bool CheckAdjacent(const vector<EdgeId>& edges, EdgeId possible_ec) const {
//        FOREACH (EdgeId e, edges) {
//            if (!StrongNeighbourCondition(e, possible_ec))
//                return false;
//        }
//        return true;
//    }
//
// public:
//
//    RelativeCoverageCondition(const Graph& g, double min_coverage_gap)
//            : base(g),
//              min_coverage_gap_(min_coverage_gap) {
//
//    }
//
//    bool Check(EdgeId e) const {
//        const Graph& g = this->g();
//        return CheckAdjacent(g.IncidentEdges(g.EdgeStart(e)), e)
//                && CheckAdjacent(g.IncidentEdges(g.EdgeEnd(e)), e);
//    }
//
// private:
//    DECL_LOGGER("RelativeCoverageCondition")
//    ;
//
//};

//todo refactor
template<class Graph>
class ThornCondition : public EdgeCondition<Graph> {
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    typedef EdgeCondition<Graph> base;

    size_t uniqueness_length_;
    size_t dijkstra_depth_;

    bool Unique(const vector<EdgeId>& edges, bool forward) const {
        return edges.size() == 1 && CheckUniqueness(*edges.begin(), forward);
    }

    bool CheckUnique(EdgeId e) const {
        TRACE("Checking conditions for edge start");
        return Unique(vector<EdgeId>(this->g().in_begin(this->g().EdgeStart(e)), this->g().in_end(this->g().EdgeStart(e))), false)
                || Unique(vector<EdgeId>(this->g().out_begin(this->g().EdgeEnd(e)), this->g().out_end(this->g().EdgeEnd(e))), true);
    }

    bool CheckThorn(EdgeId e) const {
        if (this->g().EdgeStart(e) == this->g().EdgeEnd(e))
            return false;
        if (this->g().RelatedVertices(this->g().EdgeStart(e),
                                      this->g().EdgeEnd(e))) {
            return true;
        }
        if (this->g().OutgoingEdgeCount(this->g().EdgeStart(e)) != 2)
            return false;
        if (this->g().IncomingEdgeCount(this->g().EdgeStart(e)) != 1)
            return false;
        if (this->g().OutgoingEdgeCount(this->g().EdgeEnd(e)) != 1)
            return false;
        if (this->g().IncomingEdgeCount(this->g().EdgeEnd(e)) != 2)
            return false;

        auto dij = DijkstraHelper<Graph>::CreateBoundedDijkstra(this->g(), dijkstra_depth_);
        dij.Run(this->g().EdgeStart(e));
        vector<VertexId> reached = dij.ReachedVertices();
        for (auto it = reached.begin(); it != reached.end(); ++it) {
            if (*it != this->g().EdgeEnd(e)
                    && this->g().RelatedVertices(*it, this->g().EdgeEnd(e))) {
                return true;
            }
        }
        return false;
    }

    template<class EdgeContainer>
    bool CheckAlternativeCoverage(const EdgeContainer& edges, EdgeId base) const {
        for (EdgeId e: edges) {
            if (e != base && this->g().length(e) < 400
                    && this->g().coverage(e) < 15 * this->g().coverage(base)) {
                return false;
            }
        }
        return true;
    }

    bool CheckCoverageAround(EdgeId e) const {
        return CheckAlternativeCoverage(
                this->g().IncidentEdges(this->g().EdgeStart(e)), e)
                && CheckAlternativeCoverage(
                        this->g().IncidentEdges(this->g().EdgeEnd(e)), e);
    }

    bool CheckUniqueness(EdgeId e, bool /*forward*/) const {
        return this->g().length(e) >= uniqueness_length_;
    }

 public:

    ThornCondition(Graph& g, size_t uniqueness_length, size_t dijkstra_depth)
            : base(g),
              uniqueness_length_(uniqueness_length),
              dijkstra_depth_(dijkstra_depth) {
    }

    bool Check(EdgeId e) const {
        bool tmp = (CheckUnique(e) || CheckCoverageAround(e));
        if (tmp)
            tmp &= CheckThorn(e);
        return tmp;
    }

 private:
    DECL_LOGGER("ThornCondition")
    ;

};

template<class Graph>
class MultiplicityCountingCondition : public UniquenessPlausabilityCondition<Graph> {
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    typedef pred::TypedPredicate<EdgeId> EdgePredicate;
    typedef UniquenessPlausabilityCondition<Graph> base;

    MultiplicityCounter<Graph> multiplicity_counter_;
    EdgePredicate plausiblity_condition_;

public:
    bool CheckUniqueness(EdgeId e, bool forward) const {
        TRACE( "Checking " << this->g().int_id(e) << " for uniqueness in " << (forward ? "forward" : "backward") << " direction");
        VertexId start =
                forward ? this->g().EdgeEnd(e) : this->g().EdgeStart(e);
        bool result = multiplicity_counter_.count(e, start) <= 1;
        TRACE( "Edge " << this->g().int_id(e) << " is" << (result ? "" : " not") << " unique");
        return result;
    }

    bool CheckPlausibility(EdgeId e, bool) const {
        return plausiblity_condition_(e);
    }

    MultiplicityCountingCondition(const Graph& g, size_t uniqueness_length,
                                  EdgePredicate plausiblity_condition)
            :
              //todo why 8???
              base(g),
              multiplicity_counter_(g, uniqueness_length, 8),
              plausiblity_condition_(plausiblity_condition) {

    }

 private:

    DECL_LOGGER("MultiplicityCountingCondition")
    ;
};


template<class Graph>
class ECLoopRemover : public EdgeProcessingAlgorithm<Graph> {
    typedef std::less<typename Graph::EdgeId> Comparator;
    typedef EdgeProcessingAlgorithm<Graph> base;
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;

    double ec_threshold_;
    double relative_threshold_;
    const AbstractFlankingCoverage<Graph> &flanking_coverage_;
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
    ECLoopRemover(Graph &g, const AbstractFlankingCoverage<Graph> &flanking_coverage, double ec_threshold, double relative_threshold,
                  HandlerF<Graph> removal_handler = 0): base(g),ec_threshold_(ec_threshold),
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
class HiddenECRemover: public EdgeProcessingAlgorithm<Graph> {
	typedef EdgeProcessingAlgorithm<Graph> base;
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
private:
	size_t uniqueness_length_;
	double unreliability_threshold_;
	double ec_threshold_;
	double relative_threshold_;
	const AbstractFlankingCoverage<Graph> &flanking_coverage_;
	EdgeRemover<Graph> edge_remover_;
	MultiplicityCountingCondition<Graph> condition_;
private:
	void RemoveHiddenEC(EdgeId edge) {
		if (this->g().length(edge) <= this->g().k() || (edge == this->g().conjugate(edge) && this->g().length(edge) <= 2 * this->g().k()))
			edge_remover_.DeleteEdge(edge);
		else {
			auto split_result = this->g().SplitEdge(edge, this->g().k());
			edge_remover_.DeleteEdge(split_result.first);
		}
	}

	void RemoveHiddenECWithNoCompression(EdgeId edge) {
		if (this->g().length(edge) <= this->g().k() || (edge == this->g().conjugate(edge) && this->g().length(edge) <= 2 * this->g().k())) {
			edge_remover_.DeleteEdgeWithNoCompression(edge);
		} else {
			auto split_result = this->g().SplitEdge(edge, this->g().k());
			edge_remover_.DeleteEdgeWithNoCompression(split_result.first);
		}
	}

	void DisconnectEdges(VertexId v) {
		while(!this->g().IsDeadEnd(v)) {
			RemoveHiddenECWithNoCompression(*(this->g().out_begin(v)));
		}
	}

	bool FindHiddenEC(VertexId v) {
		vector<EdgeId> edges(this->g().out_begin(v), this->g().out_end(v));
		if(flanking_coverage_.GetInCov(edges[0]) > flanking_coverage_.GetInCov(edges[1])) {
			auto tmp = edges[0];
			edges[0] = edges[1];
			edges[1] = tmp;
		}
//		cout << flanking_coverage_.GetInCov(edges[0]) << " " << flanking_coverage_.GetInCov(edges[1]) << endl;
		if(flanking_coverage_.GetInCov(edges[1]) < unreliability_threshold_) {
			DisconnectEdges(v);
//			cout << "disconnected" << endl;
			return true;
		}
		if(flanking_coverage_.GetInCov(edges[0]) * relative_threshold_ < flanking_coverage_.GetInCov(edges[1]) && flanking_coverage_.GetInCov(edges[0]) < ec_threshold_) {
			RemoveHiddenEC(edges[0]);
//			cout << "success" << endl;
			return true;
		}
		return false;
	}

	bool CheckSuspicious(VertexId v) {
		if (this->g().IncomingEdgeCount(v) != 1 || this->g().OutgoingEdgeCount(v) != 2) {
			return false;
		}
		vector<EdgeId> edges(this->g().out_begin(v), this->g().out_end(v));
		return (edges.size() == 2 && this->g().conjugate(edges[0]) == edges[1] && condition_.CheckUniqueness(this->g().GetUniqueIncomingEdge(v), false)) || this->g().length(this->g().GetUniqueIncomingEdge(v)) >= uniqueness_length_;
	}

	bool ProcessEdge(EdgeId e) {
		VertexId v = this->g().EdgeEnd(e);
		if(CheckSuspicious(v)) {
//			cout << "client: " << this->g().int_id(v) << endl;
			return FindHiddenEC(v);
		}
		return false;
	}

public:
    HiddenECRemover(Graph& g, size_t uniqueness_length,
                    const AbstractFlankingCoverage<Graph> &flanking_coverage,
                    double unreliability_threshold, double ec_threshold,
                    double relative_threshold,
                    std::function<void(EdgeId)> removal_handler = 0)
            : base(g), uniqueness_length_(uniqueness_length),
              unreliability_threshold_(unreliability_threshold * ec_threshold), ec_threshold_(ec_threshold),
              relative_threshold_(relative_threshold), flanking_coverage_(flanking_coverage),
              edge_remover_(g, removal_handler),
              condition_(g, uniqueness_length, pred::AlwaysTrue<EdgeId>()) {
    }

private:
	DECL_LOGGER("HiddenECRemover");
};

template<class Graph>
class SelfConjugateDisruptor: public EdgeProcessingAlgorithm<Graph> {
    typedef EdgeProcessingAlgorithm<Graph> base;
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    EdgeRemover<Graph> edge_remover_;
protected:

    bool ProcessEdge(EdgeId e) override {
        if (e == this->g().conjugate(e)) {
            TRACE("Disrupting self-conjugate edge " << this->g().str(e));
            EdgeId to_del = e;
            size_t len = this->g().length(e);
            if (len > 1) {
                to_del = this->g().SplitEdge(e, len / 2).second;
            }
            edge_remover_.DeleteEdge(to_del);
            return true;
        }
        return false;
    }

public:
    SelfConjugateDisruptor(Graph& g,
                           std::function<void(EdgeId)> removal_handler = 0)
            : base(g, true), edge_remover_(g, removal_handler) {
    }

private:
    DECL_LOGGER("SelfConjugateDisruptor");
};
}
