//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

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

template<class Graph, class Comparator = std::less<typename Graph::EdgeId>>
class ChimericEdgeRemovingAlgorithm : public EdgeRemovingAlgorithm<Graph,
        Comparator> {
    typedef EdgeRemovingAlgorithm<Graph, Comparator> base;
    typedef typename Graph::EdgeId EdgeId;

//    shared_ptr<func::Predicate<EdgeId>> remove_condition_;
//    boost::function<void(EdgeId)> removal_handler_;

 public:

    ChimericEdgeRemovingAlgorithm(
            Graph& g,
            shared_ptr<func::Predicate<EdgeId>> remove_condition,
            boost::function<void(EdgeId)> removal_handler = boost::none,
            const Comparator& c = Comparator(),
            shared_ptr<func::Predicate<EdgeId>> proceed_condition = make_shared<
                    func::AlwaysTrue<EdgeId>>())
            : base(g,
                   func::And<EdgeId>(
                           make_shared<AlternativesPresenceCondition<Graph>>(g),
                           remove_condition),
                   removal_handler, c, proceed_condition) {
    }

 private:
    DECL_LOGGER("EdgeRemovingAlgorithm")
    ;
};

//todo refactor
template<class Graph>
class IterativeLowCoverageEdgeRemover : public ChimericEdgeRemovingAlgorithm<
        Graph, CoverageComparator<Graph>> {
 private:
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    typedef ChimericEdgeRemovingAlgorithm<Graph, CoverageComparator<Graph>> base;

 public:
    IterativeLowCoverageEdgeRemover(
            Graph &g, double max_coverage,
            shared_ptr<Predicate<EdgeId>> condition,
            boost::function<void(EdgeId)> removal_handler)
            : base(g, condition, removal_handler, CoverageComparator<Graph>(g),
                   make_shared<CoverageUpperBound<Graph>>(g, max_coverage)) {
    }
};

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
//        return CheckAdjacent(g.AdjacentEdges(g.EdgeStart(e)), e)
//                && CheckAdjacent(g.AdjacentEdges(g.EdgeEnd(e)), e);
//    }
//
// private:
//    DECL_LOGGER("RelativeCoverageCondition")
//    ;
//
//};

//template<class Graph>
//class RelativeLowCoverageEdgeRemover : public ChimericEdgeRemovingAlgorithm<
//        Graph, CoverageComparator<Graph>> {
// private:
//    typedef typename Graph::EdgeId EdgeId;
//    typedef typename Graph::VertexId VertexId;
//    typedef ChimericEdgeRemovingAlgorithm<Graph, CoverageComparator<Graph>> base;
//
// public:
//    RelativeLowCoverageEdgeRemover(
//            Graph& g, size_t max_length, double max_coverage,
//            double coverage_gap, boost::function<void(EdgeId)> removal_handler)
//            : base(g,
//                   func::And<EdgeId>(
//                           make_shared<RelativeCoverageCondition<Graph>>(
//                                   g, coverage_gap),
//                           make_shared<LengthUpperBound<Graph>>(g, max_length)),
//                   removal_handler, CoverageComparator<Graph>(g),
//                   make_shared<CoverageUpperBound<Graph>>(g, max_coverage)) {
//    }
//};

template<class Graph>
class TopologyAndReliablityBasedChimericEdgeRemover :
        public ChimericEdgeRemovingAlgorithm<Graph, LengthComparator<Graph>> {
 private:
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    typedef ChimericEdgeRemovingAlgorithm<Graph, LengthComparator<Graph>> base;

 public:
    TopologyAndReliablityBasedChimericEdgeRemover(
            Graph& g, size_t max_length, size_t uniqueness_length,
            double max_coverage, boost::function<void(EdgeId)> removal_handler)
            : base(g,
                   func::And<EdgeId>(
                           make_shared<CoverageUpperBound<Graph>>(g,
                                                                  max_coverage),
                           make_shared<
                                   PredicateUniquenessPlausabilityCondition<
                                           Graph>>(
                                   g,
                                   /*uniqueness*/
                                   MakePathLengthLowerBound(
                                           g, UniquePathFinder<Graph>(g),
                                           uniqueness_length),
                                   /*plausibility*/make_shared<
                                           func::AlwaysTrue<EdgeId>>())),
                   removal_handler, LengthComparator<Graph>(g),
                   make_shared<LengthUpperBound<Graph>>(g, max_length)) {
    }
};

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
        dij.run(this->g().EdgeStart(e));
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
        FOREACH (EdgeId e, edges) {
            if (e != base && this->g().length(e) < 400
                    && this->g().coverage(e) < 15 * this->g().coverage(base)) {
                return false;
            }
        }
        return true;
    }

    bool CheckCoverageAround(EdgeId e) const {
        return CheckAlternativeCoverage(
                this->g().AdjacentEdges(this->g().EdgeStart(e)), e)
                && CheckAlternativeCoverage(
                        this->g().AdjacentEdges(this->g().EdgeEnd(e)), e);
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
class ThornRemover : public ChimericEdgeRemovingAlgorithm<Graph,
    CoverageComparator<Graph>> {
 private:
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    typedef ChimericEdgeRemovingAlgorithm<Graph, CoverageComparator<Graph>> base;

 public:
    ThornRemover(Graph& g, size_t max_length, size_t uniqueness_length,
                 size_t dijkstra_depth,
                 boost::function<void(EdgeId)> removal_handler)
            : base(g,
                   func::And<EdgeId>(make_shared<LengthUpperBound<Graph>>(g, max_length),
                             make_shared<ThornCondition<Graph>>(g, uniqueness_length,
                             dijkstra_depth)),
                   removal_handler, CoverageComparator<Graph>(g)) {
    }
};

//todo rename
template<class Graph>
class TopologyChimericEdgeRemover : public ChimericEdgeRemovingAlgorithm<Graph,
        LengthComparator<Graph>> {
 private:
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    typedef ChimericEdgeRemovingAlgorithm<Graph, LengthComparator<Graph>> base;

 public:
    TopologyChimericEdgeRemover(Graph& g, size_t max_length,
                                size_t uniqueness_length,
                                size_t plausibility_length,
                                boost::function<void(EdgeId)> removal_handler)
            : base(g,
                   make_shared<DefaultUniquenessPlausabilityCondition<Graph>>(
                           g, uniqueness_length, plausibility_length),
                   removal_handler, LengthComparator<Graph>(g),
                   make_shared<LengthUpperBound<Graph>>(g, max_length)) {
    }
};

template<class Graph>
class MultiplicityCountingCondition : public UniquenessPlausabilityCondition<
        Graph> {
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    typedef shared_ptr<Predicate<EdgeId>> EdgePredicate;
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
        return plausiblity_condition_->Check(e);
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
class SimpleMultiplicityCountingChimericEdgeRemover :
        public ChimericEdgeRemovingAlgorithm<Graph, LengthComparator<Graph>> {
 private:
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    typedef ChimericEdgeRemovingAlgorithm<Graph, LengthComparator<Graph>> base;

 public:
    SimpleMultiplicityCountingChimericEdgeRemover(
            Graph& g, size_t max_length, size_t uniqueness_length,
            size_t plausibility_length,
            boost::function<void(EdgeId)> removal_handler)
            : base(g,
                   make_shared<MultiplicityCountingCondition<Graph>>(
                           g,
                           uniqueness_length,
                           /*plausibility*/MakePathLengthLowerBound(
                                   g,
                                   PlausiblePathFinder<Graph>(
                                           g, 2 * plausibility_length),
                                   plausibility_length)),
                   removal_handler, LengthComparator<Graph>(g),
                   make_shared<LengthUpperBound<Graph>>(g, max_length)) {
    }
};

template<class Graph>
class HiddenECRemover: public EdgeProcessingAlgorithm<Graph,
		std::less<typename Graph::EdgeId>> {
	typedef std::less<typename Graph::EdgeId> Comparator;
	typedef EdgeProcessingAlgorithm<Graph, Comparator> base;
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
private:
	size_t uniqueness_length_;
	double unreliability_threshold_;
	double ec_threshold_;
	double relative_threshold_;
	const AbstractFlankingCoverage<Graph> &flanking_coverage_;
	EdgeRemover<Graph> edge_remover_;
	shared_ptr<MultiplicityCountingCondition<Graph>> condition_;
private:
	void RemoveHiddenEC(EdgeId edge) {
		if (this->g().length(edge) <= this->g().k())
			edge_remover_.DeleteEdge(edge);
		else {
			auto split_result = this->g().SplitEdge(edge, this->g().k());
			edge_remover_.DeleteEdge(split_result.first);
		}
	}

	void RemoveHiddenECWithNoCompression(EdgeId edge) {
		if (this->g().length(edge) <= this->g().k())
			edge_remover_.DeleteEdgeWithNoCompression(edge);
		else {
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
		return (edges.size() == 2 && this->g().conjugate(edges[0]) == edges[1] && condition_->CheckUniqueness(this->g().GetUniqueIncomingEdge(v), false)) || this->g().length(this->g().GetUniqueIncomingEdge(v)) >= uniqueness_length_;
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
                    boost::function<void(EdgeId)> removal_handler = 0)
            : base(g, Comparator(), make_shared<func::AlwaysTrue<EdgeId>>()), uniqueness_length_(uniqueness_length),
              unreliability_threshold_(unreliability_threshold * ec_threshold), ec_threshold_(ec_threshold),
              relative_threshold_(relative_threshold), flanking_coverage_(flanking_coverage),
              edge_remover_(g, removal_handler),
              condition_(new MultiplicityCountingCondition<Graph>(g, uniqueness_length,
                              make_shared<func::AlwaysTrue<EdgeId>>())) {

    }

private:
	DECL_LOGGER("HiddenECRemover");
};

template<class Graph>
class LowCoveredSelfConjEdgeRemovingAlgorithm : public EdgeRemovingAlgorithm<Graph,
        CoverageComparator<Graph>> {
    typedef EdgeRemovingAlgorithm<Graph, CoverageComparator<Graph>> base;
    typedef typename Graph::EdgeId EdgeId;

 public:

    LowCoveredSelfConjEdgeRemovingAlgorithm(
            Graph &g, size_t max_length, double max_coverage,
            boost::function<void(EdgeId)> removal_handler)
            : base(g, func::And<EdgeId>(make_shared<SelfConjugateCondition<Graph>>(g), make_shared<LengthUpperBound<Graph>>(g, max_length)),
                   removal_handler, CoverageComparator<Graph>(g),
                   make_shared<CoverageUpperBound<Graph>>(g, max_coverage)) {
    }

 private:
    DECL_LOGGER("LowCoveredSelfConjEdgeRemovingAlgorithm");
};

}
