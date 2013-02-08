//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
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

namespace omnigraph {

template<class Graph, class PathFinder>
class PathLengthLowerBound: public EdgeCondition<Graph> {
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
	typedef EdgeCondition<Graph> base;

	PathFinder path_finder_;
	size_t min_length_;

	ForwardDirection<Graph> forward_;
	BackwardDirection<Graph> backward_;

	size_t CummulativePathLength(EdgeId e,
			const AbstractDirection<Graph>& direction) const {
		return CummulativeLength(this->g(), path_finder_(e, direction));
	}

public:
	PathLengthLowerBound(const Graph& g, const PathFinder& path_finder,
			size_t min_length) :
			base(g), path_finder_(path_finder), min_length_(min_length), forward_(
					g), backward_(g) {

	}

	bool Check(EdgeId e) const {
		return std::max(CummulativePathLength(e, forward_),
				CummulativePathLength(e, backward_)) >= min_length_;
	}
};

template<class Graph, class PathFinder>
shared_ptr<PathLengthLowerBound<Graph, PathFinder>> MakePathLengthLowerBound(
		const Graph& g, const PathFinder& path_finder, size_t min_length) {
	return make_shared<PathLengthLowerBound<Graph, PathFinder>>(g, path_finder,
			min_length);
}

template<class Graph>
class UniquenessPlausabilityCondition: public EdgeCondition<Graph> {
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
	typedef EdgeCondition<Graph> base;

	virtual bool CheckUniqueness(EdgeId e, bool forward) const = 0;

	virtual bool CheckPlausibility(EdgeId e, bool forward) const = 0;

	bool SingleUnique(const vector<EdgeId>& edges, bool forward) const {
		return edges.size() == 1 && CheckUniqueness(*edges.begin(), forward);
	}

	bool ExistPlausible(EdgeId init_e, const vector<EdgeId>& edges,
			bool forward) const {
		FOREACH(EdgeId e, edges) {
			if (e == init_e)
				continue;
			if (CheckPlausibility(e, forward)) {
				return true;
			}
		}
		return false;
	}

	bool Check(EdgeId e, const AbstractDirection<Graph>& direction) const {
		return SingleUnique(direction.IncomingEdges(direction.EdgeStart(e)),
				!direction.IsForward())
				&& ExistPlausible(e,
						direction.OutgoingEdges(direction.EdgeStart(e)),
						direction.IsForward());
	}

public:

	UniquenessPlausabilityCondition(const Graph& g) :
			base(g) {

	}

	bool Check(EdgeId e) const {
		return Check(e, ForwardDirection<Graph>(this->g()))
				|| Check(e, BackwardDirection<Graph>(this->g()));
	}

};

template<class Graph>
class PredicateUniquenessPlausabilityCondition: public UniquenessPlausabilityCondition<
		Graph> {
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
	typedef shared_ptr<Predicate<EdgeId>> EdgePredicate;
	typedef UniquenessPlausabilityCondition<Graph> base;

	EdgePredicate uniqueness_condition_;
	EdgePredicate plausiblity_condition_;

	bool CheckUniqueness(EdgeId e, bool) const {
		return uniqueness_condition_->Check(e);
	}

	bool CheckPlausibility(EdgeId e, bool) const {
		return plausiblity_condition_->Check(e);
	}

public:

	PredicateUniquenessPlausabilityCondition(const Graph& g,
			EdgePredicate uniqueness_condition,
			EdgePredicate plausiblity_condition) :
			base(g), uniqueness_condition_(uniqueness_condition), plausiblity_condition_(
					plausiblity_condition) {
	}

};

template<class Graph>
class AlternativesPresenceCondition: public EdgeCondition<Graph> {
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
	typedef EdgeCondition<Graph> base;

public:

	AlternativesPresenceCondition(const Graph& g) :
			base(g) {

	}

	bool Check(EdgeId e) const {
		return this->g().OutgoingEdgeCount(this->g().EdgeStart(e)) > 1
				&& this->g().IncomingEdgeCount(this->g().EdgeEnd(e)) > 1;
	}

};

template<class Graph, class UniquePF, class PlausiblePF>
class NotRelatedVerticesCondition: public EdgeCondition<Graph> {
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
	typedef EdgeCondition<Graph> base;

public:

	NotRelatedVerticesCondition(const Graph& g) :
			base(g) {

	}

	bool Check(EdgeId e) const {
		return !(this->g().RelatedVertices(this->g().EdgeStart(e),
				this->g().EdgeEnd(e)));
	}

};

template<class Graph, class Comparator = std::less<typename Graph::EdgeId>>
class ChimericEdgeRemovingAlgorithm: public EdgeRemovingAlgorithm<Graph,
		Comparator> {
	typedef EdgeRemovingAlgorithm<Graph, Comparator> base;
	typedef typename Graph::EdgeId EdgeId;

	shared_ptr<func::Predicate<EdgeId>> remove_condition_;
	boost::function<void(EdgeId)> removal_handler_;

public:

	ChimericEdgeRemovingAlgorithm(Graph& g,
			shared_ptr<func::Predicate<EdgeId>> remove_condition,
			boost::function<void(EdgeId)> removal_handler = boost::none,
			const Comparator& c = Comparator(),
			shared_ptr<func::Predicate<EdgeId>> proceed_condition = make_shared<
					func::AlwaysTrue<EdgeId>>()) :
			base(g,
					func::And<EdgeId>(remove_condition,
							make_shared<AlternativesPresenceCondition<Graph>>(g)),
					removal_handler, c, proceed_condition) {
	}

private:
	DECL_LOGGER("EdgeRemovingAlgorithm")
	;
};

////todo isn't this one cheating?!!!
//default comparator used
template<class Graph>
class ChimericEdgeCondition: public EdgeCondition<Graph> {
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
	typedef EdgeCondition<Graph> base;

	size_t max_overlap_;

	bool CheckEnd(VertexId v) const {
		return this->g().OutgoingEdgeCount(v) == 1;
	}

	bool CheckStart(VertexId v) const {
		return this->g().IncomingEdgeCount(v) == 1;
	}

	bool CheckLength(EdgeId e) const {
		return this->g().length(e) >= this->g().k() - max_overlap_;
	}

public:

	ChimericEdgeCondition(const Graph& g, size_t max_overlap) :
			base(g), max_overlap_(max_overlap) {

	}

	bool Check(EdgeId e) const {
		return CheckEnd(this->g().EdgeEnd(e))
				&& CheckStart(this->g().EdgeEnd(e)) && CheckLength(e);
	}

};

template<class Graph>
class ChimericEdgesRemover: public ChimericEdgeRemovingAlgorithm<Graph> {
private:
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
	typedef ChimericEdgeRemovingAlgorithm<Graph> base;

public:
	ChimericEdgesRemover(Graph &g, size_t max_overlap,
			boost::function<void(EdgeId)> removal_handler) :
			base(g,
					func::And<EdgeId>(make_shared<LengthUpperBound<Graph>>(g, g.k()),
							make_shared<ChimericEdgeCondition<Graph>>(g,
									max_overlap)), removal_handler) {
	}
};

//todo refactor
template<class Graph>
class IterativeLowCoverageEdgeRemover: public ChimericEdgeRemovingAlgorithm<
		Graph, CoverageComparator<Graph>> {
private:
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
	typedef ChimericEdgeRemovingAlgorithm<Graph, CoverageComparator<Graph>> base;

public:
	IterativeLowCoverageEdgeRemover(Graph &g, double max_coverage,
			shared_ptr<Predicate<EdgeId>> condition,
			boost::function<void(EdgeId)> removal_handler) :
			base(g, condition, removal_handler, CoverageComparator<Graph>(g),
					make_shared<CoverageUpperBound<Graph>>(g, max_coverage)) {
	}
};

//coverage comparator
template<class Graph>
class RelativeCoverageCondition: public EdgeCondition<Graph> {
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
	typedef EdgeCondition<Graph> base;

	double min_coverage_gap_;

	bool CheckAlternativeCoverage(const vector<EdgeId> &edges, EdgeId e) const {
		for (auto it = edges.begin(); it != edges.end(); ++it) {
			//todo wtf?!!! remove 400 condition
			if (*it != e && this->g().length(*it) < 400
					&& this->g().coverage(*it)
							< min_coverage_gap_ * this->g().coverage(e)) {
				return false;
			}
		}

		return true;
	}

public:

	RelativeCoverageCondition(Graph& g, double min_coverage_gap) :
			base(g), min_coverage_gap_(min_coverage_gap) {

	}

	bool Check(EdgeId e) const {
		return CheckAlternativeCoverage(
				this->g().IncomingEdges(this->g().EdgeStart(e)), e)
				&& CheckAlternativeCoverage(
						this->g().OutgoingEdges(this->g().EdgeStart(e)), e)
				&& CheckAlternativeCoverage(
						this->g().IncomingEdges(this->g().EdgeEnd(e)), e)
				&& CheckAlternativeCoverage(
						this->g().OutgoingEdges(this->g().EdgeEnd(e)), e);
	}

private:
	DECL_LOGGER("RelativeCoverageCondition")
	;

};

template<class Graph>
class RelativeLowCoverageEdgeRemover: public ChimericEdgeRemovingAlgorithm<
		Graph, CoverageComparator<Graph>> {
private:
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
	typedef ChimericEdgeRemovingAlgorithm<Graph, CoverageComparator<Graph>> base;

public:
	RelativeLowCoverageEdgeRemover(Graph& g, size_t max_length,
			double max_coverage, double coverage_gap,
			boost::function<void(EdgeId)> removal_handler) :
			base(g,
					func::And<EdgeId>(
							make_shared<RelativeCoverageCondition<Graph>>(g,
									coverage_gap),
							make_shared<LengthUpperBound<Graph>>(g,
									max_length)), removal_handler,
					CoverageComparator<Graph>(g),
					make_shared<CoverageUpperBound<Graph>>(g, max_coverage)) {
	}
};

template<class T>
void Append(vector<T>& current, const vector<T>& to_append) {
	Append(current, to_append.begin(), to_append.end());
}

template<class T>
void Append(vector<T>& current, const typename vector<T>::const_iterator begin,
		const typename vector<T>::const_iterator end) {
	current.insert(current.end(), begin, end);
}

template<class T, class It>
void Append(vector<T>& current, It begin, It end) {
	current.reserve(current.size() + (end - begin));
	while (begin != end)
		current.push_back(*begin++);
}

template<class Graph>
class CheatingChimericEdgeCondition: public EdgeCondition<Graph> {
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
	typedef EdgeCondition<Graph> base;

	double coverage_gap_;
	size_t neighbour_length_threshold_;

	bool StrongNeighbourCondition(EdgeId neighbour_edge, EdgeId possible_ec) const {
		return neighbour_edge == possible_ec
				|| math::gr(this->g().coverage(neighbour_edge),
						this->g().coverage(possible_ec) * coverage_gap_)
				|| this->g().length(neighbour_edge)
						>= neighbour_length_threshold_;
	}

	bool CheckAdjacent(const vector<EdgeId>& edges, EdgeId possible_ec) const {
		FOREACH (EdgeId e, edges) {
			if (!StrongNeighbourCondition(e, possible_ec))
				return false;
		}
		return true;
	}

public:

	CheatingChimericEdgeCondition(Graph& g, double coverage_gap,
			size_t neighbour_length_threshold) :
			base(g), coverage_gap_(coverage_gap), neighbour_length_threshold_(
					neighbour_length_threshold) {

	}

	bool Check(EdgeId e) const {
		vector<EdgeId> adjacent_edges;
		VertexId start = this->g().EdgeStart(e), end = this->g().EdgeEnd(e);
		Append(adjacent_edges, this->g().out_begin(start), this->g().out_end(start));
		Append(adjacent_edges, this->g().in_begin(start), this->g().in_end(start));
		Append(adjacent_edges, this->g().out_begin(end), this->g().out_end(end));
		Append(adjacent_edges, this->g().in_begin(end), this->g().in_end(end));
		return CheckAdjacent(adjacent_edges, e);
	}

private:
	DECL_LOGGER("CheatingChimericEdgeCondition")
	;
};

////todo isn't it the same as relative coverage?!!!
template<class Graph>
class CheatingChimericEdgeRemover: public ChimericEdgeRemovingAlgorithm<Graph,
		LengthComparator<Graph>> {
private:
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
	typedef ChimericEdgeRemovingAlgorithm<Graph, LengthComparator<Graph>> base;

public:
	CheatingChimericEdgeRemover(Graph& g, size_t max_length,
			double coverage_gap, size_t neighbour_length_threshold,
			boost::function<void(EdgeId)> removal_handler) :
			base(g,
					make_shared<CheatingChimericEdgeCondition<Graph>>(g,
							coverage_gap, neighbour_length_threshold),
					removal_handler, LengthComparator<Graph>(g),
					make_shared<LengthUpperBound<Graph>>(g, max_length)) {
	}
};


template<class Graph>
class TopologyChimericEdgeRemover: public ChimericEdgeRemovingAlgorithm<Graph,
		LengthComparator<Graph>> {
private:
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
	typedef ChimericEdgeRemovingAlgorithm<Graph, LengthComparator<Graph>> base;

public:
	TopologyChimericEdgeRemover(Graph& g, size_t max_length,
			size_t uniqueness_length, size_t plausibility_length,
			boost::function<void(EdgeId)> removal_handler) :
			base(g,
					make_shared<PredicateUniquenessPlausabilityCondition<Graph>>(
							g,
							/*uniqueness*/MakePathLengthLowerBound(g,
									TrivialPathFinder<Graph>(g),
									uniqueness_length),
							/*plausibility*/MakePathLengthLowerBound(g,
									TrivialPathFinder<Graph>(g),
									plausibility_length)), removal_handler,
					LengthComparator<Graph>(g),
					make_shared<LengthUpperBound<Graph>>(g, max_length)) {
	}
};

template<class Graph>
class TopologyAndReliablityBasedChimericEdgeRemover: public ChimericEdgeRemovingAlgorithm<
		Graph, LengthComparator<Graph>> {
private:
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
	typedef ChimericEdgeRemovingAlgorithm<Graph, LengthComparator<Graph>> base;

public:
	TopologyAndReliablityBasedChimericEdgeRemover(Graph& g, size_t max_length,
			size_t uniqueness_length, double max_coverage,
			boost::function<void(EdgeId)> removal_handler) :
			base(g,
					func::And<EdgeId>(
							make_shared<CoverageUpperBound<Graph>>(g,
									max_coverage),
							make_shared<
									PredicateUniquenessPlausabilityCondition<
											Graph>>(g,
									/*uniqueness*/
									MakePathLengthLowerBound(g,
											TrivialPathFinder<Graph>(g),
											uniqueness_length),
									/*plausibility*/make_shared<
											func::AlwaysTrue<EdgeId>>())),
					removal_handler, LengthComparator<Graph>(g),
					make_shared<LengthUpperBound<Graph>>(g, max_length)) {
	}
};

//todo refactor
template<class Graph>
class ThornCondition: public EdgeCondition<Graph> {
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
		return Unique(this->g().IncomingEdges(this->g().EdgeStart(e)),
				false)
				|| Unique(this->g().OutgoingEdges(this->g().EdgeEnd(e)),
						true);
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

		BoundedDijkstra<Graph> dij(this->g(), dijkstra_depth_);
		dij.run(this->g().EdgeStart(e));
		vector<VertexId> reached = dij.ReachedVertices();
		for (auto it = reached.begin(); it != reached.end(); ++it) {
			if (*it != this->g().EdgeEnd(e)
					&& this->g().RelatedVertices(*it,
							this->g().EdgeEnd(e))) {
				return true;
			}
		}
		return false;
	}

	bool CheckAlternativeCoverage(const vector<EdgeId>& edges, EdgeId e) const {
		for (auto it = edges.begin(); it != edges.end(); ++it) {
			if (*it != e && this->g().length(*it) < 400
					&& this->g().coverage(*it)
							< 15 * this->g().coverage(e)) {
				return false;
			}
		}
		return true;
	}

	bool CheckCoverageAround(EdgeId e) const {
		return CheckAlternativeCoverage(
				this->g().IncomingEdges(this->g().EdgeStart(e)), e)
				&& CheckAlternativeCoverage(
						this->g().OutgoingEdges(this->g().EdgeStart(e)),
						e)
				&& CheckAlternativeCoverage(
						this->g().IncomingEdges(this->g().EdgeEnd(e)),
						e)
				&& CheckAlternativeCoverage(
						this->g().OutgoingEdges(this->g().EdgeEnd(e)),
						e);
	}

	bool CheckUniqueness(EdgeId e, bool forward) const {
		return this->g().length(e) >= uniqueness_length_;
	}

public:

	ThornCondition(Graph& g, size_t uniqueness_length, size_t dijkstra_depth) :
			base(g), uniqueness_length_(uniqueness_length), dijkstra_depth_(
					dijkstra_depth) {
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
class ThornRemover: public ChimericEdgeRemovingAlgorithm<Graph,
		LengthComparator<Graph>> {
private:
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
	typedef ChimericEdgeRemovingAlgorithm<Graph, LengthComparator<Graph>> base;

public:
	ThornRemover(Graph& g, size_t max_length, size_t uniqueness_length,
			size_t dijkstra_depth,
			boost::function<void(EdgeId)> removal_handler) :
			base(g,
					make_shared<ThornCondition<Graph>>(g, uniqueness_length,
							dijkstra_depth), removal_handler,
					LengthComparator<Graph>(g),
					make_shared<LengthUpperBound<Graph>>(g, max_length)) {
	}
};

template<class Graph>
class AdvancedTopologyChimericEdgeRemover: public ChimericEdgeRemovingAlgorithm<
		Graph, LengthComparator<Graph>> {
private:
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
	typedef ChimericEdgeRemovingAlgorithm<Graph, LengthComparator<Graph>> base;

public:
	AdvancedTopologyChimericEdgeRemover(Graph& g, size_t max_length,
			size_t uniqueness_length, size_t plausibility_length,
			boost::function<void(EdgeId)> removal_handler) :
			base(g,
					make_shared<PredicateUniquenessPlausabilityCondition<Graph>>(
							g,
							/*uniqueness*/MakePathLengthLowerBound(g,
									UniquePathFinder<Graph>(g),
									uniqueness_length),
							/*plausibility*/MakePathLengthLowerBound(g,
									PlausiblePathFinder<Graph>(g,
											2 * plausibility_length),
									plausibility_length)), removal_handler,
					LengthComparator<Graph>(g),
					make_shared<LengthUpperBound<Graph>>(g, max_length)) {
	}
};

template<class Graph>
class MultiplicityCountingCondition: public UniquenessPlausabilityCondition<
		Graph> {
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
	typedef shared_ptr<Predicate<EdgeId>> EdgePredicate;
	typedef UniquenessPlausabilityCondition<Graph> base;

	MultiplicityCounter<Graph> multiplicity_counter_;
	EdgePredicate plausiblity_condition_;

	bool CheckUniqueness(EdgeId e, bool forward) const {
		TRACE(
				"Checking " << this->g().int_id(e) << " for uniqueness in " << (forward ? "forward" : "backward") << " direction");
		VertexId start =
				forward ? this->g().EdgeEnd(e) : this->g().EdgeStart(e);
		bool result = multiplicity_counter_.count(e, start) <= 1;
		TRACE(
				"Edge " << this->g().int_id(e) << " is" << (result ? "" : " not") << " unique");
		return result;
	}

	bool CheckPlausibility(EdgeId e, bool) const {
		return plausiblity_condition_->Check(e);
	}

public:

	MultiplicityCountingCondition(const Graph& g, size_t uniqueness_length,
			EdgePredicate plausiblity_condition) :
			//todo why 8???
			base(g), multiplicity_counter_(g, uniqueness_length, 8), plausiblity_condition_(
					plausiblity_condition) {

	}

private:

	DECL_LOGGER("MultiplicityCountingCondition")
	;
};

template<class Graph>
class SimpleMultiplicityCountingChimericEdgeRemover: public ChimericEdgeRemovingAlgorithm<
		Graph, LengthComparator<Graph>> {
private:
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
	typedef ChimericEdgeRemovingAlgorithm<Graph, LengthComparator<Graph>> base;

public:
	SimpleMultiplicityCountingChimericEdgeRemover(Graph& g, size_t max_length,
			size_t uniqueness_length, size_t plausibility_length,
			boost::function<void(EdgeId)> removal_handler) :
			base(g,
					make_shared<MultiplicityCountingCondition<Graph>>(g,
							uniqueness_length,
							/*plausibility*/MakePathLengthLowerBound(g,
									PlausiblePathFinder<Graph>(g,
											2 * plausibility_length),
									plausibility_length)), removal_handler,
					LengthComparator<Graph>(g),
					make_shared<LengthUpperBound<Graph>>(g, max_length)) {
	}
};

template<class Graph>
class PairInfoAwareErroneousCondition: public EdgeCondition<Graph> {
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
	typedef EdgeCondition<Graph> base;

	const PairedInfoIndexT<Graph>& paired_index_;
	size_t min_neighbour_length_;
	size_t insert_size_;
	size_t read_length_;
	size_t gap_;

	bool ShouldContainInfo(EdgeId e1, EdgeId e2, size_t gap_length) const {
		//todo discuss addition of negative delta
		//todo second condition may be included into the constructor warn/assert
		TRACE(
				"Checking whether should be pair info between e1 " << PrintEdge(e1) << " and e2 " << PrintEdge(e2) << " with gap " << gap_length);
		bool should_contain = gap_length
				>= PairInfoPathLengthLowerBound(this->g().k(),
						this->g().length(e1), this->g().length(e2),
						gap_, 0.)
				&& gap_length
						<= PairInfoPathLengthUpperBound(this->g().k(),
								insert_size_, 0.);
		TRACE("Result: " << should_contain);
		return should_contain;
	}

	bool ContainsInfo(EdgeId e1, EdgeId e2, size_t ec_length) const {
		TRACE(
				"Looking for pair info between e1 " << PrintEdge(e1) << " and e2 " << PrintEdge(e2));
		const set<Point>& infos = paired_index_.GetEdgePairInfo(e1, e2);
		for (auto it = infos.begin(); it != infos.end(); ++it) {
			const Point& point = *it;
			size_t distance = this->g().length(e1) + ec_length;
			if (math::ge(distance + point.var, point.d)
					&& math::le(double(distance), point.d + point.var)) {
				TRACE("Pair info found");
				return true;
			}
		}
		TRACE("Pair info not found");
		return false;
	}

	bool CheckAnyPairInfoAbsense(EdgeId possible_ec) const {
		TRACE("Checking pair info absense");
		VertexId start = this->g().EdgeStart(possible_ec);
		for (auto I1 = this->g().in_begin(start), E1 = this->g().in_end(start); I1 != E1;
				++I1) {
			VertexId end = this->g().EdgeEnd(possible_ec);
			for (auto I2 = this->g().out_begin(end), E2 = this->g().out_end(end); I2 != E2;
					++I2)
				if (!ShouldContainInfo(*I1, *I2, this->g().length(possible_ec))
						|| ContainsInfo(*I1, *I2, this->g().length(possible_ec))) {
					TRACE("Check absense: fail");
					return false;
				}
			TRACE("Check absense: ok");
		}
		return true;
	}

	bool CheckAdjacentLengths(const vector<EdgeId>& edges, EdgeId possible_ec) const {
		TRACE("Checking adjacent lengths");
		TRACE("min_neighbour_length = " << min_neighbour_length_);
		for (auto it = edges.begin(); it != edges.end(); ++it)
			if (min_neighbour_length_ > this->g().length(*it)) {
				TRACE(
						"Check fail: edge " << PrintEdge(*it) << " was too short");
				return false;
			}
		TRACE("Check ok");
		return true;
	}

	//todo remove
	string PrintEdge(EdgeId e) const {
		stringstream ss;
		ss << this->g().int_ids().ReturnIntId(e) << "(" << e << ") "
				<< this->g().length(e) << "(" << this->g().coverage(e)
				<< ")";
		return ss.str();
	}

public:

	PairInfoAwareErroneousCondition(Graph& g,
			const PairedInfoIndexT<Graph>& paired_index,
			size_t min_neighbour_length, size_t insert_size, size_t read_length) :
			base(g), paired_index_(paired_index), min_neighbour_length_(
					min_neighbour_length), insert_size_(insert_size), read_length_(
					read_length), gap_(insert_size_ - 2 * read_length_) {
		VERIFY(insert_size_ >= 2 * read_length_);
	}

	bool Check(EdgeId e) const {
		vector<EdgeId> adjacent_edges;
		VertexId start = this->g().EdgeStart(e), end = this->g().EdgeEnd(e);
		Append(adjacent_edges, this->g().in_begin(start), this->g().in_end(start));
		Append(adjacent_edges, this->g().out_begin(end), this->g().out_end(end));
		return CheckAdjacentLengths(adjacent_edges, e)
				&& CheckAnyPairInfoAbsense(e);
	}

private:

	DECL_LOGGER("PairInfoAwareErroneousCondition")
	;
};

template<class Graph>
class PairInfoAwareErroneousEdgeRemover: public ChimericEdgeRemovingAlgorithm<
		Graph, LengthComparator<Graph>> {
private:
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
	typedef ChimericEdgeRemovingAlgorithm<Graph, LengthComparator<Graph>> base;

public:
	PairInfoAwareErroneousEdgeRemover(Graph& g,
			const PairedInfoIndexT<Graph>& paired_index, size_t max_length,
			size_t min_neighbour_length, size_t insert_size, size_t read_length,
			boost::function<void(EdgeId)> removal_handler) :
			base(g,
					make_shared<PairInfoAwareErroneousCondition<Graph>>(g,
							paired_index, min_neighbour_length, insert_size,
							read_length), removal_handler,
					LengthComparator<Graph>(g),
					make_shared<LengthUpperBound<Graph>>(g, max_length)) {
	}
};

}
