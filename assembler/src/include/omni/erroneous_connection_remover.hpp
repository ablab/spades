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

template<class Graph>
class UniquenessPlausabilityCondition: public EdgeCondition<Graph> {
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
	typedef EdgeCondition<Graph> base;

	size_t uniqueness_length_;
	size_t plausibility_length_;
//	shared_ptr<Predicate<EdgeId>> uniqueness_condition_;
//	shared_ptr<Predicate<EdgeId>> plausiblity_condition_;

	bool CheckUniqueness(EdgeId e,
			const AbstractDirection<Graph>& direction) const {
		size_t path_length = CummulativeLength(this->g(),
				unique_path_finder_(e, direction));
		return path_length >= uniqueness_length_;
	}

	bool CheckPlausibility(EdgeId e,
			const AbstractDirection<Graph>& direction) const {
		size_t path_length = CummulativeLength(this->graph(),
				plausible_path_finder_(e, direction));
		return path_length >= plausibility_length_;
	}

	bool SingleUnique(const vector<EdgeId>& edges,
			const AbstractDirection<Graph>& direction) const {
		return edges.size() == 1 && CheckUniqueness(*edges.begin(), direction);
	}

	bool ExistPlausible(EdgeId init_e, const vector<EdgeId>& edges,
			const AbstractDirection<Graph>& direction) const {
		FOREACH(EdgeId e, edges) {
			if (e == init_e)
				continue;
			if (CheckPlausibility(e, direction)) {
				return true;
			}
		}
		return false;
	}

	bool Check(EdgeId e, const AbstractDirection<Graph>& forward,
			const AbstractDirection<Graph>& backward) const {
		return SingleUnique(forward.IncomingEdges(forward.EdgeStart(e)),
				backward)
				&& ExistPlausible(e,
						forward.OutgoingEdges(forward.EdgeStart(e)), forward);
	}

public:

	UniquenessPlausabilityCondition(const Graph& g,
			const UniquePF& unique_path_finder, size_t uniqueness_length,
			const PlausiblePF& plausible_path_finder,
			size_t plausibility_length) :
			base(g), unique_path_finder_(unique_path_finder), uniqueness_length_(
					uniqueness_length), plausible_path_finder_(
					plausible_path_finder), plausibility_length_(
					plausibility_length) {

	}

	bool Check(EdgeId e) const {
		ForwardDirection<Graph> forward(this->g());
		BackwardDirection<Graph> backward(this->g());
		return Check(e, forward, backward) || Check(e, backward, forward);
	}

};

template<class Graph, class UniquePF, class PlausiblePF>
shared_ptr<UniquenessPlausabilityCondition<Graph, UniquePF, PlausiblePF>> MakeUniquenessPlausabilityCondition(
		const Graph& g, const UniquePF& unique_path_finder,
		size_t uniqueness_length, const PlausiblePF& plausible_path_finder,
		size_t plausibility_length) {
	return make_shared<
			UniquenessPlausabilityCondition<Graph, UniquePF, PlausiblePF>>(g,
			unique_path_finder, uniqueness_length, plausible_path_finder,
			plausibility_length);
}

template<class Graph, class UniquePF, class PlausiblePF>
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

template<class Graph>
class ErroneousEdgeRemover {
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
	Graph& g_;
	AbstractEdgeRemover<Graph>& edge_remover_;
	bool graph_changed_;
public:
	ErroneousEdgeRemover(Graph& g, AbstractEdgeRemover<Graph>& edge_remover) :
			g_(g), edge_remover_(edge_remover), graph_changed_(false) {

	}

	virtual ~ErroneousEdgeRemover() {

	}

	bool RemoveEdges() {
		InnerRemoveEdges();
		omnigraph::Cleaner<Graph> cleaner(g_);
		cleaner.Clean();
		return graph_changed_;
	}

protected:
	void DeleteEdge(EdgeId edge) {
		TRACE("Transferring deletion task to edge remover");
		graph_changed_ = edge_remover_.DeleteEdge(edge) || graph_changed_;
	}

	virtual void InnerRemoveEdges() = 0;

	const Graph& graph() const {
		return g_;
	}

private:
	DECL_LOGGER("ErroneousEdgeRemover")
	;
};

////todo isn't this one cheating?!!!
//template<class Graph>
//class ChimericEdgesRemover: public ErroneousEdgeRemover<Graph> {
//private:
//	typedef typename Graph::EdgeId EdgeId;
//	typedef typename Graph::VertexId VertexId;
//	typedef ErroneousEdgeRemover<Graph> base;
//
//	size_t max_overlap_;
//
//	ChimericEdgesRemover(Graph &graph, size_t max_overlap,
//			AbstractEdgeRemover<Graph>& edge_remover) :
//			base(graph, edge_remover), max_overlap_(max_overlap) {
//	}
//
//	bool CheckEnd(VertexId v) {
//		return this->graph().OutgoingEdgeCount(v) == 1
//		/*&& graph_.IncomingEdgeCount(v) >= 2*/;
////		return graph_.OutgoingEdgeCount(v) == 1
////				&& graph_.IncomingEdgeCount(v) >= 2;
//	}
//
//	bool CheckStart(VertexId v) {
//		return /*graph_.OutgoingEdgeCount(v) >= 2
//		 &&*/this->graph().IncomingEdgeCount(v) == 1;
////		return graph_.OutgoingEdgeCount(v) >= 2
////				&& graph_.IncomingEdgeCount(v) == 1;
//	}
//
//public:
//
//	void InnerRemoveEdges() {
//		for (auto it = this->graph().SmartEdgeBegin(); !it.IsEnd(); ++it) {
//			EdgeId edge = *it;
//			if (this->graph().length(edge) <= this->graph().k()
//					&& this->graph().length(edge)
//							>= this->graph().k() - max_overlap_
//					&& CheckEnd(this->graph().EdgeEnd(edge))
//					&& CheckStart(this->graph().EdgeStart(edge))) {
//				this->DeleteEdge(edge);
////				graph_.DeleteEdge(edge);
//			}
//		}
//	}
//};

//default comparator used
template<class Graph, class UniquePF, class PlausiblePF>
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
		return this->g().length(e) <= this->g().k()
				&& this->g().length(e) >= this->g().k() - max_overlap_;
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

//template<class Graph>
//class IterativeLowCoverageEdgeRemover: public ErroneousEdgeRemover<Graph> {
//	typedef typename Graph::EdgeId EdgeId;
//	typedef typename Graph::VertexId VertexId;
//	typedef ErroneousEdgeRemover<Graph> base;
//	double max_coverage_;
//	shared_ptr<Predicate<EdgeId>> condition_;
//
//public:
//	IterativeLowCoverageEdgeRemover(Graph& g, double max_coverage,
//			shared_ptr<Predicate<EdgeId>> condition,
//			AbstractEdgeRemover<Graph>& edge_remover) :
//			base(g, edge_remover), max_coverage_(max_coverage), condition_(
//					condition) {
//	}
//
//	void InnerRemoveEdges() {
//		TRACE("Removing edges")
//		CoverageComparator<Graph> comparator(this->graph());
//		for (auto it = this->graph().SmartEdgeBegin(comparator); !it.IsEnd();
//				++it) {
//			typename Graph::EdgeId e = *it;
//			TRACE("Considering edge " << this->graph().str(e));
//			if (math::gr(this->graph().coverage(e), max_coverage_)) {
//				TRACE("Max coverage " << max_coverage_ << " achieved");
//				return;
//			}
//			TRACE("Checking length");
//			if (condition_->Check(e)) {
//				TRACE("Condition ok");
//				this->DeleteEdge(e);
//			} else {
//				TRACE("Condition failed");
//			}
//			TRACE("Edge " << this->graph().str(e) << " processed");
//		}
//	}
//private:
//	DECL_LOGGER("IterativeLowCoverageEdgeRemover")
//	;
//};

//coverage comparator
template<class Graph>
class RelativeCoverageCondition: public EdgeCondition<Graph> {
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
	typedef EdgeCondition<Graph> base;

	double max_relative_coverage_;

	bool CheckAlternativeCoverage(const vector<EdgeId> &edges, EdgeId e) {
		for (auto it = edges.begin(); it != edges.end(); ++it) {
			//todo wtf?!!! remove 400 condition
			if (*it != e && this->g().length(*it) < 400
					&& this->g().coverage(*it)
							< max_relative_coverage_ * this->g().coverage(e)) {
				return false;
			}
		}

		return true;
	}

public:

	RelativeCoverageCondition(Graph& g, double max_relative_coverage) :
			base(g), max_relative_coverage_(max_relative_coverage) {

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

//template<class Graph>
//class RelativeLowCoverageEdgeRemover: public ErroneousEdgeRemover<Graph> {
//	typedef typename Graph::EdgeId EdgeId;
//	typedef typename Graph::VertexId VertexId;
//	typedef ErroneousEdgeRemover<Graph> base;
//	size_t max_length_;
//	double max_coverage_;
//	double max_relative_coverage_;
//
//public:
//	RelativeLowCoverageEdgeRemover(Graph& g, size_t max_length,
//			double max_coverage, double max_relative_coverage,
//			AbstractEdgeRemover<Graph>& edge_remover) :
//			base(g, edge_remover), max_length_(max_length), max_coverage_(
//					max_coverage), max_relative_coverage_(max_relative_coverage) {
//
//	}
//
//	void InnerRemoveEdges() {
//		TRACE("Removing edges")
//		CoverageComparator<Graph> comparator(this->graph());
//		int total_len = 0;
//		for (auto it = this->graph().SmartEdgeBegin(comparator); !it.IsEnd();
//				++it) {
//			typename Graph::EdgeId e = *it;
//			TRACE("Considering edge " << this->graph().str(e));
//
//			if (math::gr(this->graph().coverage(e), max_coverage_)) {
//				TRACE("Max coverage " << max_coverage_ << " achieved");
//				return;
//			}
//			TRACE("Checking length");
//			if (this->graph().length(e) < max_length_) {
//				TRACE("Condition ok");
//				if (CheckCoverageAround(e)) {
//					TRACE("relative recoverage condition ok");
//					total_len += this->graph().length(e);
//					this->DeleteEdge(e);
//				}
//			} else {
//				TRACE("Condition failed");
//			}
//			TRACE("Edge " << this->graph().str(e) << " processed");
//		}
//
//	}
//private:
//	bool CheckAlternativeCoverage(const vector<EdgeId> &edges, EdgeId e) {
//		for (auto it = edges.begin(); it != edges.end(); ++it) {
//			//todo wtf?!!!
//			if (*it != e && this->graph().length(*it) < 400
//					&& this->graph().coverage(*it)
//							< max_relative_coverage_
//									* this->graph().coverage(e)) {
//				return false;
//			} else if (this->graph().coverage(*it)
//					> max_relative_coverage_ * this->graph().coverage(e)) {
//				TRACE("max relatiove coverage was" << max_relative_coverage_);
//				TRACE(
//						"and edges" << this->graph().coverage(*it) <<" "<< this->graph().coverage(e));
//			}
//
//		}
//
//		return true;
//	}
//
//	bool CheckCoverageAround(EdgeId e) {
//		return CheckAlternativeCoverage(
//				this->graph().IncomingEdges(this->graph().EdgeStart(e)), e)
//				&& CheckAlternativeCoverage(
//						this->graph().OutgoingEdges(this->graph().EdgeStart(e)),
//						e)
//				&& CheckAlternativeCoverage(
//						this->graph().IncomingEdges(this->graph().EdgeEnd(e)),
//						e)
//				&& CheckAlternativeCoverage(
//						this->graph().OutgoingEdges(this->graph().EdgeEnd(e)),
//						e);
//	}
//
//	DECL_LOGGER("RelativeLowCoverageEdgeRemover")
//	;
//};

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
	size_t neighbour_length_threshold;

	bool StrongNeighbourCondition(EdgeId neighbour_edge, EdgeId possible_ec) {
		return neighbour_edge == possible_ec
				|| math::gr(this->g().coverage(neighbour_edge),
						this->g().coverage(possible_ec) * coverage_gap_)
				|| this->g().length(neighbour_edge)
						>= neighbour_length_threshold_;
	}

	bool CheckAdjacent(const vector<EdgeId>& edges, EdgeId possible_ec) {
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
		VertexId start = g.EdgeStart(e), end = g.EdgeEnd(e);
		Append(adjacent_edges, g.out_begin(start), g.out_end(start));
		Append(adjacent_edges, g.in_begin(start), g.in_end(start));
		Append(adjacent_edges, g.out_begin(end), g.out_end(end));
		Append(adjacent_edges, g.in_begin(end), g.in_end(end));
		return CheckAdjacent(adjacent_edges, e);
	}

private:
	DECL_LOGGER("CheatingChimericEdgeCondition")
	;
};

////todo isn't it the same as relative coverage?!!!
////length comparator
//template<class Graph>
//class CheatingChimericEdgeRemover: public ErroneousEdgeRemover<Graph> {
//	typedef typename Graph::EdgeId EdgeId;
//	typedef typename Graph::VertexId VertexId;
//	typedef ErroneousEdgeRemover<Graph> base;
//
//	size_t max_length_;
//	double coverage_gap_;
//	size_t neighbour_length_threshold_;
//
//public:
//	CheatingChimericEdgeRemover(Graph& g, size_t max_length,
//			double coverage_gap, size_t neighbour_length_threshold,
//			EdgeRemover<Graph>& edge_remover) :
//			base(g, edge_remover), max_length_(max_length), coverage_gap_(
//					coverage_gap), neighbour_length_threshold_(
//					neighbour_length_threshold) {
//
//	}
//
//	bool StrongNeighbourCondition(EdgeId neighbour_edge, EdgeId possible_ec) {
//		return neighbour_edge == possible_ec
//				|| math::gr(this->graph().coverage(neighbour_edge),
//						this->graph().coverage(possible_ec) * coverage_gap_)
//				|| this->graph().length(neighbour_edge)
//						>= neighbour_length_threshold_;
//	}
//
//	bool CheckAdjacent(const vector<EdgeId>& edges, EdgeId possible_ec) {
//		for (auto it = edges.begin(); it != edges.end(); ++it) {
//			if (!StrongNeighbourCondition(*it, possible_ec))
//				return false;
//		}
//		return true;
//	}
//
//	void InnerRemoveEdges() {
//		const Graph &g = this->graph();
//		LengthComparator<Graph> comparator(g);
//		for (auto it = g.SmartEdgeBegin(comparator); !it.IsEnd(); ++it) {
//			typename Graph::EdgeId e = *it;
//			if (g.length(e) > max_length_) {
//				return;
//			}
//			vector<EdgeId> adjacent_edges;
//			VertexId start = g.EdgeStart(e), end = g.EdgeEnd(e);
//			Append(adjacent_edges, g.out_begin(start), g.out_end(start));
//			Append(adjacent_edges, g.in_begin(start), g.in_end(start));
//			Append(adjacent_edges, g.out_begin(end), g.out_end(end));
//			Append(adjacent_edges, g.in_begin(end), g.in_end(end));
//
//			if (CheckAdjacent(adjacent_edges, e)) {
//				this->DeleteEdge(e);
////				VertexId start = g_.EdgeStart(e);
////				VertexId end = g_.EdgeEnd(e);
////				if (!g_.RelatedVertices(start, end)) {
////					g_.DeleteEdge(e);
////					g_.CompressVertex(start);
////					g_.CompressVertex(end);
////				}
//			}
//		}
////		omnigraph::Cleaner<Graph> cleaner(this->graph());
////		cleaner.Clean();
//	}
//};

//template<class Graph>
//class TopologyChimericEdgeRemover: public ErroneousEdgeRemover<Graph> {
//	typedef typename Graph::EdgeId EdgeId;
//	typedef typename Graph::VertexId VertexId;
//	typedef ErroneousEdgeRemover<Graph> base;
//
//	size_t max_length_;
//	size_t uniqueness_length_;
//	size_t plausibility_length_;
//
//	bool Unique(const vector<EdgeId>& edges, bool forward) const {
//		return edges.size() == 1 && CheckUniqueness(*edges.begin(), forward);
//	}
//
//	bool ExistPlausible(EdgeId edge, const vector<EdgeId>& edges,
//			bool forward) const {
//		for (auto it = edges.begin(); it != edges.end(); ++it) {
//			if (edge != *it && CheckPlausibility(*it, forward)) {
//				return true;
//			}
//		}
//		return false;
//	}
//
//	bool CheckStart(EdgeId e) const {
//		TRACE("Checking conditions for edge start");
//		bool result = Unique(
//				this->graph().IncomingEdges(this->graph().EdgeStart(e)), false)
//				&& ExistPlausible(e,
//						this->graph().OutgoingEdges(this->graph().EdgeStart(e)),
//						true);
//		TRACE("Checks for edge start " << (result ? "passed" : "failed"));
//		return result;
//	}
//
//	bool CheckEnd(EdgeId e) const {
//		TRACE("Checking conditions for edge end");
//		bool result = Unique(
//				this->graph().OutgoingEdges(this->graph().EdgeEnd(e)), true)
//				&& ExistPlausible(e,
//						this->graph().IncomingEdges(this->graph().EdgeEnd(e)),
//						false);
//		TRACE("Checks for edge end" << (result ? "passed" : "failed"));
//		return result;
//	}
//
//public:
//	TopologyChimericEdgeRemover(Graph& g, size_t max_length,
//			size_t uniqueness_length, size_t plausibility_length,
//			AbstractEdgeRemover<Graph>& edge_remover) :
//			base(g, edge_remover), max_length_(max_length), uniqueness_length_(
//					uniqueness_length), plausibility_length_(
//					plausibility_length) {
//		VERIFY(max_length < plausibility_length);
//		VERIFY(uniqueness_length > plausibility_length);
//	}
//
//protected:
//
//	size_t uniqueness_length() const {
//		return uniqueness_length_;
//	}
//
//	size_t plausibility_length() const {
//		return plausibility_length_;
//	}
//
//	virtual bool CheckUniqueness(EdgeId e, bool forward) const {
//		return this->graph().length(e) >= uniqueness_length_;
//	}
//
//	virtual bool CheckPlausibility(EdgeId e, bool forward) const {
//		return this->graph().length(e) >= plausibility_length_;
//	}
//
//	void InnerRemoveEdges() {
//		TRACE("Removing erroneous connections");
//		LengthComparator<Graph> comparator(this->graph());
//		for (auto it = this->graph().SmartEdgeBegin(comparator); !it.IsEnd();
//				++it) {
//			typename Graph::EdgeId e = *it;
//			if (this->graph().length(e) > max_length_) {
//				return;
//			}
//			TRACE("Checking edge " << this->graph().int_id(e));
//			if (CheckStart(e) || CheckEnd(e)) {
//				TRACE("Deleting edge " << this->graph().int_id(e));
//				this->DeleteEdge(e);
//				TRACE("Edge was deleted");
//			} else {
//				TRACE("Edge " << this->graph().int_id(e) << " was not deleted");
//			}
//		}
//		TRACE("Removing erroneous connections finished");
//	}
//private:
//	DECL_LOGGER("TopologyChimericEdgeRemover")
//	;
//};

//template<class Graph>
//class TopologyAndReliablityBasedChimericEdgeRemover: public ErroneousEdgeRemover<
//		Graph> {
//	typedef typename Graph::EdgeId EdgeId;
//	typedef typename Graph::VertexId VertexId;
//	typedef ErroneousEdgeRemover<Graph> base;
//
//	size_t max_length_;
//	size_t uniqueness_length_;
//	double max_coverage_;
//
//	bool Unique(const vector<EdgeId>& edges, bool forward) const {
//		return edges.size() == 1 && CheckUniqueness(*edges.begin(), forward);
//	}
//
//	bool CheckUnique(EdgeId e) const {
//		return Unique(this->graph().IncomingEdges(this->graph().EdgeStart(e)),
//				false)
//				|| Unique(this->graph().OutgoingEdges(this->graph().EdgeEnd(e)),
//						true);
//	}
//
//	bool CheckExtremelyUnreliable(EdgeId e) {
//		return this->graph().coverage(e) < max_coverage_;
//	}
//
//public:
//	TopologyAndReliablityBasedChimericEdgeRemover(Graph& g, size_t max_length,
//			size_t uniqueness_length, double max_coverage,
//			AbstractEdgeRemover<Graph>& edge_remover) :
//			base(g, edge_remover), max_length_(max_length), uniqueness_length_(
//					uniqueness_length), max_coverage_(max_coverage) {
//		VERIFY(max_length < uniqueness_length);
//	}
//
//protected:
//
//	size_t uniqueness_length() const {
//		return uniqueness_length_;
//	}
//
//	virtual bool CheckUniqueness(EdgeId e, bool forward) const {
//		return this->graph().length(e) >= uniqueness_length_;
//	}
//
//	void InnerRemoveEdges() {
//		LengthComparator<Graph> comparator(this->graph());
//		for (auto it = this->graph().SmartEdgeBegin(comparator); !it.IsEnd();
//				++it) {
//			typename Graph::EdgeId e = *it;
//			if (this->graph().length(e) > max_length_) {
//				return;
//			}
//			TRACE("Checking edge " << this->graph().length(e));
//			if (CheckUnique(e) && CheckExtremelyUnreliable(e)) {
//				TRACE("Deleting edge " << this->graph().str(e));
//				this->DeleteEdge(e);
//				TRACE("Edge was deleted");
//			} else {
//				TRACE("Edge " << this->graph().length(e) << " was not deleted");
//			}
//		}
//	}
//private:
//	DECL_LOGGER("TopologyAndReliablityBasedChimericEdgeRemover")
//	;
//};

template<class Graph>
class ThornRemover: public ErroneousEdgeRemover<Graph> {
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
	typedef ErroneousEdgeRemover<Graph> base;

	size_t max_length_;
	size_t uniqueness_length_;
	size_t dijkstra_depth_;

	bool Unique(const vector<EdgeId>& edges, bool forward) const {
		return edges.size() == 1 && CheckUniqueness(*edges.begin(), forward);
	}

	bool CheckUnique(EdgeId e) const {
		TRACE("Checking conditions for edge start");
		return Unique(this->graph().IncomingEdges(this->graph().EdgeStart(e)),
				false)
				|| Unique(this->graph().OutgoingEdges(this->graph().EdgeEnd(e)),
						true);
	}

	bool CheckThorn(EdgeId e) {
		if (this->graph().EdgeStart(e) == this->graph().EdgeEnd(e))
			return false;
		if (this->graph().RelatedVertices(this->graph().EdgeStart(e),
				this->graph().EdgeEnd(e))) {
			return true;
		}
		if (this->graph().OutgoingEdgeCount(this->graph().EdgeStart(e)) != 2)
			return false;
		if (this->graph().IncomingEdgeCount(this->graph().EdgeStart(e)) != 1)
			return false;
		if (this->graph().OutgoingEdgeCount(this->graph().EdgeEnd(e)) != 1)
			return false;
		if (this->graph().IncomingEdgeCount(this->graph().EdgeEnd(e)) != 2)
			return false;

		BoundedDijkstra<Graph> dij(this->graph(), dijkstra_depth_);
		dij.run(this->graph().EdgeStart(e));
		vector<VertexId> reached = dij.ReachedVertices();
		for (auto it = reached.begin(); it != reached.end(); ++it) {
			if (*it != this->graph().EdgeEnd(e)
					&& this->graph().RelatedVertices(*it,
							this->graph().EdgeEnd(e))) {
				return true;
			}
		}
		return false;
	}

	bool CheckAlternativeCoverage(const vector<EdgeId>& edges, EdgeId e) {
		for (auto it = edges.begin(); it != edges.end(); ++it) {
			if (*it != e && this->graph().length(*it) < 400
					&& this->graph().coverage(*it)
							< 15 * this->graph().coverage(e)) {
				return false;
			}
		}
		return true;
	}

	bool CheckCoverageAround(EdgeId e) {
		return CheckAlternativeCoverage(
				this->graph().IncomingEdges(this->graph().EdgeStart(e)), e)
				&& CheckAlternativeCoverage(
						this->graph().OutgoingEdges(this->graph().EdgeStart(e)),
						e)
				&& CheckAlternativeCoverage(
						this->graph().IncomingEdges(this->graph().EdgeEnd(e)),
						e)
				&& CheckAlternativeCoverage(
						this->graph().OutgoingEdges(this->graph().EdgeEnd(e)),
						e);
	}

	bool Check(EdgeId e) {
		bool tmp = (CheckUnique(e) || CheckCoverageAround(e));
		if (tmp)
			tmp &= CheckThorn(e);
		return tmp;
	}

public:
	ThornRemover(Graph& g, size_t max_length, size_t uniqueness_length,
			size_t dijkstra_depth, AbstractEdgeRemover<Graph>& edge_remover) :
			base(g, edge_remover), max_length_(max_length), uniqueness_length_(
					uniqueness_length), dijkstra_depth_(dijkstra_depth) {
		VERIFY(max_length < uniqueness_length);
	}

protected:

	size_t uniqueness_length() const {
		return uniqueness_length_;
	}

	virtual bool CheckUniqueness(EdgeId e, bool forward) const {
		return this->graph().length(e) >= uniqueness_length_;
	}

	void InnerRemoveEdges() {
		LengthComparator<Graph> comparator(this->graph());
		for (auto it = this->graph().SmartEdgeBegin(comparator); !it.IsEnd();
				++it) {
			typename Graph::EdgeId e = *it;
			if (this->graph().length(e) > max_length_) {
				return;
			}
			TRACE("Checking edge " << this->graph().length(e));
			if (Check(e)) {
				TRACE("Deleting edge " << this->graph().length(e));
				this->DeleteEdge(e);
				TRACE("Edge was deleted");
			} else {
				TRACE("Edge " << this->graph().length(e) << " was not deleted");
			}
		}
	}
private:
	DECL_LOGGER("ThornRemover")
	;
};

//template<class Graph>
//class AdvancedTopologyChimericEdgeRemover: public TopologyChimericEdgeRemover<
//		Graph> {
//	typedef typename Graph::EdgeId EdgeId;
//	typedef typename Graph::VertexId VertexId;
//	typedef TopologyChimericEdgeRemover<Graph> base;
//
//	UniquePathFinder<Graph> unique_path_finder_;
//
//	PlausiblePathFinder<Graph> plausible_path_finder_;
//
//public:
//	AdvancedTopologyChimericEdgeRemover(Graph& g, size_t max_length,
//			size_t uniqueness_length, size_t plausibility_length,
//			AbstractEdgeRemover<Graph>& edge_remover) :
//			base(g, max_length, uniqueness_length, plausibility_length,
//					edge_remover), unique_path_finder_(g), plausible_path_finder_(
//					g, plausibility_length * 2) {
//	}
//
//protected:
//
//	bool CheckUniqueness(EdgeId e, bool forward) const {
//		TRACE(
//				"Checking " << this->graph().length(e) << " for uniqueness in " << (forward ? "forward" : "backward") << " direction");
//		bool result = CummulativeLength(this->graph(),
//				forward ?
//						unique_path_finder_.UniquePathForward(e) :
//						unique_path_finder_.UniquePathBackward(e))
//				>= this->uniqueness_length();
//		TRACE(
//				"Edge " << this->graph().length(e) << " is" << (result ? "" : " not") << " unique");
//		return result;
//	}
//
//	bool CheckPlausibility(EdgeId e, bool forward) const {
//		TRACE(
//				"Checking " << this->graph().length(e) << " for plausibility in " << (forward ? "forward" : "backward") << " direction");
//		bool result = CummulativeLength(this->graph(),
//				forward ?
//						plausible_path_finder_(e,
//								ForwardDirection<Graph>(this->graph())) :
//						plausible_path_finder_(e,
//								BackwardDirection<Graph>(this->graph())))
//				>= this->plausibility_length();
//		TRACE(
//				"Edge " << this->graph().length(e) << " is" << (result ? "" : " not") << " plausible");
//		return result;
//	}
//private:
//	DECL_LOGGER("AdvancedTopologyChimericEdgeRemover")
//	;
//};

template<class Graph>
class SimpleMultiplicityCountingChimericEdgeRemover: public TopologyChimericEdgeRemover<
		Graph> {
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
	typedef TopologyChimericEdgeRemover<Graph> base;

	MultiplicityCounter<Graph> multiplicity_counter_;

	PlausiblePathFinder<Graph> plausible_path_finder_;

public:
	SimpleMultiplicityCountingChimericEdgeRemover(Graph& g, size_t max_length,
			size_t uniqueness_length, size_t plausibility_length,
			AbstractEdgeRemover<Graph>& edge_remover) :
			base(g, max_length, uniqueness_length, plausibility_length,
					edge_remover), multiplicity_counter_(g, uniqueness_length,
					8), plausible_path_finder_(g, plausibility_length * 2) {

	}

protected:

	bool CheckUniqueness(EdgeId e, bool forward) const {
		TRACE(
				"Checking " << this->graph().int_id(e) << " for uniqueness in " << (forward ? "forward" : "backward") << " direction");
		VertexId start =
				forward ? this->graph().EdgeEnd(e) : this->graph().EdgeStart(e);
		bool result = multiplicity_counter_.count(e, start) <= 1;
		TRACE(
				"Edge " << this->graph().int_id(e) << " is" << (result ? "" : " not") << " unique");
		return result;
	}

	bool CheckPlausibility(EdgeId e, bool forward) const {
		TRACE(
				"Checking " << this->graph().int_id(e) << " for plausibility in " << (forward ? "forward" : "backward") << " direction");
		bool result = CummulativeLength(this->graph(),
				forward ?
						plausible_path_finder_(e,
								ForwardDirection<Graph>(this->graph())) :
						plausible_path_finder_(e,
								BackwardDirection<Graph>(this->graph())))
				>= this->plausibility_length();
		TRACE(
				"Edge " << this->graph().int_id(e) << " is" << (result ? "" : " not") << " plausible");
		return result;
	}
private:
	DECL_LOGGER("SimpleMultiplicityCountingChimericEdgeRemover")
	;
};

template<class Graph>
class PairInfoAwareErroneousEdgeRemover: public ErroneousEdgeRemover<Graph> {
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
	typedef ErroneousEdgeRemover<Graph> base;

	const PairedInfoIndexT<Graph>& paired_index_;
	size_t max_length_;
	size_t min_neighbour_length_;
	size_t insert_size_;
	size_t read_length_;
	size_t gap_;

public:
	PairInfoAwareErroneousEdgeRemover(Graph& g,
			const PairedInfoIndexT<Graph>& paired_index, size_t max_length,
			size_t min_neighbour_length, size_t insert_size, size_t read_length,
			AbstractEdgeRemover<Graph>& edge_remover) :
			base(g, edge_remover), paired_index_(paired_index), max_length_(
					max_length), min_neighbour_length_(min_neighbour_length), insert_size_(
					insert_size), read_length_(read_length), gap_(
					insert_size_ - 2 * read_length_) {
		VERIFY(insert_size_ >= 2 * read_length_);
	}

	bool ShouldContainInfo(EdgeId e1, EdgeId e2, size_t gap_length) {
		//todo discuss addition of negative delta
		//todo second condition may be included into the constructor warn/assert
		TRACE(
				"Checking whether should be pair info between e1 " << PrintEdge(e1) << " and e2 " << PrintEdge(e2) << " with gap " << gap_length);
		bool should_contain = gap_length
				>= PairInfoPathLengthLowerBound(this->graph().k(),
						this->graph().length(e1), this->graph().length(e2),
						gap_, 0.)
				&& gap_length
						<= PairInfoPathLengthUpperBound(this->graph().k(),
								insert_size_, 0.);
		TRACE("Result: " << should_contain);
		return should_contain;
	}

	bool ContainsInfo(EdgeId e1, EdgeId e2, size_t ec_length) {
		TRACE(
				"Looking for pair info between e1 " << PrintEdge(e1) << " and e2 " << PrintEdge(e2));
		const set<Point>& infos = paired_index_.GetEdgePairInfo(e1, e2);
		for (auto it = infos.begin(); it != infos.end(); ++it) {
			const Point& point = *it;
			size_t distance = this->graph().length(e1) + ec_length;
			if (math::ge(distance + point.var, point.d)
					&& math::le(double(distance), point.d + point.var)) {
				TRACE("Pair info found");
				return true;
			}
		}
		TRACE("Pair info not found");
		return false;
	}

	bool CheckAnyPairInfoAbsense(EdgeId possible_ec) {
		const Graph &g = this->graph();
		TRACE("Checking pair info absense");
		VertexId start = g.EdgeStart(possible_ec);
		for (auto I1 = g.in_begin(start), E1 = g.in_end(start); I1 != E1;
				++I1) {
			VertexId end = g.EdgeEnd(possible_ec);
			for (auto I2 = g.out_begin(end), E2 = g.out_end(end); I2 != E2;
					++I2)
				if (!ShouldContainInfo(*I1, *I2, g.length(possible_ec))
						|| ContainsInfo(*I1, *I2, g.length(possible_ec))) {
					TRACE("Check absense: fail");
					return false;
				}
			TRACE("Check absense: ok");
		}
		return true;
	}

	bool CheckAdjacentLengths(const vector<EdgeId>& edges, EdgeId possible_ec) {
		TRACE("Checking adjacent lengths");
		TRACE("min_neighbour_length = " << min_neighbour_length_);
		for (auto it = edges.begin(); it != edges.end(); ++it)
			if (min_neighbour_length_ > this->graph().length(*it)) {
				TRACE(
						"Check fail: edge " << PrintEdge(*it) << " was too short");
				return false;
			}
		TRACE("Check ok");
		return true;
	}

	string PrintEdge(EdgeId e) {
		stringstream ss;
		ss << this->graph().int_ids().ReturnIntId(e) << "(" << e << ") "
				<< this->graph().length(e) << "(" << this->graph().coverage(e)
				<< ")";
		return ss.str();
	}

	void InnerRemoveEdges() {
		const Graph &g = this->graph();

		TRACE("Removing erroneous edges based on pair info");
		LengthComparator<Graph> comparator(g);
		for (auto it = g.SmartEdgeBegin(comparator); !it.IsEnd(); ++it) {
			typename Graph::EdgeId e = *it;
			TRACE("Considering edge " << PrintEdge(e));
			if (g.length(e) > max_length_) {
				TRACE("Max length bound = " << max_length_ << " was exceeded");
				return;
			}
			vector<EdgeId> adjacent_edges;
			VertexId start = g.EdgeStart(e), end = g.EdgeEnd(e);
			Append(adjacent_edges, g.in_begin(start), g.in_end(start));
			Append(adjacent_edges, g.out_begin(end), g.out_end(end));

			if (CheckAdjacentLengths(adjacent_edges, e)
					&& CheckAnyPairInfoAbsense(e)) {
//				VertexId start = g_.EdgeStart(e);
//				VertexId end = g_.EdgeEnd(e);
//				TRACE("Try deleting edge " << PrintEdge(e));
//				if (!g_.RelatedVertices(start, end)) {
//					TRACE("Vertices not related");
//					TRACE("Deleting edge " << PrintEdge(e));
//					g_.DeleteEdge(e);
//					TRACE("Compressing start");
//					g_.CompressVertex(start);
//					TRACE("Compressing end");
//					g_.CompressVertex(end);
//				} else {
//					TRACE("Vertices are related");
//				}
				this->DeleteEdge(e);
			}
		}
	}
private:
	DECL_LOGGER("PairInfoAwareErroneousEdgeRemover")
	;
};

}
