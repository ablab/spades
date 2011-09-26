/*
 * erroneous_connection_remover.hpp
 *
 *  Created on: May 31, 2011
 *      Author: sergey
 */

#ifndef ERRONEOUS_CONNECTION_REMOVER_HPP_
#define ERRONEOUS_CONNECTION_REMOVER_HPP_

namespace omnigraph {

#include "omni_tools.hpp"
#include "omni_utils.hpp"
#include "xmath.h"

template<class Graph>
class LowCoverageEdgeRemover {
	Graph& g_;
	size_t max_length_;
	double max_coverage_;
public:
	LowCoverageEdgeRemover(Graph& g, size_t max_length, double max_coverage) :
			g_(g), max_length_(max_length), max_coverage_(max_coverage) {

	}

	void RemoveEdges() {
		for (auto it = g_.SmartEdgeBegin(); !it.IsEnd(); ++it) {
			typename Graph::EdgeId e = *it;
			if (g_.length(e) < max_length_ && g_.coverage(e) < max_coverage_) {
				g_.DeleteEdge(e);
			}
		}
		omnigraph::Compressor<Graph> compressor(g_);
		compressor.CompressAllVertices();
		omnigraph::Cleaner<Graph> cleaner(g_);
		cleaner.Clean();
	}
};

template<class Graph>
class IterativeLowCoverageEdgeRemover {
	Graph& g_;
	size_t max_length_;
	double max_coverage_;
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;

public:
	IterativeLowCoverageEdgeRemover(Graph& g, size_t max_length,
			double max_coverage) :
			g_(g), max_length_(max_length), max_coverage_(max_coverage) {

	}

	void RemoveEdges() {
		TRACE("Removing edges")
		CoverageComparator<Graph> comparator(g_);
		for (auto it = g_.SmartEdgeBegin(comparator); !it.IsEnd(); ++it) {
			typename Graph::EdgeId e = *it;
			TRACE("Considering edge " << e);
			if (math::gr(g_.coverage(e), max_coverage_)) {
				TRACE("Max coverage " << max_coverage_ << " achieved");
				return;
			}
			TRACE("Checking length");
			if (g_.length(e) < max_length_) {
				TRACE("Condition ok");
				VertexId start = g_.EdgeStart(e);
				VertexId end = g_.EdgeEnd(e);
				TRACE("Start " << start);
				TRACE("End " << end);
				TRACE("Deleting edge");
				g_.DeleteEdge(e);
				TRACE("Compressing locality");
				if (!g_.RelatedVertices(start, end)) {
					TRACE("Vertices not related");
					TRACE("Compressing end");
					g_.CompressVertex(end);
				}
				TRACE("Compressing start");
				g_.CompressVertex(start);
			} else {
				TRACE("Condition failed");
			}
			TRACE("Edge " << e << " processed");
		}
		TRACE("Cleaning graph");
		omnigraph::Cleaner<Graph> cleaner(g_);
		cleaner.Clean();
		TRACE("Graph cleaned");
	}
private:
	DECL_LOGGER("IterativeLowCoverageEdgeRemover");
};

template <class T>
void Append(vector<T>& current, const vector<T>& to_append) {
	current.insert(current.end(), to_append.begin(), to_append.end());
}

template<class Graph>
class RelativelyLowCoverageEdgeRemover {
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;

	Graph& g_;
	size_t max_length_;
	double coverage_gap_;
	size_t neighbour_length_threshold_;

public:
	RelativelyLowCoverageEdgeRemover(Graph& g, size_t max_length,
			double coverage_gap, size_t neighbour_length_threshold) :
			g_(g), max_length_(max_length), coverage_gap_(coverage_gap), neighbour_length_threshold_(
					neighbour_length_threshold) {

	}

	bool StrongNeighbourCondition(EdgeId neighbour_edge,
			EdgeId possible_ec) {
		return neighbour_edge == possible_ec || math::gr(g_.coverage(neighbour_edge),
				g_.coverage(possible_ec) * coverage_gap_)
				|| g_.length(neighbour_edge) >= neighbour_length_threshold_;
	}

	bool CheckAdjacent(const vector<EdgeId>& edges, EdgeId possible_ec) {
		for (auto it = edges.begin(); it != edges.end(); ++it) {
			if (!StrongNeighbourCondition(*it, possible_ec))
				return false;
		}
		return true;
	}

	void RemoveEdges() {
		LengthComparator<Graph> comparator(g_);
		for (auto it = g_.SmartEdgeBegin(comparator); !it.IsEnd(); ++it) {
			typename Graph::EdgeId e = *it;
			if (g_.length(e) > max_length_) {
				return;
			}
			vector<EdgeId> adjacent_edges;
			Append(adjacent_edges, g_.OutgoingEdges(g_.EdgeStart(e)));
			Append(adjacent_edges, g_.IncomingEdges(g_.EdgeStart(e)));
			Append(adjacent_edges, g_.OutgoingEdges(g_.EdgeEnd(e)));
			Append(adjacent_edges, g_.IncomingEdges(g_.EdgeEnd(e)));

			if (CheckAdjacent(adjacent_edges, e)) {
				VertexId start = g_.EdgeStart(e);
				VertexId end = g_.EdgeEnd(e);
				if (!g_.RelatedVertices(start, end)) {
					g_.DeleteEdge(e);
					g_.CompressVertex(start);
					g_.CompressVertex(end);
				}
			}
		}
		omnigraph::Cleaner<Graph> cleaner(g_);
		cleaner.Clean();
	}
};

template<class Graph>
class PairInfoAwareErroneousEdgeRemover {
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;

	Graph& g_;
	PairedInfoIndex<Graph>& paired_index_;
	size_t max_length_;
	size_t min_neighbour_length_;
	size_t insert_size_;
	size_t read_length_;
	size_t gap_;

public:
	PairInfoAwareErroneousEdgeRemover(Graph& g,
			size_t max_length,
			size_t min_neighbour_length,
			size_t insert_size,
			size_t read_length) :
			g_(g),
			max_length_(max_length),
			min_neighbour_length_(min_neighbour_length),
			insert_size_(insert_size),
			read_length_(read_length),
			gap_(insert_size_ - 2 * read_length_){
		assert(insert_size_ >= 2 * read_length_);
	}

	bool ShouldContainInfo(EdgeId e1, EdgeId e2, size_t gap_length) {
		//todo discuss addition of negative delta
		//todo second condition may be included into the constructor warn/assert
		return gap_length >= PairInfoPathLengthLowerBound(g_.k, g_.length(e1), g_.length(e2), gap_, 0.)
				&& gap_length <= PairInfoPathLengthUpperBound(g_.k, insert_size_, 0.);
	}

	bool ContainsInfo(EdgeId e1, EdgeId e2, size_t ec_length) {
		vector<PairInfo<EdgeId>> infos = paired_index_.GetEdgePairInfo(e1, e2);
		for (auto it = infos.begin(); it != infos.end(); ++it) {
			PairInfo<EdgeId> info = *it;
			size_t distance = g_.length(e1) + ec_length;
			if (math::ge(distance + info.variance, info.d) && math::le(distance, info.d + info.variance)) {
				return true;
			}
		}
		return false;
	}

	bool CheckAnyPairInfoAbsense(EdgeId possible_ec) {
		vector<EdgeId> incoming = g_.IncomingEdges(g_.EdgeStart(possible_ec));
		vector<EdgeId> outgoing = g_.OutgoingEdges(g_.EdgeEnd(possible_ec));
		for (auto it1 = incoming.begin(); it1 != incoming.end(); ++it1)
			for (auto it2 = outgoing.begin(); it2 != outgoing.end(); ++it2)
				if (!ShouldContainInfo(*it1, *it2, g_.length(possible_ec)) || ContainsInfo(*it1, *it2, g_.length(possible_ec)))
					return false;
		return true;
	}

	bool CheckAdjacentLengths(const vector<EdgeId>& edges, EdgeId possible_ec) {
		for (auto it = edges.begin(); it != edges.end(); ++it)
			if (min_neighbour_length_ > g_.length(*it))
				return false;
		return true;
	}

	void RemoveEdges() {
		LengthComparator<Graph> comparator(g_);
		for (auto it = g_.SmartEdgeBegin(comparator); !it.IsEnd(); ++it) {
			typename Graph::EdgeId e = *it;
			if (g_.length(e) > max_length_) {
				return;
			}
			vector<EdgeId> adjacent_edges;
			Append(adjacent_edges, g_.IncomingEdges(g_.EdgeStart(e)));
			Append(adjacent_edges, g_.OutgoingEdges(g_.EdgeEnd(e)));

			if (CheckAdjacentLengths(adjacent_edges, e) && CheckAnyPairInfoAbsense(e)) {
				VertexId start = g_.EdgeStart(e);
				VertexId end = g_.EdgeEnd(e);
				if (!g_.RelatedVertices(start, end)) {
					g_.DeleteEdge(e);
					g_.CompressVertex(start);
					g_.CompressVertex(end);
				}
			}
		}
		omnigraph::Cleaner<Graph> cleaner(g_);
		cleaner.Clean();
	}
};

}

#endif /* ERRONEOUS_CONNECTION_REMOVER_HPP_ */
