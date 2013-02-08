#pragma once

#include "graph_processing_algorithm.hpp"
#include "basic_edge_conditions.hpp"
#include "omni_tools.hpp"
#include "omni_utils.hpp"
#include "func.hpp"
#include "xmath.h"

namespace omnigraph {

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
					func::And<EdgeId>(
							make_shared<LengthUpperBound<Graph>>(g, g.k()),
							make_shared<ChimericEdgeCondition<Graph>>(g,
									max_overlap)), removal_handler) {
	}
};

template<class Graph>
class SimpleTopologyChimericEdgeRemover: public ChimericEdgeRemovingAlgorithm<
		Graph, LengthComparator<Graph>> {
private:
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
	typedef ChimericEdgeRemovingAlgorithm<Graph, LengthComparator<Graph>> base;

public:
	SimpleTopologyChimericEdgeRemover(Graph& g, size_t max_length,
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
						this->g().length(e1), this->g().length(e2), gap_, 0.)
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
		for (auto I1 = this->g().in_begin(start), E1 = this->g().in_end(start);
				I1 != E1; ++I1) {
			VertexId end = this->g().EdgeEnd(possible_ec);
			for (auto I2 = this->g().out_begin(end), E2 = this->g().out_end(
					end); I2 != E2; ++I2)
				if (!ShouldContainInfo(*I1, *I2, this->g().length(possible_ec))
						|| ContainsInfo(*I1, *I2,
								this->g().length(possible_ec))) {
					TRACE("Check absense: fail");
					return false;
				}
			TRACE("Check absense: ok");
		}
		return true;
	}

	bool CheckAdjacentLengths(const vector<EdgeId>& edges,
			EdgeId possible_ec) const {
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
				<< this->g().length(e) << "(" << this->g().coverage(e) << ")";
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
		Append(adjacent_edges, this->g().in_begin(start),
				this->g().in_end(start));
		Append(adjacent_edges, this->g().out_begin(end),
				this->g().out_end(end));
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

template<class Graph>
bool CheatingRemoveErroneousEdges(Graph &g,
		const debruijn_config::simplification::cheating_erroneous_connections_remover& cec_config,
		boost::function<void(typename Graph::EdgeId)> removal_handler) {
	INFO("Cheating removal of erroneous edges started");
	size_t max_length = LengthThresholdFinder::MaxErroneousConnectionLength(
			g.k(), cec_config.max_ec_length_coefficient);
	double coverage_gap = cec_config.coverage_gap;
	size_t sufficient_neighbour_length = cec_config.sufficient_neighbour_length;
	return omnigraph::CheatingChimericEdgeRemover < Graph
			> (g, max_length, coverage_gap, sufficient_neighbour_length, removal_handler).Process();
}

template<class Graph>
bool ChimericRemoveErroneousEdges(Graph &g,
		boost::function<void(typename Graph::EdgeId)> removal_handler) {
	INFO("Simple removal of chimeric edges based only on length started");
	ChimericEdgesRemover<Graph> remover(g, 10, removal_handler);
	bool changed = remover.Process();
	DEBUG("Removal of chimeric edges finished");
	return changed;
}

template<class Graph>
void RemoveEroneousEdgesUsingPairedInfo(Graph& g,
		const PairedInfoIndexT<Graph>& paired_index,
		boost::function<void(typename Graph::EdgeId)> removal_handler) {
	INFO("Removing erroneous edges using paired info");
	size_t max_length = LengthThresholdFinder::MaxErroneousConnectionLength(
			g.k(), cfg::get().simp.piec.max_ec_length_coefficient);
	size_t min_neighbour_length = cfg::get().simp.piec.min_neighbour_length;
	omnigraph::PairInfoAwareErroneousEdgeRemover<Graph> erroneous_edge_remover(
			g, paired_index, max_length, min_neighbour_length,
			*cfg::get().ds.IS, *cfg::get().ds.RL, removal_handler);
	erroneous_edge_remover.Process();

	DEBUG("Erroneous edges using paired info removed");
}

template<class graph_pack>
void RemoveErroneousEdgesWithPI(graph_pack& gp,
		const PairedInfoIndexT<typename graph_pack::graph_t>& paired_index) {
	typedef typename graph_pack::graph_t Graph;
	typedef typename Graph::EdgeId EdgeId;
	EdgeQuality<Graph> quality_handler(gp.g, gp.index, gp.kmer_mapper,
			gp.genome);
	QualityLoggingRemovalHandler<Graph> qual_removal_handler(gp.g,
			quality_handler);
	boost::function<void(EdgeId)> removal_handler_f = boost::bind(
			&QualityLoggingRemovalHandler < Graph > ::HandleDelete,
			&qual_removal_handler, _1);
	INFO("Pair info aware ErroneousConnectionsRemoval");
	RemoveEroneousEdgesUsingPairedInfo(gp.g, paired_index, removal_handler_f);
	INFO("Pair info aware ErroneousConnectionsRemoval stats");
	CountStats(gp.g, gp.index, gp.genome, gp.k_value);
}

}
