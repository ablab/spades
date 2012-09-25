#pragma once

#include "standard.hpp"
#include "omni/id_track_handler.hpp"

namespace debruijn_graph {

template<class Graph>
class IdNormalizer {
	typedef IdTrackHandler<Graph> IdHandler;
	typedef typename IdHandler::realIdType id_t;
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::SmartVertexIt SmartVertexIt;
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::SmartEdgeIt SmartEdgeIt;

	const Graph& g_;
	IdTrackHandler<Graph>& int_ids_;
public:
	IdNormalizer(const Graph& g, IdTrackHandler<Graph>& int_ids)
		: g_(g), int_ids_(int_ids) {

	}

	template<class SmartIt, class IdElemF>
	void Normalize(SmartIt& it, IdElemF id_elem_f) {
		vector<pair<Sequence, id_t>> sequence_id_pairs;
		while (!it.IsEnd()) {
			sequence_id_pairs.push_back(make_pair(g_.nucls(*it)
					, int_ids_.ReturnIntId(*it)));
			++it;
		}
		sort(sequence_id_pairs.begin(), sequence_id_pairs.end());
		for (size_t i = 0, n = sequence_id_pairs.size(); i < n; ++i) {
			int_ids_.add(id_elem_f(sequence_id_pairs[i].second), (id_t) i);
		}
	}

	void operator() () {
		auto smart_vertex_it = g_.SmartVertexBegin();
		Normalize(smart_vertex_it
				, bind(&IdHandler::ReturnVertexId, &int_ids_, _1));
		auto smart_edge_it = g_.SmartEdgeBegin();
		Normalize(smart_edge_it
				, bind(&IdHandler::ReturnEdgeId, &int_ids_, _1));
	}
};

}
