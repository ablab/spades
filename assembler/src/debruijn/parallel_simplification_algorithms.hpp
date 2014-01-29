#pragma once

#include "standard_base.hpp"
#include "omni/graph_processing_algorithm.hpp"
#include "omni/basic_edge_conditions.hpp"
#include "omni/bulge_remover.hpp"

namespace debruijn {

namespace simplification {

template<class Graph>
class ParallelTipClippingFunctor {
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
	typedef boost::function<void(EdgeId)> HandlerF;

	Graph& g_;
	size_t length_bound_;
	HandlerF handler_f_;

	bool IsIncomingTip(EdgeId e) const {
		return g_.length() <= length_bound_ &&
				g_.IncomingEdgeCount(g_.EdgeStart(e)) == 0;
	}

public:

	ParallelTipClippingFunctor(Graph& g, size_t length_bound, HandlerF handler_f) :
	g_(g), length_bound_(length_bound), handler_f_(handler_f) {

	}

	bool operator()(VertexId v) const {
		vector<EdgeId> tips;
		for (EdgeId e : g_.IncomingEdges(v)) {
			if (IsIncomingTip(e)) {
				tips.push_back(e);
			}
		}

		//if all of edges are tips, leave the longest one
		if (tips.size() == g_.IncomingEdgeCount(v)) {
			sort(tips.begin(), tips.end(), omnigraph::LengthComparator<Graph>(g_));
			tips.pop_back();
		}

		for (EdgeId e : tips) {
		    if (handler_f_) {
		        handler_f_(e);
		    }
			g_.RemoveEdge(e);
		}
		return false;
	}
};

template<class Graph>
class ParallelSimpleBRFunctor {
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;

	Graph& g_;
    size_t max_length_;
    double max_coverage_;
    double max_relative_coverage_;
    size_t max_delta_;
    double max_relative_delta_;

public:

	ParallelSimpleBRFunctor(Graph& g,
	                        size_t max_length,
	                        double max_coverage,
	                        double max_relative_coverage,
	                        size_t max_delta,
	                        double max_relative_delta,
	                        boost::function<void(EdgeId)> handler_f)
        : g_(g),
          max_length_(max_length),
          max_coverage_(max_coverage),
          max_relative_coverage_(max_relative_coverage),
          max_delta_(max_delta),
          max_relative_delta_(max_relative_delta)
    {

	}

	bool operator()(VertexId v) const {

	}
};

}

}
