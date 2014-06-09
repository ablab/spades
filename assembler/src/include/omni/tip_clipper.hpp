//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/*
 * tip_clipper.hpp
 *
 *  Created on: Mar 25, 2011
 *      Author: sergey
 */

#pragma once

#include <set>

#include "omni_utils.hpp"
#include "xmath.h"
#include "basic_edge_conditions.hpp"
#include "graph_processing_algorithm.hpp"

namespace omnigraph {

template<class Graph>
class RelativeCoverageTipCondition: public EdgeCondition<Graph> {
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
	typedef EdgeCondition<Graph> base;

	const double max_relative_coverage_;

	template<class IteratorType>
	double MaxCompetitorCoverage(EdgeId tip, IteratorType begin, IteratorType end) const {
		double result = 0;
		for (auto it = begin; it != end; ++it) {
			if (*it != tip)
				result = std::max(result, this->g().coverage(*it));
		}
		return result;
	}

	double MaxCompetitorCoverage(EdgeId tip) const {
		const Graph &g = this->g();
		VertexId start = g.EdgeStart(tip), end = g.EdgeEnd(tip);
		auto out = g.OutgoingEdges(start);
		auto in = g.IncomingEdges(end);
		return std::max(
						MaxCompetitorCoverage(tip, out.begin(),	out.end()),
						MaxCompetitorCoverage(tip, in.begin(), in.end()));
//		return std::max(
//				MaxCompetitorCoverage(tip, g.out_begin(start),
//						g.out_end(start)),
//				MaxCompetitorCoverage(tip, g.in_begin(end), g.in_end(end)));
	}

public:

	RelativeCoverageTipCondition(const Graph& g, double max_relative_coverage) :
			base(g), max_relative_coverage_(max_relative_coverage) {
	}

	bool Check(EdgeId e) const {
		//+1 is a trick to deal with edges of 0 coverage from iterative run
		double max_coverage = MaxCompetitorCoverage(e) + 1;
		return math::le(this->g().coverage(e),
				max_relative_coverage_ * max_coverage);
	}
};

template<class Graph>
class TipCondition : public EdgeCondition<Graph> {
    typedef EdgeCondition<Graph> base;

    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;

    /**
     * This method checks if given vertex topologically looks like end of tip
     * @param v vertex to be checked
     * @return true if vertex judged to be tip and false otherwise.
     */
    bool IsTip(VertexId v) const {
        return this->g().IncomingEdgeCount(v) + this->g().OutgoingEdgeCount(v) == 1;
    }

public:
    TipCondition(const Graph& g) : base(g) {
    }

    /**
     * This method checks if given edge topologically looks like a tip.
     * @param edge edge vertex to be checked
     * @return true if edge judged to be tip and false otherwise.
     */
    /*virtual*/ bool Check(EdgeId e) const {
        return (IsTip(this->g().EdgeEnd(e)) || IsTip(this->g().EdgeStart(e)))
                && (this->g().OutgoingEdgeCount(this->g().EdgeStart(e))
                        + this->g().IncomingEdgeCount(this->g().EdgeEnd(e)) > 2);
    }

};

template<class Graph>
class TipClipper: public EdgeRemovingAlgorithm<Graph, LengthComparator<Graph>> {
    typedef EdgeRemovingAlgorithm<Graph, LengthComparator<Graph>> base;

    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;

public:

    TipClipper(Graph& g, size_t max_tip_length,
            const shared_ptr<Predicate<EdgeId>>& condition = make_shared<func::AlwaysTrue<EdgeId>>(),
            boost::function<void(EdgeId)> removal_handler = 0) :
            base(g,
                 And<EdgeId>(make_shared<TipCondition<Graph>>(g), condition),
                 removal_handler, LengthComparator<Graph>(g),
                 make_shared<LengthUpperBound<Graph>>(g, max_tip_length)) {
    }

private:
    DECL_LOGGER("TipClipper")
};

template<class Graph>
class DefaultTipClipper: public TipClipper<Graph> {
	typedef TipClipper<Graph> base;

	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;

public:

	DefaultTipClipper(Graph& g, size_t max_tip_length, size_t max_coverage,
			double max_relative_coverage,
			boost::function<void(EdgeId)> removal_handler = 0) :
			base(g, max_tip_length,
			     And<EdgeId>(make_shared<CoverageUpperBound<Graph>>(g,
									max_coverage),
							make_shared<RelativeCoverageTipCondition<Graph>>(
									g, max_relative_coverage)),
					removal_handler) {
	}

private:
	DECL_LOGGER("DefaultTipClipper")
};

template<class Graph>
class TopologyTipClipper: public TipClipper<Graph> {
    typedef TipClipper<Graph> base;

    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;

public:

    TopologyTipClipper(Graph& g, size_t max_tip_length,
                       size_t uniqueness_length, size_t plausibility_length,
                       boost::function<void(EdgeId)> removal_handler = 0) :
            base(g, max_tip_length,
                 make_shared<DefaultUniquenessPlausabilityCondition<Graph>>(g, uniqueness_length,
                 plausibility_length),
                 removal_handler) {
    }

private:
    DECL_LOGGER("TopologyTipClipper")
};

template<class Graph>
bool ClipTips(
        Graph& graph,
        size_t max_tip_length,
        shared_ptr<Predicate<typename Graph::EdgeId>> condition 
            = make_shared<func::AlwaysTrue<typename Graph::EdgeId>>(),
        boost::function<void(typename Graph::EdgeId)> raw_removal_handler = 0) {

    DEBUG("Max tip length: " << max_tip_length);

    omnigraph::TipClipper<Graph> tc(graph, max_tip_length, condition,
                                    raw_removal_handler);

    return tc.Process();
}

} // namespace omnigraph
