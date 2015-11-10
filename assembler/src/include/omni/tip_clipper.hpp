//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

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
#include "func.hpp"
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
class MismatchTipCondition : public EdgeCondition<Graph> {
    typedef EdgeCondition<Graph> base;

    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;

    size_t max_diff_;
    size_t Hamming(EdgeId edge1, EdgeId edge2) const {
        size_t len = std::min(this->g().length(edge1), this->g().length(edge2));
        size_t cnt = 0;
        Sequence seq1 = this->g().EdgeNucls(edge1);
        Sequence seq2 = this->g().EdgeNucls(edge2);
        for(size_t i = 0; i < len; i++) {
            if(seq1[i] != seq2[i])
                cnt++;
        }
        return cnt;
    }

public:
    static const size_t INF = size_t(-1);
    MismatchTipCondition(const Graph& g, size_t max_diff) : base(g), max_diff_(max_diff) {
    }

    /**
     * This method checks if given edge topologically looks like a tip.
     * @param edge edge vertex to be checked
     * @return true if edge judged to be tip and false otherwise.
     */
    /*virtual*/ bool Check(EdgeId e) const {
        if(max_diff_ == INF) {
            return true;
        }
        auto alternatives = this->g().OutgoingEdges(this->g().EdgeStart(e));
        for(auto it = alternatives.begin(); it != alternatives.end(); ++it) {
            if(e != *it && this->g().length(e) < this->g().length(*it) && Hamming(e, *it) <= max_diff_) {
                return true;
            }
        }
        return false;
    }

};

template<class Graph>
shared_ptr<func::Predicate<typename Graph::EdgeId>> AddTipCondition(const Graph& g,
                                                                  shared_ptr<func::Predicate<typename Graph::EdgeId>> condition) {
    return func::And<typename Graph::EdgeId>(
            make_shared<TipCondition<Graph>>(g),
            condition);
}

template<class Graph>
shared_ptr<func::Predicate<typename Graph::EdgeId>>
NecessaryTipCondition(const Graph& g, size_t max_length, double max_coverage) {
    return AddTipCondition(g, func::And<typename Graph::EdgeId>(std::make_shared<LengthUpperBound<Graph>>(g, max_length),
                               std::make_shared<CoverageUpperBound<Graph>>(g, max_coverage)));
}

//template<class Graph>
//bool ClipTips(
//        Graph& g,
//        size_t max_length,
//        shared_ptr<Predicate<typename Graph::EdgeId>> condition
//            = make_shared<func::AlwaysTrue<typename Graph::EdgeId>>(),
//        std::function<void(typename Graph::EdgeId)> removal_handler = 0) {
//
//    omnigraph::EdgeRemovingAlgorithm<Graph> tc(g,
//                                               AddTipCondition(g, condition),
//                                               removal_handler);
//
//    return tc.Run(LengthComparator<Graph>(g),
//                      make_shared<LengthUpperBound<Graph>>(g, max_length));
//}

} // namespace omnigraph
