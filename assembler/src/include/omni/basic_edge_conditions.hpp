#pragma once

#include "func.hpp"

namespace omnigraph {

using namespace func;

template<class Graph>
class EdgeCondition: public Predicate<typename Graph::EdgeId> {
	typedef typename Graph::EdgeId EdgeId;

	const Graph& g_;
protected:

	EdgeCondition(const Graph& g) :
			g_(g) {
	}

	const Graph& g() const {
		return g_;
	}

};

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
class CoverageUpperBound: public EdgeCondition<Graph> {
	typedef typename Graph::EdgeId EdgeId;
	typedef EdgeCondition<Graph> base;
	//todo why size_t???
	const size_t max_coverage_;

public:

	CoverageUpperBound(const Graph& g, size_t max_coverage) :
			base(g), max_coverage_(max_coverage) {
	}

	bool Check(EdgeId e) const {
		return this->g().coverage(e) < max_coverage_;
	}

};

template<class Graph>
class LengthUpperBound: public EdgeCondition<Graph> {
	typedef typename Graph::EdgeId EdgeId;
	typedef EdgeCondition<Graph> base;

	const size_t max_length_;

public:

	LengthUpperBound(const Graph& g, size_t max_length) :
			base(g), max_length_(max_length) {
	}

	bool Check(EdgeId e) const {
		return this->g().length(e) <= max_length_;
	}

};

}
