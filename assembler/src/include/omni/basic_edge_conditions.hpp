//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

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
class CoverageUpperBound: public EdgeCondition<Graph> {
	typedef typename Graph::EdgeId EdgeId;
	typedef EdgeCondition<Graph> base;
	const double max_coverage_;

public:

	CoverageUpperBound(const Graph& g, double max_coverage) :
			base(g), max_coverage_(max_coverage) {
	}

	bool Check(EdgeId e) const {
		return math::le(this->g().coverage(e), max_coverage_);
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
