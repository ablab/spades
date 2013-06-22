//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include "omni/abstract_conjugate_graph.hpp"
#include "omni/abstract_nonconjugate_graph.hpp"
#include "omni/coverage.hpp"
#include "omni/id_track_handler.hpp"
#include "sequence/sequence_tools.hpp"

namespace debruijn_graph {
using omnigraph::CoverageIndex;

class DeBruijnMaster;

class DeBruijnVertexData {
	friend class NewDeBruijnGraph;
	friend class DeBruinMaster;
public:
	DeBruijnVertexData() {

	}
};

class DeBruijnEdgeData {
	friend class NewDeBruijnGraph;
	friend class DeBruinMaster;
	Sequence nucls_;
public:

	DeBruijnEdgeData(const Sequence &nucls) :
			nucls_(nucls) {
	}

	const Sequence& nucls() const {
		return nucls_;
	}
};

class DeBruijnMaster {
private:
	const size_t k_;

public:
	typedef DeBruijnVertexData VertexData;
	typedef DeBruijnEdgeData EdgeData;

	DeBruijnMaster(size_t k) :
			k_(k) {
	}

	const EdgeData MergeData(const vector<const EdgeData*>& to_merge) const {
		vector<Sequence> ss;
		ss.reserve(to_merge.size());
		for (auto it = to_merge.begin(); it != to_merge.end(); ++it) {
			ss.push_back((*it)->nucls());
		}
		return EdgeData(MergeOverlappingSequences(ss, k_));
	}

	pair<VertexData, pair<EdgeData, EdgeData> > SplitData(const EdgeData &edge,
			size_t position) const {
		return make_pair(
				VertexData(),
				make_pair(EdgeData(edge.nucls().Subseq(0, position + k_)),
						EdgeData(edge.nucls().Subseq(position))));
	}

	EdgeData GlueData(const EdgeData & /*data1*/, const EdgeData &data2) const {
		return data2;
	}

	bool isSelfConjugate(const EdgeData &data) const {
		return data.nucls() == !(data.nucls());
	}

	EdgeData conjugate(const EdgeData &data) const {
		return EdgeData(!(data.nucls()));
	}

	VertexData conjugate(const VertexData & /*data*/) const {
		return VertexData();
	}

	size_t length(const EdgeData& data) const {
		return data.nucls().size() - k_;
	}

};

template<class T>
class DeBruijnGraph: public T {
public:
	typedef T base;
	typedef typename base::VertexId VertexId;
	typedef typename base::EdgeId EdgeId;
	typedef typename base::VertexData VertexData;
	typedef typename base::EdgeData EdgeData;
	typedef typename base::VertexIterator VertexIterator;
	typedef VertexIterator iterator; // for for_each
	typedef const VertexIterator const_iterator; // for for_each
private:
	const size_t k_;
	CoverageIndex<DeBruijnGraph> coverage_index_;

public:
	DeBruijnGraph(size_t k) :
			base(DeBruijnMaster(k)), k_(k), coverage_index_(*this) {
	}

	CoverageIndex<DeBruijnGraph>& coverage_index() {
		return coverage_index_;
	}

	/**
	 * Method returns average coverage of the edge
	 */
	double coverage(EdgeId edge) const {
		return coverage_index_.coverage(edge);
	}

	using base::AddVertex;
	using base::AddEdge;

	VertexId AddVertex() {
		return AddVertex(VertexData());
	}

	EdgeId AddEdge(VertexId from, VertexId to, const Sequence &nucls) {
		VERIFY(nucls.size() > k_);
		return AddEdge(from, to, EdgeData(nucls));
	}

	size_t k() const {
		return k_;
	}

	/**
	 * Method returns Sequence stored in the edge
	 */
	const Sequence& EdgeNucls(EdgeId edge) const {
		return this->data(edge).nucls();
	}

	const Sequence VertexNucls(VertexId v) const {
		if (this->OutgoingEdges(v).size() > 0) {
			return EdgeNucls(this->OutgoingEdges(v)[0]).Subseq(0, k_);
		} else if (this->IncomingEdges(v).size() > 0) {
			EdgeId inc = this->IncomingEdges(v)[0];
			size_t length = EdgeNucls(inc).size();
			return EdgeNucls(inc).Subseq(length - k_, length);
		}
		VERIFY(false);
		return Sequence();
	}

private:
	DECL_LOGGER("DeBruijnGraph")
};

typedef DeBruijnGraph<AbstractConjugateGraph<DeBruijnMaster>> ConjugateDeBruijnGraph;
typedef DeBruijnGraph<AbstractNonconjugateGraph<DeBruijnMaster>> NonconjugateDeBruijnGraph;

typedef ConjugateDeBruijnGraph Graph;
typedef Graph::EdgeId EdgeId;
typedef Graph::VertexId VertexId;
typedef NonconjugateDeBruijnGraph NCGraph;
}
