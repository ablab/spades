//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
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
	CoverageData coverage_;
	CoverageData flanking_cov_;
	Sequence nucls_;
public:

	DeBruijnEdgeData(const Sequence &nucls) :
			nucls_(nucls) {
	}

	const Sequence& nucls() const {
		return nucls_;
	}

    void inc_raw_coverage(int value) {
        coverage_.inc_coverage(value);
    }

    void set_raw_coverage(unsigned coverage) {
        coverage_.set_coverage(coverage);
    }

    unsigned raw_coverage() const {
        return coverage_.coverage();
    }

    void inc_flanking_coverage(int value) {
        flanking_cov_.inc_coverage(value);
    }

    void set_flanking_coverage(unsigned flanking_coverage) {
        flanking_cov_.set_coverage(flanking_coverage);
    }

    //not length normalized
    unsigned flanking_coverage() const {
        return flanking_cov_.coverage();
    }

    size_t size() const {
		return nucls_.size();
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

	const EdgeData MergeData(const vector<const EdgeData*>& to_merge,
			bool safe_merging = true) const {
		vector<Sequence> ss;
		ss.reserve(to_merge.size());
		for (auto it = to_merge.begin(); it != to_merge.end(); ++it) {
			ss.push_back((*it)->nucls());
		}
		return EdgeData(MergeOverlappingSequences(ss, k_,safe_merging));
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

	const CoverageIndex<DeBruijnGraph>& coverage_index() const {
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
	    //todo add verify on vertex nucls consistency
		if (this->OutgoingEdgeCount(v) > 0) {
			return EdgeNucls(*(this->out_begin(v))).Subseq(0, k_);
		} else if (this->IncomingEdgeCount(v) > 0) {
			EdgeId inc = *(this->in_begin(v));
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
