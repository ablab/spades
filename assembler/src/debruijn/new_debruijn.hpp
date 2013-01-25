//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once
//#ifndef NEW_DEBRUIJN_HPP_
//#define NEW_DEBRUIJN_HPP_

#include "omni/abstract_conjugate_graph.hpp"
#include "omni/abstract_nonconjugate_graph.hpp"
#include "omni/coverage.hpp"
#include "omni/id_track_handler.hpp"
#include "sequence/sequence_tools.hpp"
#include "omni/concurrent_graph_component.hpp"

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

	EdgeData GlueData(const EdgeData &data1, const EdgeData &data2) const {
		return data2;
	}

	bool isSelfConjugate(const EdgeData &data) const {
		return data.nucls() == !(data.nucls());
	}

	EdgeData conjugate(const EdgeData &data) const {
		return EdgeData(!(data.nucls()));
	}

	VertexData conjugate(const VertexData &data) const {
		return VertexData();
	}

	const size_t length(const EdgeData& data) const {
		return data.nucls().size() - k_;
	}

};

template<class T>
class DeBruijnGraph: public T {

	friend class omnigraph::ConcurrentGraphComponent<DeBruijnGraph>;

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

//class ConjugateDeBruijnGraph: public AbstractConjugateGraph<DeBruijnMaster> {
//	typedef AbstractConjugateGraph<DeBruijnMaster> base;
//public:
//
//	typedef base::VertexId VertexId;
//	typedef base::EdgeId EdgeId;
//	typedef base::VertexData VertexData;
//	typedef base::EdgeData EdgeData;
//	typedef base::VertexIterator VertexIterator;
//
////	typedef typename super::SmartVertexIt SmartVertexIt;
////	typedef typename super::SmartEdgeIt SmartEdgeIt;
//private:
//	const size_t k_;
//	CoverageIndex<ConjugateDeBruijnGraph> coverage_index_;
//
//public:
//	ConjugateDeBruijnGraph(size_t k) :
//			base(DeBruijnMaster(k)), k_(k), coverage_index_(*this) {
//
//	}
//
//	virtual ~ConjugateDeBruijnGraph() {
//		DEBUG("~ConjugateDeBruijnGraph()");
//	}
//
//	CoverageIndex<ConjugateDeBruijnGraph>& coverage_index() {
//		return coverage_index_;
//	}
//
////	/**
////	 * Method sets coverage value for the edge
////	 */
////	void SetCoverage(EdgeId edge, size_t cov) {
////		coverage_index_->SetCoverage(edge, cov);
////	}
//
//	/**
//	 * Method returns average coverage of the edge
//	 */
//	double coverage(EdgeId edge) const {
//		return coverage_index_.coverage(edge);
//	}
//
////	/**
////	 * Method increases coverage value
////	 */
////	void IncCoverage(EdgeId edge, int toAdd) {
////		coverage_index_->IncCoverage(edge, toAdd);
////	}
////
////	/**
////	 * Method increases coverage value by 1
////	 */
////	void IncCoverage(EdgeId edge) {
////		coverage_index_->IncCoverage(edge);
////	}
//
//	/**
//	 * Method returns Sequence stored in the edge
//	 */
//	const Sequence& EdgeNucls(EdgeId edge) const {
//		return data(edge).nucls();
//	}
//
//	using base::AddVertex;
//	using base::AddEdge;
//
//	VertexId AddVertex() {
//		return AddVertex(VertexData());
//	}
//
//	EdgeId AddEdge(VertexId from, VertexId to, const Sequence &nucls) {
//		return AddEdge(from, to, EdgeData(nucls));
//	}
//
//	size_t k() const {
//		return k_;
//	}
//
//	Sequence VertexNucls(VertexId v) const {
//		if (OutgoingEdges(v).size() > 0) {
//			return EdgeNucls(OutgoingEdges(v)[0]).Subseq(0, k_);
//		} else if (OutgoingEdges(conjugate(v)).size() > 0) {
//			return !VertexNucls(conjugate(v));
//		}
//		VERIFY(false);
//		return Sequence();
//	}
//	//TODO:
//	/* It seems, that these two functions must be called  default_label, and str must sign lower two. */
//	std::string str(EdgeId edge) const {
//		//		return " ";
//
//		stringstream ss;
////		ss << /*edge << " " << */length(edge) << "(" << coverage(edge) << ")";
//		ss << edge << " " << length(edge) << "(" << coverage(edge) << ")";
//		return ss.str();
//
//	}
//
//	std::string str(VertexId v) const {
//		return " ";
//	}
//private:
//	DECL_LOGGER("ConjugateDeBruijnGraph")
//};
//
//class NonconjugateDeBruijnGraph: public AbstractNonconjugateGraph<DeBruijnMaster> {
//private:
//	typedef omnigraph::AbstractNonconjugateGraph<DeBruijnMaster> base;
//	const size_t k_;
//	CoverageIndex<NonconjugateDeBruijnGraph> coverage_index_;
//
//public:
//	NonconjugateDeBruijnGraph(size_t k) :
//			base(DeBruijnMaster(k)), k_(k), coverage_index_(*this) {
//	}
//
//	virtual ~NonconjugateDeBruijnGraph() {
//		DEBUG("~NonconjugateDeBruijnGraph()");
//	}
//
//	CoverageIndex<NonconjugateDeBruijnGraph>& coverage_index() {
//		return coverage_index_;
//	}
//
//	/**
//	 * Method returns Sequence stored in the edge
//	 */
//	const Sequence& EdgeNucls(EdgeId edge) const {
//		return data(edge).nucls();
//	}
//
//	size_t k() const {
//		return k_;
//	}
//
//	Sequence VertexNucls(VertexId v) const {
//		if (OutgoingEdges(v).size() > 0) {
//			return EdgeNucls(OutgoingEdges(v)[0]).Subseq(0, k_);
//		} else if (IncomingEdges(v).size() > 0) {
//			EdgeId inc = IncomingEdges(v)[0];
//			size_t length = EdgeNucls(inc).size();
//			return EdgeNucls(inc).Subseq(length - k_, length);
//		}
//		VERIFY(false);
//		return Sequence();
//	}
//
////	/**
////	 * Method sets coverage value for the edge
////	 */
////	void SetCoverage(EdgeId edge, size_t cov) {
////		coverage_index_->SetCoverage(edge, cov);
////	}
//
//	/**
//	 * Method returns average coverage of the edge
//	 */
//	double coverage(EdgeId edge) const {
//		return coverage_index_->coverage(edge);
//	}
//
////	/**
////	 * Method increases coverage value
////	 */
////	void IncCoverage(EdgeId edge, int toAdd) {
////		coverage_index_->IncCoverage(edge, toAdd);
////	}
////
////	/**
////	 * Method increases coverage value by 1
////	 */
////	void IncCoverage(EdgeId edge) {
////		coverage_index_->IncCoverage(edge);
////	}
//
//	using base::AddVertex;
//	using base::AddEdge;
//
//	virtual VertexId AddVertex() {
//		return AddVertex(VertexData());
//	}
//
//	virtual EdgeId AddEdge(VertexId from, VertexId to, const Sequence &nucls) {
//		return AddEdge(from, to, EdgeData(nucls));
//	}
//
//	//todo refactor
//	std::string str(EdgeId edge) const {
//		//		return " ";
//		stringstream ss;
//		ss << length(edge) << "(" << coverage(edge) << ")";
//		return ss.str();
//
//	}
//
//	std::string str(VertexId v) const {
//		return " ";
//	}
//
//private:
//	DECL_LOGGER("NonconjugateDeBruijnGraph")
//};
//
//typedef ConjugateDeBruijnGraph Graph;
//typedef Graph::EdgeId EdgeId;
//typedef Graph::VertexId VertexId;
//typedef NonconjugateDeBruijnGraph NCGraph;
//
//}

//#endif /* NEW_DEBRUIJN_HPP_ */

