#ifndef NEW_DEBRUIJN_HPP_
#define NEW_DEBRUIJN_HPP_

#include "omni/abstract_conjugate_graph.hpp"
#include "omni/abstract_nonconjugate_graph.hpp"
#include "omni/coverage.hpp"
#include "omni/ID_track_handler.hpp"

using omnigraph::CoverageIndex;
namespace debruijn_graph {

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
public:
	const Sequence nucls_;

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

	const EdgeData MergeData(const vector<const EdgeData *>& toMerge) const {
		SequenceBuilder sb;
		sb.append(toMerge[0]->nucls_.Subseq(0, k_));
		for (size_t i = 0; i < toMerge.size(); i++) {
			sb.append(toMerge[i]->nucls_.Subseq(k_));
		}
		return EdgeData(sb.BuildSequence());
	}

	pair<VertexData, pair<EdgeData, EdgeData> > SplitData(const EdgeData &edge,
			size_t position) const {
		return make_pair(
				VertexData(),
				make_pair(EdgeData(edge.nucls_.Subseq(0, position + k_)),
						EdgeData(edge.nucls_.Subseq(position))));
	}

	EdgeData GlueData(const EdgeData &data1, const EdgeData &data2) const {
		return data2;
	}

	bool isSelfConjugate(const EdgeData &data) const {
		return data.nucls_ == !(data.nucls_);
	}

	EdgeData conjugate(const EdgeData &data) const {
		return EdgeData(!(data.nucls_));
	}

	VertexData conjugate(const VertexData &data) const {
		return VertexData();
	}

	std::string str(const EdgeData &edge) const {
		return " ";
	}

	std::string str(const VertexData &v) const {
		return "";
	}

	const size_t length(EdgeData data) const {
		return data.nucls_.size() - k_;
	}

};

class ConjugateDeBruijnGraph: public AbstractConjugateGraph<DeBruijnMaster> {
	typedef AbstractConjugateGraph<DeBruijnMaster> base;
public:

private:
	const size_t k_;
	CoverageIndex<ConjugateDeBruijnGraph>* coverage_index_;
	DECL_LOGGER("ConjugateDeBruijnGraph")

public:
	ConjugateDeBruijnGraph(size_t k) :
		base(DeBruijnMaster(k)), k_(k) {
		coverage_index_ = new CoverageIndex<ConjugateDeBruijnGraph> (*this);
	}

	virtual ~ConjugateDeBruijnGraph() {
		DEBUG("~ConjugateDeBruijnGraph() start");
		delete coverage_index_;
		DEBUG("~ConjugateDeBruijnGraph() finished");
	}

	CoverageIndex<ConjugateDeBruijnGraph>& coverage_index() {
		return *coverage_index_;
	}

//	/**
//	 * Method sets coverage value for the edge
//	 */
//	void SetCoverage(EdgeId edge, size_t cov) {
//		coverage_index_->SetCoverage(edge, cov);
//	}

	/**
	 * Method returns average coverage of the edge
	 */
	double coverage(EdgeId edge) const {
		return coverage_index_->coverage(edge);
	}

//	/**
//	 * Method increases coverage value
//	 */
//	void IncCoverage(EdgeId edge, int toAdd) {
//		coverage_index_->IncCoverage(edge, toAdd);
//	}
//
//	/**
//	 * Method increases coverage value by 1
//	 */
//	void IncCoverage(EdgeId edge) {
//		coverage_index_->IncCoverage(edge);
//	}

	/**
	 * Method returns Sequence stored in the edge
	 */
	const Sequence& EdgeNucls(EdgeId edge) const {
		return data(edge).nucls();
	}

	VertexId AddVertex() {
		return base::AddVertex(VertexData());
	}

	EdgeId AddEdge(VertexId from, VertexId to, const Sequence &nucls) {
		return base::AddEdge(from, to, EdgeData(nucls));
	}

	size_t k() const {
		return k_;
	}

	Sequence VertexNucls(VertexId v) const {
		if (OutgoingEdges(v).size() > 0) {
			return EdgeNucls(OutgoingEdges(v)[0]).Subseq(0, k_);
		} else if (OutgoingEdges(conjugate(v)).size() > 0) {
			return !VertexNucls(conjugate(v));
		}
		assert(false);
	}
	//TODO:
	/* It seems, that these two functions must be called  default_label, and str must sign lower two. */
	std::string str(EdgeId edge) const {
		//		return " ";

		stringstream ss;
//		ss << /*edge << " " << */length(edge) << "(" << coverage(edge) << ")";
		ss << edge << " " << length(edge) << "(" << coverage(edge) << ")";
		return ss.str();

	}

	std::string str(VertexId v) const {
		return " ";
	}

	//todo extract from here!!!
	std::string toPrint(VertexId v,
			IdTrackHandler<ConjugateDeBruijnGraph> &id_handler) const {
		stringstream ss;
		ss << "Vertex " << id_handler.ReturnIntId(v) << " ~ "
				<< id_handler.ReturnIntId(conjugate(v)) << " .";
		return ss.str();
	}

	//todo extract from here!!!
	std::string toPrint(EdgeId e,
			IdTrackHandler<ConjugateDeBruijnGraph> &id_handler) const {
		stringstream ss;
		ss << "Edge " << id_handler.ReturnIntId(e) << " : "
				<< id_handler.ReturnIntId(EdgeStart(e)) << " -> "
				<< id_handler.ReturnIntId(EdgeEnd(e)) << ", l = " << length(e)
				<< " ~ " << id_handler.ReturnIntId(conjugate(e)) << " .";
		return ss.str();
	}

};

class NonconjugateDeBruijnGraph: public AbstractNonconjugateGraph<DeBruijnMaster> {
private:
	typedef omnigraph::AbstractNonconjugateGraph<DeBruijnMaster> base;
	const size_t k_;
	CoverageIndex<NonconjugateDeBruijnGraph>* coverage_index_;DECL_LOGGER("NonconjugateDeBruijnGraph")

public:
	NonconjugateDeBruijnGraph(size_t k) :
		base(DeBruijnMaster(k)), k_(k) {
		coverage_index_ = new CoverageIndex<NonconjugateDeBruijnGraph> (*this);
	}

	virtual ~NonconjugateDeBruijnGraph() {
		DEBUG("~NonconjugateDeBruijnGraph()");
	}

	CoverageIndex<NonconjugateDeBruijnGraph>& coverage_index() {
		return *coverage_index_;
	}

	/**
	 * Method returns Sequence stored in the edge
	 */
	const Sequence& EdgeNucls(EdgeId edge) const {
		return data(edge).nucls();
	}

	size_t k() const {
		return k_;
	}

	Sequence VertexNucls(VertexId v) const {
		if (OutgoingEdges(v).size() > 0) {
			return EdgeNucls(OutgoingEdges(v)[0]).Subseq(0, k_);
		} else if (IncomingEdges(v).size() > 0) {
			EdgeId inc = IncomingEdges(v)[0];
			size_t length = EdgeNucls(inc).size();
			return EdgeNucls(inc).Subseq(length - k_, length);
		}
		assert(false);
	}

//	/**
//	 * Method sets coverage value for the edge
//	 */
//	void SetCoverage(EdgeId edge, size_t cov) {
//		coverage_index_->SetCoverage(edge, cov);
//	}

	/**
	 * Method returns average coverage of the edge
	 */
	double coverage(EdgeId edge) const {
		return coverage_index_->coverage(edge);
	}

//	/**
//	 * Method increases coverage value
//	 */
//	void IncCoverage(EdgeId edge, int toAdd) {
//		coverage_index_->IncCoverage(edge, toAdd);
//	}
//
//	/**
//	 * Method increases coverage value by 1
//	 */
//	void IncCoverage(EdgeId edge) {
//		coverage_index_->IncCoverage(edge);
//	}

	virtual VertexId AddVertex() {
		return base::AddVertex(VertexData());
	}

	virtual EdgeId AddEdge(VertexId from, VertexId to, const Sequence &nucls) {
		return base::AddEdge(from, to, EdgeData(nucls));
	}

	std::string str(EdgeId edge) const {
		//		return " ";
		stringstream ss;
		ss << length(edge) << "(" << coverage(edge) << ")";
		return ss.str();

	}

	std::string str(VertexId v) const {
		return " ";
	}

	//todo extract from here!!!
	std::string toPrint(VertexId v,
			IdTrackHandler<NonconjugateDeBruijnGraph> &id_handler) const {
		stringstream ss;
		ss << "Vertex " << id_handler.ReturnIntId(v) << " .";
		return ss.str();
	}

	//todo extract from here!!!
	std::string toPrint(EdgeId e,
			IdTrackHandler<NonconjugateDeBruijnGraph> &id_handler) const {
		stringstream ss;
		ss << "Edge " << id_handler.ReturnIntId(e) <<" : " << id_handler.ReturnIntId(EdgeStart(e))
				<< " -> " << id_handler.ReturnIntId(EdgeEnd(e))<<", l = "<< length(e) <<" .";
		return ss.str();
	}

};

typedef ConjugateDeBruijnGraph Graph;
typedef Graph::EdgeId EdgeId;
typedef Graph::VertexId VertexId;
typedef NonconjugateDeBruijnGraph NCGraph;

}

#endif /* NEW_DEBRUIJN_HPP_ */

