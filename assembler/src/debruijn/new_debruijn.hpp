#ifndef NEW_DEBRUIJN_HPP_
#define NEW_DEBRUIJN_HPP_
#include "abstract_conjugate_graph.hpp"
#include "abstract_nonconjugate_graph.hpp"
#include "coverage.hpp"
#include "ID_track_handler.hpp"

using omnigraph::CoverageIndex;
namespace debruijn_graph {

class DeBruijnMaster;

class VertexData {
	friend class NewDeBruijnGraph;
	friend class DeBruinMaster;
public:
	VertexData() {

	}
};

class EdgeData {
	friend class NewDeBruijnGraph;
	friend class DeBruinMaster;
public:
	const Sequence nucls_;

	EdgeData(const Sequence &nucls) :
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
	DeBruijnMaster(size_t k) :
		k_(k) {
	}

	const EdgeData MergeData(const vector<EdgeData *> &toMerge) const {
		SequenceBuilder sb;
		sb.append(toMerge[0]->nucls_.Subseq(0, k_));
		for (size_t i = 0; i < toMerge.size(); i++) {
			sb.append(toMerge[i]->nucls_.Subseq(k_));
		}
		return EdgeData(sb.BuildSequence());
	}

	pair<VertexData, pair<EdgeData, EdgeData> > SplitData(EdgeData &edge,
			size_t position) {
		return make_pair(
				VertexData(),
				make_pair(EdgeData(edge.nucls_.Subseq(0, position + k_)),
						EdgeData(edge.nucls_.Subseq(position))));
	}

	EdgeData GlueData(EdgeData &data1, EdgeData &data2) {
		return data2;
	}

	bool isSelfConjugate(const EdgeData &data) {
		return data.nucls_ == !(data.nucls_);
	}

	EdgeData conjugate(const EdgeData &data) {
		return EdgeData(!(data.nucls_));
	}

	VertexData conjugate(const VertexData &data) {
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

class ConjugateDeBruijnGraph: public AbstractConjugateGraph<VertexData,
		EdgeData, DeBruijnMaster> {
public:
	typedef SmartVertexIterator<ObservableGraph<VertexId, EdgeId> >
			SmartVertexItarator;
	typedef SmartEdgeIterator<ObservableGraph<VertexId, EdgeId>>
			SmartEdgeItarator;
private:
	typedef AbstractConjugateGraph<VertexData, EdgeData, DeBruijnMaster> super;
	const size_t k_;
	CoverageIndex<ConjugateDeBruijnGraph>* coverage_index_;
	DECL_LOGGER("ConjugateDeBruijnGraph")

public:
	ConjugateDeBruijnGraph(size_t k) :
		super(DeBruijnMaster(k)), k_(k) {
		coverage_index_ = new CoverageIndex<ConjugateDeBruijnGraph> (*this);
	}

	virtual ~ConjugateDeBruijnGraph() {
		DEBUG("~ConjugateDeBruijnGraph() start");
		delete coverage_index_;
		DEBUG("~ConjugateDeBruijnGraph() finished");
	}

	template<class Stream, class ReadThreader>
	void FillCoverage(Stream& stream, const ReadThreader& threader) {
		coverage_index_->FillIndex(stream, threader);
	}

	/**
	 * Method sets coverage value for the edge
	 */
	void SetCoverage(EdgeId edge, size_t cov) {
		coverage_index_->SetCoverage(edge, cov);
	}

	/**
	 * Method returns average coverage of the edge
	 */
	double coverage(EdgeId edge) const {
		return coverage_index_->coverage(edge);
	}

	/**
	 * Method increases coverage value
	 */
	void IncCoverage(EdgeId edge, int toAdd) {
		coverage_index_->IncCoverage(edge, toAdd);
	}

	/**
	 * Method increases coverage value by 1
	 */
	void IncCoverage(EdgeId edge) {
		coverage_index_->IncCoverage(edge);
	}

	/**
	 * Method returns Sequence stored in the edge
	 */
	const Sequence& EdgeNucls(EdgeId edge) const {
		return data(edge).nucls();
	}

	using super::AddVertex;
	using super::AddEdge;

	VertexId AddVertex() {
		return super::AddVertex(VertexData());
	}

	EdgeId AddEdge(VertexId from, VertexId to, const Sequence &nucls) {
		return super::AddEdge(from, to, EdgeData(nucls));
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
		ss << /*edge << " " << */length(edge) << "(" << coverage(edge) << ")";
		return ss.str();

	}

	std::string str(VertexId v) const {
		return " ";
	}

	std::string toPrint(VertexId v,
			IdTrackHandler<ConjugateDeBruijnGraph> &id_handler) const {
		stringstream ss;
		ss << "Vertex " << id_handler.ReturnIntId(v) << " ~ "
				<< id_handler.ReturnIntId(conjugate(v)) << " .";
		return ss.str();
	}

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

class NonconjugateDeBruijnGraph: public AbstractNonconjugateGraph<VertexData,
		EdgeData, DeBruijnMaster> {
private:
	typedef AbstractNonconjugateGraph<VertexData, EdgeData, DeBruijnMaster>
			super;
	const size_t k_;
	CoverageIndex<NonconjugateDeBruijnGraph>* coverage_index_;DECL_LOGGER("NonconjugateDeBruijnGraph")

public:
	NonconjugateDeBruijnGraph(size_t k) :
		super(DeBruijnMaster(k)), k_(k) {
		coverage_index_ = new CoverageIndex<NonconjugateDeBruijnGraph> (*this);
	}

	virtual ~NonconjugateDeBruijnGraph() {
		DEBUG("~NonconjugateDeBruijnGraph()");
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

	template<class Stream, class ReadThreader>
	void FillCoverage(Stream& stream, const ReadThreader& threader) {
		coverage_index_->FillIndex(stream, threader);
	}

	/**
	 * Method sets coverage value for the edge
	 */
	void SetCoverage(EdgeId edge, size_t cov) {
		coverage_index_->SetCoverage(edge, cov);
	}

	/**
	 * Method returns average coverage of the edge
	 */
	double coverage(EdgeId edge) const {
		return coverage_index_->coverage(edge);
	}

	/**
	 * Method increases coverage value
	 */
	void IncCoverage(EdgeId edge, int toAdd) {
		coverage_index_->IncCoverage(edge, toAdd);
	}

	/**
	 * Method increases coverage value by 1
	 */
	void IncCoverage(EdgeId edge) {
		coverage_index_->IncCoverage(edge);
	}

	using super::AddVertex;
	using super::AddEdge;

	virtual VertexId AddVertex() {
		return super::AddVertex(VertexData());
	}

	virtual EdgeId AddEdge(VertexId from, VertexId to, const Sequence &nucls) {
		return super::AddEdge(from, to, EdgeData(nucls));
	}

	std::string str(EdgeId edge) const {
		//		return " ";
		stringstream ss;
		ss << master_.length(data(edge)) << "(" << coverage(edge) << ")";
		return ss.str();

	}

	std::string str(VertexId v) const {
		return " ";
	}

	std::string toPrint(VertexId v,
			IdTrackHandler<NonconjugateDeBruijnGraph> &id_handler) const {
		stringstream ss;
		ss << "Vertex " << id_handler.ReturnIntId(v) << " .";
		return ss.str();
	}

	std::string toPrint(EdgeId e,
			IdTrackHandler<NonconjugateDeBruijnGraph> &id_handler) const {
		stringstream ss;
		ss << "Edge " << id_handler.ReturnIntId(e) <<" : " << id_handler.ReturnIntId(EdgeStart(e)) << " -> " << id_handler.ReturnIntId(EdgeEnd(e))<<", l = "<< master_.length(data(e)) <<" .";
		return ss.str();
	}

};
}

#endif /* NEW_DEBRUIJN_HPP_ */

