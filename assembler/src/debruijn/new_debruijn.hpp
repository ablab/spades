#ifndef NEW_DEBRUIJN_HPP_
#define NEW_DEBRUIJN_HPP_
#include "abstract_conjugate_graph.hpp"
#include "abstract_nonconjugate_graph.hpp"
#include "coverage.hpp"

using omnigraph::CoverageIndex;
namespace debruijn_graph {

class DeBruijnMaster;

class VertexData {
	friend class NewDeBruijnGraph;
	friend class DeBruinMaster;
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

};

class NewConjugateDeBruijnGraph: public AbstractConjugateGraph<VertexData,
		EdgeData, DeBruijnMaster> {
	friend class CoverageIndex<NewConjugateDeBruijnGraph> ;
public:
	typedef SmartVertexIterator<ObservableGraph<VertexId, EdgeId> >
			SmartVertexItarator;
	typedef SmartEdgeIterator<ObservableGraph<VertexId, EdgeId>>
			SmartEdgeItarator;
private:
	typedef AbstractConjugateGraph<VertexData, EdgeData, DeBruijnMaster> super;
	const size_t k_;
	CoverageIndex<NewConjugateDeBruijnGraph>* coverage_index_;

public:
	NewConjugateDeBruijnGraph(size_t k) :
		super(DeBruijnMaster(k)), k_(k) {
		coverage_index_ = new CoverageIndex<NewConjugateDeBruijnGraph> (*this);
		AddActionHandler(coverage_index_);
	}

	virtual ~NewConjugateDeBruijnGraph() {
		delete coverage_index_;
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

	const size_t length(EdgeId edge) const {
		return data(edge).nucls_.size() - k_;
	}

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

	std::string str(EdgeId edge) {
		//		return " ";

		stringstream ss;
		ss << length(edge) << "(" << coverage(edge) << ")";
		return ss.str();

	}

	std::string str(VertexId v) {
		return " ";
		//
		//		stringstream ss;
		//		ss << length(edge) << "(" << coverage(edge) << ")";
		//		return ss.str();

	}
};

class NewNonconjugateDeBruijnGraph: public AbstractNonconjugateGraph<
		VertexData, EdgeData, DeBruijnMaster> {
private:
	typedef AbstractNonconjugateGraph<VertexData, EdgeData, DeBruijnMaster>
			super;
	const size_t k_;
	CoverageIndex<NewNonconjugateDeBruijnGraph>* coverage_index_;

public:
	NewNonconjugateDeBruijnGraph(size_t k) :
		super(DeBruijnMaster(k)), k_(k) {
		coverage_index_ = new CoverageIndex<NewNonconjugateDeBruijnGraph> (*this);
		AddActionHandler(coverage_index_);
	}

	virtual ~NewNonconjugateDeBruijnGraph() {
	}

	/**
	 * Method returns Sequence stored in the edge
	 */
	const Sequence& EdgeNucls(EdgeId edge) const {
		return data(edge).nucls();
	}

	const size_t length(EdgeId edge) const {
		return data(edge).nucls_.size() - k_;
	}

	size_t k() const {
		return k_;
	}

	Sequence VertexNucls(VertexId v) const {
		if (OutgoingEdges(v).size() > 0) {
			return EdgeNucls(OutgoingEdges(v)[0]).Subseq(0, k_);
		} else if (IncomingEdges(v).size() > 0) {
			EdgeId inc = IncomingEdges(v)[0];
			return EdgeNucls(inc).Subseq(length(inc) - k_, length(inc));
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

	VertexId AddVertex() {
		return super::AddVertex(VertexData());
	}

	EdgeId AddEdge(VertexId from, VertexId to, const Sequence &nucls) {
		return super::AddEdge(from, to, EdgeData(nucls));
	}

};
}

#endif /* NEW_DEBRUIJN_HPP_ */

