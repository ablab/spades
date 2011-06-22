#ifndef NEW_DEBRUIJN_HPP_
#define NEW_DEBRUIJN_HPP_
#include "abstract_conjugate_graph.hpp"
#include "abstract_nonconjugate_graph.hpp"

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
				make_pair(
						EdgeData(edge.nucls_.Subseq(0, position + k_)),
						EdgeData(edge.nucls_.Subseq(position))));
	}

	void GlueData(EdgeData &data1, EdgeData &data2) {
	}

	bool isSelfConjugate(EdgeData &data) {
		return data.nucls_ == !(data.nucls_);
	}

	EdgeData conjugate(EdgeData &data) {
		return EdgeData(!(data.nucls_));
	}
};

class NewConjugateDeBruijnGraph: public AbstractConjugateGraph<VertexData, EdgeData,
		DeBruijnMaster> {
private:
	typedef AbstractConjugateGraph<VertexData, EdgeData, DeBruijnMaster> super;
	const size_t k_;
public:
	NewConjugateDeBruijnGraph(size_t k) :
		super(DeBruijnMaster(k)), k_(k) {
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

};

class NewNonconjugateDeBruijnGraph: public AbstractNonconjugateGraph<VertexData, EdgeData,
		DeBruijnMaster> {
private:
	typedef AbstractNonconjugateGraph<VertexData, EdgeData, DeBruijnMaster> super;
	const size_t k_;
public:
	NewNonconjugateDeBruijnGraph(size_t k) :
		super(DeBruijnMaster(k)), k_(k) {
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

};
}

#endif /* NEW_DEBRUIJN_HPP_ */
