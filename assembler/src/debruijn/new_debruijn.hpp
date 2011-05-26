#ifndef NEW_DEBRUIJN_HPP_
#define NEW_DEBRUIJN_HPP_
#include "abstract_conjugate_graph.hpp"

namespace debruijn_graph {

class DeBruijnMaster;

class VertexData {
	friend class NewDeBruijnGraph;
	friend class DeBruinMaster;
};

class EdgeData {
	friend class NewDeBruijnGraph;
	friend class DeBruinMaster;
private:
	const Sequence nucls_;
	size_t coverage_;

	EdgeData(const Sequence &nucls, size_t coverage) :
		nucls_(nucls), coverage_(coverage) {
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

	const EdgeData MergeData(const vector<EdgeData&> &toMerge) const {
		SequenceBuilder sb;
		sb.append(toMerge[0].nucls_.Subseq(0, k_));
		for (size_t i = 0; i < toMerge.size(); i++) {
			sb.append(toMerge[i].nucls_.Subseq(k_));
		}
	}
};

class NewDeBruijnGraph: public AbstractConjugateGraph<VertexData, EdgeData,
		DeBruijnMaster> {
private:
	typedef AbstractConjugateGraph<VertexData, EdgeData, DeBruijnMaster> super;
	const size_t k_;
public:
	NewDeBruijnGraph(size_t k) :
		super(DeBruijnMaster(k)), k_(k) {
	}

	/**
	 * Method returns Sequence stored in the edge
	 */
	const Sequence& EdgeNucls(EdgeId edge) const {
		return data(edge).nucls();
	}

};

}

#endif /* NEW_DEBRUIJN_HPP_ */
