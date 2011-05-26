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
public:
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

	const EdgeData MergeData(const vector<EdgeData *> &toMerge) const {
		SequenceBuilder sb;
		size_t coverage = 0;
		sb.append(toMerge[0]->nucls_.Subseq(0, k_));
		for (size_t i = 0; i < toMerge.size(); i++) {
			sb.append(toMerge[i]->nucls_.Subseq(k_));
			coverage += toMerge[i]->coverage_;
		}
		return EdgeData(sb.BuildSequence(), coverage);
	}

	pair<VertexData, pair<EdgeData, EdgeData> > SplitData(EdgeData &edge,
			size_t position) {
		size_t length = edge.nucls_.size() - k_;
		size_t coverage1 = edge.coverage_ * position / length;
		size_t coverage2 = edge.coverage_ * (length - position) / length;
		if (coverage1 == 0)
			coverage1 = 1;
		if (coverage2 == 0)
			coverage2 = 1;
		return make_pair(
				VertexData(),
				make_pair(
						EdgeData(edge.nucls_.Subseq(0, position + k_),
								coverage1),
						EdgeData(edge.nucls_.Subseq(position), coverage2)));
	}

	void GlueData(EdgeData &data1, EdgeData &data2) {
		data2.coverage_ += data1.coverage_;
	}

	bool isSelfConjugate(EdgeData &data) {
		return data.nucls_ == !(data.nucls_);
	}

	EdgeData conjugate(EdgeData &data) {
		return EdgeData(!(data.nucls_), data.coverage_);
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

	const size_t length(EdgeId edge) const {
		return data(edge).nucls_.size() - k_;
	}

};

}

#endif /* NEW_DEBRUIJN_HPP_ */
