#ifndef OMNIGRAPH_HPP_
#define OMNIGRAPH_HPP_

#include "abstract_conjugate_graph.hpp"

namespace omnigraph {

class OmniVertex {
	size_t length_;
public:
	OmniVertex(size_t length) : length_(length) {}
	OmniVertex(Sequence sequence) : length_(sequence.size()) {}

	size_t length() const {
		return length_;
	}

	std::string str() const {
		return ToString(length_);
	}
};

class OmniEdge {
	size_t length_;
public:
	OmniEdge(size_t length) : length_(length) {}
	OmniEdge(size_t length, Sequence sequence) : length_(length) {}

	size_t length() const {
		return length_;
	}

	std::string str() const {
		return ToString(length_);
	}

	bool operator==(const OmniEdge &data) const {
		return length_ == data.length_;
	}
};

class OmniDataMaster {
public:
	OmniDataMaster() {}

	bool isSelfConjugate(const OmniEdge &data) {
		return false;
	}

	OmniVertex conjugate(const OmniVertex &data) {
		return data;
	}

	OmniEdge conjugate(const OmniEdge &data) {
		return data;
	}

	std::string str(const OmniEdge &edge) const {
		return edge.str();
	}

	std::string str(const OmniVertex &v) const {
		return v.str();
	}

	size_t length(const OmniEdge &edge) const {
		return edge.length();
	}

	size_t length(const OmniVertex &v) const {
		return v.length();
	}

	bool equals(const OmniEdge &data1, const OmniEdge &data2) {
		return data1 == data2;
	}

	OmniEdge MergeData(vector<OmniEdge*> toMerge) {
		size_t length = 0;
		for (auto it = toMerge.begin(); it != toMerge.end(); ++it) {
			length += (*it)->length();
		}
		return OmniEdge(length);
	}
};

class Omnigraph : public AbstractConjugateGraph<OmniVertex, OmniEdge, OmniDataMaster> {
	CoverageIndex<Omnigraph>* coverage_index_;

public:
	Omnigraph() : AbstractConjugateGraph<OmniVertex, OmniEdge, OmniDataMaster>(OmniDataMaster()) {
		coverage_index_ = new CoverageIndex<Omnigraph>(*this);
	}

	virtual ~Omnigraph() {
		delete coverage_index_;
	}

	void SetCoverage(EdgeId edge, size_t cov) {
		coverage_index_->SetCoverage(edge, cov);
	}

	double coverage(EdgeId edge) const {
		return coverage_index_->coverage(edge);
	}

	void IncCoverage(EdgeId edge, int toAdd) {
		coverage_index_->IncCoverage(edge, toAdd);
	}

	void IncCoverage(EdgeId edge) {
		coverage_index_->IncCoverage(edge);
	}
};

}
#endif /* OMNIGRAPH_HPP_ */
