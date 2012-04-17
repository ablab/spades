//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef OMNIGRAPH_NUCL_HPP_
#define OMNIGRAPH_NUCL_HPP_

#include "abstract_conjugate_graph.hpp"

namespace omnigraph {

class OmniNuclVertex {
	Sequence sequence_;
public:
	OmniNuclVertex(Sequence sequence) : sequence_(sequence) {}

	size_t length() const {
		return sequence_.size();
	}

	std::string str() const {
		size_t len = length();
		return ToString(len);
	}

	Sequence nucls() const {
		return sequence_;
	}
};

class OmniNuclEdge {
	size_t length_;
	Sequence sequence_;
public:
	OmniNuclEdge(size_t length, Sequence sequence) : length_(length), sequence_(sequence) {}

	size_t length() const {
		return length_;
	}

	std::string str() const {
		size_t len = length();
		return ToString(len);
	}

	Sequence nucls() const {
		return sequence_;
	}

	bool operator==(const OmniNuclEdge &data) const {
		return length_ == data.length_ && sequence_ == data.sequence_;
	}
};

class OmniNuclDataMaster {
public:
	OmniNuclDataMaster() {}

	bool isSelfConjugate(const OmniNuclEdge &data) const {
		return data.nucls() == !data.nucls();
	}

	OmniNuclVertex conjugate(const OmniNuclVertex &data) const {
		return OmniNuclVertex(!data.nucls());
	}

	OmniNuclEdge conjugate(const OmniNuclEdge &data) const {
		return OmniNuclEdge(data.length(), !data.nucls());
	}

	std::string str(const OmniNuclEdge &edge) const {
		return edge.str();
	}

	std::string str(const OmniNuclVertex &v) const {
		return v.str();
	}

	size_t length(const OmniNuclEdge &edge) const {
		return edge.length();
	}

	size_t length(const OmniNuclVertex &v) const {
		return v.length();
	}

	Sequence nucls(const OmniNuclEdge &edge) const {
		return edge.nucls();
	}

	Sequence nucls(const OmniNuclVertex &v) const {
		return v.nucls();
	}

	bool equals(const OmniNuclEdge &data1, const OmniNuclEdge &data2) const {
		return data1 == data2;
	}

	OmniNuclEdge MergeData(vector<OmniNuclEdge const*> toMerge) const {
		size_t length = 0;
		SequenceBuilder sb;
		for (auto it = toMerge.begin(); it != toMerge.end(); ++it) {
			if (length == 0) {
				sb.append((*it)->nucls());
			} else {
				sb.append((*it)->nucls().Subseq((*it)->nucls().size() - (*it)->length()));
			}
			length += (*it)->length();
		}
		return OmniNuclEdge(length, sb.BuildSequence());
	}

	OmniNuclEdge GlueData(const OmniNuclEdge &edge1, const OmniNuclEdge &edge2) const {
		return edge2;
	}

	pair<OmniNuclVertex, pair<OmniNuclEdge, OmniNuclEdge> > SplitData(OmniNuclEdge const& edge,
			size_t position) const {
		size_t k = edge.nucls().size() - edge.length();
		return make_pair(
				OmniNuclVertex(edge.nucls().Subseq(position, position + k)),
				make_pair(OmniNuclEdge(position, edge.nucls().Subseq(0, position + k)),
						OmniNuclEdge(edge.length() - position, edge.nucls().Subseq(position))));
	}
};

class OmnigraphNucl : public AbstractConjugateGraph<OmniNuclVertex, OmniNuclEdge, OmniNuclDataMaster> {
	CoverageIndex<OmnigraphNucl>* coverage_index_;

public:
	OmnigraphNucl() : AbstractConjugateGraph<OmniNuclVertex, OmniNuclEdge, OmniNuclDataMaster>(OmniNuclDataMaster()) {
		coverage_index_ = new CoverageIndex<OmnigraphNucl>(*this);
	}

	virtual ~OmnigraphNucl() {
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

	Sequence nucls(const EdgeId edge) const {
		return master().nucls(data(edge));
	}

	Sequence nucls(const VertexId v) const {
		return master().nucls(data(v));
	}

};

}
#endif /* OMNIGRAPH_NUCL_HPP_ */
