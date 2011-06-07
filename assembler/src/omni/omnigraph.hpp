#ifndef OMNIGRAPH_HPP_
#define OMNIGRAPH_HPP_

#include "abstract_conjugate_graph.hpp"

namespace omnigraph {

class OmniVertex {
	size_t length_;
public:
	OmniVertex(size_t length) : length_(length) {}

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

	OmniEdge conjugate(const OmniEdge &data) {
		return data;
	}

	std::string str(const OmniEdge &edge) {
		return edge.str();
	}

	std::string str(const OmniVertex &v) {
		return v.str();
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
public:
	Omnigraph() : AbstractConjugateGraph<OmniVertex, OmniEdge, OmniDataMaster>(OmniDataMaster()) {

	}
};

}
#endif /* OMNI_GRAPH_HPP_ */
