#ifndef OMNIGRAPH_HPP_
#define OMNIGRAPH_HPP_

#include "abstract_conjugate_graph.hpp"

namespace omnigraph {

class OmniVertex {
	size_t length_;
public:
	OmniVertex(size_t length) : length_(length) {}

	size_t length() {
		return length_;
	}
};

class OmniEdge {
	size_t length_;
public:
	OmniEdge(size_t length) : length_(length) {}

	size_t length() {
		return length_;
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
};

class Omnigraph : public AbstractConjugateGraph<OmniVertex, OmniEdge, OmniDataMaster> {
public:
	Omnigraph() : AbstractConjugateGraph(OmniDataMaster()) {

	}
};

}
#endif /* OMNI_GRAPH_HPP_ */
