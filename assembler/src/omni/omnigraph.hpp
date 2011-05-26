#ifndef OMNIGRAPH_HPP_
#define OMNIGRAPH_HPP_

#include "abstract_conjugate_graph.hpp"

namespace omnigraph {

class OmniVertex {
	int length_;
public:
	OmniVertex(int length) : length_(length) {}
};

class OmniEdge {
	int length_;
public:
	OmniEdge(int length) : length_(length) {}
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
