#ifndef OMNIGRAPH_HPP_
#define OMNIGRAPH_HPP_

#include "abstract_conjugate_graph.hpp"

namespace omnigraph {

class OmniDataMaster {
public:
	OmniDataMaster();
};

class OmniVertex {
	int length_;
public:
	OmniVertex(int length) : length_(length) {};
};

class OmniEdge {
	int length_;
public:
	OmniEdge(int length) : length_(length) {}
};

class Omnigraph : public AbstractConjugateGraph<OmniVertex, OmniEdge, OmniDataMaster> {
public:
	Omnigraph() : AbstractConjugateGraph(OmniDataMaster()) {

	}
};

}
#endif /* OMNI_GRAPH_HPP_ */
