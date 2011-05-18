#ifndef OMNI_TOOLS_HPP_
#define OMNI_TOOLS_HPP_

#include "omni_utils.hpp"

namespace omnigraph {

LOGGER("omg.graph");

template<class Graph>
class Compresser {
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;

	Graph &graph_;

	bool GoUniqueWay(EdgeId &e) {
		VertexId u = graph_.EdgeEnd(e);
		if (!graph_.CheckUniqueOutgiongEdge(u)
				|| !graph_.CheckUniqueIncomingEdge(u)) {
			return false;
		}
		e = graph_.GetUniqueOutgoingEdge(u);
		return true;
	}

public:
	Compresser(Graph &graph) :
		graph_(graph) {
	}

	bool CompressVertex(VertexId v) {
		if (!graph_.CheckUniqueOutgiongEdge(v) || !graph_.CheckUniqueIncomingEdge(v))
			return false;
		EdgeId e = graph_.GetUniqueOutgoingEdge(v);
		while (GoUniqueWay(e)) {
		}
		vector<EdgeId> mergeList;
		e = graph_.Complement(e);
		do {
			mergeList.push_back(e);
		} while (GoUniqueWay(e));
		graph_.MergePath(mergeList);
		return true;
	}

	void CompressAllVertices() {
		SmartVertexIterator<Graph> end = graph_.SmartVertexEnd();
		for (SmartVertexIterator<Graph> it = graph_.SmartVertexBegin(); it
				!= end; ++it) {
			VertexId v = *it;
			CompressVertex(v);
		}
	}
};

}

#endif /* OMNI_TOOLS_HPP_ */
