#ifndef GRAPH_LABELER_HPP_
#define GRAPH_LABELER_HPP_

namespace omnigraph {

/**
 * (Interface)
 * Provides string labels for vertices and edges of some graph.
 * Used with GraphPrinter to visualize graphs.
 */
template<class Graph>
class GraphLabeler {
protected:
public:
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;

	virtual ~GraphLabeler() {
	}

	virtual std::string label(VertexId vertexId) const  = 0;

	virtual std::string label(EdgeId edgeId) const = 0;
};

/**
 * Trivial implementation of GraphLabeler.
 * All labels are "".
 */
template<class Graph>
class EmptyGraphLabeler : public GraphLabeler<Graph> {
	typedef GraphLabeler<Graph> super;
	typedef typename super::EdgeId EdgeId;
	typedef typename super::VertexId VertexId;
public:
	EmptyGraphLabeler() {}

	std::string label(VertexId vertexId) const {
		return "";
	}

	std::string label(EdgeId edgeId) const {
		return "";
	}
};

/**
 * Implementation of GraphLabeler for Graphs that have methods
 * str(VertexId) and str(EdgeId), such as AbstractConjugateGraph.

 */
template<class Graph>
class StrGraphLabeler : public GraphLabeler<Graph> {
protected:
	typedef GraphLabeler<Graph> super;
	typedef typename super::EdgeId EdgeId;
	typedef typename super::VertexId VertexId;
	Graph& g_;
public:
	StrGraphLabeler(Graph& g) : g_(g) {}

	virtual std::string label(VertexId vertexId) const {
		return g_.str(vertexId);
	}

	virtual std::string label(EdgeId edgeId) const {
		return g_.str(edgeId);
	}
};

}

#endif /* GRAPH_LABELER_HPP_ */
