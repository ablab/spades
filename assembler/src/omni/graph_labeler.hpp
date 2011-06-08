#ifndef GRAPH_LABELER_HPP_
#define GRAPH_LABELER_HPP_

namespace gvis {

template<class Graph>
class GraphLabeler {
protected:
public:
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
//	GraphPrinter(const string &name, ostream &out = cout) : name_(name), out_(out) {
//	}

	virtual ~GraphLabeler() {
	}

	virtual std::string label(VertexId vertexId) = 0;

	virtual std::string label(EdgeId edgeId) = 0;
};

template<class Graph>
class EmptyGraphLabeler : public GraphLabeler<Graph> {
	typedef GraphLabeler<Graph> super;
	typedef typename super::EdgeId EdgeId;
	typedef typename super::VertexId VertexId;
public:
	EmptyGraphLabeler() {}

	std::string label(VertexId vertexId) {
		return "";
	}

	std::string label(EdgeId edgeId) {
		return "";
	}
};

template<class Graph>
class StrGraphLabeler : public GraphLabeler<Graph> {
protected:
	typedef GraphLabeler<Graph> super;
	typedef typename super::EdgeId EdgeId;
	typedef typename super::VertexId VertexId;
	Graph& g_;
public:
	StrGraphLabeler(Graph& g) : g_(g) {}

	std::string label(VertexId vertexId) {
		return g_.str(vertexId);
	}

	std::string label(EdgeId edgeId) {
		return g_.str(edgeId);
	}
};

}

#endif /* GRAPH_LABELER_HPP_ */
