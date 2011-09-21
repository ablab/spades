#ifndef GRAPH_LABELER_HPP_
#define GRAPH_LABELER_HPP_

#include "simple_tools.hpp"

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


template<class Graph>
class StrCoverageGraphLabeler : public GraphLabeler<Graph> {
protected:
	typedef GraphLabeler<Graph> super;
	typedef typename super::EdgeId EdgeId;
	typedef typename super::VertexId VertexId;
	Graph& g_;
public:
	StrCoverageGraphLabeler(Graph& g) : g_(g) {}

	virtual std::string label(VertexId vertexId) const {
		return g_.str(vertexId);
	}

	virtual std::string label(EdgeId edgeId) const {
		double coverage = g_.coverage(edgeId);
		return g_.str(edgeId) + " {" + ToString(coverage) + "}";
	}
};

template<class Graph>
class LabelerList : public GraphLabeler<Graph> {
private:
	typedef GraphLabeler<Graph> super;
protected:
	typedef typename super::EdgeId EdgeId;
	typedef typename super::VertexId VertexId;
private:
	vector<GraphLabeler<Graph>*> list_;

	template<typename ElementId>
	string ConstructLabel(ElementId id) const {
		vector<string> to_print;
		for(size_t i = 0; i < list_.size(); i++) {
			string next = list_[i]->label(id);
			if(next.size() != 0) {
				to_print.push_back(next);
			}
		}
		string result = "";
		for(size_t i = 0; i < to_print.size(); i++) {
			result += to_print[i];
			if(i + 1 < to_print.size())
				result += "\\n";
		}
		return result;
	}

public:
	LabelerList() {
	}

	LabelerList(GraphLabeler<Graph> &labeler1, GraphLabeler<Graph> &labeler2) {
		AddLabeler(labeler1);
		AddLabeler(labeler2);
	}

	virtual ~LabelerList() {
	}

	void AddLabeler(GraphLabeler<Graph> &labeler) {
		list_.push_back(&labeler);
	}

	virtual std::string label(VertexId vertexId) const {
		return ConstructLabel<VertexId>(vertexId);
	}

	virtual std::string label(EdgeId edgeId) const {
		return ConstructLabel<EdgeId>(edgeId);
	}
};

}

#endif /* GRAPH_LABELER_HPP_ */
