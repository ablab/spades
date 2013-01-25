//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef GRAPH_LABELER_HPP_
#define GRAPH_LABELER_HPP_

#include "simple_tools.hpp"
#include "edges_position_handler.hpp"

namespace omnigraph {

/**
 * (Interface)
 * Provides string labels for vertices and edges of some graph.
 * Used with GraphPrinter to visualize graphs.
 */
template<class Graph>
class GraphLabeler {
public:
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;

	virtual ~GraphLabeler() {
	}

	virtual string label(VertexId v) const = 0;

	virtual string label(EdgeId e) const = 0;

};

template<class Graph>
class MapGraphLabeler {
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
	map<EdgeId, string> edge_map_;
	map<VertexId, string> vertex_map_;

public:

	virtual string label(VertexId v) const {
		auto it = vertex_map_.find(v);
		if (it == vertex_map_.end())
			return "";
		else
			return it->second;
	}

	virtual string label(EdgeId e) const {
		auto it = edge_map_.find(e);
		if (it == edge_map_.end())
			return "";
		else
			return it->second;
	}

};

template<class Graph>
class AbstractGraphLabeler: public GraphLabeler<Graph> {
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
	const Graph& g_;
protected:
	AbstractGraphLabeler(const Graph& g): g_(g) {

	}

	const Graph& graph() const {
		return g_;
	}

public:
	/*virtual*/ std::string label(VertexId v) const {
		return "";
	}

	/*virtual*/ std::string label(EdgeId e) const {
		return "";
	}

};

/**
 * Trivial implementation of GraphLabeler.
 * All labels are "".
 */
template<class Graph>
class EmptyGraphLabeler : public GraphLabeler<Graph> {
	typedef GraphLabeler<Graph> base;
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
public:
	EmptyGraphLabeler() {}

	std::string label(VertexId v) const {
		return "";
	}

	std::string label(EdgeId e) const {
		return "";
	}
};

/**
 * Implementation of GraphLabeler for Graphs that have methods
 * str(VertexId) and str(EdgeId), such as AbstractConjugateGraph.
 */
template<class Graph>
class StrGraphLabeler : public AbstractGraphLabeler<Graph> {
	typedef AbstractGraphLabeler<Graph> base;
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
public:
	StrGraphLabeler(const Graph& g) : base(g) {}

	/*virtual*/ std::string label(VertexId v) const {
		return this->graph().str(v);
	}

	/*virtual*/ std::string label(EdgeId e) const {
		return this->graph().str(e);
	}

	/*virtual*/ ~StrGraphLabeler() {

	}
};

template <class Graph>
auto_ptr<GraphLabeler<Graph>> StrGraphLabelerInstance(const Graph& g) {
	return auto_ptr<GraphLabeler<Graph>>(new StrGraphLabeler<Graph>(g));
}

template<class Graph>
class LengthIdGraphLabeler : public StrGraphLabeler<Graph> {
	typedef StrGraphLabeler<Graph> base;
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
public:
	LengthIdGraphLabeler(const Graph& g) : base(g) {}

	/*virtual*/ std::string label(EdgeId e) const {
		std::stringstream ss;
		ss << this->graph().length(e) << " (id: " << this->graph().int_id(e) << ")";
		return ss.str();
	}

	/*virtual*/ ~LengthIdGraphLabeler() {

	}
};

template<class Graph>
class CoverageGraphLabeler : public AbstractGraphLabeler<Graph> {
	typedef AbstractGraphLabeler<Graph> base;
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
public:
	CoverageGraphLabeler(const Graph& g) : base(g) {}

	virtual std::string label(EdgeId e) const {
		double coverage = this->graph().coverage(e);
		return " {Cov:" + ToString(coverage) + "}";
	}
};

template<class Graph>
class CompositeLabeler : public GraphLabeler<Graph> {
private:
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
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
	CompositeLabeler() {
	}

	CompositeLabeler(GraphLabeler<Graph> &labeler1, GraphLabeler<Graph> &labeler2, GraphLabeler<Graph> &labeler3) {
		AddLabeler(labeler1);
		AddLabeler(labeler2);
		AddLabeler(labeler3);
	}

	CompositeLabeler(GraphLabeler<Graph> &labeler1, GraphLabeler<Graph> &labeler2) {
		AddLabeler(labeler1);
		AddLabeler(labeler2);
	}

	virtual ~CompositeLabeler() {
	}

	void AddLabeler(GraphLabeler<Graph> &labeler) {
		list_.push_back(&labeler);
	}

	virtual string label(VertexId vertexId) const {
		return ConstructLabel<VertexId>(vertexId);
	}

	virtual string label(EdgeId edgeId) const {
		return ConstructLabel<EdgeId>(edgeId);
	}
};

template<class Graph>
class EdgePosGraphLabeler: public AbstractGraphLabeler<Graph> {
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
public:
	const EdgesPositionHandler<Graph>& edge_pos_;

	EdgePosGraphLabeler(const Graph& g, const EdgesPositionHandler<Graph>& edge_pos) :
		AbstractGraphLabeler<Graph>(g), edge_pos_(edge_pos) {
	}

	virtual std::string label(EdgeId edgeId) const {
		return "Positions: " + edge_pos_.str(edgeId);
	}

	virtual ~EdgePosGraphLabeler() {
		TRACE("~EdgePosGraphLabeler");
	}

};

}

#endif /* GRAPH_LABELER_HPP_ */
