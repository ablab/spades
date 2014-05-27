//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

namespace omnigraph {

//todo make handler!!!
template<class Graph>
class GraphComponent {
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
	typedef typename std::set<VertexId>::const_iterator vertex_iterator;
	typedef typename std::set<EdgeId>::const_iterator edge_iterator;
	const Graph& graph_;
	std::set<VertexId> vertices_;
	std::set<EdgeId> edges_;
	string name_;


	template<class VertexIt>
	void FillVertices(VertexIt begin, VertexIt end) {
		for (auto it = begin; it != end; ++it) {
			vertices_.insert(*it);
		}
	}

	template<class VertexIt>
	void FillVertices(VertexIt begin, VertexIt end, bool add_conjugate) {
		for (auto it = begin; it != end; ++it) {
			vertices_.insert(*it);
			if (add_conjugate)
				vertices_.insert(graph_.conjugate(*it));
		}
	}

	void FillEdges() {
		for (auto v_it = vertices_.begin(); v_it != vertices_.end(); ++v_it) {
			TRACE("working with vertex " << graph_.str(*v_it));
			FOREACH (EdgeId e,  graph_.OutgoingEdges(*v_it)) {
				VertexId edge_end = graph_.EdgeEnd(e);
				TRACE(graph_.coverage(e) << " " << graph_.length(e));
				if (vertices_.count(edge_end) > 0) {
					edges_.insert(e);
					TRACE("Edge added");
				}
			}
		}
	}

	template<class VertexIt>
	void Fill(VertexIt begin, VertexIt end) {
		FillVertices(begin, end);
		FillEdges();
	}

	template<class VertexIt>
	void Fill(VertexIt begin, VertexIt end, bool add_conjugate) {
		FillVertices(begin, end, add_conjugate);
		FillEdges();
	}

public:
	template<class VertexIt>
	GraphComponent(const Graph &g, VertexIt begin, VertexIt end, const string &name = "") :
		graph_(g), name_(name) {
		Fill(begin, end);
	}

	//todo refactor and get rid of hack
	template<class VertexIt>
	GraphComponent(const Graph &g, VertexIt begin, VertexIt end,
			bool add_conjugate, const string &name = "") : graph_(g), name_(name) {
		Fill(begin, end, add_conjugate);
	}

	//Full graph component
	GraphComponent(const Graph &g, const string &name = "") : graph_(g), name_(name) {
		Fill(g.begin(), g.end());
	}

	//may be used for conjugate closure
	GraphComponent(const GraphComponent& component, bool add_conjugate, const string &name = "") : graph_(component.graph_), name_(name)
//		vertices_(component.vertices_.begin(), component.vertices_.end()),
//		edges_(component.edges_.begin(), component.edges_.end())
	{
		Fill(component.v_begin(), component.v_end(), add_conjugate);
	}

	GraphComponent<Graph> &operator=(const GraphComponent<Graph> &that) {
	    VERIFY(&this->graph_ == &that.graph_);
	    this->vertices_ = that.vertices_;
        this->edges_ = that.edges_;
	    this->name_ = that.name_;
	    return *this;
	}

	const Graph& g() const {
		return graph_;
	}

	string name() const {
	    return name_;
	}

	size_t v_size() const {
		return vertices_.size();
	}

	size_t e_size() const {
		return edges_.size();
	}

	bool contains(EdgeId e) const {
		return edges_.count(e) > 0;
	}

	bool contains(VertexId v) const {
		return vertices_.count(v) > 0;
	}

	edge_iterator e_begin() const {
		return edges_.begin();
	}
	edge_iterator e_end() const {
		return edges_.end();
	}

	const std::set<EdgeId>& edges() const {
		return edges_;
	}

	const std::set<VertexId>& vertices() const{
		return vertices_;
	}

	vertex_iterator v_begin() const {
		return vertices_.begin();
	}
	vertex_iterator v_end() const {
		return vertices_.end();
	}

	bool IsBorder(VertexId v) const {
		if(vertices_.count(v) == 0)
			return false;
		FOREACH (EdgeId e, graph_.AdjacentEdges(v)) {
			if (vertices_.count(graph_.EdgeStart(e)) == 0
					|| vertices_.count(graph_.EdgeEnd(e)) == 0) {
				return true;
			}
		}
		return false;
	}

};
}
