//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
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
	const Graph& g_;
	std::set<VertexId> vertices_;
	std::set<EdgeId> edges_;


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
				vertices_.insert(g_.conjugate(*it));
		}
	}

	void FillEdges() {
		for (auto v_it = vertices_.begin(); v_it != vertices_.end(); ++v_it) {
			const vector<EdgeId> edges = g_.OutgoingEdges(*v_it);
			TRACE("working with vertex " << *v_it);
			for (auto e_it = edges.begin(); e_it != edges.end(); ++e_it) {
				VertexId edge_end = g_.EdgeEnd(*e_it);
				TRACE(g_.coverage(*e_it) << " " << g_.length(*e_it));
				if (vertices_.count(edge_end) > 0) {
					edges_.insert(*e_it);
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
	GraphComponent(const Graph &g, VertexIt begin, VertexIt end) :
			g_(g) {
		Fill(begin, end);
	}

	//todo refactor and get rid of hack
	template<class VertexIt>
	GraphComponent(const Graph &g, VertexIt begin, VertexIt end,
			bool add_conjugate) : g_(g) {
		Fill(begin, end, add_conjugate);
	}

	//Full graph component
	GraphComponent(const Graph &g) : g_(g) {
		Fill(g.begin(), g.end());
	}

	//may be used for conjugate closure
	GraphComponent(const GraphComponent& component, bool add_conjugate) : g_(component.g_)
//		vertices_(component.vertices_.begin(), component.vertices_.end()),
//		edges_(component.edges_.begin(), component.edges_.end())
	{
		Fill(component.v_begin(), component.v_end(), add_conjugate);
	}

	const Graph& g() const {
		return g_;
	}

	size_t v_size() const {
		return vertices_.size();
	}

	size_t e_size() const {
		return edges_.size();
	}

	bool contains(EdgeId e) const {
		if (edges_.count(e) > 0)
			return true;
		return false;
	}

	bool contains(VertexId v) const {
		if (vertices_.count(v) > 0)
			return true;
		return false;
	}

	edge_iterator e_begin() const {
		return edges_.begin();
	}
	edge_iterator e_end() const {
		return edges_.end();
	}

	const std::set<EdgeId>& edges() {
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

};
}
