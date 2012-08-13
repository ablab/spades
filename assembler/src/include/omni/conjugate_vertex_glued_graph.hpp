//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************


/*
 * conjugate_vertexe_glued_graph.hpp
 *
 *  Created on: Sep 3, 2012
 *      Author: Alexander Opeykin (alexander.opeykin@gmail.com)
 */


#ifndef CONJUGATE_VERTEX_GLUED_GRAPH_HPP_
#define CONJUGATE_VERTEX_GLUED_GRAPH_HPP_

#include <boost/foreach.hpp>

#include "omni_utils.hpp"
#include "standard_base.hpp"


namespace omnigraph {

/*
 * This class is used as a graph wrapper for DevisibleTree class.
 * It decorates a couple of methods to hide differences between
 * vertex and it's conjugate.
 */

template <class Graph>
class ConjugateVertexGluedGraph {
public:

	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
	typedef typename set<VertexId>::iterator iterator;
	typedef typename set<VertexId>::const_iterator const_iterator;


	ConjugateVertexGluedGraph(Graph& graph)
			: graph_(graph) {

		BOOST_FOREACH(const VertexId& vertex, graph_) {
			vertices_.insert(GetMinWithConjugate(vertex));
		}
	}

	VertexId GetMinWithConjugate(const VertexId& vertex) const{
		VertexId conjugate = graph_.conjugate(vertex);

		return vertex < conjugate ? vertex : conjugate;
	}

	const_iterator begin() const {
		return vertices_.begin();
	}

	const_iterator end() const {
		return vertices_.end();
	}

	VertexId EdgeStart(const EdgeId& edge) const {
		return GetMinWithConjugate(graph_.EdgeStart(edge));
	}

	VertexId EdgeEnd(const EdgeId& edge) const {
		return GetMinWithConjugate(graph_.EdgeEnd(edge));
	}

	SmartEdgeIterator<Graph> SmartEdgeBegin() const {
		return SmartEdgeIterator<Graph>(graph_);
	}

private:
	Graph& graph_;
	set<VertexId> vertices_;
};

} // namespace omnigraph

#endif /* CONJUGATE_VERTEX_GLUED_GRAPH_HPP_ */
