//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

namespace omnigraph {

template<class Graph>
struct vertex_neighbour {
protected:
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
public:
	VertexId 	vertex;
	EdgeId 		edge;

	vertex_neighbour(VertexId new_vertex, EdgeId new_edge) :
		vertex(new_vertex),
		edge(new_edge) { }
};

template<class Graph>
class NeighbourIterator {
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
protected:
	const Graph &graph_;
	VertexId vertex_;
public:
	NeighbourIterator(const Graph &graph, VertexId vertex) :
		graph_(graph),
		vertex_(vertex) { }

	virtual bool HasNext() = 0;
	virtual vertex_neighbour<Graph> Next() = 0;
	virtual  ~NeighbourIterator() { }
};

template<class Graph>
class ForwardNeighbourIterator : public NeighbourIterator<Graph>{
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
	typedef typename VertexId::type::edge_const_iterator edge_const_iterator;

	pair<edge_const_iterator, edge_const_iterator> out_edges_;
public:
	ForwardNeighbourIterator(const Graph &graph, VertexId vertex) :
		NeighbourIterator<Graph>(graph, vertex),
		out_edges_(make_pair(graph.OutgoingEdges(vertex).begin(),
				graph.OutgoingEdges(vertex).end())) { }

	bool HasNext(){
		return out_edges_.first != out_edges_.second;
	}

	vertex_neighbour<Graph> Next() {
		vertex_neighbour<Graph> res(this->graph_.EdgeEnd(*out_edges_.first), *out_edges_.first);
		out_edges_.first++;
		return res;
	}
};

template<class Graph>
class BackwardNeighbourIterator : public NeighbourIterator<Graph>{
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
	typedef typename VertexId::type::edge_const_iterator edge_const_iterator;

	pair<edge_const_iterator, edge_const_iterator>  in_edges_;
public:
	BackwardNeighbourIterator(const Graph &graph, VertexId vertex) :
		NeighbourIterator<Graph>(graph, vertex),
		in_edges_(make_pair(graph.IncomingEdges(vertex).begin(),
				graph.IncomingEdges(vertex).end())) { }

	bool HasNext(){
		return in_edges_.first != in_edges_.second;
	}

	vertex_neighbour<Graph> Next() {
		vertex_neighbour<Graph> res(this->graph_.EdgeStart(*in_edges_.first), *in_edges_.first);
		in_edges_.first++;
		return res;
	}
};

template<class Graph>
class UnorientedNeighbourIterator : public NeighbourIterator<Graph>{
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
	typedef typename VertexId::type::edge_const_iterator edge_const_iterator;

	pair<edge_const_iterator, edge_const_iterator>  in_edges_;
	pair<edge_const_iterator, edge_const_iterator>  out_edges_;
public:
	UnorientedNeighbourIterator(const Graph &graph, VertexId vertex) :
		NeighbourIterator<Graph>(graph, vertex),
		in_edges_(make_pair(graph.IncomingEdges(vertex).begin(),
				graph.IncomingEdges(vertex).end())),
		out_edges_(make_pair(graph.OutgoingEdges(vertex).begin(),
				graph.OutgoingEdges(vertex).end())) { }

	bool HasNext(){
		return in_edges_.first != in_edges_.second;
	}

	// first all outgoing edges are visited
	// then all incoming
	vertex_neighbour<Graph> Next() {
		if(out_edges_.first != out_edges_.second){
			vertex_neighbour<Graph> res(this->graph_.EdgeEnd(*out_edges_.first), *out_edges_.first);
			out_edges_.first++;
			return res;
		}
		vertex_neighbour<Graph> res(this->graph_.EdgeStart(*in_edges_.first), *in_edges_.first);
		in_edges_.first++;
		return res;
	}
};

template<class Graph>
class NeighbourIteratorFactory {
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
protected:
	const Graph &graph_;
public:
	NeighbourIteratorFactory(const Graph &graph) : graph_(graph) { }
	virtual shared_ptr<NeighbourIterator<Graph> > CreateIterator(VertexId vertex) = 0;
	virtual ~NeighbourIteratorFactory() { }
};

template<class Graph>
class ForwardNeighbourIteratorFactory : public NeighbourIteratorFactory<Graph>{
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
public:
	ForwardNeighbourIteratorFactory(const Graph &graph) :
		NeighbourIteratorFactory<Graph>(graph) { }
	shared_ptr<NeighbourIterator<Graph> > CreateIterator(VertexId vertex){
		return shared_ptr<NeighbourIterator<Graph> >(new ForwardNeighbourIterator<Graph>(this->graph_, vertex));
	}
};

template<class Graph>
class BackwardNeighbourIteratorFactory : public NeighbourIteratorFactory<Graph>{
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
public:
	BackwardNeighbourIteratorFactory(const Graph &graph) :
		NeighbourIteratorFactory<Graph>(graph) { }
	shared_ptr<NeighbourIterator<Graph> > CreateIterator(VertexId vertex){
		return shared_ptr<NeighbourIterator<Graph> >(new BackwardNeighbourIterator<Graph>(this->graph_, vertex));
	}
};

template<class Graph>
class UnorientedNeighbourIteratorFactory : public NeighbourIteratorFactory<Graph>{
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
public:
	UnorientedNeighbourIteratorFactory(const Graph &graph) :
		NeighbourIteratorFactory<Graph>(graph) { }
	shared_ptr<NeighbourIterator<Graph> > CreateIterator(VertexId vertex){
		return shared_ptr<NeighbourIterator<Graph> >(new UnorientedNeighbourIterator<Graph>(this->graph_, vertex));
	}
};

}
