/**
 * @file graph.hpp
 *
 * @author vyahhi
 * @version 1.0
 *
 * @section LICENSE
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of
 * the License, or (at your option) any later version.
 *
 * @section DESCRIPTION
 *
 * Simple graph interface and implementation
 *
 */

#ifndef GRAPH_HPP_
#define GRAPH_HPP_

#include <map>

template<typename Node, typename Edge, typename NodeData, typename EdgeData>
class Graph {
public:
	void putNode(const Node &n, const NodeData &nd);
	void putEdge(const Edge &e, const EdgeData &ed);
	NodeData& getNode(const Node &n);
	EdgeData& getEdge(const Edge &e);
private:
	std::map<Node,NodeData> nodes_;
	std::map<Edge,EdgeData> edges_;
};


template<typename Node, typename Edge, typename NodeData, typename EdgeData>
void Graph<Node,Edge,NodeData,EdgeData>::putNode(const Node &n, const NodeData &nd) {
	nodes_[n] = nd;
}
template<typename Node, typename Edge, typename NodeData, typename EdgeData>
void Graph<Node,Edge,NodeData,EdgeData>::putEdge(const Edge &e, const EdgeData &ed) {
	edges_[e] = ed;
}

template<typename Node, typename Edge, typename NodeData, typename EdgeData>
NodeData& Graph<Node,Edge,NodeData,EdgeData>::getNode(const Node &n) {
	return nodes_[n];
}

template<typename Node, typename Edge, typename NodeData, typename EdgeData>
EdgeData& Graph<Node,Edge,NodeData,EdgeData>::getEdge(const Edge &e) {
	return edges_[e];
}

#endif /* GRAPH_HPP_ */
