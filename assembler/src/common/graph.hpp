/*
 * graph.hpp
 *
 *  Created on: 02.03.2011
 *      Author: vyahhi
 */

#ifndef GRAPH_HPP_
#define GRAPH_HPP_

#include <map>

// Simple graph interface and implementation
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

// Stupid Reference Implementation

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
