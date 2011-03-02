/*
 * graph.hpp
 *
 *  Created on: 02.03.2011
 *      Author: vyahhi
 */

#ifndef GRAPH_HPP_
#define GRAPH_HPP_

template<typename Node, typename NodeData, typename Edge, typename EdgeData>
class Graph {
public:
	addEdge(const Edge &e, const EdgeData &ed);
	addNode(const Node &e, const NodeData &ed);
};

#endif /* GRAPH_HPP_ */
