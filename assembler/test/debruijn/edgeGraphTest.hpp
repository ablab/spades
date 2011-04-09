#ifndef EDGEGRAPHTEST_HPP_
#define EDGEGRAPHTEST_HPP_
#include "edge_graph.hpp"

namespace edge_graph {

void EmptyGraphTest() {
	EdgeGraph g(11);
	ASSERT_EQUAL(11, g.k());
	ASSERT_EQUAL(0u, g.size());
}

void OneVertexGraphTest() {
	EdgeGraph g(11);
	g.AddVertex();
	ASSERT_EQUAL(2u, g.size());
	Vertex *v = *(g.begin());
	Vertex *rcv = g.ComplementVertex(v);
	ASSERT(v != rcv);
	ASSERT_EQUAL(v, g.ComplementVertex(rcv));
}

pair<vector<Vertex *> , vector<Edge *> > createGraph(EdgeGraph &graph,
		int edgeNumber) {
	vector<Vertex *> v;
	vector<Edge *> e;
	v.push_back(graph.AddVertex());
	for (int i = 0; i < edgeNumber; i++) {
		v.push_back(graph.AddVertex());
		e.push_back(
				graph.AddEdge(v[v.size() - 2], v[v.size() - 1],
						Sequence("AAAAAAAAAAAAAAAAA")));
	}
	return make_pair(v, e);
}

void OneEdgeGraphTest() {
	EdgeGraph g(11);
	pair<vector<Vertex *> , vector<Edge *> > data = createGraph(g, 1);
	ASSERT_EQUAL(1u, g.OutgoingEdgeCount(data.first[0]));
	ASSERT_EQUAL(0u, g.OutgoingEdgeCount(data.first[1]));
	ASSERT_EQUAL(data.second[0], g.GetUniqueOutgoingEdge(data.first[0]));
	ASSERT_EQUAL(g.ComplementEdge(data.second[0]),
			g.GetUniqueOutgoingEdge(g.ComplementVertex(data.first[1])));
	ASSERT_EQUAL(data.second[0],
			g.ComplementEdge(g.ComplementEdge(data.second[0])));
	ASSERT_EQUAL(!(g.EdgeNucls(data.second[0])),
			g.EdgeNucls(g.ComplementEdge(data.second[0])));
}

void EdgeMethodsSimpleTest() {
	EdgeGraph g(11);
	pair<vector<Vertex *> , vector<Edge *> > data = createGraph(g, 2);
	ASSERT_EQUAL(data.second[0], &g.GetData(data.second[0]));
	ASSERT_EQUAL(
			true,
			g.AreLinkable(data.first[0], data.first[1],
					Sequence("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")));
	ASSERT_EQUAL(
			false,
			g.AreLinkable(data.first[0], data.first[1],
					Sequence("CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC")));
}

void VertexMethodsSimpleTest() {
	EdgeGraph g(11);
	pair<vector<Vertex *> , vector<Edge *> > data = createGraph(g, 2);
	ASSERT_EQUAL(data.second[0], g.GetUniqueIncomingEdge(data.first[1]));
	ASSERT_EQUAL(data.second[0], g.GetUniqueOutgoingEdge(data.first[0]));
	ASSERT_EQUAL(false, g.CanCompressVertex(data.first[0]));
	ASSERT_EQUAL(true, g.CanCompressVertex(data.first[1]));
	ASSERT_EQUAL(false, g.CheckUniqueIncomingEdge(data.first[0]));
	ASSERT_EQUAL(true, g.CheckUniqueIncomingEdge(data.first[1]));
	ASSERT_EQUAL(false, g.CheckUniqueOutgiongEdge(data.first[2]));
	ASSERT_EQUAL(true, g.CheckUniqueOutgiongEdge(data.first[1]));
	ASSERT_EQUAL(true, g.IsDeadEnd(data.first[2]));
	ASSERT_EQUAL(false, g.IsDeadEnd(data.first[1]));
	ASSERT_EQUAL(true, g.IsDeadStart(data.first[0]));
	ASSERT_EQUAL(false, g.IsDeadStart(data.first[1]));
}

void GraphMethodsSimpleTest() {
	EdgeGraph g(11);
	pair<vector<Vertex *> , vector<Edge *> > data = createGraph(g, 2);
	ASSERT_EQUAL(vector<GraphActionHandler<EdgeGraph> *> (), g.GetHandlers());
	GraphActionHandler<EdgeGraph> *handler =
			new GraphActionHandler<EdgeGraph> ();
	g.AddActionHandler(handler);
	vector<GraphActionHandler<EdgeGraph> *> handlers = g.GetHandlers();
	ASSERT_EQUAL(1u, handlers.size());
	ASSERT_EQUAL(handler, handlers[0]);
	g.RemoveActionHandler(handler);
	ASSERT_EQUAL(vector<GraphActionHandler<EdgeGraph> *> (), g.GetHandlers());
}

void SmartIteratorTest() {
	EdgeGraph g(11);
	pair<vector<Vertex *> , vector<Edge *> > data = createGraph(g, 4);
	size_t num = 0;
	set<Vertex *> visited;
	SmartVertexIterator<EdgeGraph> endIt = g.SmartVertexEnd();
	for (SmartVertexIterator<EdgeGraph> it = g.SmartVertexBegin(); g.SmartVertexEnd()
			!= it; ++it) {
		num++;
		visited.insert(*it);
	}
	ASSERT_EQUAL(num, data.first.size() * 2);
	for (size_t i = 0; i < data.first.size(); i++) {
		ASSERT(visited.find(data.first[i]) != visited.end());
		ASSERT(visited.find(g.ComplementVertex(data.first[i])) != visited.end());
	}
}

}

using namespace edge_graph ;
cute::suite EdgeGraphSuite() {
	cute::suite s;
	s.push_back(CUTE(EmptyGraphTest));
	s.push_back(CUTE(OneVertexGraphTest));
	s.push_back(CUTE(OneEdgeGraphTest));
	s.push_back(CUTE(EdgeMethodsSimpleTest));
	s.push_back(CUTE(VertexMethodsSimpleTest));
	s.push_back(CUTE(GraphMethodsSimpleTest));
	s.push_back(CUTE(SmartIteratorTest));
	return s;
}

#endif /* EDGEGRAPHTEST_HPP_ */
