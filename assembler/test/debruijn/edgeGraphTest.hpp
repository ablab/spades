#ifndef EDGEGRAPHTEST_HPP_
#define EDGEGRAPHTEST_HPP_
#include "edge_graph.hpp"

using namespace edge_graph;

void EmptyGraphTest() {
	EdgeGraph g(11);
	ASSERT_EQUAL(11, g.k());
	ASSERT_EQUAL(0u, g.size());
}

void  OneVertexGraphTest() {
	EdgeGraph g(11);
	g.AddVertex();
	ASSERT_EQUAL(2u, g.size());
	Vertex *v = *(g.begin());
	Vertex *rcv = g.ComplementVertex(v);
	ASSERT(v != rcv);
	ASSERT_EQUAL(v, g.ComplementVertex(rcv));
}

void  OneEdgeGraphTest() {
	EdgeGraph g(11);
	Vertex *v1 = g.AddVertex();
	Vertex *v2 = g.AddVertex();
	Edge *e = g.AddEdge(v1, v2, Sequence("AAAAAAAAAAAAAAAAA"));
	ASSERT_EQUAL(1u, g.OutgoingEdgeCount(v1));
	ASSERT_EQUAL(0u, g.OutgoingEdgeCount(v2));
	ASSERT_EQUAL(e, g.GetUniqueOutgoingEdge(v1));
	ASSERT_EQUAL(g.ComplementEdge(e), g.GetUniqueOutgoingEdge(g.ComplementVertex(v2)));
	ASSERT_EQUAL(e, g.ComplementEdge(g.ComplementEdge(e)));
	ASSERT_EQUAL(!(g.EdgeNucls(e)), g.EdgeNucls(g.ComplementEdge(e)));
}

cute::suite EdgeGraphSuite() {
	cute::suite s;
//	s.push_back(CUTE(EmptyGraphTest));
//	s.push_back(CUTE(OneVertexGraphTest));
	s.push_back(CUTE(OneEdgeGraphTest));
	return s;
}

#endif /* EDGEGRAPHTEST_HPP_ */
