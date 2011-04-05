#ifndef EDGEGRAPHTEST_HPP_
#define EDGEGRAPHTEST_HPP_
#include "edge_graph.hpp"

using namespace edge_graph;

void EmptyGraphTest() {
	EdgeGraph g(11);
	ASSERT_EQUAL(11, g.k());
	ASSERT_EQUAL(0u, g.vertices().size());
}

void  OneVertexGraphTest() {
	EdgeGraph g(11);
	g.AddVertex();
	ASSERT_EQUAL(2u, g.vertices().size());
	Vertex *v = *(g.vertices().begin());
	Vertex *rcv = v->complement();
	ASSERT(v != rcv);
	ASSERT_EQUAL(v, rcv->complement());
}

void  OneEdgeGraphTest() {
	EdgeGraph g(11);
	Vertex *v1 = g.AddVertex();
	Vertex *v2 = g.AddVertex();
	Edge *e = g.AddEdge(v1, v2, Sequence("AAAAAAAAAAAAAAAAA"));
	ASSERT_EQUAL(1u, v1->OutgoingEdgeCount());
	ASSERT_EQUAL(1u, v2->OutgoingEdgeCount());
	ASSERT_EQUAL(e, (*(v1->begin())));
	ASSERT_EQUAL(g.complementEdge(e), *(v2->begin()));
	ASSERT_EQUAL(e, g.complementEdge(g.complementEdge(e)));
	ASSERT_EQUAL(!(e->nucls()), g.complementEdge(e)->nucls());
}

cute::suite EdgeGraphSuite() {
	cute::suite s;
	s.push_back(CUTE(EmptyGraphTest));
	s.push_back(CUTE(OneVertexGraphTest));
	s.push_back(CUTE(OneEdgeGraphTest));
	return s;
}

#endif /* EDGEGRAPHTEST_HPP_ */
