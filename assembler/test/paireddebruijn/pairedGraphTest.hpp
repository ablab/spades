#ifndef PAIREDGRAPHTEST_HPP_
#define PAIREDGRAPHTEST_HPP_
#include "cute.h"
#include "pairedGraph.hpp"

using namespace paired_assembler;

void TestEmptyGraphCreate() {
	PairedGraph g;
}

void TestEmptyGraph() {
	PairedGraph g;
	ASSERT_EQUAL(0, g.VertexCount);
	ASSERT_EQUAL(0, g.EdgeId);
	VertexIterator *it = g.vertexIterator();
	ASSERT_EQUAL(false, it->hasNext());
}

void TestOneVertexGraph() {
	PairedGraph g;
	VertexPrototype *v = new VertexPrototype(0,
			new Sequence("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"), 0);
	g.addVertex(v);
	ASSERT_EQUAL(1, g.VertexCount);
	ASSERT_EQUAL(0, g.EdgeId);
	VertexIterator *it = g.vertexIterator();
	ASSERT_EQUAL(false, it->hasNext());
}

vector<VertexPrototype *> addVertices(PairedGraph &graph, int number) {
	vector<VertexPrototype *> v;
	for (int i = 0; i < number; i++) {
		v.push_back(
				new VertexPrototype(i,
						new Sequence("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"), 0));
		graph.addVertex(v[i]);
		ASSERT_EQUAL(i + 1, graph.VertexCount);
	}
	ASSERT_EQUAL(0, graph.EdgeId);
	ASSERT_EQUAL(number, graph.VertexCount);
	return v;
}

void TestNoEdgeGraph() {
	PairedGraph g;
	vector<VertexPrototype *> vertices = addVertices(g, 10);
	VertexIterator *it = g.vertexIterator();
	ASSERT_EQUAL(false, it->hasNext());
}

void TestOneEdgeGraph() {
	PairedGraph g;
	vector<VertexPrototype *> vertices = addVertices(g, 2);
	string s = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
	Sequence *upper = new Sequence(s);
	Sequence *lower = new Sequence(s);
	Edge * e = new Edge(upper, lower, 0, 1, 15, 0, 0);
	g.addEdge(e);
	ASSERT_EQUAL(0, e->EdgeId);
	ASSERT_EQUAL(2, g.VertexCount);
	ASSERT_EQUAL(1, g.EdgeId);
	VertexIterator *it = g.vertexIterator();
	for (int i = 0; i < 2; i++) {
		ASSERT_EQUAL(true, it->hasNext());
		ASSERT_EQUAL(vertices[i], it->next());
	}
	ASSERT_EQUAL(1, g.rightDegree(vertices[0]));
	ASSERT_EQUAL(1, g.leftDegree(vertices[1]));
	ASSERT_EQUAL(0, g.rightDegree(vertices[1]));
	ASSERT_EQUAL(0, g.leftDegree(vertices[0]));
	ASSERT_EQUAL(e, g.rightEdge(vertices[0], 0));
	ASSERT_EQUAL(e, g.leftEdge(vertices[1], 0));
}

cute::suite PairedGraphSuite() {
	cute::suite s;
	s.push_back(CUTE(TestEmptyGraphCreate));
	s.push_back(CUTE(TestEmptyGraph));
	s.push_back(CUTE(TestOneVertexGraph));
	s.push_back(CUTE(TestNoEdgeGraph));
	s.push_back(CUTE(TestOneEdgeGraph));
	return s;
}

#endif /* PAIREDGRAPHTEST_HPP_ */
