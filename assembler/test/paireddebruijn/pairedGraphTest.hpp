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

pair<vector<VertexPrototype *> , vector<Edge *> > createCycleGraph(
		PairedGraph &graph, int length) {
	vector<VertexPrototype *> v = addVertices(graph, length);
	string s = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
	vector<Edge *> edges;
	for (int i = 0; i < length; i++) {
		Sequence *upper = new Sequence(s);
		Sequence *lower = new Sequence(s);
		Edge * e = new Edge(upper, lower, i, (i + 1) % length, 15, 0, 0);
		edges.push_back(e);
		graph.addEdge(e);
	}
	ASSERT_EQUAL(3, graph.VertexCount);
	ASSERT_EQUAL(3, graph.EdgeId);
	for (int i = 0; i < 3; i++) {
		ASSERT_EQUAL(i, edges[i]->EdgeId);
	}
	return make_pair(v, edges);
}

void TestThreeEdgesGraph() {
	PairedGraph g;
	pair<vector<VertexPrototype *>, vector<Edge *> > data = createCycleGraph(g, 3);
	vector<VertexPrototype *> vertices =  data.first;
	vector<Edge *> edges = data.second;
	VertexIterator *it = g.vertexIterator();
	for (int i = 0; i < 3; i++) {
		ASSERT_EQUAL(true, it->hasNext());
		ASSERT_EQUAL(vertices[i], it->next());
		ASSERT_EQUAL(1, g.rightDegree(i));
		ASSERT_EQUAL(1, g.leftDegree(i));
		ASSERT_EQUAL(edges[i], g.rightEdge(g.vertexList_[i], 0));
		ASSERT_EQUAL(edges[(i + 2) % 3], g.leftEdge(g.vertexList_[i], 0));
	}
}

void checkVertices(PairedGraph &g, vector<VertexPrototype *> &vertices) {
	VertexIterator *it = g.vertexIterator();
	for(int i = 0; i < vertices.size(); i++) {
		ASSERT_EQUAL(true, it->hasNext());
		ASSERT_EQUAL(vertices[i], it->next());
	}
	ASSERT_EQUAL(false, it->hasNext());
}

void TestRemoveEdge() {
	PairedGraph g;
	pair<vector<VertexPrototype *>, vector<Edge *> > data = createCycleGraph(g, 3);
	vector<VertexPrototype *> vertices =  data.first;
	vector<Edge *> edges = data.second;
	g.removeEdge(edges[0]);
	checkVertices(g, vertices);
	ASSERT_EQUAL(0, g.rightDegree(vertices[0]));
	ASSERT_EQUAL(0, g.leftDegree(vertices[1]));
	ASSERT_EQUAL(1, g.rightDegree(vertices[1]));
	ASSERT_EQUAL(1, g.leftDegree(vertices[2]));
	ASSERT_EQUAL(1, g.rightDegree(vertices[2]));
	ASSERT_EQUAL(1, g.leftDegree(vertices[0]));
}

void TestRemoveVertex() {
	PairedGraph g;
	pair<vector<VertexPrototype *>, vector<Edge *> > data = createCycleGraph(g, 3);
	vector<VertexPrototype *> vertices =  data.first;
	vector<Edge *> edges = data.second;
	g.removeVertex(vertices[0]);
	VertexIterator *it = g.vertexIterator();
	ASSERT_EQUAL(true, it->hasNext());
	ASSERT_EQUAL(vertices[1], it->next());
	ASSERT_EQUAL(true, it->hasNext());
	ASSERT_EQUAL(vertices[2], it->next());
	ASSERT_EQUAL(false, it->hasNext());
	ASSERT_EQUAL(0, g.leftDegree(vertices[1]));
	ASSERT_EQUAL(0, g.rightDegree(vertices[2]));
	ASSERT_EQUAL(1, g.rightDegree(vertices[1]));
	ASSERT_EQUAL(1, g.leftDegree(vertices[2]));
}

cute::suite PairedGraphSuite() {
	cute::suite s;
	s.push_back(CUTE(TestEmptyGraphCreate));
	s.push_back(CUTE(TestEmptyGraph));
	s.push_back(CUTE(TestOneVertexGraph));
	s.push_back(CUTE(TestNoEdgeGraph));
	s.push_back(CUTE(TestOneEdgeGraph));
	s.push_back(CUTE(TestThreeEdgesGraph));
	s.push_back(CUTE(TestRemoveEdge));
	s.push_back(CUTE(TestRemoveVertex));
	return s;
}

#endif /* PAIREDGRAPHTEST_HPP_ */
