#include "graphVisualizer.hpp"
#include "iostream"
#include "string"
#include "cute.h"

using namespace gvis;

void testOfflineEmpty() {
	stringstream ss;
	GraphScheme<int> g1("", ss);
	IGraphPrinter<int> *g = &g1;
	g->output();
	string tmp;
	tmp.assign((std::istreambuf_iterator<char>(ss)), std::istreambuf_iterator<
					char>());
	ASSERT_EQUAL("digraph  {\n}\n", tmp);
}

void testOfflineEmptyWithName() {
	stringstream ss;
	GraphScheme<int> g1("myName", ss);
	IGraphPrinter<int> *g = &g1;
	g->output();
	string tmp;
	tmp.assign((std::istreambuf_iterator<char>(ss)), std::istreambuf_iterator<
					char>());
	ASSERT_EQUAL("digraph myName {\n}\n", tmp);
}

void testOfflineSingleVertex() {
	stringstream ss;
	GraphScheme<int> g1("myName", ss);
	IGraphPrinter<int> *g = &g1;
	g->addVertex(0, "oppa");
	g->output();
	string tmp;
	tmp.assign((std::istreambuf_iterator<char>(ss)), std::istreambuf_iterator<
					char>());
	ASSERT_EQUAL("digraph myName {\n0 [label=\"oppa\",fillcolor=\"white\"]\n}\n", tmp);
}

void testOfflineSingleEdge() {
	stringstream ss;
	GraphScheme<int> g1("myName", ss);
	IGraphPrinter<int> *g = &g1;
	g->addEdge(0, 1, "oppa");
	g->output();
	string tmp;
	tmp.assign((std::istreambuf_iterator<char>(ss)), std::istreambuf_iterator<
					char>());
	ASSERT_EQUAL("digraph myName {\n0->1[label=\"oppa\",color=\"black\"]\n}\n", tmp);
}

cute::suite offlineGraphVisualizerSuite() {
	cute::suite s;
	s.push_back(CUTE(testOfflineEmpty));
	s.push_back(CUTE(testOfflineEmptyWithName));
	s.push_back(CUTE(testOfflineSingleVertex));
	s.push_back(CUTE(testOfflineSingleEdge));
	return s;
}
