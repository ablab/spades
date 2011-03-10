#include "graphVisualizer.hpp"
#include "iostream"
#include "string"
#include "cute.h"

using namespace gvis;

void testOnline() {
	stringstream ss;
	GraphPrinter<int> g1("", ss);
	IGraphPrinter<int> *g = &g1;
	string tmp;
	tmp.assign((std::istreambuf_iterator<char>(ss)),
			std::istreambuf_iterator<char>());
	ASSERT_EQUAL("digraph  {\n", tmp);
}

void testOnlineEmpty() {
	stringstream ss;
	GraphPrinter<int> g1("", ss);
	IGraphPrinter<int> *g = &g1;
	g->output();
	string tmp;
	tmp.assign((std::istreambuf_iterator<char>(ss)),
			std::istreambuf_iterator<char>());
	ASSERT_EQUAL("digraph  {\n}\n", tmp);
}

void testOnlineEmptyWithName() {
	stringstream ss;
	GraphPrinter<int> g1("myName", ss);
	IGraphPrinter<int> *g = &g1;
	g->output();
	string tmp;
	tmp.assign((std::istreambuf_iterator<char>(ss)),
			std::istreambuf_iterator<char>());
	ASSERT_EQUAL("digraph myName {\n}\n", tmp);
}

void testOnlineSingleVertex() {
	stringstream ss;
	GraphPrinter<int> g1("myName", ss);
	IGraphPrinter<int> *g = &g1;
	g->addVertex(0, "oppa");
	g->output();
	string tmp;
	//	cout <<
	tmp.assign((std::istreambuf_iterator<char>(ss)),
			std::istreambuf_iterator<char>());
	ASSERT_EQUAL("digraph myName {\n0 [label=\"oppa\",fillcolor=\"white\"]\n}\n", tmp);
}

void testOnlineSingleEdge() {
	stringstream ss;
	GraphPrinter<int> g1("myName", ss);
	IGraphPrinter<int> *g = &g1;
	g->addEdge(0, 1, "oppa");
	g->output();
	string tmp;
	tmp.assign((std::istreambuf_iterator<char>(ss)),
			std::istreambuf_iterator<char>());
	ASSERT_EQUAL("digraph myName {\n0->1[label=\"oppa\",color=\"black\"]\n}\n", tmp);
}

cute::suite onlineGraphVisualizerSuite() {
	cute::suite s;
	s.push_back(CUTE(testOnline));
	s.push_back(CUTE(testOnlineEmpty));
	s.push_back(CUTE(testOnlineEmptyWithName));
	s.push_back(CUTE(testOnlineSingleVertex));
	s.push_back(CUTE(testOnlineSingleEdge));
	return s;
}
