#include "graphVisualizer.hpp"
#include "iostream"
#include "string"
#include "cute.h"

using namespace gvis;

void testOfflineEmpty() {
	string s;
	stringstream ss(s);
	GraphScheme<int> g("", ss);
	g.output();
	string tmp;
	tmp.assign((std::istreambuf_iterator<char>(ss)), std::istreambuf_iterator<
					char>());
	ASSERT_EQUAL("digraph  {\n}\n", tmp);
}

void testOfflineEmptyWithName() {
	string s;
	stringstream ss(s);
	GraphScheme<int> g("myName", ss);
	g.output();
	string tmp;
	tmp.assign((std::istreambuf_iterator<char>(ss)), std::istreambuf_iterator<
					char>());
	ASSERT_EQUAL("digraph myName {\n}\n", tmp);
}

void testOfflineSingleVertex() {
	string s;
	stringstream ss(s);
	GraphScheme<int> g("myName", ss);
	g.addVertex(0, "oppa");
	g.output();
	string tmp;
	tmp.assign((std::istreambuf_iterator<char>(ss)), std::istreambuf_iterator<
					char>());
	ASSERT_EQUAL("digraph myName {\n0 [label=oppa]\n}\n", tmp);
}

void testOfflineSingleEdge() {
	string s;
	stringstream ss(s);
	OnlineGraphPrinter<int> g("myName", ss);
	g.addEdge(0, 1, "oppa");
	g.output();
	string tmp;
	tmp.assign((std::istreambuf_iterator<char>(ss)), std::istreambuf_iterator<
					char>());
	ASSERT_EQUAL("digraph myName {\n0->1 [label=oppa]\n}\n", tmp);
}

cute::suite offlineGraphVisualizerSuite() {
	cute::suite s;
	s.push_back(CUTE(testOfflineEmpty));
	s.push_back(CUTE(testOfflineEmptyWithName));
	s.push_back(CUTE(testOfflineSingleVertex));
	s.push_back(CUTE(testOfflineSingleEdge));
	return s;
}
