#include "graphVisualizer.hpp"
#include "iostream"
#include "string"
#include "cute.h"

using namespace gvis;

void testOnline() {
	string s;
	stringstream ss(s);
	OnlineGraphPrinter<int> g("", ss);
	string tmp;
	tmp.assign((std::istreambuf_iterator<char>(ss)), std::istreambuf_iterator<
					char>());
	ASSERT_EQUAL("digraph  {\n", tmp);
}

void testOnlineEmpty() {
	string s;
	stringstream ss(s);
	OnlineGraphPrinter<int> g("", ss);
	g.output();
	string tmp;
	tmp.assign((std::istreambuf_iterator<char>(ss)), std::istreambuf_iterator<
					char>());
	ASSERT_EQUAL("digraph  {\n}\n", tmp);
}

void testOnlineEmptyWithName() {
	string s;
	stringstream ss(s);
	OnlineGraphPrinter<int> g("myName", ss);
	g.output();
	string tmp;
	tmp.assign((std::istreambuf_iterator<char>(ss)), std::istreambuf_iterator<
					char>());
	ASSERT_EQUAL("digraph myName {\n}\n", tmp);
}

void testOnlineSingleVertex() {
	string s;
	stringstream ss(s);
	OnlineGraphPrinter<int> g("myName", ss);
	g.addVertex(0, "oppa");
	g.output();
	string tmp;
//	cout <<
	tmp.assign((std::istreambuf_iterator<char>(ss)), std::istreambuf_iterator<
					char>());
	ASSERT_EQUAL("digraph myName {\n0 [label=oppa]\n}\n", tmp);
}

void testOnlineSingleEdge() {
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

cute::suite onlineGraphVisualizerSuite() {
	cute::suite s;
	s.push_back(CUTE(testOnline));
	s.push_back(CUTE(testOnlineEmpty));
	s.push_back(CUTE(testOnlineEmptyWithName));
	s.push_back(CUTE(testOnlineSingleVertex));
	s.push_back(CUTE(testOnlineSingleEdge));
	return s;
}
