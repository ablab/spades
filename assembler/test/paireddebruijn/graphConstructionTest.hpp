#include "cute.h"
#include "graphConstruction.hpp"
#include "iostream"

using namespace gvis;

string repeatLetter(char c, int num) {
	string res;
	for (int i = 0; i < num; i++)
		res += c;
	return res;
}

string streamToString(istream &s) {
	string res;
	res.assign((std::istreambuf_iterator<char>(s)),
			std::istreambuf_iterator<char>());
	return res;
}

void TestStoreVertexEmpty() {
	resetVertexCount();
	string s;
	stringstream ss(s);
	GraphScheme<int> g("oppa", ss);
	verticesMap map;
	Sequence seq(repeatLetter('A', k));
	ASSERT_EQUAL(0, storeVertex(g, map, 0, &seq));
	g.output();
	ASSERT_EQUAL("digraph oppa {\n0 [label=" + repeatLetter('A', k - 1) + "]\n}\n", streamToString(ss));
}

void TestStoreVertexOneEntry() {
	resetVertexCount();
	string s;
	stringstream ss(s);
	GraphScheme<int> g("oppa", ss);
	verticesMap map;
	Sequence seq(repeatLetter('A', l));
	ASSERT_EQUAL(0, storeVertex(g, map, 0, &seq));
	ASSERT_EQUAL(0, storeVertex(g, map, 0, &seq));
	ASSERT_EQUAL(0, storeVertex(g, map, 0, &seq));
	ASSERT_EQUAL(0, storeVertex(g, map, 0, &seq));
	g.output();
	ASSERT_EQUAL("digraph oppa {\n0 [label=" + repeatLetter('A', k - 1) + "]\n}\n", streamToString(ss));
}

void TestStoreVertexTwoEntries() {
	resetVertexCount();
	string s;
	stringstream ss(s);
	GraphScheme<int> g("oppa", ss);
	verticesMap map;
	Sequence seqA(repeatLetter('A', l));
	Sequence seqB(repeatLetter('C', l));
	ASSERT_EQUAL(0, storeVertex(g, map, 0, &seqA));
	ASSERT_EQUAL(1, storeVertex(g, map, 1, &seqB));
	g.output();
	string result = "digraph oppa {\n0 [label=" + repeatLetter('A', k - 1)
			+ "]\n1 [label=" + repeatLetter('A', k - 2) + "C" + "]\n}\n";
	ASSERT_EQUAL(result, streamToString(ss));
}

void TestStoreVertexTwoEntriesWithSameKmer() {
	resetVertexCount();
	string s;
	stringstream ss(s);
	GraphScheme<int> g("oppa", ss);
	verticesMap map;
	Sequence seqA(repeatLetter('A', l));
	Sequence seqB(repeatLetter('C', l));
	ASSERT_EQUAL(0, storeVertex(g, map, 0, &seqA));
	ASSERT_EQUAL(1, storeVertex(g, map, 0, &seqB));
	g.output();
	string result = "digraph oppa {\n0 [label=" + repeatLetter('A', k - 1)
			+ "]\n1 [label=" + repeatLetter('A', k - 1) + "]\n}\n";
	ASSERT_EQUAL(result, streamToString(ss));
}

cute::suite GraphConstructionSuite() {
	cute::suite s;
	s.push_back(CUTE(TestStoreVertexEmpty));
	s.push_back(CUTE(TestStoreVertexOneEntry));
	s.push_back(CUTE(TestStoreVertexTwoEntries));
	s.push_back(CUTE(TestStoreVertexTwoEntriesWithSameKmer));
	return s;
}

