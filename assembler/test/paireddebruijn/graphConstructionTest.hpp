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

void TestCheckUniqueWayRightEmpty() {
	edgesMap edges;
	Sequence seqA(repeatLetter('A', l));
	int res = checkUniqueWayRight(edges, 0, &seqA);
	ASSERT_EQUAL(0, res);
}

void TestCheckUniqueWayRightSingle() {
	edgesMap edges;
	Sequence seqA(repeatLetter('A', l));
	VertexPrototype vp(&seqA, 0);
	edges[1].push_back(&vp);
	cout << edges[0].size() << " " << edges[1].size() << " " << edges[2].size()
			<< endl;
	int res = checkUniqueWayRight(edges, 0, &seqA);
	ASSERT_EQUAL(1, res);
}

void TestCheckUniqueWayRightDoublePair() {
	edgesMap edges;
	Sequence seqA(repeatLetter('A', l));
	Sequence seqAG(repeatLetter('A', l) + "G");
	Sequence seqAC(repeatLetter('A', l) + "C");
	VertexPrototype vpAC(&seqAC, 0);
	VertexPrototype vpAG(&seqAG, 0);
	edges[1].push_back(&vpAC);
	edges[1].push_back(&vpAG);
	int res = checkUniqueWayRight(edges, 0, &seqA);
	ASSERT_EQUAL(0, res);
}

void TestCheckUniqueWayRightDoubleKmer() {
	edgesMap edges;
	Sequence seqA(repeatLetter('A', l));
	VertexPrototype vp1(&seqA, 0);
	VertexPrototype vp2(&seqA, 1);
	edges[0].push_back(&vp1);
	edges[1].push_back(&vp2);
	int res = checkUniqueWayRight(edges, 0, &seqA);
	ASSERT_EQUAL(0, res);
}

void TestCheckUniqueWayLeftEmpty() {
	edgesMap edges;
	Sequence seqA(repeatLetter('A', l));
	int res = checkUniqueWayLeft(edges, 0, &seqA);
	ASSERT_EQUAL(0, res);
}

void TestCheckUniqueWayLeftSingle() {
	edgesMap edges;
	Sequence seqA(repeatLetter('A', l));
	VertexPrototype vp(&seqA, 0);
	edges[(ll) 1 << (2 * (k - 1))].push_back(&vp);
	int res = checkUniqueWayLeft(edges, 0, &seqA);
	ASSERT_EQUAL(1, res);
}

void TestCheckUniqueWayLeftDoublePair() {
	edgesMap edges;
	Sequence seqA(repeatLetter('A', l));
	Sequence seqAG(repeatLetter('A', l) + "G");
	Sequence seqAC(repeatLetter('A', l) + "C");
	VertexPrototype vpAC(&seqAC, 0);
	VertexPrototype vpAG(&seqAG, 0);
	edges[(ll) 1 << (2 * (k - 1))].push_back(&vpAC);
	edges[(ll) 1 << (2 * (k - 1))].push_back(&vpAG);
	int res = checkUniqueWayLeft(edges, 0, &seqA);
	ASSERT_EQUAL(0, res);
}

void TestCheckUniqueWayLeftDoubleKmer() {
	edgesMap edges;
	Sequence seqA(repeatLetter('A', l));
	VertexPrototype vp1(&seqA, 0);
	VertexPrototype vp2(&seqA, 1);
	edges[(ll) 1 << (2 * (k - 1))].push_back(&vp1);
	edges[(ll) 2 << (2 * (k - 1))].push_back(&vp2);
	int res = checkUniqueWayLeft(edges, 0, &seqA);
	ASSERT_EQUAL(0, res);
}

cute::suite GraphConstructionSuite() {
	cute::suite s;
	s.push_back(CUTE(TestStoreVertexEmpty));
	s.push_back(CUTE(TestStoreVertexOneEntry));
	s.push_back(CUTE(TestStoreVertexTwoEntries));
	s.push_back(CUTE(TestStoreVertexTwoEntriesWithSameKmer));
	s.push_back(CUTE(TestCheckUniqueWayRightEmpty));
	s.push_back(CUTE(TestCheckUniqueWayRightSingle));
	s.push_back(CUTE(TestCheckUniqueWayRightDoublePair));
	s.push_back(CUTE(TestCheckUniqueWayRightDoubleKmer));
	s.push_back(CUTE(TestCheckUniqueWayLeftEmpty));
	s.push_back(CUTE(TestCheckUniqueWayLeftSingle));
	s.push_back(CUTE(TestCheckUniqueWayLeftDoublePair));
	s.push_back(CUTE(TestCheckUniqueWayLeftDoubleKmer));
	return s;
}

